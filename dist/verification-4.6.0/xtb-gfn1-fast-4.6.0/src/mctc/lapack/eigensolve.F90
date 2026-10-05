! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
! Copyright (C) 2020, NVIDIA CORPORATION. All rights reserved.
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Wrapper for eigensolver routines
module xtb_mctc_lapack_eigensolve
   use xtb_mctc_accuracy, only : sp, dp
   use, intrinsic :: iso_fortran_env, only : int64
   use xtb_mctc_blas_level3, only : blas_trsm
   use xtb_mctc_lapack_geneigval, only : lapack_sygvd
   use xtb_mctc_lapack_stdeigval, only : lapack_syevd, lapack_syevr
   use xtb_mctc_lapack_gst, only : lapack_sygst
   use xtb_mctc_lapack_trf, only : mctc_potrf
   use xtb_type_environment, only : TEnvironment
#ifdef USE_CUSOLVER
   use xtb_mctc_global
   use cusolverDn
#endif
   implicit none
   private

   public :: TEigenSolver, init
   integer, parameter, public :: full_eig_backend_auto = 0
   integer, parameter, public :: full_eig_backend_syevd = 1
   integer, parameter, public :: full_eig_backend_syevr = 2


   type :: TEigenSolver
      private
      integer :: n
      integer, allocatable :: iwork(:)
      real(sp), allocatable :: swork(:)
      real(sp), allocatable :: sbmat(:, :)
      real(dp), allocatable :: dwork(:)
      real(dp), allocatable :: dbmat(:, :)
      ! GFN1-fast 2.2.4: reusable indexed-eigensolver workspaces. Avoid
      ! allocate/free on every SCC partial solve and use LAPACK *SYEVR.
      integer, allocatable :: partial_iwork(:), partial_isuppz(:)
      real(sp), allocatable :: spartial_work(:), sz(:, :)
      real(dp), allocatable :: dpartial_work(:), dz(:, :)
      ! GFN1-fast 2.2.4: phase-resolved timing of the cached-factor full solve.
      real(dp) :: full_reduce_time = 0.0_dp
      real(dp) :: full_diag_time = 0.0_dp
      real(dp) :: full_back_time = 0.0_dp
      integer :: full_profile_calls = 0
      ! GFN1-fast 2.2.6: full-spectrum backend selection.  AUTO benchmarks
      ! SYEVD and SYEVR on the same reduced matrix once, validates eigenvalues,
      ! then keeps the faster backend for the remainder of this solver lifetime.
      integer :: full_backend_requested = full_eig_backend_auto
      integer :: full_backend_selected = full_eig_backend_syevd
      logical :: full_backend_autotune = .true.
      logical :: full_backend_tuned = .false.
      real(dp) :: full_backend_max_ratio = 0.98_dp
      integer :: full_backend_autotune_min_n = 384
      real(dp) :: full_syevd_time = 0.0_dp
      real(dp) :: full_syevr_time = 0.0_dp
      integer :: full_syevd_calls = 0
      integer :: full_syevr_calls = 0
      integer :: full_backend_trials = 0
      integer :: full_backend_fallbacks = 0
      real(sp), allocatable :: sfull_alt(:, :), sfull_eval(:)
      real(dp), allocatable :: dfull_alt(:, :), dfull_eval(:)
#ifdef USE_CUSOLVER
      integer :: lwork
#endif
   contains
      generic :: solve => sgen_solve, dgen_solve
      procedure :: sgen_solve => mctc_ssygvd
      procedure :: dgen_solve => mctc_dsygvd
      ! GFN1-fast 0.6.0: reuse the Cholesky factorisation of the overlap
      ! matrix prepared by init(). This is mathematically equivalent to
      ! solving the generalized symmetric eigenproblem from scratch, but
      ! avoids refactorizing the invariant overlap matrix every SCC step.
      generic :: fact_solve => sfact_solve, dfact_solve
      procedure :: sfact_solve => mctc_ssygvd_factorized
      procedure :: dfact_solve => mctc_dsygvd_factorized
      ! GFN1-fast 0.8.0: exact indexed partial eigensolver for the lowest
      ! eigenpairs, with one extra boundary-probe eigenvalue. The same cached
      ! overlap Cholesky factor is reused.
      generic :: partial_fact_solve => spartial_fact_solve, dpartial_fact_solve
      procedure :: spartial_fact_solve => mctc_ssyevr_factorized
      procedure :: dpartial_fact_solve => mctc_dsyevr_factorized
      procedure :: reset_full_profile => resetFullEigenProfile
      procedure :: get_full_profile => getFullEigenProfile
      procedure :: configure_full_backend => configureFullEigenBackend
      procedure :: get_full_backend_profile => getFullEigenBackendProfile
   end type TEigenSolver


   interface init
      module procedure :: initSEigenSolver
      module procedure :: initDEigenSolver
   end interface init


contains


subroutine initSEigenSolver(self, env, bmat)
   character(len=*), parameter :: source = 'mctc_lapack_sygvd'
   class(TEigenSolver), intent(out) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: bmat(:, :)

   self%n = size(bmat, 1)

   allocate(self%swork(1 + 6*self%n + 2*self%n**2))
   allocate(self%iwork(3 + 5*self%n))
   allocate(self%spartial_work(max(1,26*self%n)))
   allocate(self%partial_iwork(max(1,10*self%n)))
   allocate(self%partial_isuppz(max(2,2*self%n)))

   self%sbmat = bmat
   ! Check for Cholesky factorisation
   call mctc_potrf(env, self%sbmat)

end subroutine initSEigenSolver


subroutine initDEigenSolver(self, env, bmat)
   character(len=*), parameter :: source = 'mctc_lapack_sygvd'
   class(TEigenSolver), intent(out) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: bmat(:, :)
#ifdef USE_CUSOLVER
   integer :: istat, lwork
   ! dummy is only a dummy argument used to query the workspace size needed
   ! for cuSolverDnDsygvd -- it is okay to pass an empty array to cuSolverDnDsygvd_bufferSize
   real(dp) :: dummy(:) 
#endif

   self%n = size(bmat, 1)

#ifdef USE_CUSOLVER
   istat = cusolverDnDsygvd_bufferSize(cusolverDnH, CUSOLVER_EIG_TYPE_1, &
     CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, self%n, dummy,    &
     self%n, dummy, self%n, dummy, lwork)
   if (istat /= 0) then
      call env%error("failed to get dygvd buffer size", source)
   end if

   self%lwork = lwork
   allocate(self%dwork(lwork))
#else
   allocate(self%dwork(1 + 6*self%n + 2*self%n**2))
   allocate(self%iwork(3 + 5*self%n))
#endif
   allocate(self%dpartial_work(max(1,26*self%n)))
   allocate(self%partial_iwork(max(1,10*self%n)))
   allocate(self%partial_isuppz(max(2,2*self%n)))

   self%dbmat = bmat
   ! Check for Cholesky factorisation
   call mctc_potrf(env, self%dbmat)

end subroutine initDEigenSolver


subroutine mctc_ssygvd(self, env, amat, bmat, eval)
   character(len=*), parameter :: source = 'mctc_lapack_sygvd'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(out) :: eval(:)
   integer :: info, lswork, liwork

   self%sbmat(:, :) = bmat

   lswork = size(self%swork)
   liwork = size(self%iwork)
   call lapack_sygvd(1, 'v', 'u', self%n, amat, self%n, self%sbmat, self%n, eval, &
      & self%swork, lswork, self%iwork, liwork, info)

   if (info /= 0) then
      call env%error("Failed to solve eigenvalue problem", source)
   end if

end subroutine mctc_ssygvd


subroutine mctc_dsygvd(self, env, amat, bmat, eval)
   character(len=*), parameter :: source = 'mctc_lapack_sygvd'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(out) :: eval(:)
   integer :: info, ldwork, liwork
#ifdef USE_CUSOLVER
   integer :: istat
#endif

   self%dbmat(:, :) = bmat

#ifdef USE_CUSOLVER
   !$acc enter data copyin(amat, self%dbmat, eval, info) create(self%dwork)

   !$acc host_data use_device(amat, self%dbmat, eval, self%dwork, info)
   istat = cusolverDnDsygvd(cusolverDnH, CUSOLVER_EIG_TYPE_1, &
     CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, self%n, amat, self%n, &
     self%dbmat, self%n, eval, self%dwork, self%lwork, info)
   !$acc end host_data

   !$acc exit data copyout(amat, self%dbmat, eval, info) delete(self%dwork)

   if (istat /= 0) then
      call env%error("cuSovlerDnDsygvd failed", source)
   end if
#else
   ldwork = size(self%dwork)
   liwork = size(self%iwork)
   call lapack_sygvd(1, 'v', 'u', self%n, amat, self%n, self%dbmat, self%n, eval, &
      & self%dwork, ldwork, self%iwork, liwork, info)
#endif

   if (info /= 0) then
      call env%error("Failed to solve eigenvalue problem", source)
   end if

end subroutine mctc_dsygvd



subroutine mctc_ssygvd_factorized(self, env, amat, eval)
   character(len=*), parameter :: source = 'mctc_lapack_ssygvd_factorized'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   integer :: info

   ! Reduce A x = lambda B x to a standard symmetric problem using the
   ! Cholesky factor of B cached in self%sbmat by initSEigenSolver().
   block
      real(dp) :: t0
      t0 = eigensolverWallTime()
      call lapack_sygst(1, 'u', self%n, amat, self%n, self%sbmat, self%n, info)
      self%full_reduce_time = self%full_reduce_time + eigensolverWallTime() - t0
   end block
   if (info /= 0) then
      call env%error("Failed to reduce eigenvalue problem", source)
      return
   end if

   block
      real(dp) :: t0
      t0 = eigensolverWallTime()
      call diagonalizeFullSP(self, amat, eval, info)
      self%full_diag_time = self%full_diag_time + eigensolverWallTime() - t0
   end block
   if (info /= 0) then
      call env%error("Failed to compute eigenvalues and eigenvectors", source)
      return
   end if

   ! Backtransform eigenvectors to the original non-orthogonal AO basis.
   block
      real(dp) :: t0
      t0 = eigensolverWallTime()
      call blas_trsm('l', 'u', 'n', 'n', self%n, self%n, 1.0_sp, self%sbmat, &
         & self%n, amat, self%n)
      self%full_back_time = self%full_back_time + eigensolverWallTime() - t0
   end block
   self%full_profile_calls = self%full_profile_calls + 1
end subroutine mctc_ssygvd_factorized


subroutine mctc_dsygvd_factorized(self, env, amat, eval)
   character(len=*), parameter :: source = 'mctc_lapack_dsygvd_factorized'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   integer :: info

   block
      real(dp) :: t0
      t0 = eigensolverWallTime()
      call lapack_sygst(1, 'u', self%n, amat, self%n, self%dbmat, self%n, info)
      self%full_reduce_time = self%full_reduce_time + eigensolverWallTime() - t0
   end block
   if (info /= 0) then
      call env%error("Failed to reduce eigenvalue problem", source)
      return
   end if

   block
      real(dp) :: t0
      t0 = eigensolverWallTime()
      call diagonalizeFullDP(self, amat, eval, info)
      self%full_diag_time = self%full_diag_time + eigensolverWallTime() - t0
   end block
   if (info /= 0) then
      call env%error("Failed to compute eigenvalues and eigenvectors", source)
      return
   end if

   block
      real(dp) :: t0
      t0 = eigensolverWallTime()
      call blas_trsm('l', 'u', 'n', 'n', self%n, self%n, 1.0_dp, self%dbmat, &
         & self%n, amat, self%n)
      self%full_back_time = self%full_back_time + eigensolverWallTime() - t0
   end block
   self%full_profile_calls = self%full_profile_calls + 1
end subroutine mctc_dsygvd_factorized


! GFN1-fast 2.2.6: full-spectrum backend dispatcher. AUTO performs one
! same-matrix SYEVD/SYEVR audit.  The reference SYEVD result is retained unless
! SYEVR succeeds, agrees in all eigenvalues to LAPACK-scale tolerance, and is
! faster by the configured margin.  The one-time audit therefore cannot change
! the SCC state when SYEVR is rejected.
subroutine diagonalizeFullSP(self, amat, eval, info)
   class(TEigenSolver), intent(inout) :: self
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   integer, intent(out) :: info
   integer :: info_alt, mcalc
   real(dp) :: t_syevd, t_syevr
   real(sp) :: scale, tol
   logical :: eig_ok

   select case (self%full_backend_requested)
   case (full_eig_backend_syevd)
      call runFullSYEVDSP(self, amat, eval, info, t_syevd)
      self%full_backend_selected = full_eig_backend_syevd
      self%full_backend_tuned = .true.
   case (full_eig_backend_syevr)
      call runFullSYEVRSP(self, amat, eval, info, t_syevr)
      if (info /= 0) then
         self%full_backend_fallbacks = self%full_backend_fallbacks + 1
         call runFullSYEVDSP(self, amat, eval, info, t_syevd)
         self%full_backend_selected = full_eig_backend_syevd
      else
         self%full_backend_selected = full_eig_backend_syevr
      end if
      self%full_backend_tuned = .true.
   case default
      if (.not.self%full_backend_autotune) then
         call runFullSYEVDSP(self, amat, eval, info, t_syevd)
         self%full_backend_selected = full_eig_backend_syevd
         self%full_backend_tuned = .true.
      else if (.not.self%full_backend_tuned) then
         call ensureFullWorkspaceSP(self)
         self%sfull_alt(:, :) = amat(:, :)
         call runFullSYEVDSP(self, amat, eval, info, t_syevd)
         if (info /= 0) return
         call runFullSYEVRRawSP(self, self%sfull_alt, self%sfull_eval, info_alt, mcalc, t_syevr)
         self%full_backend_trials = self%full_backend_trials + 1
         scale = max(1.0_sp, maxval(abs(eval)))
         tol = 8192.0_sp*epsilon(1.0_sp)*scale
         eig_ok = .false.
         if (info_alt == 0 .and. mcalc == self%n) then
            eig_ok = maxval(abs(self%sfull_eval-eval)) <= tol
         end if
         if (eig_ok .and. t_syevr <= self%full_backend_max_ratio*t_syevd) then
            eval(:) = self%sfull_eval(:)
            amat(:, :) = self%sz(:, :)
            self%full_backend_selected = full_eig_backend_syevr
         else
            if (.not.eig_ok) self%full_backend_fallbacks = self%full_backend_fallbacks + 1
            self%full_backend_selected = full_eig_backend_syevd
         end if
         self%full_backend_tuned = .true.
      else if (self%full_backend_selected == full_eig_backend_syevr) then
         call runFullSYEVRSP(self, amat, eval, info, t_syevr)
         if (info /= 0) then
            self%full_backend_fallbacks = self%full_backend_fallbacks + 1
            call runFullSYEVDSP(self, amat, eval, info, t_syevd)
            self%full_backend_selected = full_eig_backend_syevd
         end if
      else
         call runFullSYEVDSP(self, amat, eval, info, t_syevd)
      end if
   end select
end subroutine diagonalizeFullSP


subroutine diagonalizeFullDP(self, amat, eval, info)
   class(TEigenSolver), intent(inout) :: self
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   integer, intent(out) :: info
   integer :: info_alt, mcalc
   real(dp) :: t_syevd, t_syevr, scale, tol
   logical :: eig_ok

   select case (self%full_backend_requested)
   case (full_eig_backend_syevd)
      call runFullSYEVDDP(self, amat, eval, info, t_syevd)
      self%full_backend_selected = full_eig_backend_syevd
      self%full_backend_tuned = .true.
   case (full_eig_backend_syevr)
      call runFullSYEVRDP(self, amat, eval, info, t_syevr)
      if (info /= 0) then
         self%full_backend_fallbacks = self%full_backend_fallbacks + 1
         call runFullSYEVDDP(self, amat, eval, info, t_syevd)
         self%full_backend_selected = full_eig_backend_syevd
      else
         self%full_backend_selected = full_eig_backend_syevr
      end if
      self%full_backend_tuned = .true.
   case default
      if (.not.self%full_backend_autotune) then
         call runFullSYEVDDP(self, amat, eval, info, t_syevd)
         self%full_backend_selected = full_eig_backend_syevd
         self%full_backend_tuned = .true.
      else if (.not.self%full_backend_tuned) then
         call ensureFullWorkspaceDP(self)
         self%dfull_alt(:, :) = amat(:, :)
         call runFullSYEVDDP(self, amat, eval, info, t_syevd)
         if (info /= 0) return
         call runFullSYEVRRawDP(self, self%dfull_alt, self%dfull_eval, info_alt, mcalc, t_syevr)
         self%full_backend_trials = self%full_backend_trials + 1
         scale = max(1.0_dp, maxval(abs(eval)))
         tol = 8192.0_dp*epsilon(1.0_dp)*scale
         eig_ok = .false.
         if (info_alt == 0 .and. mcalc == self%n) then
            eig_ok = maxval(abs(self%dfull_eval-eval)) <= tol
         end if
         if (eig_ok .and. t_syevr <= self%full_backend_max_ratio*t_syevd) then
            eval(:) = self%dfull_eval(:)
            amat(:, :) = self%dz(:, :)
            self%full_backend_selected = full_eig_backend_syevr
         else
            if (.not.eig_ok) self%full_backend_fallbacks = self%full_backend_fallbacks + 1
            self%full_backend_selected = full_eig_backend_syevd
         end if
         self%full_backend_tuned = .true.
      else if (self%full_backend_selected == full_eig_backend_syevr) then
         call runFullSYEVRDP(self, amat, eval, info, t_syevr)
         if (info /= 0) then
            self%full_backend_fallbacks = self%full_backend_fallbacks + 1
            call runFullSYEVDDP(self, amat, eval, info, t_syevd)
            self%full_backend_selected = full_eig_backend_syevd
         end if
      else
         call runFullSYEVDDP(self, amat, eval, info, t_syevd)
      end if
   end select
end subroutine diagonalizeFullDP


subroutine runFullSYEVDSP(self, amat, eval, info, elapsed)
   class(TEigenSolver), intent(inout) :: self
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   integer, intent(out) :: info
   real(dp), intent(out) :: elapsed
   integer :: lswork, liwork
   real(dp) :: t0
   lswork = size(self%swork); liwork = size(self%iwork)
   t0 = eigensolverWallTime()
   call lapack_syevd('v','u',self%n,amat,self%n,eval,self%swork,lswork,self%iwork,liwork,info)
   elapsed = eigensolverWallTime()-t0
   self%full_syevd_time = self%full_syevd_time + elapsed
   self%full_syevd_calls = self%full_syevd_calls + 1
end subroutine runFullSYEVDSP

subroutine runFullSYEVDDP(self, amat, eval, info, elapsed)
   class(TEigenSolver), intent(inout) :: self
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   integer, intent(out) :: info
   real(dp), intent(out) :: elapsed
   integer :: ldwork, liwork
   real(dp) :: t0
   ldwork = size(self%dwork); liwork = size(self%iwork)
   t0 = eigensolverWallTime()
   call lapack_syevd('v','u',self%n,amat,self%n,eval,self%dwork,ldwork,self%iwork,liwork,info)
   elapsed = eigensolverWallTime()-t0
   self%full_syevd_time = self%full_syevd_time + elapsed
   self%full_syevd_calls = self%full_syevd_calls + 1
end subroutine runFullSYEVDDP

subroutine ensureFullWorkspaceSP(self)
   class(TEigenSolver), intent(inout) :: self
   if (.not.allocated(self%sfull_alt)) allocate(self%sfull_alt(self%n,self%n))
   if (.not.allocated(self%sfull_eval)) allocate(self%sfull_eval(self%n))
   if (.not.allocated(self%sz) .or. size(self%sz,1) /= self%n .or. size(self%sz,2) < self%n) then
      if (allocated(self%sz)) deallocate(self%sz)
      allocate(self%sz(self%n,self%n))
   end if
end subroutine ensureFullWorkspaceSP

subroutine ensureFullWorkspaceDP(self)
   class(TEigenSolver), intent(inout) :: self
   if (.not.allocated(self%dfull_alt)) allocate(self%dfull_alt(self%n,self%n))
   if (.not.allocated(self%dfull_eval)) allocate(self%dfull_eval(self%n))
   if (.not.allocated(self%dz) .or. size(self%dz,1) /= self%n .or. size(self%dz,2) < self%n) then
      if (allocated(self%dz)) deallocate(self%dz)
      allocate(self%dz(self%n,self%n))
   end if
end subroutine ensureFullWorkspaceDP

subroutine runFullSYEVRRawSP(self, awork, eval, info, mcalc, elapsed)
   class(TEigenSolver), intent(inout) :: self
   real(sp), intent(inout) :: awork(:, :)
   real(sp), intent(out) :: eval(:)
   integer, intent(out) :: info, mcalc
   real(dp), intent(out) :: elapsed
   integer :: lwork, liwork
   real(dp) :: t0
   call ensureFullWorkspaceSP(self)
   lwork=size(self%spartial_work); liwork=size(self%partial_iwork)
   t0=eigensolverWallTime()
   call lapack_syevr('v','a','u',self%n,awork,self%n,0.0_sp,0.0_sp,0,0,0.0_sp, &
      & mcalc,eval,self%sz,self%n,self%partial_isuppz,self%spartial_work,lwork, &
      & self%partial_iwork,liwork,info)
   elapsed=eigensolverWallTime()-t0
   self%full_syevr_time=self%full_syevr_time+elapsed
   self%full_syevr_calls=self%full_syevr_calls+1
end subroutine runFullSYEVRRawSP

subroutine runFullSYEVRRawDP(self, awork, eval, info, mcalc, elapsed)
   class(TEigenSolver), intent(inout) :: self
   real(dp), intent(inout) :: awork(:, :)
   real(dp), intent(out) :: eval(:)
   integer, intent(out) :: info, mcalc
   real(dp), intent(out) :: elapsed
   integer :: lwork, liwork
   real(dp) :: t0
   call ensureFullWorkspaceDP(self)
   lwork=size(self%dpartial_work); liwork=size(self%partial_iwork)
   t0=eigensolverWallTime()
   call lapack_syevr('v','a','u',self%n,awork,self%n,0.0_dp,0.0_dp,0,0,0.0_dp, &
      & mcalc,eval,self%dz,self%n,self%partial_isuppz,self%dpartial_work,lwork, &
      & self%partial_iwork,liwork,info)
   elapsed=eigensolverWallTime()-t0
   self%full_syevr_time=self%full_syevr_time+elapsed
   self%full_syevr_calls=self%full_syevr_calls+1
end subroutine runFullSYEVRRawDP

subroutine runFullSYEVRSP(self, amat, eval, info, elapsed)
   class(TEigenSolver), intent(inout) :: self
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   integer, intent(out) :: info
   real(dp), intent(out) :: elapsed
   integer :: mcalc
   call ensureFullWorkspaceSP(self)
   self%sfull_alt(:,:)=amat(:,:)
   call runFullSYEVRRawSP(self,self%sfull_alt,eval,info,mcalc,elapsed)
   if (info == 0 .and. mcalc == self%n) then
      amat(:,:)=self%sz(:,:)
   else if (info == 0) then
      info = -1001
   end if
end subroutine runFullSYEVRSP

subroutine runFullSYEVRDP(self, amat, eval, info, elapsed)
   class(TEigenSolver), intent(inout) :: self
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   integer, intent(out) :: info
   real(dp), intent(out) :: elapsed
   integer :: mcalc
   call ensureFullWorkspaceDP(self)
   self%dfull_alt(:,:)=amat(:,:)
   call runFullSYEVRRawDP(self,self%dfull_alt,eval,info,mcalc,elapsed)
   if (info == 0 .and. mcalc == self%n) then
      amat(:,:)=self%dz(:,:)
   else if (info == 0) then
      info = -1001
   end if
end subroutine runFullSYEVRDP


subroutine mctc_ssyevr_factorized(self, env, amat, nroots, eval, nfound, success, probeEval)
   character(len=*), parameter :: source = 'mctc_lapack_ssyevr_factorized'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(in) :: nroots
   real(sp), intent(out) :: eval(:)
   integer, intent(out) :: nfound
   logical, intent(out) :: success
   real(sp), intent(out), optional :: probeEval
   integer :: info, ncalc, mcalc, lwork, liwork

   success = .false.
   nfound = 0
   if (present(probeEval)) probeEval = huge(1.0_sp)
   if (nroots < 1 .or. nroots > self%n) return
   ncalc = nroots
   if (present(probeEval)) ncalc = min(self%n, nroots + 1)

   call lapack_sygst(1, 'u', self%n, amat, self%n, self%sbmat, self%n, info)
   if (info /= 0) return

   if (.not.allocated(self%sz) .or. size(self%sz,1) /= self%n .or. size(self%sz,2) < ncalc) then
      if (allocated(self%sz)) deallocate(self%sz)
      allocate(self%sz(self%n,ncalc))
   end if
   lwork = size(self%spartial_work)
   liwork = size(self%partial_iwork)
   call lapack_syevr('v', 'i', 'u', self%n, amat, self%n, 0.0_sp, 0.0_sp, &
      & 1, ncalc, 0.0_sp, mcalc, eval, self%sz, self%n, self%partial_isuppz, &
      & self%spartial_work, lwork, self%partial_iwork, liwork, info)
   if (info == 0 .and. mcalc == ncalc) then
      if (present(probeEval) .and. ncalc > nroots) probeEval = eval(nroots+1)
      amat(:, :) = 0.0_sp
      amat(:, 1:nroots) = self%sz(:, 1:nroots)
      call blas_trsm('l', 'u', 'n', 'n', self%n, nroots, 1.0_sp, self%sbmat, &
         & self%n, amat, self%n)
      nfound = nroots
      success = .true.
   end if
end subroutine mctc_ssyevr_factorized


subroutine mctc_dsyevr_factorized(self, env, amat, nroots, eval, nfound, success, probeEval)
   character(len=*), parameter :: source = 'mctc_lapack_dsyevr_factorized'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(in) :: nroots
   real(dp), intent(out) :: eval(:)
   integer, intent(out) :: nfound
   logical, intent(out) :: success
   real(dp), intent(out), optional :: probeEval
   integer :: info, ncalc, mcalc, lwork, liwork

   success = .false.
   nfound = 0
   if (present(probeEval)) probeEval = huge(1.0_dp)
   if (nroots < 1 .or. nroots > self%n) return
   ncalc = nroots
   if (present(probeEval)) ncalc = min(self%n, nroots + 1)

   call lapack_sygst(1, 'u', self%n, amat, self%n, self%dbmat, self%n, info)
   if (info /= 0) return

   if (.not.allocated(self%dz) .or. size(self%dz,1) /= self%n .or. size(self%dz,2) < ncalc) then
      if (allocated(self%dz)) deallocate(self%dz)
      allocate(self%dz(self%n,ncalc))
   end if
   lwork = size(self%dpartial_work)
   liwork = size(self%partial_iwork)
   call lapack_syevr('v', 'i', 'u', self%n, amat, self%n, 0.0_dp, 0.0_dp, &
      & 1, ncalc, 0.0_dp, mcalc, eval, self%dz, self%n, self%partial_isuppz, &
      & self%dpartial_work, lwork, self%partial_iwork, liwork, info)
   if (info == 0 .and. mcalc == ncalc) then
      if (present(probeEval) .and. ncalc > nroots) probeEval = eval(nroots+1)
      amat(:, :) = 0.0_dp
      amat(:, 1:nroots) = self%dz(:, 1:nroots)
      call blas_trsm('l', 'u', 'n', 'n', self%n, nroots, 1.0_dp, self%dbmat, &
         & self%n, amat, self%n)
      nfound = nroots
      success = .true.
   end if
end subroutine mctc_dsyevr_factorized


subroutine configureFullEigenBackend(self, requested, autotune, maxRatio, minN)
   class(TEigenSolver), intent(inout) :: self
   integer, intent(in) :: requested
   logical, intent(in) :: autotune
   real(dp), intent(in) :: maxRatio
   integer, intent(in) :: minN
   integer :: req

   req = requested
   if (req < full_eig_backend_auto .or. req > full_eig_backend_syevr) req = full_eig_backend_auto
   if (self%full_backend_requested /= req .or. &
      & self%full_backend_autotune .neqv. autotune .or. &
      & abs(self%full_backend_max_ratio-maxRatio) > 8.0_dp*epsilon(1.0_dp)) then
      self%full_backend_tuned = .false.
   end if
   self%full_backend_requested = req
   self%full_backend_autotune_min_n = max(1,minN)
   self%full_backend_autotune = autotune .and. self%n >= self%full_backend_autotune_min_n
   self%full_backend_max_ratio = min(2.0_dp,max(0.10_dp,maxRatio))
   if (req == full_eig_backend_syevd) then
      self%full_backend_selected = full_eig_backend_syevd
      self%full_backend_tuned = .true.
   else if (req == full_eig_backend_syevr) then
      self%full_backend_selected = full_eig_backend_syevr
      self%full_backend_tuned = .true.
   else if (.not.self%full_backend_autotune) then
      self%full_backend_selected = full_eig_backend_syevd
      self%full_backend_tuned = .true.
   end if
end subroutine configureFullEigenBackend

subroutine getFullEigenBackendProfile(self, requested, selected, autotune, maxRatio, &
      & syevdTime, syevrTime, syevdCalls, syevrCalls, trials, fallbacks)
   class(TEigenSolver), intent(in) :: self
   integer, intent(out) :: requested, selected, syevdCalls, syevrCalls, trials, fallbacks
   logical, intent(out) :: autotune
   real(dp), intent(out) :: maxRatio, syevdTime, syevrTime
   requested = self%full_backend_requested
   selected = self%full_backend_selected
   autotune = self%full_backend_autotune
   maxRatio = self%full_backend_max_ratio
   syevdTime = self%full_syevd_time
   syevrTime = self%full_syevr_time
   syevdCalls = self%full_syevd_calls
   syevrCalls = self%full_syevr_calls
   trials = self%full_backend_trials
   fallbacks = self%full_backend_fallbacks
end subroutine getFullEigenBackendProfile

subroutine resetFullEigenProfile(self)
   class(TEigenSolver), intent(inout) :: self
   self%full_reduce_time = 0.0_dp
   self%full_diag_time = 0.0_dp
   self%full_back_time = 0.0_dp
   self%full_profile_calls = 0
   self%full_syevd_time = 0.0_dp
   self%full_syevr_time = 0.0_dp
   self%full_syevd_calls = 0
   self%full_syevr_calls = 0
   self%full_backend_trials = 0
   self%full_backend_fallbacks = 0
end subroutine resetFullEigenProfile

subroutine getFullEigenProfile(self, reduceTime, diagTime, backTime, calls)
   class(TEigenSolver), intent(in) :: self
   real(dp), intent(out) :: reduceTime, diagTime, backTime
   integer, intent(out) :: calls
   reduceTime = self%full_reduce_time
   diagTime = self%full_diag_time
   backTime = self%full_back_time
   calls = self%full_profile_calls
end subroutine getFullEigenProfile

real(dp) function eigensolverWallTime() result(t)
   integer(int64) :: count, rate
   call system_clock(count=count, count_rate=rate)
   if (rate > 0_int64) then
      t = real(count,dp)/real(rate,dp)
   else
      call cpu_time(t)
   end if
end function eigensolverWallTime

end module xtb_mctc_lapack_eigensolve
