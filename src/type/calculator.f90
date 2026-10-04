! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

!> abstract calculator that hides implementation details from calling codes
module xtb_type_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_model, only : TSolvModel
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy, getGFN1FastPolicy, &
      & useGFN1StaticHessianSchedule, useGFN1ParallelHessian
   implicit none

   public :: TCalculator
   private


   !> Base calculator
   type, abstract :: TCalculator

      real(wp) :: accuracy
      logical :: lSolv = .false.
      type(TSolvModel), allocatable :: solvation
      logical :: threadsafe = .true.

   contains

      !> Perform single point calculation
      procedure(singlepoint), deferred :: singlepoint

      !> Perform hessian calculation
      procedure :: hessian

      !> Write informative printout
      procedure(writeInfo), deferred :: writeInfo

   end type TCalculator


   abstract interface
      subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
            & energy, gradient, sigma, hlgap, results)
         import :: TCalculator, TEnvironment, TMolecule, TRestart, wp
         import :: scc_results

         !> Calculator instance
         class(TCalculator), intent(inout) :: self

         !> Computational environment
         type(TEnvironment), intent(inout) :: env

         !> Molecular structure data
         type(TMolecule), intent(inout) :: mol

         !> Wavefunction data
         type(TRestart), intent(inout) :: chk

         !> Print level for IO
         integer, intent(in) :: printlevel

         !> Restart from previous results
         logical, intent(in) :: restart

         !> Total energy
         real(wp), intent(out) :: energy

         !> Molecular gradient
         real(wp), intent(out) :: gradient(:, :)

         !> Strain derivatives
         real(wp), intent(out) :: sigma(:, :)

         !> HOMO-LUMO gap
         real(wp), intent(out) :: hlgap

         !> Detailed results
         type(scc_results), intent(out) :: results

      end subroutine singlepoint


      subroutine writeInfo(self, unit, mol)
         import :: TCalculator, TMolecule

         !> Calculator instance
         class(TCalculator), intent(in) :: self

         !> Unit for I/O
         integer, intent(in) :: unit

         !> Molecular structure data
         type(TMolecule), intent(in) :: mol

      end subroutine writeInfo
   end interface


contains


!> Evaluate hessian by finite difference for all atoms
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad)
   !$ use omp_lib, only : omp_sched_kind, omp_sched_static, omp_sched_dynamic, &
   !$& omp_get_schedule, omp_set_schedule, omp_get_max_threads
   character(len=*), parameter :: source = "hessian_numdiff_numdiff2"
   !> Single point calculator
   class(TCalculator), intent(inout) :: self
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol0
   !> Restart data
   type(TRestart), intent(in) :: chk0
   !> List of atoms to displace
   integer, intent(in) :: list(:)
   !> Step size for numerical differentiation
   real(wp), intent(in) :: step
   !> Array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)
   !> Array to add dipole gradient to
   real(wp), intent(inout) :: dipgrad(:, :)
   !> Array to add polarizability gradient to
   real(wp), intent(inout), optional :: polgrad(:, :)

   integer :: iat, jat, kat, ic, jc, ii, jj
   type(TMolecule) :: mol
   type(TRestart) :: chk
   real(wp) :: er, el, dr(3), dl(3), sr(3, 3), sl(3, 3), egap_local, step2
   real(wp) :: alphal(3, 3), alphar(3, 3)
   real(wp) :: t0, t1, w0, w1
   real(wp), allocatable :: gr(:, :), gl(:, :)
   type(scc_results) :: rr, rl
   type(TGFN1FastPolicy) :: fastPolicy
   logical :: hessianStatic, hessianParallel, fastProfile
   integer :: hessianThreads
   !$ integer(kind=omp_sched_kind) :: oldScheduleKind
   !$ integer :: oldScheduleChunk

   call getGFN1FastPolicy(fastPolicy)
   hessianStatic = useGFN1StaticHessianSchedule(fastPolicy, chk0%wfn%nao)
   fastProfile = fastPolicy%profile
   hessianThreads = 1
   !$ hessianThreads = omp_get_max_threads()
   hessianParallel = self%threadsafe .and. useGFN1ParallelHessian(fastPolicy, &
      & chk0%wfn%nao, 3*size(list), hessianThreads)
   !$ call omp_get_schedule(oldScheduleKind, oldScheduleChunk)
   !$ if (hessianStatic) then
   !$    call omp_set_schedule(omp_sched_static, 1)
   !$ else
   !$    call omp_set_schedule(omp_sched_dynamic, 1)
   !$ end if

   call timing(t0, w0)
   step2 = 0.5_wp / step

   ! GFN1-fast 1.5.0: numerical-Hessian displacement-level parallelism.
   ! The stock implementation already exposed an OpenMP work-sharing loop,
   ! but its private allocatable scratch arrays and shared HOMO-LUMO-gap
   ! scalar made the execution model unnecessarily fragile.  Use an explicit
   ! parallel region so each worker allocates its gradient scratch only once,
   ! keeps all scalar results thread-local, and dynamically balances Cartesian
   ! displacements whose SCC iteration counts can differ substantially.
   !
   ! Each (atom,Cartesian) work item owns one Hessian column and evaluates the
   ! +step/-step pair serially.  Therefore no floating-point reduction or
   ! accumulation order is changed relative to the serial algorithm.
   !$omp parallel if(hessianParallel) default(none) &
   !$omp shared(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad, step2, t0, w0, hessianParallel) &
   !$omp private(kat, ic, iat, jat, jc, jj, ii, er, el, gr, gl, sr, sl, rr, rl, dr, dl, &
   !$omp& alphar, alphal, mol, chk, t1, w1, egap_local)
   allocate(gr(3, mol0%n), gl(3, mol0%n))

   ! GFN1-fast 1.6.0: initialize large thread-local objects once.  A Hessian
   ! displacement changes only one Cartesian coordinate; singlepoint() calls
   ! mol%update itself, so rebuilding/reallocating the complete TMolecule for
   ! every +/- point is unnecessary.  TRestart%copy now also reuses the
   ! already-allocated wavefunction buffers.
   call mol%copy(mol0)
   call chk%copy(chk0)

   !$omp do schedule(runtime) collapse(2)
   do kat = 1, size(list)
      do ic = 1, 3
         iat = list(kat)
         ii = 3*(iat - 1) + ic
         er = 0.0_wp
         el = 0.0_wp
         gr = 0.0_wp
         gl = 0.0_wp

         mol%xyz(ic, iat) = mol0%xyz(ic, iat) + step
         call chk%copy(chk0)
         call self%singlepoint(env, mol, chk, -1, .true., er, gr, sr, egap_local, rr)
         dr = rr%dipole
         if (present(polgrad)) then
            alphar(:, :) = rr%alpha
         endif

         mol%xyz(ic, iat) = mol0%xyz(ic, iat) - step
         call chk%copy(chk0)
         call self%singlepoint(env, mol, chk, -1, .true., el, gl, sl, egap_local, rl)
         dl = rl%dipole
         ! Restore the worker molecule to the reference geometry for its next
         ! dynamically scheduled displacement. Derived geometry data are
         ! refreshed by singlepoint() before use.
         mol%xyz(ic, iat) = mol0%xyz(ic, iat)

         if (present(polgrad)) then
            alphal(:, :) = rl%alpha
            polgrad(1, ii) = (alphar(1, 1) - alphal(1, 1)) * step2
            polgrad(2, ii) = (alphar(1, 2) - alphal(1, 2)) * step2
            polgrad(3, ii) = (alphar(2, 2) - alphal(2, 2)) * step2
            polgrad(4, ii) = (alphar(1, 3) - alphal(1, 3)) * step2
            polgrad(5, ii) = (alphar(2, 3) - alphal(2, 3)) * step2
            polgrad(6, ii) = (alphar(3, 3) - alphal(3, 3)) * step2
         endif
         dipgrad(:, ii) = (dr - dl) * step2
         do jat = 1, mol0%n
            do jc = 1, 3
               jj = 3*(jat - 1) + jc
               hess(jj, ii) = hess(jj, ii) &
                  & + (gr(jc, jat) - gl(jc, jat)) * step2
            end do
         end do

         if (kat == 3 .and. ic == 3) then
            !$omp critical(xtb_numdiff2)
            call timing(t1, w1)
            write(*,'("estimated CPU  time",F10.2," min")') &
               & 0.3333333_wp*size(list)*(t1-t0)/60.0_wp
            write(*,'("estimated wall time",F10.2," min")') &
               & 0.3333333_wp*size(list)*(w1-w0)/60.0_wp
            !$omp end critical(xtb_numdiff2)
         endif

      end do
   end do
   !$omp end do

   deallocate(gr, gl)
   !$omp end parallel
   !$ call omp_set_schedule(oldScheduleKind, oldScheduleChunk)

   if (fastProfile) then
      call timing(t1, w1)
      write(env%unit,'(/,1x,a)') 'GFN1-fast 2.0.3 Hessian profile'
      write(env%unit,'(3x,a,i0)') 'NAO                 = ',chk0%wfn%nao
      write(env%unit,'(3x,a,i0)') 'Cartesian jobs       = ',3*size(list)
      write(env%unit,'(3x,a,l1,a,i0)') 'outer OpenMP         = ',hessianParallel, &
         & ', max threads=',hessianThreads
      if (hessianStatic) then
         write(env%unit,'(3x,a)') 'OpenMP schedule       = static,1'
      else
         write(env%unit,'(3x,a)') 'OpenMP schedule       = dynamic,1'
      end if
      write(env%unit,'(3x,a,i0)') 'static NAO threshold = ',fastPolicy%hessianStaticMaxNao
      write(env%unit,'(3x,a,f10.6,a)') 'finite-diff wall     = ',w1-w0,' s'
   end if
end subroutine hessian


end module xtb_type_calculator
