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

!> general functions for core functionalities of the SCC
module xtb_scc_core
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_la, only : contract
   use xtb_mctc_lapack, only : lapack_sygvd
   use xtb_mctc_blas, only : blas_gemm, blas_syrk, mctc_symv, mctc_gemm
   use xtb_mctc_lapack_eigensolve, only : TEigenSolver
   use xtb_type_environment, only : TEnvironment
   use xtb_type_solvation, only : TSolvation
   use xtb_xtb_data
   use xtb_xtb_coulomb
   use xtb_xtb_dispersion
   use xtb_xtb_multipole
   use xtb_broyden
   use xtb_gfn1_davidson, only : gfn1BlockDavidsonMF, buildGFN1SparsePairs, &
      & updateGFN1PairHamiltonian, applyGFN1PairHamiltonianBlock, &
      & applyGFN1PairMetricBlock, gfn1DavidsonGuardCount, &
      & validateGFN1DavidsonEigensystem
   implicit none
   private

   public :: build_h0, scc, electro, solve, solve4
   public :: fermismear, occ, occu, dmat, get_unrestricted_wiberg
   public :: get_wiberg, mpopall, mpop0, mpopao, mpop, mpopsh, mpopsh_h0_upper, qsh2qat, lpop
   public :: iniqshell, setzshell
   public :: shellPoly, h0scal


   integer, private, parameter :: mmm(20)=(/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/)


contains

!! ========================================================================
!  build GFN2 core Hamiltonian
!! ========================================================================
subroutine build_h0(hData,H0,n,at,ndim,nmat,matlist, &
   &                xyz,selfEnergy,S,aoat2,lao2,valao2,aoexp,ao2sh)
   type(THamiltonianData), intent(in) :: hData
   real(wp),intent(out) :: H0(ndim*(ndim+1)/2)
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nmat
   integer, intent(in)  :: matlist(2,nmat)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: selfEnergy(:)
   real(wp),intent(in)  :: S(ndim,ndim)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: lao2(ndim)
   integer, intent(in)  :: valao2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   real(wp),intent(in)  :: aoexp(ndim)

   integer  :: i,j,k,m
   integer  :: iat,jat,ish,jsh,il,jl,iZp,jZp
   real(wp) :: hdii,hdjj,hav
   real(wp) :: km

   H0=0.0_wp

   do m = 1, nmat
      i = matlist(1,m)
      j = matlist(2,m)
      k = j+i*(i-1)/2
      iat = aoat2(i)
      jat = aoat2(j)
      ish = ao2sh(i)
      jsh = ao2sh(j)
      iZp = at(iat)
      jZp = at(jat)
      il = mmm(lao2(i))
      jl = mmm(lao2(j))
      hdii = selfEnergy(ish)
      hdjj = selfEnergy(jsh)
      call h0scal(hData,il,jl,izp,jzp,valao2(i).ne.0,valao2(j).ne.0, &
      &           km)
      km = km*(2*sqrt(aoexp(i)*aoexp(j))/(aoexp(i)+aoexp(j)))**hData%wExp
      hav = 0.5d0*(hdii+hdjj)* &
      &      shellPoly(hData%shellPoly(il, iZp), hData%shellPoly(jl, jZp), &
      &                hData%atomicRad(iZp), hData%atomicRad(jZp),xyz(:,iat),xyz(:,jat))
      H0(k) = S(j,i)*km*hav
   enddo
!  diagonal
   k=0
   do i=1,ndim
      k=k+i
      iat = aoat2(i)
      ish = ao2sh(i)
      il = mmm(lao2(i))
      H0(k) = selfEnergy(ish)
   enddo

end subroutine build_h0

!> build isotropic H1/Fockian
subroutine buildIsotropicH1(n, at, ndim, nshell, nmat, matlist, H, &
      & H0, S, shellShift, aoat2, ao2sh)
   use xtb_mctc_convert, only : autoev,evtoau
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nshell
   integer, intent(in)  :: nmat
   integer, intent(in)  :: matlist(2,nmat)
   real(wp),intent(in)  :: H0(ndim*(1+ndim)/2)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(in)  :: shellShift(nshell)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   real(wp),intent(out) :: H(ndim,ndim)

   integer  :: m,i,j,k
   integer  :: ishell,jshell
   integer  :: ii,jj,kk
   real(wp) :: dum
   real(wp) :: eh1,t8,t9,tgb,h1

   H = 0.0_wp

   do m = 1, nmat
      i = matlist(1,m)
      j = matlist(2,m)
      k = j+i*(i-1)/2
      ishell = ao2sh(i)
      jshell = ao2sh(j)
      ! SCC terms
      eh1 = autoev*(shellShift(ishell) + shellShift(jshell))
      H1 = -S(j,i)*eh1*0.5_wp
      H(j,i) = H0(k) + H1
      H(i,j) = H(j,i)
   enddo

end subroutine buildIsotropicH1

!> build isotropic H1/Fockian using precomputed AO/shell metadata.
!> GFN1-fast 0.5.0: removes invariant index arithmetic and AO->shell
!> lookups from every SCC iteration without changing the Hamiltonian formula.
subroutine buildIsotropicH1Cached(ndim, nmat, matlist, h0idx, matISh, matJSh, H, &
      & H0, S, shellShift, useParallel)
   use xtb_mctc_convert, only : autoev
   integer, intent(in) :: ndim, nmat
   integer, intent(in) :: matlist(2,nmat), h0idx(nmat), matISh(nmat), matJSh(nmat)
   real(wp), intent(in) :: H0(ndim*(1+ndim)/2), S(ndim,ndim), shellShift(:)
   logical, intent(in) :: useParallel
   real(wp), intent(out) :: H(ndim,ndim)
   integer :: m, i, j
   real(wp) :: eh1, h1

   H = 0.0_wp
   ! GFN1-fast 1.0.0: each packed AO pair is independent. Parallelize only
   ! this write-once loop, avoiding any floating-point reductions so the
   ! numerical result is unchanged by the OpenMP thread count. Keep a true
   ! serial branch for small matrices / one-thread runs to avoid OpenMP
   ! runtime overhead.
   if (useParallel) then
      !$omp parallel do default(none) schedule(static) &
      !$omp shared(nmat,matlist,h0idx,matISh,matJSh,H,H0,S,shellShift) &
      !$omp private(m,i,j,eh1,h1)
      do m = 1, nmat
         i = matlist(1,m)
         j = matlist(2,m)
         eh1 = autoev*(shellShift(matISh(m)) + shellShift(matJSh(m)))
         h1 = -S(j,i)*eh1*0.5_wp
         H(j,i) = H0(h0idx(m)) + h1
         H(i,j) = H(j,i)
      enddo
      !$omp end parallel do
   else
      do m = 1, nmat
         i = matlist(1,m)
         j = matlist(2,m)
         eh1 = autoev*(shellShift(matISh(m)) + shellShift(matJSh(m)))
         h1 = -S(j,i)*eh1*0.5_wp
         H(j,i) = H0(h0idx(m)) + h1
         H(i,j) = H(j,i)
      enddo
   end if
end subroutine buildIsotropicH1Cached

!> build anisotropic H1/Fockian
subroutine addAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp,matlist,mdlst,mqlst,&
                         H,S,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
   use xtb_mctc_convert, only : autoev,evtoau
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nshell
   integer, intent(in)  :: nmat
   integer, intent(in)  :: ndp
   integer, intent(in)  :: nqp
   integer, intent(in)  :: matlist(2,nmat)
   integer, intent(in)  :: mdlst(2,ndp)
   integer, intent(in)  :: mqlst(2,nqp)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(in)  :: dpint(3,ndim,ndim)
   real(wp),intent(in)  :: qpint(6,ndim,ndim)
   real(wp),intent(in)  :: vs(n)
   real(wp),intent(in)  :: vd(3,n)
   real(wp),intent(in)  :: vq(6,n)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   real(wp),intent(inout) :: H(ndim,ndim)

   integer, external :: lin
   integer  :: m,i,j,k,l
   integer  :: ii,jj,kk
   integer  :: ishell,jshell
   real(wp) :: dum,eh1,t8,t9,tgb

   !> overlap dependent terms
   do m=1,nmat
      i=matlist(1,m)
      j=matlist(2,m)
      k=j+i*(i-1)/2
      ii=aoat2(i)
      jj=aoat2(j)
      dum=S(j,i)
      ! CAMM potential
      eh1=0.50d0*dum*(vs(ii)+vs(jj))*autoev
      H(j,i)=H(j,i)+eh1
      H(i,j)=H(j,i)
   enddo
   !> dipolar terms
   do m=1,ndp
      i=mdlst(1,m)
      j=mdlst(2,m)
      k=lin(j,i)
      ii=aoat2(i)
      jj=aoat2(j)
      eh1=0.0d0
      do l=1,3
         eh1=eh1+dpint(l,i,j)*(vd(l,ii)+vd(l,jj))
      enddo
      eh1=0.50d0*eh1*autoev
      H(i,j)=H(i,j)+eh1
      H(j,i)=H(i,j)
   enddo
   !> quadrupole-dependent terms
   do m=1,nqp
      i=mqlst(1,m)
      j=mqlst(2,m)
      ii=aoat2(i)
      jj=aoat2(j)
      k=lin(j,i)
      eh1=0.0d0
      ! note: these come in the following order
      ! xx, yy, zz, xy, xz, yz
      do l=1,6
         eh1=eh1+qpint(l,i,j)*(vq(l,ii)+vq(l,jj))
      enddo
      eh1=0.50d0*eh1*autoev
      H(i,j)=H(i,j)+eh1
      H(j,i)=H(i,j)
   enddo

end subroutine addAnisotropicH1


!> self consistent charge iterator
subroutine scc(env,xtbData,solver,n,nel,nopen,ndim,ndp,nqp,nmat,nshell, &
      &        at,matlist,mdlst,mqlst,aoat2,ao2sh,ash, &
      &        q,dipm,qp,qq,qlmom,qsh,zsh, &
      &        xyz,aes, &
      &        cm5,cm5a,gborn,solvation, &
      &        scD4, &
      &        broy,broydamp,damp0, &
      &        pcem,shellShift,externShift, &
      &        et,focc,focca,foccb,efa,efb, &
      &        eel,ees,eaes,epol,ed,epcem,egap,emo,ihomo,ihomoa,ihomob, &
      &        H0,H,S,dpint,qpint,P,ies, &
      &        maxiter,startpdiag,scfconv,qconv,restart, &
      &        minpr,pr, &
      &        fail,jter)
   use xtb_mctc_convert, only : autoev,evtoau
   use xtb_mctc_constants, only : kB
   !$ use omp_lib, only : omp_get_max_threads, omp_in_parallel

   use xtb_disp_dftd4,  only: disppot,edisp_scc
   use xtb_aespot, only : gfn2broyden_diff,gfn2broyden_out,gfn2broyden_save, &
   &                  mmompop,aniso_electro,setvsdq
   use xtb_embedding, only : electro_pcem
   use xtb_solv_gbsa, only : TBorn
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy, getGFN1FastPolicy, &
      & useGFN1PartialEigensolver, useGFN1Davidson, useGFN1SparseDavidson, &
      & gfn1SpectralPartialUnsafe, gfn1PartialSpectrumReady, &
      & gfn1PartialSpectrumReadyThermal, gfn1ThermalRootGuard, gfn1PartialAutotuneKeep
   character(len=*), parameter :: source = 'scc_core'

   type(TEnvironment), intent(inout) :: env

   type(TxTBData), intent(in) :: xtbData
   type(TEigenSolver), intent(inout) :: solver

   integer, intent(in)  :: n
   integer, intent(in)  :: nel
   integer, intent(in)  :: nopen
   integer, intent(in)  :: ndim
   integer, intent(in)  :: ndp
   integer, intent(in)  :: nqp
   integer, intent(in)  :: nmat
   integer, intent(in)  :: nshell
!! ------------------------------------------------------------------------
!  general options for the iterator
   integer, intent(in)  :: maxiter
   integer, intent(in)  :: startpdiag
   real(wp),intent(in)  :: scfconv
   real(wp),intent(in)  :: qconv
   ! GFN1-fast 1.4.0: restart means that the incoming wavefunction/charges
   ! originate from a converged nearby geometry.  Numerical Hessian calls
   ! already provide exactly such a reference checkpoint for every +/-
   ! displacement.  We use it only to seed the exact partial-eigensolver
   ! root count; all existing residual/occupation/degeneracy gates remain.
   logical, intent(in)  :: restart
   logical, intent(in)  :: minpr
   logical, intent(in)  :: pr
   logical, intent(out) :: fail
!! ------------------------------------------------------------------------
   integer, intent(in)  :: at(n)
!! ------------------------------------------------------------------------
   integer, intent(in)  :: matlist(2,nmat)
   integer, intent(in)  :: mdlst(2,ndp)
   integer, intent(in)  :: mqlst(2,nqp)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   integer, intent(in)  :: ash(:)
!! ------------------------------------------------------------------------
!  a bunch of charges and CAMMs
   real(wp),intent(inout) :: q(n)
   real(wp),intent(inout) :: dipm(3,n)
   real(wp),intent(inout) :: qp(6,n)
   real(wp),intent(inout) :: qq(n)
   real(wp),intent(inout) :: qlmom(3,n)
   real(wp),intent(inout) :: qsh(nshell)
   real(wp),intent(in)    :: zsh(nshell)
!! ------------------------------------------------------------------------
!  anisotropic electrostatic
   type(TxTBMultipole), intent(in), optional :: aes
   real(wp),intent(in)    :: xyz(3,n)
   real(wp), allocatable :: vs(:)
   real(wp), allocatable :: vd(:, :)
   real(wp), allocatable :: vq(:, :)
!! ------------------------------------------------------------------------
!  continuum solvation model GBSA
   class(TSolvation), allocatable, intent(inout) :: solvation
   real(wp),intent(in)    :: cm5a(n)
   real(wp),intent(inout) :: cm5(n)
   real(wp),intent(inout) :: gborn
!! ------------------------------------------------------------------------
!  selfconsistent DFT-D4 dispersion correction
   type(TxTBDispersion), intent(inout), optional :: scD4
!! ------------------------------------------------------------------------
!  point charge embedding potentials
   logical, intent(in)    :: pcem
   real(wp),intent(inout) :: shellShift(nshell)
   real(wp),intent(inout) :: externShift(nshell)
   real(wp), allocatable :: atomicShift(:)
   integer, allocatable :: matH0Idx(:), matISh(:), matJSh(:)
!! ------------------------------------------------------------------------
!  Fermi-smearing
   real(wp),intent(in)    :: et
   real(wp),intent(inout) :: focc(ndim)
   real(wp),intent(inout) :: foccb(ndim),focca(ndim)
   real(wp),intent(inout) :: efa,efb
!! ------------------------------------------------------------------------
!  Convergence accelerators, a simple damping as well as a Broyden mixing
!  are available. The Broyden mixing is used by default seems reliable.
   real(wp),intent(in)    :: damp0
   real(wp)               :: damp
!  Broyden
   integer                :: nbr
   logical, intent(in)    :: broy
   real(wp),intent(inout) :: broydamp
   real(wp)               :: omegap
   real(wp),allocatable   :: df(:,:)
   real(wp),allocatable   :: u(:,:)
   real(wp),allocatable   :: a(:,:)
   real(wp),allocatable   :: q_in(:)
   real(wp),allocatable   :: dq(:)
   real(wp),allocatable   :: qlast_in(:)
   real(wp),allocatable   :: dqlast(:)
   real(wp),allocatable   :: omega(:)
   ! GFN1-fast 0.3.0: reusable density-matrix scratch space
   real(wp),allocatable   :: dmatWork(:,:), eigWork(:,:)
   ! GFN1-fast 2.1.1: compact previous-SCC eigenspace and Davidson output.
   ! These are NAO x Nactive rather than full NAO x NAO work matrices.
   real(wp),allocatable   :: davidsonSeed(:,:), davidsonC(:,:), davidsonCarry(:,:), &
      & pairH(:), pairS(:), pairH0(:), auditEval(:)
   ! GFN1-fast 1.2.0: exact linear propagation of the ALPB/GBSA atomic
   ! potential through the SCC mixer.  This removes the Born-matrix SYMV
   ! from addShift after the first SCC step; the final verification step is
   ! deliberately refreshed from the matrix.
   real(wp),allocatable   :: solvPotIn(:), solvPotOut(:), solvPotAux(:)
   real(wp),allocatable   :: solvPotAuxLast(:), solvDPot(:), solvDPotLast(:)
   real(wp),allocatable   :: solvUAux(:,:)
!! ------------------------------------------------------------------------
!  results of the SCC iterator
   real(wp),intent(out)   :: eel
   real(wp),intent(out)   :: epcem
   real(wp),intent(out)   :: ees
   real(wp),intent(out)   :: eaes
   real(wp),intent(out)   :: epol
   real(wp),intent(out)   :: ed
   real(wp),intent(out)   :: egap
   real(wp),intent(out)   :: emo(ndim)
   integer, intent(inout) :: ihomo
   integer, intent(inout) :: ihomoa
   integer, intent(inout) :: ihomob
!! ------------------------------------------------------------------------
   real(wp),intent(in)    :: H0(ndim*(ndim+1)/2)
   real(wp),intent(out)   :: H(ndim,ndim)
   real(wp),intent(inout) :: P(ndim,ndim)
   real(wp),intent(in)    :: S(ndim,ndim)
   real(wp),intent(in)    :: dpint(:,:,:)
   real(wp),intent(in)    :: qpint(:,:,:)
   type(TxTBCoulomb), intent(inout) :: ies

   integer, intent(inout) :: jter
!! ------------------------------------------------------------------------
!  local variables
   integer,external :: lin
   integer  :: i,ii,j,jj,k,kk,l,m
   integer  :: ishell,jshell
   real(wp) :: t8,t9
   real(wp) :: eh1,dum,tgb
   real(wp) :: eold
   real(wp) :: h0EnergyCompact
   real(wp) :: ga,gb
   real(wp) :: rmsq
   real(wp) :: nfoda,nfodb
   logical  :: fulldiag
   logical  :: lastdiag
   integer  :: iter
   integer  :: thisiter
   integer  :: nsub, nroots, nfound, last_occ, required_roots, thermalGuardRoots
   integer  :: davidsonIter, davidsonRestart, davidsonSeedCols, davidsonCarryCols
   integer  :: davidsonGuard, seedCols
   real(wp) :: maxres, maxnormerr, maxadjovlp, maxprojerr, bkt, efmax
   real(wp) :: davidsonMaxRes, davidsonStart, davidsonElapsed, lastFullSolveTime
   real(wp) :: partialPathStart, partialPathElapsed
   real(wp) :: boundary_scale, boundary_tol, occdiff
   real(wp) :: maxFractionality
   real(wp) :: auditEigMax, auditSubspaceError, auditOverlap2
   real(wp) :: dampUsed
   logical  :: partial_used, partial_ok, fermi_covered, occupation_covered
   logical  :: matrixFreePartial, hamiltonianBuilt, fullSolvedEarly
   logical  :: davidsonDisabled
   logical  :: spectralPartialDisabled, spectralUnsafe
   logical  :: partialPathDisabled, partialSpectrumReady, fullSpectrumThisIter
   logical  :: partialAutotuneDisabled
   logical  :: occupation_degeneracy_safe
   logical  :: parallelH1
   logical  :: bornResponseCache, solvPotValid
   logical  :: converged
   logical  :: compactDensityActive, compactThisIter
   logical  :: econverged
   logical  :: qconverged
   ! GFN1-fast 2.0.2 adaptive execution/profiling.
   type(TGFN1FastPolicy) :: fastPolicy
   logical :: fastProfile, partial_attempted
   integer :: profPartialAttempts, profPartialAccepted, profPartialFallbacks
   integer :: profPartialDeferred, certifiedFullSolves, partialFallbackReason
   integer :: profThermalRetries, profThermalDisables, thermalFermiFallbacks
   integer :: profFallbackSolver, profFallbackSpectral, profFallbackOccupation
   integer :: profFallbackDegeneracy, profFallbackFermi, profFallbackAudit
   integer :: profDavidsonAttempts, profDavidsonAccepted, profDavidsonIters, profDavidsonRestarts
   integer :: profFullSolves, profDavidsonAudits, profDavidsonAuditFailures
   integer :: profSpectralDisables, profAutotuneDisables, profPartialTimed
   integer :: profCompactDensityIterations
   real(wp) :: profStart
   real(wp) :: profPotential, profHamiltonian, profEigensolver
   real(wp) :: profSubspace, profDensity, profPopulation, profEnergy
   real(wp) :: profConvergence, profMixing
   real(wp) :: profFullCoreTime, profPartialPathTime
   real(wp) :: profFullReduceTime, profFullDiagTime, profFullBackTime
   real(wp) :: profFullSYEVDTime, profFullSYEVRTime, profFullBackendMaxRatio
   integer :: profFullPhaseCalls
   integer :: profFullBackendRequested, profFullBackendSelected
   integer :: profFullSYEVDCalls, profFullSYEVRCalls, profFullBackendTrials, profFullBackendFallbacks
   logical :: profFullBackendAutotune
   character(len=8) :: profFullBackendRequestedName, profFullBackendSelectedName

   converged = .false.
   lastdiag = .false.
   nsub = 0
   call getGFN1FastPolicy(fastPolicy)
   call solver%configure_full_backend(fastPolicy%fullEigBackend, fastPolicy%fullEigAutotune, &
      & real(fastPolicy%fullEigAutotuneMaxRatio,wp), fastPolicy%fullEigAutotuneMinNao)
   compactDensityActive = fastPolicy%compactDensitySCC .and. xtbData%level == 1 .and. &
      & .not.present(aes) .and. ndim >= fastPolicy%compactDensityMinNao
   ! Only the user-visible SCC call emits detailed timings. Numerical Hessian
   ! displacement calls use a negative/minimal print level and would otherwise
   ! flood a shared output stream from multiple OpenMP workers.
   fastProfile = fastPolicy%profile .and. minpr
   profPotential = 0.0_wp
   profHamiltonian = 0.0_wp
   profEigensolver = 0.0_wp
   profSubspace = 0.0_wp
   profDensity = 0.0_wp
   profPopulation = 0.0_wp
   profEnergy = 0.0_wp
   profConvergence = 0.0_wp
   profMixing = 0.0_wp
   profFullCoreTime = 0.0_wp
   profPartialPathTime = 0.0_wp
   profFullReduceTime = 0.0_wp
   profFullDiagTime = 0.0_wp
   profFullBackTime = 0.0_wp
   profFullSYEVDTime = 0.0_wp
   profFullSYEVRTime = 0.0_wp
   profFullBackendMaxRatio = real(fastPolicy%fullEigAutotuneMaxRatio,wp)
   profFullPhaseCalls = 0
   profFullBackendRequested = fastPolicy%fullEigBackend
   profFullBackendSelected = 1
   profFullSYEVDCalls = 0
   profFullSYEVRCalls = 0
   profFullBackendTrials = 0
   profFullBackendFallbacks = 0
   profFullBackendAutotune = fastPolicy%fullEigAutotune
   profFullBackendRequestedName = 'auto'
   profFullBackendSelectedName = 'syevd'
   profPartialAttempts = 0
   profPartialAccepted = 0
   profPartialFallbacks = 0
   profPartialDeferred = 0
   profThermalRetries = 0
   profThermalDisables = 0
   thermalFermiFallbacks = 0
   profFallbackSolver = 0
   profFallbackSpectral = 0
   profFallbackOccupation = 0
   profFallbackDegeneracy = 0
   profFallbackFermi = 0
   profFallbackAudit = 0
   profDavidsonAudits = 0
   profDavidsonAuditFailures = 0
   auditEigMax = 0.0_wp
   auditSubspaceError = 0.0_wp
   profDavidsonAttempts = 0
   profDavidsonAccepted = 0
   profDavidsonIters = 0
   profDavidsonRestarts = 0
   profFullSolves = 0
   profSpectralDisables = 0
   profAutotuneDisables = 0
   profPartialTimed = 0
   profCompactDensityIterations = 0
   davidsonDisabled = .false.
   spectralPartialDisabled = .false.
   partialPathDisabled = .false.
   partialAutotuneDisabled = .false.
   certifiedFullSolves = 0
   davidsonSeedCols = 0
   davidsonCarryCols = 0
   davidsonRestart = 0
   lastFullSolveTime = 0.0_wp
   thermalGuardRoots = fastPolicy%thermalGuardRoots
   if (fastProfile) call solver%reset_full_profile()

   ! GFN1-fast 2.1.5: a restarted finite-temperature wavefunction can already
   ! reveal that the Fermi edge is strongly fractional before the first SCC
   ! eigensolve. In that case skip all partial eigensolvers immediately. The
   ! full solver remains the exact reference path.
   if (xtbData%level == 1 .and. restart .and. et > 0.1_wp .and. .not.fastPolicy%thermalPartial) then
      maxFractionality = 0.0_wp
      do m = 1, ndim
         maxFractionality = max(maxFractionality, &
            & abs(focc(m)-real(nint(focc(m)),wp)))
      end do
      if (maxFractionality > real(fastPolicy%partialMaxFractionality,wp)) then
         spectralPartialDisabled = .true.
         profSpectralDisables = profSpectralDisables + 1
      end if
   end if

   ! GFN1-fast 1.4.0/2.1.7 restart seed.  Preserve the converged occupation
   ! pattern as a root-count/eigenspace seed, but 2.1.7 no longer allows it to
   ! trigger a speculative partial solve before this SCC has completed at
   ! least one exact full diagonalization.  This keeps Hessian/geometry warm
   ! starts while making the first spectrum certification exact.
   if (xtbData%level == 1 .and. restart .and. .not.spectralPartialDisabled .and. &
      & ndim >= fastPolicy%partialMinNao) then
      last_occ = max(ihomo, max(ihomoa, ihomob))
      do m = ndim, 1, -1
         if (focc(m) /= 0.0_wp) then
            last_occ = max(last_occ, m)
            exit
         end if
      end do
      if (last_occ > 0) nsub = min(ndim, last_occ + 1)
   end if

   bornResponseCache = .false.
   solvPotValid = .false.
   parallelH1 = .false.
   ! GFN1-fast 1.5.0: numerical Hessians are parallelized across independent
   ! Cartesian displacements.  Do not create a nested Hamiltonian OpenMP team
   ! inside such an outer region; keep all cores available to displacement jobs.
   !$ parallelH1 = nmat >= fastPolicy%h1OMPMinNmat .and. omp_get_max_threads() > 1 .and. .not. omp_in_parallel()
   ! number of iterations for this iterator
   thisiter = maxiter - jter

   damp = damp0
   if (present(aes)) then
      nbr = nshell + 9*n
      allocate(vs(n), vd(3, n), vq(6, n))
   else
      nbr = nshell
   end if
!  broyden data storage and init
   allocate( df(thisiter,nbr),u(thisiter,nbr),a(thisiter,thisiter), &
   &         dq(nbr),dqlast(nbr),qlast_in(nbr),omega(thisiter), &
   &         q_in(nbr),atomicShift(n), source = 0.0_wp )

   if (xtbData%level == 1 .and. allocated(solvation) .and. &
      & ndim >= fastPolicy%solventCacheMinNao) then
      select type(solvation)
      type is (TBorn)
         bornResponseCache = .true.
         allocate(solvPotIn(n), solvPotOut(n), solvPotAux(n), &
            & solvPotAuxLast(n), solvDPot(n), solvDPotLast(n), &
            & solvUAux(thisiter,n), source = 0.0_wp)
      end select
   end if

   ! GFN1-fast 0.5.0: geometry/basis-invariant Hamiltonian metadata.
   ! matlist and ao2sh do not change during SCC, so compute these once.
   allocate(matH0Idx(nmat), matISh(nmat), matJSh(nmat))
   if (xtbData%level == 1 .and. ndim >= fastPolicy%davidsonMinNao .and. &
      & useGFN1SparseDavidson(fastPolicy,ndim,nmat)) then
      allocate(pairH(nmat),pairS(nmat),pairH0(nmat))
   end if
   do m = 1, nmat
      i = matlist(1,m)
      j = matlist(2,m)
      matH0Idx(m) = j + i*(i-1)/2
      matISh(m) = ao2sh(i)
      matJSh(m) = ao2sh(j)
   enddo
   if (allocated(pairH)) then
      call buildGFN1SparsePairs(H0,S,nmat,matlist,matH0Idx,pairS,pairH0)
   end if

!! ------------------------------------------------------------------------
!  Iteration entry point
   scc_iterator: do iter = 1, thisiter

   partial_attempted = .false.
   partialFallbackReason = 0
   fullSpectrumThisIter = .false.
   partialSpectrumReady = gfn1PartialSpectrumReadyThermal(fastPolicy,certifiedFullSolves, &
      & et,spectralPartialDisabled,partialPathDisabled)
   if (.not.partialSpectrumReady .and. xtbData%level == 1 .and. .not.lastdiag .and. &
      & nsub > 0 .and. nsub < ndim .and. ndim >= fastPolicy%partialMinNao) then
      profPartialDeferred = profPartialDeferred + 1
   end if
   if (fastProfile) profStart = gfn1FastWallTime()
   ! set up ES potential
   atomicShift(:) = 0.0_wp
   shellShift(:) = externShift
   call ies%addShift(q, qsh, atomicShift, shellShift)
   ! compute potential intermediates
   if (present(aes)) then
      call setvsdq(aes,n,at,xyz,q,dipm,qp,aes%gab3,aes%gab5,vs,vd,vq)
   end if
   ! Solvation contributions
   if (allocated(solvation)) then
      cm5(:) = q + cm5a
      ! The final SCC verification iteration always recomputes A*q directly
      ! from the Born matrix, bounding any accumulated roundoff from propagated
      ! intermediate potentials without changing the convergence criteria.
      if (bornResponseCache .and. lastdiag) then
         select type(solvation)
         type is (TBorn)
            call solvation%clearSCCShiftCache()
         end select
      end if
      call solvation%addShift(env, cm5, qsh, atomicShift, shellShift)
   end if
   ! self consistent dispersion contributions
   if (present(scD4)) then
      call scD4%addShift(at, q, atomicShift)
   end if
   ! expand all atomic potentials to shell resolved potentials
   call addToShellShift(ash, atomicShift, shellShift)
   if (fastProfile) then
      profPotential = profPotential + gfn1FastWallTime() - profStart
      profStart = gfn1FastWallTime()
   end if

   ! GFN1-fast 2.1.1: defer dense Hamiltonian construction.  The Davidson
   ! path applies H directly from H0/S/shellShift and AO-pair metadata.  Only
   ! LAPACK fallback/medium-size paths materialize the dense H matrix.
   hamiltonianBuilt = .false.
   fullSolvedEarly = .false.
   matrixFreePartial = .false.

   ! ------------------------------------------------------------------------
   ! solve HC=SCemo
   ! ------------------------------------------------------------------------
   fulldiag=.false.
   if(iter.lt.startpdiag) fulldiag=.true.
   if(lastdiag)           fulldiag=.true.

   partial_used = .false.
   partial_ok = .false.
   nroots = ndim

   ! True matrix-free generalized block-Davidson.  Successful iterations never
   ! construct H(:,:) at all.  The existing parity gates remain authoritative.
   if (xtbData%level == 1 .and. .not.lastdiag .and. partialSpectrumReady .and. &
      & allocated(davidsonSeed) .and. &
      & nsub > 0 .and. nsub < ndim .and. ndim >= fastPolicy%partialMinNao) then
      required_roots = max(1,max(ihomo,max(ihomoa,ihomob))+1)
      ! GFN1-fast 2.1.4: solve one additional spectral certification root.
      ! The extra unoccupied root is never used in occupations/density; it
      ! provides a residual-checked guard above the HOMO/LUMO coverage edge.
      thermalGuardRoots = gfn1ThermalRootGuard(fastPolicy,et,thermalGuardRoots)
      nroots = min(ndim,max(nsub + thermalGuardRoots,required_roots)+1)
      if (allocated(pairH) .and. .not.davidsonDisabled .and. &
         & useGFN1Davidson(fastPolicy,ndim,nroots) .and. &
         & useGFN1SparseDavidson(fastPolicy,ndim,nmat) .and. &
         & size(davidsonSeed,2) >= nroots) then
         partial_attempted = .true.
         profPartialAttempts = profPartialAttempts + 1
         profDavidsonAttempts = profDavidsonAttempts + 1
         call ensureMatrixWorkspace(davidsonC,ndim,nroots)
         davidsonGuard = gfn1DavidsonGuardCount(nroots,ndim)
         call ensureMatrixWorkspace(davidsonCarry,ndim,nroots+davidsonGuard)
         seedCols = min(size(davidsonSeed,2),max(nroots,davidsonSeedCols))
         davidsonCarryCols = 0
         davidsonRestart = 0
         davidsonStart = gfn1FastWallTime()
         call updateGFN1PairHamiltonian(pairH0,pairS,nmat,matISh,matJSh, &
            & shellShift,pairH)
         call gfn1BlockDavidsonMF(pairH,pairS,S,nmat,matlist,davidsonSeed(:,1:seedCols), &
            & nroots,emo(1:nroots),davidsonC(:,1:nroots),partial_ok, &
            & davidsonIter,davidsonMaxRes,davidsonCarry,davidsonCarryCols,davidsonRestart)
         davidsonElapsed = gfn1FastWallTime()-davidsonStart
         profDavidsonIters = profDavidsonIters + davidsonIter
         profDavidsonRestarts = profDavidsonRestarts + davidsonRestart
         if (.not.partial_ok) davidsonDisabled = .true.
         if (lastFullSolveTime > 1.0e-6_wp .and. &
            & davidsonElapsed > 1.05_wp*lastFullSolveTime) davidsonDisabled = .true.
         nfound = merge(nroots,0,partial_ok)
         if (partial_ok) then
            call ensureMatrixWorkspace(dmatWork,ndim,nroots)
            call ensureMatrixWorkspace(eigWork,ndim,nroots)
            call validateGFN1DavidsonEigensystem(pairH,pairS,nmat,matlist, &
               & davidsonC(:,1:nroots),emo(1:nroots),dmatWork,eigWork,maxres,maxnormerr, &
               & maxadjovlp,maxprojerr,partial_ok)
         end if
         if (.not.partial_ok) partialFallbackReason = 1
         ! Optional 2.1.4 parity audit.  This deliberately pays for one full
         ! LAPACK solve and therefore is disabled by default.  It is intended
         ! for validation/CI and difficult systems: compare the Davidson roots
         ! and invariant subspace at the *same SCC Hamiltonian*, then continue
         ! from the full solution so audit mode itself cannot change results.
         if (partial_ok .and. fastPolicy%davidsonAudit) then
            profDavidsonAudits = profDavidsonAudits + 1
            if (allocated(auditEval)) deallocate(auditEval)
            allocate(auditEval(nroots))
            auditEval(:)=emo(1:nroots)
            call buildIsotropicH1Cached(ndim,nmat,matlist,matH0Idx,matISh,matJSh, &
               & H,H0,S,shellShift,parallelH1)
            if (present(aes)) call addAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp, &
               & matlist,mdlst,mqlst,H,S,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
            hamiltonianBuilt=.true.
            profFullSolves=profFullSolves+1
            fullSpectrumThisIter = .true.
            davidsonStart = gfn1FastWallTime()
            call solver%fact_solve(env,H,emo)
            lastFullSolveTime = gfn1FastWallTime() - davidsonStart
            profFullCoreTime = profFullCoreTime + lastFullSolveTime
            call env%check(fail)
            if(fail)then
               call env%error("Audit diagonalization of Hamiltonian failed",source)
               return
            endif
            auditEigMax=maxval(abs(emo(1:nroots)-auditEval(:)))
            call ensureMatrixWorkspace(eigWork,ndim,nroots)
            call applyGFN1PairMetricBlock(pairS,nmat,matlist,H(:,1:nroots),eigWork(:,1:nroots))
            auditOverlap2=0.0_wp
            do j=1,nroots
               do i=1,nroots
                  auditOverlap2=auditOverlap2 + &
                     & dot_product(davidsonC(:,i),eigWork(:,j))**2
               enddo
            enddo
            auditSubspaceError=sqrt(max(0.0_wp,1.0_wp-auditOverlap2/real(nroots,wp)))
            if (auditEigMax > 5.0e-10_wp*max(1.0_wp,maxval(abs(emo(1:nroots)))) .or. &
               & auditSubspaceError > 1.0e-7_wp) then
               davidsonDisabled=.true.
               profDavidsonAuditFailures=profDavidsonAuditFailures+1
            endif
            deallocate(auditEval)
            partialFallbackReason = 6
            partial_ok=.false.
            partial_used=.false.
            matrixFreePartial=.false.
            fullSolvedEarly=.true.
            nroots=ndim
         end if
         if (partial_ok) then
            profDavidsonAccepted = profDavidsonAccepted + 1
            partial_used = .true.
            matrixFreePartial = .true.
            H(:, :) = 0.0_wp
            H(:,1:nroots) = davidsonC(:,1:nroots)
            if (nroots < ndim) emo(nroots+1:ndim)=1.0e6_wp
         end if
      end if
   end if

   ! Medium-size exact indexed LAPACK path.  Materialize H only if Davidson was
   ! not attempted; a failed Davidson goes directly to the robust full solve.
   if (.not.partial_used .and. .not.partial_attempted .and. xtbData%level == 1 .and. &
      & .not.lastdiag .and. partialSpectrumReady .and. nsub > 0 .and. nsub < ndim) then
      required_roots = max(1,max(ihomo,max(ihomoa,ihomob))+1)
      ! GFN1-fast 2.1.4: solve one additional spectral certification root.
      ! The extra unoccupied root is never used in occupations/density; it
      ! provides a residual-checked guard above the HOMO/LUMO coverage edge.
      thermalGuardRoots = gfn1ThermalRootGuard(fastPolicy,et,thermalGuardRoots)
      nroots = min(ndim,max(nsub + thermalGuardRoots,required_roots)+1)
      if (useGFN1PartialEigensolver(fastPolicy,ndim,nroots)) then
         if (fastProfile) profStart = gfn1FastWallTime()
         call buildIsotropicH1Cached(ndim,nmat,matlist,matH0Idx,matISh,matJSh, &
            & H,H0,S,shellShift,parallelH1)
         if (present(aes)) call addAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp, &
            & matlist,mdlst,mqlst,H,S,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
         hamiltonianBuilt = .true.
         if (fastProfile) then
            profHamiltonian = profHamiltonian + gfn1FastWallTime()-profStart
            profStart = gfn1FastWallTime()
         end if
         partial_attempted = .true.
         profPartialAttempts = profPartialAttempts + 1
         P(:, :) = H(:, :)
         partialPathStart = gfn1FastWallTime()
         call solver%partial_fact_solve(env,H,nroots,emo,nfound,partial_ok)
         if (partial_ok .and. nfound == nroots) then
            call ensureMatrixWorkspace(dmatWork,ndim,nroots)
            call ensureMatrixWorkspace(eigWork,ndim,nroots)
            call validatePartialEigensystem(ndim,nroots,P,S,H,emo,dmatWork,eigWork, &
               & maxres,maxnormerr,maxadjovlp,partial_ok)
         end if
         partialPathElapsed = gfn1FastWallTime() - partialPathStart
         profPartialPathTime = profPartialPathTime + partialPathElapsed
         profPartialTimed = profPartialTimed + 1
         if (partial_ok .and. fastPolicy%partialAutotune .and. lastFullSolveTime > 1.0e-9_wp) then
            if (.not.gfn1PartialAutotuneKeep(fastPolicy,partialPathElapsed,lastFullSolveTime)) then
               partialAutotuneDisabled = .true.
            end if
         end if
         if (partial_ok) then
            partial_used=.true.
            if (nroots<ndim) then
               emo(nroots+1:ndim)=1.0e6_wp
               H(:,nroots+1:ndim)=0.0_wp
            end if
         else
            partialFallbackReason = 1
            H(:, :)=P(:, :)
         end if
      end if
   end if

   if (.not.partial_used .and. .not.fullSolvedEarly) then
      if (.not.hamiltonianBuilt) then
         if (fastProfile) profStart = gfn1FastWallTime()
         call buildIsotropicH1Cached(ndim,nmat,matlist,matH0Idx,matISh,matJSh, &
            & H,H0,S,shellShift,parallelH1)
         if (present(aes)) call addAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp, &
            & matlist,mdlst,mqlst,H,S,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
         hamiltonianBuilt=.true.
         if (fastProfile) then
            profHamiltonian=profHamiltonian+gfn1FastWallTime()-profStart
            profStart=gfn1FastWallTime()
         end if
      end if
      profFullSolves=profFullSolves+1
      fullSpectrumThisIter = .true.
      davidsonStart=gfn1FastWallTime()
      call solver%fact_solve(env,H,emo)
      lastFullSolveTime=gfn1FastWallTime()-davidsonStart
      profFullCoreTime = profFullCoreTime + lastFullSolveTime
      call env%check(fail)
      if(fail)then
         call env%error("Diagonalization of Hamiltonian failed",source)
         return
      endif
      nroots=ndim
      matrixFreePartial=.false.
   end if

   if(ihomo+1.le.ndim.and.ihomo.ge.1)egap=emo(ihomo+1)-emo(ihomo)
   ! automatic reset to small value
   if(egap.lt.0.1.and.iter.eq.0) broydamp=0.03

   call update_occupations(ndim,nel,nopen,et,emo,focc,focca,foccb, &
      & ihomoa,ihomob,nfoda,nfodb,efa,efb,ga,gb)

   ! GFN1-fast 2.1.5 spectral safety dispatcher. Partial diagonalization is
   ! counterproductive near the Fermi edge: the spectrum is sensitive, Fermi
   ! occupations become fractional, and validation/fallback overhead can
   ! exceed the factorized full solve. Use the first trustworthy spectrum to
   ! disable partial/Davidson for the remainder of this SCC. If a partial solve
   ! itself reveals the unsafe regime, redo *this same iteration* with the full
   ! factorized solver so energy/gradient never consume the partial result.
   maxFractionality = 0.0_wp
   do m = 1, ndim
      maxFractionality = max(maxFractionality, &
         & abs(focc(m)-real(nint(focc(m)),wp)))
   end do
   spectralUnsafe = et > 0.1_wp .and. .not.fastPolicy%thermalPartial .and. &
      & maxFractionality > real(fastPolicy%partialMaxFractionality,wp)
   if (ihomo >= 1 .and. ihomo < ndim) spectralUnsafe = gfn1SpectralPartialUnsafe( &
      & fastPolicy, egap, et, maxFractionality)

   if (xtbData%level == 1 .and. spectralUnsafe) then
      if (.not.spectralPartialDisabled) profSpectralDisables = profSpectralDisables + 1
      spectralPartialDisabled = .true.
      davidsonDisabled = .true.
      davidsonSeedCols = 0
      davidsonCarryCols = 0
      nsub = 0
      if (partial_used) then
         partialFallbackReason = 2
         if (matrixFreePartial) then
            call buildIsotropicH1Cached(ndim,nmat,matlist,matH0Idx,matISh,matJSh, &
               & H,H0,S,shellShift,parallelH1)
            if (present(aes)) call addAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp, &
               & matlist,mdlst,mqlst,H,S,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
            hamiltonianBuilt = .true.
         else
            H(:, :) = P(:, :)
         end if
         profFullSolves = profFullSolves + 1
         fullSpectrumThisIter = .true.
         davidsonStart = gfn1FastWallTime()
         call solver%fact_solve(env,H,emo)
         lastFullSolveTime = gfn1FastWallTime() - davidsonStart
         profFullCoreTime = profFullCoreTime + lastFullSolveTime
         call env%check(fail)
         if(fail)then
            call env%error("Low-gap fallback diagonalization failed",source)
            return
         endif
         partial_used = .false.
         matrixFreePartial = .false.
         nroots = ndim
         if(ihomo+1.le.ndim.and.ihomo.ge.1) egap=emo(ihomo+1)-emo(ihomo)
         call update_occupations(ndim,nel,nopen,et,emo,focc,focca,foccb, &
            & ihomoa,ihomob,nfoda,nfodb,efa,efb,ga,gb)
      end if
   end if

   ! GFN1-fast 0.8.0 occupation/gap coverage gate. At T=0 the number of
   ! occupied orbitals is fixed, but this check makes the assumption explicit.
   ! It also guarantees that the HOMO-LUMO gap never uses a synthetic omitted
   ! eigenvalue. Any failure restores the original Hamiltonian and reruns the
   ! exact full eigensolver.
   if (partial_used .and. nroots < ndim) then
      last_occ = 0
      do m = ndim, 1, -1
         if (focc(m) /= 0.0_wp) then
            last_occ = m
            exit
         end if
      end do
      occupation_covered = last_occ <= nroots
      if (ihomo > 0 .and. ihomo < ndim) occupation_covered = occupation_covered .and. ihomo+1 <= nroots
      if (.not.occupation_covered) then
         partialFallbackReason = 3
         if (matrixFreePartial) then
            call buildIsotropicH1Cached(ndim,nmat,matlist,matH0Idx,matISh,matJSh, &
               & H,H0,S,shellShift,parallelH1)
            if (present(aes)) call addAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp, &
               & matlist,mdlst,mqlst,H,S,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
            hamiltonianBuilt = .true.
         else
            H(:, :) = P(:, :)
         end if
         profFullSolves = profFullSolves + 1
         fullSpectrumThisIter = .true.
         davidsonStart = gfn1FastWallTime()
         call solver%fact_solve(env, H, emo)
         lastFullSolveTime = gfn1FastWallTime() - davidsonStart
         profFullCoreTime = profFullCoreTime + lastFullSolveTime
         call env%check(fail)
         if(fail)then
            call env%error("Diagonalization of Hamiltonian failed", source)
            return
         endif
         partial_used = .false.
         matrixFreePartial = .false.
         nroots = ndim
         if(ihomo+1.le.ndim.and.ihomo.ge.1)egap=emo(ihomo+1)-emo(ihomo)
         call update_occupations(ndim,nel,nopen,et,emo,focc,focca,foccb, &
            & ihomoa,ihomob,nfoda,nfodb,efa,efb,ga,gb)
      end if
   end if

   ! GFN1-fast 0.9.0 occupied-projector degeneracy gate.
   ! Degeneracies are only dangerous when they cross two different occupation
   ! weights (for example 2->1 or 1->0).  Rotations inside an equal-occupation
   ! degenerate block leave the density projector exactly invariant.  If a
   ! degenerate boundary crosses different occupations, reproduce the reference
   ! full eigensolver instead of accepting an arbitrary partial-basis rotation.
   if (partial_used) then
      occupation_degeneracy_safe = .true.
      do m = 1, min(nroots-1,ndim-1)
         boundary_scale = max(1.0_wp,abs(emo(m)),abs(emo(m+1)))
         boundary_tol = 256.0_wp*epsilon(1.0_wp)*boundary_scale
         occdiff = abs(focc(m+1)-focc(m))
         if (abs(emo(m+1)-emo(m)) <= boundary_tol .and. &
            & occdiff > 64.0_wp*epsilon(1.0_wp)) then
            occupation_degeneracy_safe = .false.
            exit
         end if
      end do
      if (.not.occupation_degeneracy_safe) then
         partialFallbackReason = 4
         if (matrixFreePartial) then
            call buildIsotropicH1Cached(ndim,nmat,matlist,matH0Idx,matISh,matJSh, &
               & H,H0,S,shellShift,parallelH1)
            if (present(aes)) call addAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp, &
               & matlist,mdlst,mqlst,H,S,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
            hamiltonianBuilt = .true.
         else
            H(:, :) = P(:, :)
         end if
         profFullSolves = profFullSolves + 1
         fullSpectrumThisIter = .true.
         davidsonStart = gfn1FastWallTime()
         call solver%fact_solve(env, H, emo)
         lastFullSolveTime = gfn1FastWallTime() - davidsonStart
         profFullCoreTime = profFullCoreTime + lastFullSolveTime
         call env%check(fail)
         if(fail)then
            call env%error("Diagonalization of Hamiltonian failed", source)
            return
         endif
         partial_used = .false.
         matrixFreePartial = .false.
         nroots = ndim
         if(ihomo+1.le.ndim.and.ihomo.ge.1)egap=emo(ihomo+1)-emo(ihomo)
         call update_occupations(ndim,nel,nopen,et,emo,focc,focca,foccb, &
            & ihomoa,ihomob,nfoda,nfodb,efa,efb,ga,gb)
      end if
   end if

   ! With finite-temperature filling, omitted eigenpairs are safe only if the
   ! highest computed root is beyond the exact zero-occupation cutoff used by
   ! fermismear. 2.2.2 keeps this at >=50 kBT (or a stricter user value). Since the omitted roots are higher
   ! than root nroots, this condition proves that all omitted occupations are
   ! exactly zero in the original algorithm.
   if (partial_used .and. et > 0.1_wp .and. nroots < ndim) then
      bkt = kB*autoev*et
      efmax = -huge(1.0_wp)
      if (ihomoa > 0) efmax = max(efmax, efa)
      if (ihomob > 0) efmax = max(efmax, efb)
      fermi_covered = efmax > -0.5_wp*huge(1.0_wp)
      if (fermi_covered) fermi_covered = emo(nroots)-efmax >= real(fastPolicy%thermalTailKBT,wp)*bkt
      if (.not.fermi_covered) then
         partialFallbackReason = 5
         if (matrixFreePartial) then
            call buildIsotropicH1Cached(ndim,nmat,matlist,matH0Idx,matISh,matJSh, &
               & H,H0,S,shellShift,parallelH1)
            if (present(aes)) call addAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp, &
               & matlist,mdlst,mqlst,H,S,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
            hamiltonianBuilt = .true.
         else
            H(:, :) = P(:, :)
         end if
         profFullSolves = profFullSolves + 1
         fullSpectrumThisIter = .true.
         davidsonStart = gfn1FastWallTime()
         call solver%fact_solve(env, H, emo)
         lastFullSolveTime = gfn1FastWallTime() - davidsonStart
         profFullCoreTime = profFullCoreTime + lastFullSolveTime
         call env%check(fail)
         if(fail)then
            call env%error("Diagonalization of Hamiltonian failed", source)
            return
         endif
         partial_used = .false.
         matrixFreePartial = .false.
         nroots = ndim
         if(ihomo+1.le.ndim.and.ihomo.ge.1)egap=emo(ihomo+1)-emo(ihomo)
         call update_occupations(ndim,nel,nopen,et,emo,focc,focca,foccb, &
            & ihomoa,ihomob,nfoda,nfodb,efa,efb,ga,gb)
      end if
   end if

   ! GFN1-fast 2.1.7: only a full solve that was reached without a
   ! speculative partial attempt certifies the spectrum for later iterations.
   ! Any real partial fallback disables further partial attempts in this SCC;
   ! retrying the same unsafe/invalid path was a measurable net loss in 2.1.6.
   if (xtbData%level == 1 .and. fullSpectrumThisIter .and. .not.partial_attempted .and. &
      & .not.spectralUnsafe .and. .not.spectralPartialDisabled) then
      certifiedFullSolves = certifiedFullSolves + 1
   end if

   fulldiag = .not.partial_used
   if (partial_attempted) then
      if (partial_used) then
         profPartialAccepted = profPartialAccepted + 1
         if (partialAutotuneDisabled .and. .not.partialPathDisabled) then
            partialPathDisabled = .true.
            profAutotuneDisables = profAutotuneDisables + 1
         end if
      else
         profPartialFallbacks = profPartialFallbacks + 1
         select case(partialFallbackReason)
         case(1)
            profFallbackSolver = profFallbackSolver + 1
            partialPathDisabled = .true.
         case(2)
            profFallbackSpectral = profFallbackSpectral + 1
            partialPathDisabled = .true.
         case(3)
            profFallbackOccupation = profFallbackOccupation + 1
            partialPathDisabled = .true.
         case(4)
            profFallbackDegeneracy = profFallbackDegeneracy + 1
            partialPathDisabled = .true.
         case(5)
            profFallbackFermi = profFallbackFermi + 1
            if (fastPolicy%thermalPartial .and. et > 0.1_wp .and. thermalFermiFallbacks < 1) then
               ! The previous full spectrum underestimated how far the thermal
               ! occupation tail moved. Grow the guard window once and retry
               ! after the exact fallback has refreshed nsub. A second miss
               ! disables partial solving for the remainder of this SCC.
               thermalFermiFallbacks = thermalFermiFallbacks + 1
               thermalGuardRoots = min(ndim,max(thermalGuardRoots + 8,2*thermalGuardRoots))
               profThermalRetries = profThermalRetries + 1
            else
               thermalFermiFallbacks = thermalFermiFallbacks + 1
               partialPathDisabled = .true.
               profThermalDisables = profThermalDisables + 1
            end if
         case(6)
            profFallbackAudit = profFallbackAudit + 1
         case default
            profFallbackSolver = profFallbackSolver + 1
            partialPathDisabled = .true.
         end select
      end if
   end if
   if (fastProfile) then
      profEigensolver = profEigensolver + gfn1FastWallTime() - profStart
      profStart = gfn1FastWallTime()
   end if

   ! GFN1-fast 0.9.0: track the occupied invariant subspace rather than an
   ! arbitrary eight-orbital virtual guard band.  One zero-occupation root is
   ! retained for the HOMO-LUMO gap and as a spectral guard.  For finite-T
   ! filling, "last_occ" already includes every orbital whose occupation is
   ! exactly non-zero under the original 50 kBT rule.
   last_occ = max(ihomo, max(ihomoa, ihomob))
   ! GFN1-fast 2.0.3: below the partial-eigensolver threshold the SCC uses a
   ! lean dense path. At zero electronic temperature the occupied range is
   ! already known exactly from ihomoa/ihomob, so avoid rescanning the full
   ! occupation vector and do not maintain partial-subspace state at all.
   if (et > 0.1_wp .or. &
      & (.not.spectralPartialDisabled .and. .not.partialPathDisabled .and. ndim >= fastPolicy%partialMinNao)) then
      do m = ndim, 1, -1
         if (focc(m) /= 0.0_wp) then
            last_occ = max(last_occ, m)
            exit
         end if
      end do
   end if
   if (.not.spectralPartialDisabled .and. .not.partialPathDisabled .and. ndim >= fastPolicy%partialMinNao) then
      nsub = min(ndim, max(1, last_occ + 1))
   else
      nsub = 0
   end if
   ! GFN1-fast 2.1.3: preserve a bounded low-energy guard subspace across SCC
   ! iterations.  Only the first nsub columns are accepted eigenvectors; the
   ! additional carried Ritz directions are seed-only and never enter the
   ! occupation, density, energy, or convergence tests.
   if (xtbData%level == 1 .and. nsub > 0 .and. nsub < ndim) then
      davidsonGuard = gfn1DavidsonGuardCount(nsub,ndim)
      if (partial_used .and. matrixFreePartial .and. davidsonCarryCols >= nsub) then
         davidsonSeedCols = min(davidsonCarryCols,nsub+davidsonGuard)
         call ensureMatrixWorkspace(davidsonSeed,ndim,davidsonSeedCols)
         davidsonSeed(:,1:davidsonSeedCols) = davidsonCarry(:,1:davidsonSeedCols)
      else if (.not.partial_used) then
         davidsonSeedCols = min(ndim,nsub+davidsonGuard)
         call ensureMatrixWorkspace(davidsonSeed,ndim,davidsonSeedCols)
         davidsonSeed(:,1:davidsonSeedCols) = H(:,1:davidsonSeedCols)
      else
         davidsonSeedCols = nsub
         call ensureMatrixWorkspace(davidsonSeed,ndim,davidsonSeedCols)
         davidsonSeed(:,1:davidsonSeedCols) = H(:,1:davidsonSeedCols)
      end if
   else
      davidsonSeedCols = 0
   end if
   ! GFN1-fast 1.6.0: the old X(:,1:last_occ) MO copy was retained
   ! from an abandoned subspace-seeding path. SYEVX does not consume
   ! it, so removing the copy saves NAO*Nocc memory traffic exactly.

   ! save q
   q_in(1:nshell)=qsh(1:nshell)
   if (present(aes)) then
      k=nshell
      call gfn2broyden_save(n,k,nbr,dipm,qp,q_in)
   end if

   if (fastProfile) then
      profSubspace = profSubspace + gfn1FastWallTime() - profStart
      profStart = gfn1FastWallTime()
   end if

   ! GFN1-fast 2.2.6: intermediate native-GFN1 iterations need only the
   ! upper density triangle because both Mulliken populations and Tr(P H0)
   ! consume P(j,i), j<=i.  At finite electronic temperature build that
   ! triangle directly as P = (C sqrt(f))(C sqrt(f))^T with SYRK.  The final
   ! verification iteration keeps the legacy full-density path so the
   ! converged Energy/Gradient reference remains unchanged.
   compactThisIter = compactDensityActive .and. .not.lastdiag
   call ensureMatrixWorkspace(dmatWork,ndim,max(1,last_occ))
   call dmat(ndim,focc,H,P,dmatWork,compactThisIter)
   if (compactThisIter) profCompactDensityIterations = profCompactDensityIterations + 1
   if (fastProfile) then
      profDensity = profDensity + gfn1FastWallTime() - profStart
      profStart = gfn1FastWallTime()
   end if

   ! Fuse the two identical upper-triangle traversals previously performed by
   ! mpopsh() and electro(): shell populations and the fixed-H0 trace are
   ! accumulated in the original summation order in one pass.
   call mpopsh_h0_upper(ndim,nshell,ao2sh,S,P,H0,qsh,h0EnergyCompact)
   qsh = zsh - qsh
   call qsh2qat(ash,qsh,q)
   if (fastProfile) then
      profPopulation = profPopulation + gfn1FastWallTime() - profStart
      profStart = gfn1FastWallTime()
   end if

   eold=eel
   call ies%getEnergy(q,qsh,ees)
   eel = ees + h0EnergyCompact
   ! multipole electrostatic
   if (present(aes)) then
      call mmompop(n,ndim,aoat2,xyz,p,s,dpint,qpint,dipm,qp)
      ! evaluate energy
      call aniso_electro(aes,n,at,xyz,q,dipm,qp,aes%gab3,aes%gab5,eaes,epol)
      eel=eel+eaes+epol
   end if

   ! Self consistent dispersion
   if (present(scD4)) then
      call scD4%getEnergy(at, q, ed)
      eel = eel + ed
   endif

   ! point charge contribution
   if(pcem) then
      call electro_pcem(nshell,qsh,externShift,epcem,eel)
   end if

   ! new cm5 charges and gborn energy
   if (allocated(solvation)) then
      cm5=q+cm5a
      call solvation%getEnergy(env, cm5, qsh, gborn)
      eel = eel + gborn
      solvPotValid = .false.
      if (bornResponseCache) then
         select type(solvation)
         type is (TBorn)
            call solvation%getSCCPotentials(solvPotIn,solvPotOut,solvPotValid)
         end select
      end if
   end if

   ! add el. entropies*T
   eel=eel+ga+gb
   if (fastProfile) then
      profEnergy = profEnergy + gfn1FastWallTime() - profStart
      profStart = gfn1FastWallTime()
   end if

   ! ------------------------------------------------------------------------
   ! check for energy convergence
   econverged = abs(eel - eold) < scfconv
   ! ------------------------------------------------------------------------

   dq(1:nshell)=qsh(1:nshell)-q_in(1:nshell)
   if (present(aes)) then
      k=nshell
      call gfn2broyden_diff(n,k,nbr,dipm,qp,q_in,dq) ! CAMM case
   end if
   rmsq=sum(dq(1:nbr)**2)/dble(n)
   rmsq=sqrt(rmsq)

   ! ------------------------------------------------------------------------
   ! end of SCC convergence part
   qconverged = rmsq < qconv
   ! ------------------------------------------------------------------------
   if (fastProfile) then
      profConvergence = profConvergence + gfn1FastWallTime() - profStart
      profStart = gfn1FastWallTime()
   end if

   ! SCC convergence acceleration
   if(.not.broy)then

      ! simple damp
      if(iter.gt.0) then
         omegap=egap
         dampUsed = damp
         ! monopoles only
         do i=1,nshell
            qsh(i)=dampUsed*qsh(i)+(1.0d0-dampUsed)*q_in(i)
         enddo
         if (solvPotValid) then
            solvPotAux(:) = dampUsed*solvPotOut(:) + &
               & (1.0_wp-dampUsed)*solvPotIn(:)
         end if
         if (present(aes)) then
            ! CAMM
            k=nshell
            do i=1,n
               do j=1,3
                  k=k+1
                  dipm(j,i)=damp*dipm(j,i)+(1.0d0-damp)*q_in(k)
               enddo
               do j=1,6
                  k=k+1
                  qp(j,i)=damp*qp(j,i)+(1.0d0-damp)*q_in(k)
               enddo
            enddo
         end if
         if(eel-eold.lt.0) then
            damp=damp*1.15
         else
            damp=damp0
         endif
         damp=min(damp,1.0)
         if(egap.lt.1.0)damp=min(damp,0.5)
      endif

   else

      ! Broyden mixing
      omegap=0.0d0
      if (solvPotValid) then
         solvPotAux(:) = solvPotIn(:)
         solvDPot(:) = solvPotOut(:)-solvPotIn(:)
         call broyden(nbr,q_in,qlast_in,dq,dqlast,iter,thisiter,broydamp,omega,df,u,a, &
            & solvPotAux,solvPotAuxLast,solvDPot,solvDPotLast,solvUAux)
      else
         call broyden(nbr,q_in,qlast_in,dq,dqlast,iter,thisiter,broydamp,omega,df,u,a)
      end if
      qsh(1:nshell)=q_in(1:nshell)
      if (present(aes)) then
         k=nshell
         call gfn2broyden_out(n,k,nbr,q_in,dipm,qp) ! CAMM case
      end if
      if(iter.gt.1) omegap=omega(iter-1)
   endif ! Broyden?

   call qsh2qat(ash, qsh, q) !new qat

   if(allocated(solvation)) then
      cm5 = q+cm5a
      if (bornResponseCache .and. solvPotValid) then
         select type(solvation)
         type is (TBorn)
            call solvation%cacheSCCShift(cm5,solvPotAux)
         end select
      end if
   end if

   if (fastProfile) profMixing = profMixing + gfn1FastWallTime() - profStart

   if(minpr)write(env%unit,'(i4,F15.7,E14.6,E11.3,f8.2,2x,f8.1,l3)') &
   &  iter+jter,eel,eel-eold,rmsq,egap,omegap,fulldiag
   qq=q

!  end of SCC convergence part

!! ------------------------------------------------------------------------
   if (econverged.and.qconverged) then
      converged = .true.
      if (lastdiag) exit scc_iterator
      lastdiag = .true.
   endif
!! ------------------------------------------------------------------------

   enddo scc_iterator

   ! If SCC exits without the normal full-density verification iteration,
   ! materialize a complete P for downstream diagnostics/properties.
   if (compactDensityActive .and. .not.converged) then
      last_occ = 0
      do m = ndim, 1, -1
         if (focc(m) /= 0.0_wp) then
            last_occ = m
            exit
         end if
      end do
      call ensureMatrixWorkspace(dmatWork,ndim,max(1,last_occ))
      call dmat(ndim,focc,H,P,dmatWork,.false.)
   end if

   jter = jter + min(iter,thisiter)
   fail = .not.converged

   if (fastProfile) then
      call solver%get_full_profile(profFullReduceTime,profFullDiagTime,profFullBackTime,profFullPhaseCalls)
      call solver%get_full_backend_profile(profFullBackendRequested,profFullBackendSelected, &
         & profFullBackendAutotune,profFullBackendMaxRatio,profFullSYEVDTime,profFullSYEVRTime, &
         & profFullSYEVDCalls,profFullSYEVRCalls,profFullBackendTrials,profFullBackendFallbacks)
      select case (profFullBackendRequested)
      case (1)
         profFullBackendRequestedName = 'syevd'
      case (2)
         profFullBackendRequestedName = 'syevr'
      case default
         profFullBackendRequestedName = 'auto'
      end select
      select case (profFullBackendSelected)
      case (2)
         profFullBackendSelectedName = 'syevr'
      case default
         profFullBackendSelectedName = 'syevd'
      end select
      write(env%unit,'(/,1x,a)') 'GFN1-fast 2.2.6 SCC profile'
      write(env%unit,'(3x,a,i0,a,i0,a,i0)') 'NAO=',ndim,', partial threshold=',fastPolicy%partialMinNao, &
         & ', Davidson threshold=',fastPolicy%davidsonMinNao
      write(env%unit,'(3x,a,f8.4,a,f6.2,a,es10.3,a,l1)') 'partial min gap=', &
         & real(fastPolicy%partialMinGapEV,wp),' eV, min gap/kBT=', &
         & real(fastPolicy%partialMinGapKBT,wp),', max fractionality=', &
         & real(fastPolicy%partialMaxFractionality,wp),', spectral disabled=',spectralPartialDisabled
      write(env%unit,'(3x,a,l1,a,i0,a,i0,a,f6.1,a,i0)') 'thermal partial=',fastPolicy%thermalPartial, &
         & ', warmup full=',fastPolicy%thermalWarmupFullSolves,', guard roots=',thermalGuardRoots, &
         & ', tail cutoff=',real(fastPolicy%thermalTailKBT,wp),' kBT, Fermi fallbacks=',thermalFermiFallbacks
      write(env%unit,'(3x,a,l1,a,f6.3,a,l1,a,i0)') 'partial autotune=',fastPolicy%partialAutotune, &
         & ', max ratio=',real(fastPolicy%partialAutotuneMaxRatio,wp),', disabled=',partialAutotuneDisabled, &
         & ', disables=',profAutotuneDisables
      write(env%unit,'(3x,a,i0)') 'partial probe min NAO=',fastPolicy%partialProbeMinNao
      write(env%unit,'(3x,a,l1,a,i0,a,i0)') 'compact density SCC=',compactDensityActive, &
         & ', min NAO=',fastPolicy%compactDensityMinNao,', iterations=',profCompactDensityIterations
      write(env%unit,'(3x,a,l1,a,i0)') 'H1 OpenMP=',parallelH1,', nmat threshold=',fastPolicy%h1OMPMinNmat
      write(env%unit,'(3x,a,i0,a,i0,a,l1)') 'AO pairs=',nmat,', sparse max %=', &
         & fastPolicy%sparseMaxPairPercent,', sparse Davidson=', &
         & useGFN1SparseDavidson(fastPolicy,ndim,nmat)
      write(env%unit,'(3x,a,l1,a,i0)') 'solvent SCC cache=',bornResponseCache, &
         & ', NAO threshold=',fastPolicy%solventCacheMinNao
      write(env%unit,'(3x,a,f10.6,a)') 'potential/shift : ',profPotential,' s'
      write(env%unit,'(3x,a,f10.6,a)') 'Hamiltonian     : ',profHamiltonian,' s'
      write(env%unit,'(3x,a,f10.6,a)') 'eigensolver     : ',profEigensolver,' s'
      write(env%unit,'(3x,a,f10.6,a,f10.6,a)') 'full solve core : ',profFullCoreTime,' s, avg=', &
         & profFullCoreTime/real(max(1,profFullSolves),wp),' s'
      write(env%unit,'(3x,a,f10.6,a,f10.6,a)') '  SYGST reduce  : ',profFullReduceTime,' s, avg=', &
         & profFullReduceTime/real(max(1,profFullPhaseCalls),wp),' s'
      write(env%unit,'(3x,a,f10.6,a,f10.6,a)') '  SYEVD diagonal: ',profFullDiagTime,' s, avg=', &
         & profFullDiagTime/real(max(1,profFullPhaseCalls),wp),' s'
      write(env%unit,'(3x,a,f10.6,a,f10.6,a)') '  TRSM backtrans: ',profFullBackTime,' s, avg=', &
         & profFullBackTime/real(max(1,profFullPhaseCalls),wp),' s'
      write(env%unit,'(3x,a,a,a,a,a,l1,a,f6.3)') 'full eig backend: requested=', &
         & trim(profFullBackendRequestedName),', selected=',trim(profFullBackendSelectedName), &
         & ', autotune=',profFullBackendAutotune,', max ratio=',profFullBackendMaxRatio
      write(env%unit,'(3x,a,i0)') 'full eig autotune min NAO=',fastPolicy%fullEigAutotuneMinNao
      write(env%unit,'(3x,a,f10.6,a,i0,a,f10.6,a,i0)') '  SYEVD total   : ',profFullSYEVDTime, &
         & ' s, calls=',profFullSYEVDCalls,', SYEVR total=',profFullSYEVRTime, &
         & ' s, calls=',profFullSYEVRCalls
      write(env%unit,'(3x,a,i0,a,i0)') '  backend trials=',profFullBackendTrials, &
         & ', fallbacks=',profFullBackendFallbacks
      write(env%unit,'(3x,a,f10.6,a,f10.6,a,i0)') 'partial path    : ',profPartialPathTime,' s, avg=', &
         & profPartialPathTime/real(max(1,profPartialTimed),wp),' s, timed=',profPartialTimed
      write(env%unit,'(3x,a,f10.6,a)') 'subspace/setup  : ',profSubspace,' s'
      write(env%unit,'(3x,a,f10.6,a)') 'density matrix  : ',profDensity,' s'
      write(env%unit,'(3x,a,f10.6,a)') 'Mulliken/pop.   : ',profPopulation,' s'
      write(env%unit,'(3x,a,f10.6,a)') 'energy terms    : ',profEnergy,' s'
      write(env%unit,'(3x,a,f10.6,a)') 'convergence     : ',profConvergence,' s'
      write(env%unit,'(3x,a,f10.6,a)') 'mixing          : ',profMixing,' s'
      write(env%unit,'(3x,a,i0,a,i0,a,i0,a,i0)') 'partial attempts=',profPartialAttempts, &
         & ', accepted=',profPartialAccepted,', fallback=',profPartialFallbacks, &
         & ', full solves=',profFullSolves
      write(env%unit,'(3x,a,i0,a,i0,a,l1)') 'partial warmup full solves=', &
         & fastPolicy%partialWarmupFullSolves,', certified=',certifiedFullSolves, &
         & ', path disabled=',partialPathDisabled
      write(env%unit,'(3x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'partial deferred=',profPartialDeferred, &
         & ', fail solver/validation=',profFallbackSolver,', spectral=',profFallbackSpectral, &
         & ', occupation=',profFallbackOccupation,', degeneracy=',profFallbackDegeneracy, &
         & ', Fermi=',profFallbackFermi,', audit=',profFallbackAudit
      write(env%unit,'(3x,a,i0)') 'spectral safety disables=',profSpectralDisables
      write(env%unit,'(3x,a,i0,a,i0)') 'thermal retries=',profThermalRetries,', thermal disables=',profThermalDisables
      write(env%unit,'(3x,a,i0,a,i0,a,i0,a,i0,a,l1)') 'Davidson attempts=',profDavidsonAttempts, &
         & ', accepted=',profDavidsonAccepted,', total iterations=',profDavidsonIters, &
         & ', thick restarts=',profDavidsonRestarts,', auto-disabled=',davidsonDisabled
      if (profDavidsonAudits > 0) write(env%unit,'(3x,a,i0,a,i0,a,es10.3,a,es10.3)') &
         & 'Davidson audits=',profDavidsonAudits,', failures=',profDavidsonAuditFailures, &
         & ', max eig diff=',auditEigMax,', subspace err=',auditSubspaceError
   end if

end subroutine scc


real(wp) function gfn1FastWallTime() result(t)
   integer(kind=8) :: count, rate
   call system_clock(count=count, count_rate=rate)
   if (rate > 0) then
      t = real(count,wp)/real(rate,wp)
   else
      call cpu_time(t)
   end if
end function gfn1FastWallTime


subroutine update_occupations(ndim,nel,nopen,et,emo,focc,focca,foccb, &
      & ihomoa,ihomob,nfoda,nfodb,efa,efb,ga,gb)
   integer, intent(in) :: ndim, nel, nopen
   real(wp), intent(in) :: et, emo(ndim)
   real(wp), intent(inout) :: focc(ndim), focca(ndim), foccb(ndim)
   integer, intent(inout) :: ihomoa, ihomob
   real(wp), intent(out) :: nfoda, nfodb, ga, gb
   real(wp), intent(inout) :: efa, efb

   nfoda = 0.0_wp
   nfodb = 0.0_wp
   if(et.gt.0.1_wp)then
      if(nel.gt.0) then
         call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      else
         focca=0.0_wp
         foccb=0.0_wp
         ihomoa=0
         ihomob=0
      endif
      if (ihomoa+1.le.ndim) then
         call fermismear(.false.,ndim,ihomoa,et,emo,focca,nfoda,efa,ga)
      else
         ga = 0.0_wp
      endif
      if (ihomob+1.le.ndim) then
         call fermismear(.false.,ndim,ihomob,et,emo,foccb,nfodb,efb,gb)
      else
         gb = 0.0_wp
      endif
      focc = focca + foccb
   else
      ga = 0.0_wp
      gb = 0.0_wp
   endif
end subroutine update_occupations


!> Ensure a compact NAO x Nactive scratch matrix exists.
!> The contents are disposable; only capacity matters.
subroutine ensureMatrixWorkspace(work,nrow,ncol)
   real(wp), allocatable, intent(inout) :: work(:,:)
   integer, intent(in) :: nrow, ncol
   integer :: needcol

   needcol = max(1,ncol)
   if (.not.allocated(work)) then
      allocate(work(nrow,needcol))
   else if (size(work,1) /= nrow .or. size(work,2) < needcol) then
      deallocate(work)
      allocate(work(nrow,needcol))
   end if
end subroutine ensureMatrixWorkspace


!> Validate a partial solution of H C = S C e without approximations.
!> The result is accepted only when all requested roots have a small relative
!> residual, full S-orthogonality, projected-energy consistency, and finite ordered eigenvalues.
subroutine validatePartialEigensystem(ndim,nroots,Horig,S,C,eval,workHC,workSC, &
      & maxRelResidual,maxNormError,maxOrthOverlap,success)
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   integer, intent(in) :: ndim, nroots
   real(wp), intent(in) :: Horig(ndim,ndim), S(ndim,ndim)
   real(wp), intent(in) :: C(:,:), eval(ndim)
   real(wp), intent(inout) :: workHC(:,:), workSC(:,:)
   real(wp), intent(out) :: maxRelResidual, maxNormError, maxOrthOverlap
   logical, intent(out) :: success
   integer :: i,j
   real(wp) :: r2,hc2,sc2,denom,relres,escale,hscale,maxProjectedError
   real(wp), allocatable :: gram(:,:)
   real(wp), parameter :: residualTol=1.0e-11_wp,normTol=1.0e-10_wp
   real(wp), parameter :: orthTol=1.0e-9_wp,projectedTol=2.0e-10_wp

   success=.false.; maxRelResidual=huge(1.0_wp); maxNormError=huge(1.0_wp)
   maxOrthOverlap=huge(1.0_wp)
   if(nroots<1 .or. nroots>ndim)return

   call mctc_gemm(Horig,C(:,1:nroots),workHC(:,1:nroots))
   ! 2.1.4 bug fix: the medium-size SYEVX path does not allocate the sparse
   ! Davidson pair metric.  Validate it against its actual dense S matrix.
   call mctc_gemm(S,C(:,1:nroots),workSC(:,1:nroots))

   maxRelResidual=0.0_wp
   do j=1,nroots
      if(.not.ieee_is_finite(eval(j)))return
      if(j>1)then
         escale=max(1.0_wp,abs(eval(j)),abs(eval(j-1)))
         if(eval(j)<eval(j-1)-64.0_wp*epsilon(1.0_wp)*escale)return
      endif
      r2=sum((workHC(:,j)-eval(j)*workSC(:,j))**2)
      hc2=sum(workHC(:,j)**2);sc2=sum(workSC(:,j)**2)
      if(.not.ieee_is_finite(r2).or..not.ieee_is_finite(hc2).or..not.ieee_is_finite(sc2))return
      denom=sqrt(max(0.0_wp,hc2))+abs(eval(j))*sqrt(max(0.0_wp,sc2))+tiny(1.0_wp)
      relres=sqrt(max(0.0_wp,r2))/denom
      if(.not.ieee_is_finite(relres))return
      maxRelResidual=max(maxRelResidual,relres)
   enddo

   allocate(gram(nroots,nroots))
   call mctc_gemm(C(:,1:nroots),workSC(:,1:nroots),gram,transa='t')
   maxNormError=0.0_wp;maxOrthOverlap=0.0_wp;maxProjectedError=0.0_wp
   do j=1,nroots
      if(.not.ieee_is_finite(gram(j,j)))then
         deallocate(gram);return
      endif
      maxNormError=max(maxNormError,abs(gram(j,j)-1.0_wp))
      hscale=max(1.0_wp,abs(eval(j)))
      maxProjectedError=max(maxProjectedError, &
         & abs(dot_product(C(:,j),workHC(:,j))-eval(j))/hscale)
      do i=1,nroots
         if(.not.ieee_is_finite(gram(i,j)))then
            deallocate(gram);return
         endif
         if(i/=j) maxOrthOverlap=max(maxOrthOverlap,abs(gram(i,j)))
      enddo
   enddo
   success=maxRelResidual<=residualTol.and.maxNormError<=normTol.and. &
      & maxOrthOverlap<=orthTol.and.maxProjectedError<=projectedTol
   deallocate(gram)
end subroutine validatePartialEigensystem




!> H0 off-diag scaling
subroutine h0scal(hData,il,jl,izp,jzp,valaoi,valaoj,km)
  !$acc routine seq
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in)  :: il
   integer, intent(in)  :: jl
   integer, intent(in)  :: izp
   integer, intent(in)  :: jzp
   logical, intent(in)  :: valaoi
   logical, intent(in)  :: valaoj
   real(wp),intent(out) :: km
   real(wp) :: den, enpoly

   km = 0.0_wp

!  valence
   if(valaoi.and.valaoj) then
      den=(hData%electronegativity(izp)-hData%electronegativity(jzp))**2
      enpoly = (1.0_wp+hData%enScale(jl-1,il-1)*den*(1.0_wp+hData%enScale4*den))
      km=hData%kScale(jl-1,il-1)*enpoly*hData%pairParam(izp,jzp)
      return
   endif

!  "DZ" functions (on H for GFN or 3S for EA calc on all atoms)
   if((.not.valaoi).and.(.not.valaoj)) then
      km=hData%kDiff
      return
   endif
   if(.not.valaoi.and.valaoj) then
      km=0.5*(hData%kScale(jl-1,jl-1)+hData%kDiff)
      return
   endif
   if(.not.valaoj.and.valaoi) then
      km=0.5*(hData%kScale(il-1,il-1)+hData%kDiff)
   endif


end subroutine h0scal


!! ========================================================================
!  total energy for GFN1
!! ========================================================================
pure subroutine electro(n,at,nbf,nshell,ies,H0,P,dq,dqsh,es,scc)
   use xtb_mctc_convert, only : evtoau
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   integer, intent(in) :: nbf
   integer, intent(in) :: nshell
   real(wp),intent(in)  :: H0(nbf*(nbf+1)/2)
   real(wp),intent(in)  :: P (nbf,nbf)
   type(TxTBCoulomb), intent(inout) :: ies
   real(wp),intent(in)  :: dq(n)
   real(wp),intent(in)  :: dqsh(nshell)
   real(wp),intent(out) :: es
   real(wp),intent(out) :: scc
   real(wp) :: ehb ! not used

   integer  :: i,j,k
   real(wp) :: h,t

   ! ES energy in Eh
   call ies%getEnergy(dq, dqsh, es)

!  H0 part
   k=0
   h=0.0_wp
   do i=1,nbf
      do j=1,i-1
         k=k+1
         h=h+P(j,i)*H0(k)
      enddo
      k=k+1
      h=h+P(i,i)*H0(k)*0.5_wp
   enddo

!  Etotal in Eh
   scc = es + 2.0_wp*h*evtoau

end subroutine electro


!! ========================================================================
!  S(R) enhancement factor
!! ========================================================================
pure function shellPoly(iPoly,jPoly,iRad,jRad,xyz1,xyz2)
   use xtb_mctc_convert, only : aatoau
   !$acc routine seq
   real(wp), intent(in) :: iPoly,jPoly
   real(wp), intent(in) :: iRad,jRad
   real(wp), intent(in) :: xyz1(3),xyz2(3)
   real(wp) :: shellPoly
   real(wp) :: rab,k1,rr,r,rf1,rf2,dx,dy,dz,a

   a=0.5           ! R^a dependence 0.5 in GFN1

   dx=xyz1(1)-xyz2(1)
   dy=xyz1(2)-xyz2(2)
   dz=xyz1(3)-xyz2(3)

   rab=sqrt(dx**2+dy**2+dz**2)

   ! this sloppy conv. factor has been used in development, keep it
   rr=jRad+iRad

   r=rab/rr

   rf1=1.0d0+0.01*iPoly*r**a
   rf2=1.0d0+0.01*jPoly*r**a

   shellPoly= rf1*rf2

end function shellPoly

!! ========================================================================
!  set up Coulomb potential due to 2nd order fluctuation
!! ========================================================================
pure subroutine setespot(nshell, qsh, jmat, shellShift)

   !> Dimension of the Coulomb matrix
   integer, intent(in) :: nshell

   !> Shell resolved partial charges
   real(wp),intent(in) :: qsh(:)

   !> Coulomb matrix
   real(wp),intent(in) :: jmat(:, :)

   !> shell-resolved potential shift
   real(wp),intent(inout) :: shellShift(:)

   !call blas_symv('l', nshell, 1.0_wp, jmat, nshell, qsh, 1, 1.0_wp, shellShift, 1)
   call mctc_symv(jmat, qsh, shellShift, beta=1.0_wp)

end subroutine setespot


pure subroutine addToShellShift(shellAtom, atomicShift, shellShift)

   !> Atom each shell is centered on
   integer, intent(in) :: shellAtom(:)

   !> Atomic potential shift
   real(wp), intent(in) :: atomicShift(:)

   !> Shell-resolved potential shift
   real(wp), intent(inout) :: shellShift(:)

   integer :: ii

   do ii = 1, size(shellShift, dim=1)
      shellShift(ii) = shellShift(ii) + atomicShift(shellAtom(ii))
   end do

end subroutine addToShellShift


!! ========================================================================
!  eigenvalue solver single-precision
!! ========================================================================
subroutine solve4(full,ndim,ihomo,acc,H,S,X,P,e,fail)
   use xtb_mctc_accuracy, only : sp
   integer, intent(in)   :: ndim
   logical, intent(in)   :: full
   integer, intent(in)   :: ihomo
   real(wp),intent(inout):: H(ndim,ndim)
   real(wp),intent(in)   :: S(ndim,ndim)
   real(wp),intent(out)  :: X(ndim,ndim)
   real(wp),intent(out)  :: P(ndim,ndim)
   real(wp),intent(out)  :: e(ndim)
   real(wp),intent(in)   :: acc
   logical, intent(out)  :: fail

   integer i,j,info,lwork,liwork,nfound,iu,nbf
   integer, allocatable :: iwork(:),ifail(:)
   real(wp),allocatable :: aux  (:)
   real(wp) w0,w1,t0,t1

   real(sp),allocatable :: H4(:,:)
   real(sp),allocatable :: S4(:,:)
   real(sp),allocatable :: X4(:,:)
   real(sp),allocatable :: P4(:,:)
   real(sp),allocatable :: e4(:)
   real(sp),allocatable :: aux4(:)


   allocate(H4(ndim,ndim),S4(ndim,ndim))
   allocate(X4(ndim,ndim),P4(ndim,ndim),e4(ndim))

   H4 = H
   S4 = S

   fail =.false.
!  standard first full diag call
   if(full) then
!                                                     call timing(t0,w0)
!     if(ndim.gt.0)then
!     USE DIAG IN NON-ORTHORGONAL BASIS
      allocate (aux4(1),iwork(1),ifail(ndim))
      P4 = s4
      call lapack_sygvd(1,'v','u',ndim,h4,ndim,p4,ndim,e4,aux4, &!workspace query
     &           -1,iwork,liwork,info)
      lwork=int(aux4(1))
      liwork=iwork(1)
      deallocate(aux4,iwork)
      allocate (aux4(lwork),iwork(liwork))              !do it
      call lapack_sygvd(1,'v','u',ndim,h4,ndim,p4,ndim,e4,aux4, &
     &           lwork,iwork,liwork,info)
      if(info.ne.0) then
         fail=.true.
         return
      endif
      X4 = H4 ! save
      deallocate(aux4,iwork,ifail)

!     else
!        USE DIAG IN ORTHOGONAL BASIS WITH X=S^-1/2 TRAFO
!        nbf = ndim
!        lwork  = 1 + 6*nbf + 2*nbf**2
!        allocate (aux(lwork))
!        call blas_gemm('N','N',nbf,nbf,nbf,1.0d0,H,nbf,X,nbf,0.0d0,P,nbf)
!        call blas_gemm('T','N',nbf,nbf,nbf,1.0d0,X,nbf,P,nbf,0.0d0,H,nbf)
!        call SYEV('V','U',nbf,H,nbf,e,aux,lwork,info)
!        if(info.ne.0) error stop 'diag error'
!        call blas_gemm('N','N',nbf,nbf,nbf,1.0d0,X,nbf,H,nbf,0.0d0,P,nbf)
!        H = P
!        deallocate(aux)
!     endif
!                                                     call timing(t1,w1)
!                                    call prtime(6,t1-t0,w1-w0,'dsygvd')

   else
!                                                     call timing(t0,w0)
!     go to MO basis using trafo(X) from first iteration (=full diag)
!      call blas_gemm('N','N',ndim,ndim,ndim,1.d0,H4,ndim,X4,ndim,0.d0,P4,ndim)
!      call blas_gemm('T','N',ndim,ndim,ndim,1.d0,X4,ndim,P4,ndim,0.d0,H4,ndim)
!                                                     call timing(t1,w1)
!                       call prtime(6,1.5*(t1-t0),1.5*(w1-w0),'3xdgemm')
!                                                     call timing(t0,w0)
!      call pseudodiag(ndim,ihomo,H4,e4)
!                                                     call timing(t1,w1)
!                                call prtime(6,t1-t0,w1-w0,'pseudodiag')

!     C = X C', P=scratch
!      call blas_gemm('N','N',ndim,ndim,ndim,1.d0,X4,ndim,H4,ndim,0.d0,P4,ndim)
!     save and output MO matrix in AO basis
!      H4 = P4
   endif

   H = H4
   P = P4
   X = X4
   e = e4

   deallocate(e4,P4,X4,S4,H4)

end subroutine solve4


!! ========================================================================
!  eigenvalue solver
!! ========================================================================
subroutine solve(full,ndim,ihomo,acc,H,S,X,P,e,fail)
   integer, intent(in)   :: ndim
   logical, intent(in)   :: full
   integer, intent(in)   :: ihomo
   real(wp),intent(inout):: H(ndim,ndim)
   real(wp),intent(in)   :: S(ndim,ndim)
   real(wp),intent(out)  :: X(ndim,ndim)
   real(wp),intent(out)  :: P(ndim,ndim)
   real(wp),intent(out)  :: e(ndim)
   real(wp),intent(in)   :: acc
   logical, intent(out)  :: fail

   integer i,j,info,lwork,liwork,nfound,iu,nbf
   integer, allocatable :: iwork(:),ifail(:)
   real(wp),allocatable :: aux  (:)
   real(wp) w0,w1,t0,t1

   fail =.false.

!  standard first full diag call
   if(full) then
!                                                     call timing(t0,w0)
!     if(ndim.gt.0)then
!     USE DIAG IN NON-ORTHORGONAL BASIS
      allocate (aux(1),iwork(1),ifail(ndim))
      P = s
      call lapack_sygvd(1,'v','u',ndim,h,ndim,p,ndim,e,aux, &!workspace query
     &           -1,iwork,liwork,info)
      lwork=int(aux(1))
      liwork=iwork(1)
      deallocate(aux,iwork)
      allocate (aux(lwork),iwork(liwork))              !do it
      call lapack_sygvd(1,'v','u',ndim,h,ndim,p,ndim,e,aux, &
     &           lwork,iwork,liwork,info)
      !write(*,*)'SYGVD INFO', info
      if(info.ne.0) then
         fail=.true.
         return
      endif
      X = H ! save
      deallocate(aux,iwork,ifail)

!     else
!        USE DIAG IN ORTHOGONAL BASIS WITH X=S^-1/2 TRAFO
!        nbf = ndim
!        lwork  = 1 + 6*nbf + 2*nbf**2
!        allocate (aux(lwork))
!        call blas_gemm('N','N',nbf,nbf,nbf,1.0d0,H,nbf,X,nbf,0.0d0,P,nbf)
!        call blas_gemm('T','N',nbf,nbf,nbf,1.0d0,X,nbf,P,nbf,0.0d0,H,nbf)
!        call SYEV('V','U',nbf,H,nbf,e,aux,lwork,info)
!        if(info.ne.0) error stop 'diag error'
!        call blas_gemm('N','N',nbf,nbf,nbf,1.0d0,X,nbf,H,nbf,0.0d0,P,nbf)
!        H = P
!        deallocate(aux)
!     endif
!                                                     call timing(t1,w1)
!                                    call prtime(6,t1-t0,w1-w0,'dsygvd')

   else
!                                                     call timing(t0,w0)
!     go to MO basis using trafo(X) from first iteration (=full diag)
      call blas_gemm('N','N',ndim,ndim,ndim,1.d0,H,ndim,X,ndim,0.d0,P,ndim)
      call blas_gemm('T','N',ndim,ndim,ndim,1.d0,X,ndim,P,ndim,0.d0,H,ndim)
!                                                     call timing(t1,w1)
!                       call prtime(6,1.5*(t1-t0),1.5*(w1-w0),'3xdgemm')
!                                                     call timing(t0,w0)
      call pseudodiag(ndim,ihomo,H,e)
!                                                     call timing(t1,w1)
!                                call prtime(6,t1-t0,w1-w0,'pseudodiag')

!     C = X C', P=scratch
      call blas_gemm('N','N',ndim,ndim,ndim,1.d0,X,ndim,H,ndim,0.d0,P,ndim)
!     save and output MO matrix in AO basis
      H = P
   endif

end subroutine solve


subroutine fermismear(prt,norbs,nel,t,eig,occ,fod,e_fermi,s)
   use xtb_mctc_convert, only : autoev
   use xtb_mctc_constants, only : kB
   integer, intent(in)  :: norbs
   integer, intent(in)  :: nel
   real(wp),intent(in)  :: eig(norbs)
   real(wp),intent(out) :: occ(norbs)
   real(wp),intent(in)  :: t
   real(wp),intent(out) :: fod
   real(wp),intent(out) :: e_fermi
   logical, intent(in)  :: prt

   real(wp) :: boltz,bkt,occt,total_number,thr
   real(wp) :: total_dfermi,dfermifunct,fermifunct,s,change_fermi

   parameter (boltz = kB*autoev)
   parameter (thr   = 1.d-9)
   integer :: ncycle,i,j,m,k,i1,i2

   bkt = boltz*t

   e_fermi = 0.5*(eig(nel)+eig(nel+1))
   occt=nel

   do ncycle = 1, 200  ! this loop would be possible instead of gotos
      total_number = 0.0
      total_dfermi = 0.0
      do i = 1, norbs
         fermifunct = 0.0
         if((eig(i)-e_fermi)/bkt.lt.50) then
            fermifunct = 1.0/(exp((eig(i)-e_fermi)/bkt)+1.0)
            dfermifunct = exp((eig(i)-e_fermi)/bkt) / &
            &       (bkt*(exp((eig(i)-e_fermi)/bkt)+1.0)**2)
         else
            dfermifunct = 0.0
         end if
         occ(i) = fermifunct
         total_number = total_number + fermifunct
         total_dfermi = total_dfermi + dfermifunct
      end do
      change_fermi = (occt-total_number)/total_dfermi
      e_fermi = e_fermi+change_fermi
      if (abs(occt-total_number).le.thr) exit
   enddo

   fod=0
   s  =0
   do i=1,norbs
      if(occ(i).gt.thr.and.1.0d00-occ(i).gt.thr) &
      &   s=s+occ(i)*log(occ(i))+(1.0d0-occ(i))*log(1.0d00-occ(i))
      if (eig(i).lt.e_fermi) then
         fod=fod+1.0d0-occ(i)
      else
         fod=fod+      occ(i)
      endif
   enddo
   s=s*kB*t

   if (prt) then
      write(*,'('' t,e(fermi),nfod : '',2f10.3,f10.6)') t,e_fermi,fod
   endif

end subroutine fermismear


subroutine occ(ndim,nel,nopen,ihomo,focc)
   integer  :: nel
   integer  :: nopen
   integer  :: ndim
   integer  :: ihomo
   real(wp) :: focc(ndim)
   integer  :: i,na,nb

   focc=0
!  even nel
   if(mod(nel,2).eq.0)then
      ihomo=nel/2
      do i=1,ihomo
         focc(i)=2.0d0
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1.0d0
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1.0
            focc(ihomo+i)=focc(ihomo+i)+1.0
         enddo
      endif
!  odd nel
   else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na
         focc(i)=focc(i)+1.
      enddo
      do i=1,nb
         focc(i)=focc(i)+1.
      enddo
   endif

   do i=1,ndim
      if(focc(i).gt.0.99) ihomo=i
   enddo

end subroutine occ


subroutine occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
   integer  :: nel
   integer  :: nopen
   integer  :: ndim
   integer  :: ihomoa
   integer  :: ihomob
   real(wp) :: focca(ndim)
   real(wp) :: foccb(ndim)
   integer  :: focc(ndim)
   integer  :: i,na,nb,ihomo

   focc=0
   focca=0
   foccb=0
!  even nel
   if(mod(nel,2).eq.0)then
      ihomo=nel/2
      do i=1,ihomo
         focc(i)=2
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1
            focc(ihomo+i)=focc(ihomo+i)+1
         enddo
      endif
!  odd nel
   else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na
         focc(i)=focc(i)+1
      enddo
      do i=1,nb
         focc(i)=focc(i)+1
      enddo
   endif

   do i=1,ndim
      if(focc(i).eq.2)then
         focca(i)=1.0d0
         foccb(i)=1.0d0
      endif
      if(focc(i).eq.1)focca(i)=1.0d0
   enddo

   ihomoa=0
   ihomob=0
   do i=1,ndim
      if(focca(i).gt.0.99) ihomoa=i
      if(foccb(i).gt.0.99) ihomob=i
   enddo

end subroutine occu


!> density matrix
! C: MO coefficient
! X: scratch
! P  dmat
subroutine dmat(ndim,focc,C,P,work,compactUpper)
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: focc(:)
   real(wp),intent(in)  :: C(:,:)
   real(wp),intent(out) :: P(:,:)
   real(wp),intent(inout),optional :: work(:,:)
   logical, intent(in), optional :: compactUpper
   integer :: i,j,m,nocc,ndoub,nsing,icol
   real(wp),allocatable :: Ptmp(:,:)
   logical :: ownWork, integerOccupation, useCompactUpper

   ! GFN1-fast 0.2.0: only exactly non-zero occupied columns contribute.
   nocc = 0
   do m = 1, ndim
      if (focc(m) /= 0.0_wp) nocc = m
   end do

   if (nocc == 0) then
      P = 0.0_wp
      return
   end if

   ownWork = .not.present(work)
   if (ownWork) then
      allocate(Ptmp(ndim,nocc))
   else
      if (size(work,1) < ndim .or. size(work,2) < nocc) then
         error stop 'dmat: workspace too small'
      end if
   end if

   useCompactUpper = .false.
   if (present(compactUpper)) useCompactUpper = compactUpper
   if (useCompactUpper) then
      ! Exact low-rank density factorization for arbitrary (including Fermi)
      ! occupations.  SYRK computes only the upper triangle required by the
      ! SCC contractions, halving the dense GEMM output/work of the legacy
      ! fractional-occupation path.
      do m = 1, nocc
         do i = 1, ndim
            if (ownWork) then
               Ptmp(i,m) = C(i,m)*sqrt(max(0.0_wp,focc(m)))
            else
               work(i,m) = C(i,m)*sqrt(max(0.0_wp,focc(m)))
            end if
         end do
      end do
      if (ownWork) then
         call blas_syrk('U','N',ndim,nocc,1.0_wp,Ptmp,ndim,0.0_wp,P,ndim)
      else
         call blas_syrk('U','N',ndim,nocc,1.0_wp,work,ndim,0.0_wp,P,ndim)
      end if
      if (ownWork) deallocate(Ptmp)
      return
   end if

   ! GFN1-fast 0.9.0: for the usual zero-temperature xTB occupations the
   ! density is the sum of invariant projectors for the doubly and singly
   ! occupied subspaces:
   !
   !   P = 2 C_d C_d^T + C_s C_s^T .
   !
   ! Build those projectors with SYRK, which evaluates only one triangle.
   ! This is mathematically identical to the GEMM path and independent of any
   ! rotation inside an equal-occupation occupied subspace.  Fractional/Fermi
   ! occupations retain the general GEMM path below without approximation.
   integerOccupation = .true.
   ndoub = 0
   nsing = 0
   do m = 1, nocc
      if (focc(m) == 2.0_wp) then
         ndoub = ndoub + 1
      else if (focc(m) == 1.0_wp) then
         nsing = nsing + 1
      else if (focc(m) /= 0.0_wp) then
         integerOccupation = .false.
         exit
      end if
   end do

   if (integerOccupation) then
      P = 0.0_wp

      if (ndoub > 0) then
         icol = 0
         do m = 1, nocc
            if (focc(m) == 2.0_wp) then
               icol = icol + 1
               if (ownWork) then
                  Ptmp(:,icol) = C(:,m)
               else
                  work(:,icol) = C(:,m)
               end if
            end if
         end do
         if (ownWork) then
            call blas_syrk('U','N',ndim,ndoub,2.0_wp,Ptmp,ndim,0.0_wp,P,ndim)
         else
            call blas_syrk('U','N',ndim,ndoub,2.0_wp,work,ndim,0.0_wp,P,ndim)
         end if
      end if

      if (nsing > 0) then
         icol = 0
         do m = 1, nocc
            if (focc(m) == 1.0_wp) then
               icol = icol + 1
               if (ownWork) then
                  Ptmp(:,icol) = C(:,m)
               else
                  work(:,icol) = C(:,m)
               end if
            end if
         end do
         if (ownWork) then
            call blas_syrk('U','N',ndim,nsing,1.0_wp,Ptmp,ndim, &
               & merge(1.0_wp,0.0_wp,ndoub>0),P,ndim)
         else
            call blas_syrk('U','N',ndim,nsing,1.0_wp,work,ndim, &
               & merge(1.0_wp,0.0_wp,ndoub>0),P,ndim)
         end if
      end if

      ! SYRK writes the upper triangle used by the original SCC population
      ! routines. Mirror it once so downstream gradient code sees a full
      ! symmetric density matrix exactly as before.
      do j = 1, ndim
         do i = j+1, ndim
            P(i,j) = P(j,i)
         end do
      end do

      if (ownWork) deallocate(Ptmp)
      return
   end if

   ! General exact occupation path (including finite-temperature smearing).
   do m=1,nocc
      do i=1,ndim
         if (ownWork) then
            Ptmp(i,m)=C(i,m)*focc(m)
         else
            work(i,m)=C(i,m)*focc(m)
         end if
      enddo
   enddo

   if (ownWork) then
      call blas_gemm('N','T',ndim,ndim,nocc,1.0_wp,C,ndim,Ptmp,ndim, &
         & 0.0_wp,P,ndim)
   else
      call blas_gemm('N','T',ndim,ndim,nocc,1.0_wp,C,ndim,work,ndim, &
         & 0.0_wp,P,ndim)
   end if

   if (ownWork) deallocate(Ptmp)

end subroutine dmat

! Reference: I. Mayer, "Simple Theorems, Proofs, and Derivations in Quantum Chemistry", formula (7.35)
subroutine get_wiberg(n,ndim,at,xyz,P,S,wb,fila2)
   integer, intent(in)  :: n,ndim,at(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: P(ndim,ndim)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(out) :: wb (n,n)
   integer, intent(in)  :: fila2(:,:)

   real(wp),allocatable :: Ptmp(:,:)
   real(wp) xsum,rab
   integer i,j,k,m

   allocate(Ptmp(ndim,ndim))
   call blas_gemm('N','N',ndim,ndim,ndim,1.0d0,P,ndim,S,ndim,0.0d0,Ptmp,ndim)
   wb = 0
   do i = 1, n
      do j = 1, i-1
         xsum = 0.0_wp
         rab = sum((xyz(:,i) - xyz(:,j))**2)
         if(rab < 100.0_wp)then
            do k = fila2(1,i), fila2(2,i) ! AOs on atom i
               do m = fila2(1,j), fila2(2,j) ! AOs on atom j
                  xsum = xsum + Ptmp(k,m)*Ptmp(m,k)
               enddo
            enddo
         endif
         wb(i,j) = xsum
         wb(j,i) = xsum
      enddo
   enddo
   deallocate(Ptmp)

end subroutine get_wiberg

! Reference: I. Mayer, "Simple Theorems, Proofs, and Derivations in Quantum Chemistry", formula (7.36)
subroutine get_unrestricted_wiberg(n,ndim,at,xyz,Pa,Pb,S,wb,fila2)
   integer, intent(in)  :: n,ndim,at(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: Pa(ndim,ndim)
   real(wp),intent(in)  :: Pb(ndim,ndim)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(out) :: wb (n,n)
   integer, intent(in)  :: fila2(:,:)

   real(wp),allocatable :: Ptmp_a(:,:)
   real(wp),allocatable :: Ptmp_b(:,:)
   real(wp) xsum,rab
   integer i,j,k,m

   allocate(Ptmp_a(ndim,ndim))
   allocate(Ptmp_b(ndim,ndim))

   ! P^(alpha) * S !
   call blas_gemm('N','N',ndim,ndim,ndim,1.0d0,Pa,ndim,S,ndim,0.0d0,Ptmp_a,ndim)
   
   ! P^(beta) * S !
   call blas_gemm('N','N',ndim,ndim,ndim,1.0d0,Pb,ndim,S,ndim,0.0d0,Ptmp_b,ndim)
   
   wb = 0
   do i = 1, n
      do j = 1, i-1
         xsum = 0.0_wp
         rab = sum((xyz(:,i) - xyz(:,j))**2)
         if(rab < 100.0_wp)then
            do k = fila2(1,i), fila2(2,i) ! AOs on atom i
               do m = fila2(1,j), fila2(2,j) ! AOs on atom j
                  xsum = xsum + Ptmp_a(k,m)*Ptmp_a(m,k) + Ptmp_b(k,m)*Ptmp_b(m,k) 
               enddo
            enddo
         endif
         wb(i,j) = 2*xsum
         wb(j,i) = 2*xsum
      enddo
   enddo
   deallocate(Ptmp_a)
   deallocate(Ptmp_b)

end subroutine get_unrestricted_wiberg

!> Mulliken pop + AO pop
subroutine mpopall(n,nao,aoat,S,P,qao,q)
   integer nao,n,aoat(nao)
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  qao(nao),q(n),ps

   integer i,j,ii,jj,ij,is,js

   q  = 0
   qao= 0
   do i=1,nao
      ii=aoat(i)
      do j=1,i-1
         jj=aoat(j)
         ps=p(j,i)*s(j,i)
         q(ii)=q(ii)+ps
         q(jj)=q(jj)+ps
         qao(i)=qao(i)+ps
         qao(j)=qao(j)+ps
      enddo
      ps=p(i,i)*s(i,i)
      q(ii)=q(ii)+ps
      qao(i)=qao(i)+ps
   enddo

end subroutine mpopall


!> Mulliken pop
subroutine mpop0(n,nao,aoat,S,P,q)
   integer nao,n,aoat(nao)
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  q(n),ps

   integer i,j,ii,jj,ij,is,js

   q = 0
   do i=1,nao
      ii=aoat(i)
      do j=1,i-1
         jj=aoat(j)
         ps=p(j,i)*s(j,i)
         q(ii)=q(ii)+ps
         q(jj)=q(jj)+ps
      enddo
      ps=p(i,i)*s(i,i)
      q(ii)=q(ii)+ps
   enddo

end subroutine mpop0


!> Mulliken AO pop
subroutine mpopao(n,nao,S,P,qao)
   integer nao,n
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  qao(nao),ps

   integer i,j

   qao = 0
   do i=1,nao
      do j=1,i-1
         ps=p(j,i)*s(j,i)
         qao(i)=qao(i)+ps
         qao(j)=qao(j)+ps
      enddo
      ps=p(i,i)*s(i,i)
      qao(i)=qao(i)+ps
   enddo

end subroutine mpopao


!> Mulliken pop
subroutine mpop(n,nao,aoat,lao,S,P,q,ql)
   integer nao,n,aoat(nao),lao(nao)
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  q(n),ps
   real(wp)  ql(3,n)

   integer i,j,ii,jj,ij,is,js,mmm(20)
   data    mmm/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/

   ql= 0
   q = 0
   do i=1,nao
      ii=aoat(i)
      is=mmm(lao(i))
      do j=1,i-1
         jj=aoat(j)
         js=mmm(lao(j))
         ps=p(j,i)*s(j,i)
         q(ii)=q(ii)+ps
         q(jj)=q(jj)+ps
         ql(is,ii)=ql(is,ii)+ps
         ql(js,jj)=ql(js,jj)+ps
      enddo
      ps=p(i,i)*s(i,i)
      q(ii)=q(ii)+ps
      ql(is,ii)=ql(is,ii)+ps
   enddo

end subroutine mpop


!> Fused shell-Mulliken and fixed-H0 contraction using the upper P triangle.
!> GFN1-fast 2.2.6: mpopsh() and electro() traversed the same AO triangle
!> separately.  This routine preserves each accumulator's original order while
!> loading P/S/H0 only once.
subroutine mpopsh_h0_upper(nao,nshell,ao2sh,S,P,H0,popsh,h0EnergyEh)
   use xtb_mctc_convert, only : evtoau
   integer, intent(in) :: nao, nshell, ao2sh(nao)
   real(wp), intent(in) :: S(nao,nao), P(nao,nao)
   real(wp), intent(in) :: H0(nao*(nao+1)/2)
   real(wp), intent(out) :: popsh(nshell), h0EnergyEh
   integer :: i,j,ii,jj,k
   real(wp) :: ps,h

   popsh = 0.0_wp
   h = 0.0_wp
   k = 0
   do i = 1, nao
      ii = ao2sh(i)
      do j = 1, i-1
         k = k + 1
         jj = ao2sh(j)
         ps = P(j,i)*S(j,i)
         popsh(ii) = popsh(ii) + ps
         popsh(jj) = popsh(jj) + ps
         h = h + P(j,i)*H0(k)
      end do
      k = k + 1
      ps = P(i,i)*S(i,i)
      popsh(ii) = popsh(ii) + ps
      h = h + 0.5_wp*P(i,i)*H0(k)
   end do
   h0EnergyEh = 2.0_wp*h*evtoau
end subroutine mpopsh_h0_upper


!> Mulliken pop shell wise
subroutine mpopsh(n,nao,nshell,ao2sh,S,P,qsh)
   integer nao,n,nshell,ao2sh(nao)
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  qsh(nshell),ps

   integer i,j,ii,jj,ij

   qsh=0
   do i=1,nao
      ii =ao2sh(i)
      do j=1,i-1
         jj =ao2sh(j)
         ps=p(j,i)*s(j,i)
         qsh(ii)=qsh(ii)+ps
         qsh(jj)=qsh(jj)+ps
      enddo
      ps=p(i,i)*s(i,i)
      qsh(ii)=qsh(ii)+ps
   enddo

end subroutine mpopsh


subroutine qsh2qat(ash,qsh,qat)
   integer, intent(in) :: ash(:)
   real(wp), intent(in) :: qsh(:)
   real(wp), intent(out) :: qat(:)

   integer :: iSh

   qat(:) = 0.0_wp
   do iSh = 1, size(qsh)
      qat(ash(iSh)) = qat(ash(iSh)) + qsh(iSh)
   enddo

end subroutine qsh2qat


!> Loewdin pop
subroutine lpop(n,nao,aoat,lao,occ,C,f,q,ql)
   integer nao,n,aoat(nao),lao(nao)
   real(wp)  C (nao,nao)
   real(wp)  occ(nao)
   real(wp)  q(n)
   real(wp)  ql(3,n)
   real(wp)  f

   integer i,j,ii,jj,js,mmm(20)
   data    mmm/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/
   real(wp)  cc

   do i=1,nao
      if(occ(i).lt.1.d-8) cycle
      do j=1,nao
         cc=f*C(j,i)*C(j,i)*occ(i)
         jj=aoat(j)
         js=mmm(lao(j))
         q(jj)=q(jj)+cc
         ql(js,jj)=ql(js,jj)+cc
      enddo
   enddo

end subroutine lpop


!> atomic valence shell pops and total atomic energy
subroutine iniqshell(xtbData,n,at,z,nshell,q,qsh,gfn_method)
   type(TxTBData), intent(in) :: xtbData
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: nshell
   integer, intent(in)  :: gfn_method
   real(wp),intent(in)  :: z(n)
   real(wp),intent(in)  :: q(n)
   real(wp),intent(out) :: qsh(nshell)
   real(wp) :: zshell
   real(wp) :: ntot,fracz
   integer  :: i,j,k,m,l,ll(0:3),iat,lll,iver
   data ll /1,3,5,7/

   qsh = 0.0_wp

   k=0
   do i=1,n
      iat=at(i)
      ntot=-1.d-6
      do m=1,xtbData%nShell(iat)
         l=xtbData%hamiltonian%angShell(m,iat)
         k=k+1
         zshell=xtbData%hamiltonian%referenceOcc(m,iat)
         ntot=ntot+zshell
         if(ntot.gt.z(i)) zshell=0
         fracz=zshell/z(i)
         qsh(k)=fracz*q(i)
      enddo
   enddo

end subroutine iniqshell


subroutine setzshell(xtbData,n,at,nshell,z,zsh,e,gfn_method)
   type(TxTBData), intent(in) :: xtbData
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: nshell
   integer, intent(in)  :: gfn_method
   real(wp),intent(in)  :: z(n)
   real(wp),intent(out) :: zsh(nshell)
!   integer, intent(out) :: ash(nshell)
!   integer, intent(out) :: lsh(nshell)
   real(wp),intent(out) :: e

   real(wp)  ntot,fracz
   integer i,j,k,m,l,ll(0:3),iat,lll,iver
   data ll /1,3,5,7/

   k=0
   e=0.0_wp
   do i=1,n
      iat=at(i)
      ntot=-1.d-6
      do m=1,xtbData%nShell(iat)
         l=xtbData%hamiltonian%angShell(m,iat)
         k=k+1
         zsh(k)=xtbData%hamiltonian%referenceOcc(m,iat)
!         lsh(k)=l
!         ash(k)=i
         ntot=ntot+zsh(k)
         if(ntot.gt.z(i)) zsh(k)=0
         e=e+xtbData%hamiltonian%selfEnergy(m,iat)*zsh(k)
      enddo
   enddo

end subroutine setzshell


end module xtb_scc_core
