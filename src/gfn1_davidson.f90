! This file is part of the GFN1-fast development branch of xtb.
!
! GFN1-fast 2.1.4: generalized block-Davidson eigensolver for the lowest
! occupied/guard subspace.  The routine works directly in the non-orthogonal
! AO basis and therefore avoids the O(N^3) dense tridiagonalization used by
! SYEVD/SYEVX when a high-quality SCC-to-SCC starting subspace is available.
module xtb_gfn1_davidson
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_gemm
   use xtb_mctc_lapack, only : lapack_syevd
   implicit none
   private

   public :: gfn1BlockDavidson, gfn1BlockDavidsonMF, applyGFN1HamiltonianBlock
   public :: buildGFN1PairHamiltonian, applyGFN1PairHamiltonianBlock
   public :: buildGFN1SparsePairs, updateGFN1PairHamiltonian
   public :: applyGFN1PairMetricBlock, gfn1DavidsonGuardCount
   public :: validateGFN1DavidsonEigensystem

   integer, parameter :: maxDavidsonIter = 20
   integer, parameter :: maxCorrectionBlock = 64
   real(wp), parameter :: davidsonResidualTol = 5.0e-12_wp
   real(wp), parameter :: minMetricNorm = 1.0e-12_wp
   real(wp), parameter :: preconditionerFloor = 5.0e-2_wp

contains


!> Extract the screened overlap and core-Hamiltonian values into the same
!> symmetric AO-pair list used by the SCC Hamiltonian.  This is exact with
!> respect to the existing xtb screening: entries absent from matlist have
!> already been set to zero in the dense overlap matrix.
subroutine buildGFN1SparsePairs(H0,S,nmat,matlist,h0idx,pairS,pairH0)
   real(wp), intent(in) :: H0(:),S(:,:)
   integer, intent(in) :: nmat,matlist(2,nmat),h0idx(nmat)
   real(wp), intent(out) :: pairS(nmat),pairH0(nmat)
   integer :: m,i,j
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(H0,S,nmat,matlist,h0idx,pairS,pairH0) private(m,i,j) if(nmat>=8192)
   do m=1,nmat
      i=matlist(1,m); j=matlist(2,m)
      pairS(m)=S(j,i)
      pairH0(m)=H0(h0idx(m))
   end do
   !$omp end parallel do
end subroutine buildGFN1SparsePairs

!> Update only the charge-dependent pair Hamiltonian values.  H0 and S are
!> already compressed, eliminating indirect dense-matrix reads in every SCC.
subroutine updateGFN1PairHamiltonian(pairH0,pairS,nmat,matISh,matJSh, &
      & shellShift,pairH)
   use xtb_mctc_convert, only : autoev
   real(wp), intent(in) :: pairH0(nmat),pairS(nmat),shellShift(:)
   integer, intent(in) :: nmat,matISh(nmat),matJSh(nmat)
   real(wp), intent(out) :: pairH(nmat)
   integer :: m
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(pairH0,pairS,nmat,matISh,matJSh,shellShift,pairH) private(m) if(nmat>=8192)
   do m=1,nmat
      pairH(m)=pairH0(m)-0.5_wp*autoev*pairS(m) * &
         & (shellShift(matISh(m))+shellShift(matJSh(m)))
   end do
   !$omp end parallel do
end subroutine updateGFN1PairHamiltonian

!> Apply a compressed symmetric AO-pair matrix to a dense block.
subroutine applyGFN1PairBlock(pairA,nmat,matlist,X,Y)
   real(wp), intent(in) :: pairA(nmat),X(:,:)
   integer, intent(in) :: nmat,matlist(2,nmat)
   real(wp), intent(out) :: Y(:,:)
   integer :: m,i,j,k
   Y=0.0_wp
   !$omp parallel do default(none) schedule(static) shared(pairA,nmat,matlist,X,Y) &
   !$omp private(k,m,i,j) if(size(X,2)>=4)
   do k=1,size(X,2)
      do m=1,nmat
         i=matlist(1,m); j=matlist(2,m)
         Y(i,k)=Y(i,k)+pairA(m)*X(j,k)
         if(i/=j)Y(j,k)=Y(j,k)+pairA(m)*X(i,k)
      end do
   end do
   !$omp end parallel do
end subroutine applyGFN1PairBlock

subroutine applyGFN1PairMetricBlock(pairS,nmat,matlist,X,Y)
   real(wp), intent(in) :: pairS(nmat),X(:,:)
   integer, intent(in) :: nmat,matlist(2,nmat)
   real(wp), intent(out) :: Y(:,:)
   call applyGFN1PairBlock(pairS,nmat,matlist,X,Y)
end subroutine applyGFN1PairMetricBlock

!> Build the compressed GFN1 Hamiltonian values for the existing AO-pair list.
!> This is O(Npair) once per SCC iteration and avoids recomputing shell shifts
!> inside every Davidson block product.
subroutine buildGFN1PairHamiltonian(H0,S,nmat,matlist,h0idx,matISh,matJSh, &
      & shellShift,pairH)
   use xtb_mctc_convert, only : autoev
   real(wp), intent(in) :: H0(:),S(:,:),shellShift(:)
   integer, intent(in) :: nmat,matlist(2,nmat),h0idx(nmat),matISh(nmat),matJSh(nmat)
   real(wp), intent(out) :: pairH(nmat)
   integer :: m,i,j
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(H0,S,nmat,matlist,h0idx,matISh,matJSh,shellShift,pairH) private(m,i,j) if(nmat>=8192)
   do m=1,nmat
      i=matlist(1,m); j=matlist(2,m)
      pairH(m)=H0(h0idx(m))-0.5_wp*autoev*S(j,i) * &
         & (shellShift(matISh(m))+shellShift(matJSh(m)))
   end do
   !$omp end parallel do
end subroutine buildGFN1PairHamiltonian

!> Apply a compressed symmetric AO-pair Hamiltonian to a dense block.
subroutine applyGFN1PairHamiltonianBlock(pairH,nmat,matlist,X,Y)
   real(wp), intent(in) :: pairH(nmat),X(:,:)
   integer, intent(in) :: nmat,matlist(2,nmat)
   real(wp), intent(out) :: Y(:,:)
   call applyGFN1PairBlock(pairH,nmat,matlist,X,Y)
end subroutine applyGFN1PairHamiltonianBlock

!> Apply the GFN1 isotropic Hamiltonian to a block without materializing H.
!> The packed H0 and AO-pair list are exactly the same data used by
!> buildIsotropicH1Cached, so this changes only the execution order, not the
!> model or screening threshold.  Parallelism is over block columns; each
!> column has private accumulation and is therefore race-free.
subroutine applyGFN1HamiltonianBlock(H0,S,nmat,matlist,h0idx,matISh,matJSh, &
      & shellShift,X,Y)
   real(wp), intent(in) :: H0(:),S(:,:),shellShift(:),X(:,:)
   integer, intent(in) :: nmat,matlist(2,nmat),h0idx(nmat),matISh(nmat),matJSh(nmat)
   real(wp), intent(out) :: Y(:,:)
   real(wp), allocatable :: pairH(:)
   allocate(pairH(nmat))
   call buildGFN1PairHamiltonian(H0,S,nmat,matlist,h0idx,matISh,matJSh,shellShift,pairH)
   call applyGFN1PairHamiltonianBlock(pairH,nmat,matlist,X,Y)
   deallocate(pairH)
end subroutine applyGFN1HamiltonianBlock

!> Number of extra low-energy Ritz directions carried across SCC steps.
!> The guard space is bounded independently of system size so the seed remains
!> compact while retaining the most useful correction history.
pure integer function gfn1DavidsonGuardCount(nroots,n) result(nguard)
   integer, intent(in) :: nroots,n
   if (nroots <= 0 .or. n <= nroots) then
      nguard=0
   else
      nguard=min(n-nroots,min(32,max(8,nroots/8)))
   end if
end function gfn1DavidsonGuardCount

!> Strict matrix-free validation of a candidate generalized eigensystem.
!>
!> 2.1.4 strengthens the production acceptance gate beyond per-root residuals:
!>   * all requested roots must be finite and monotonically ordered,
!>   * every generalized residual H*c-S*c*e is small,
!>   * the complete S-metric Gram matrix C^T*S*C must be the identity,
!>   * the projected Hamiltonian C^T*H*C must be diagonal and agree with e.
!>
!> Checking the full Gram/projected matrices catches non-adjacent loss of
!> orthogonality and duplicate/misrotated Ritz vectors that an adjacent-only
!> test can miss.  No model tolerance or SCC convergence threshold is changed.
subroutine validateGFN1DavidsonEigensystem(pairH,pairS,nmat,matlist,C,eval, &
      & workHC,workSC,maxRelResidual,maxNormError,maxOrthError, &
      & maxProjectedError,success)
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   real(wp), intent(in) :: pairH(nmat),pairS(nmat),C(:,:),eval(:)
   integer, intent(in) :: nmat,matlist(2,nmat)
   real(wp), intent(inout) :: workHC(:,:),workSC(:,:)
   real(wp), intent(out) :: maxRelResidual,maxNormError,maxOrthError,maxProjectedError
   logical, intent(out) :: success
   integer :: n,nroots,i,j
   real(wp) :: r2,hc2,sc2,denom,relres,escale,hscale
   real(wp), allocatable :: gram(:,:)
   real(wp), parameter :: residualTol=1.0e-11_wp
   real(wp), parameter :: normTol=1.0e-10_wp
   real(wp), parameter :: orthTol=1.0e-9_wp
   real(wp), parameter :: projectedTol=2.0e-10_wp

   success=.false.
   maxRelResidual=huge(1.0_wp)
   maxNormError=huge(1.0_wp)
   maxOrthError=huge(1.0_wp)
   maxProjectedError=huge(1.0_wp)
   n=size(C,1); nroots=size(eval)
   if(nroots<1 .or. size(C,2)<nroots) return
   if(size(workHC,1)<n .or. size(workHC,2)<nroots) return
   if(size(workSC,1)<n .or. size(workSC,2)<nroots) return

   call applyGFN1PairHamiltonianBlock(pairH,nmat,matlist,C(:,1:nroots),workHC(:,1:nroots))
   call applyGFN1PairMetricBlock(pairS,nmat,matlist,C(:,1:nroots),workSC(:,1:nroots))

   maxRelResidual=0.0_wp
   do j=1,nroots
      if(.not.ieee_is_finite(eval(j))) return
      if(j>1)then
         escale=max(1.0_wp,abs(eval(j)),abs(eval(j-1)))
         if(eval(j)<eval(j-1)-64.0_wp*epsilon(1.0_wp)*escale) return
      endif
      r2=sum((workHC(:,j)-eval(j)*workSC(:,j))**2)
      hc2=sum(workHC(:,j)**2); sc2=sum(workSC(:,j)**2)
      if(.not.ieee_is_finite(r2).or..not.ieee_is_finite(hc2).or..not.ieee_is_finite(sc2)) return
      denom=sqrt(max(0.0_wp,hc2))+abs(eval(j))*sqrt(max(0.0_wp,sc2))+tiny(1.0_wp)
      relres=sqrt(max(0.0_wp,r2))/denom
      if(.not.ieee_is_finite(relres)) return
      maxRelResidual=max(maxRelResidual,relres)
   enddo

   allocate(gram(nroots,nroots))
   call mctc_gemm(C(:,1:nroots),workSC(:,1:nroots),gram,transa='t')

   maxNormError=0.0_wp; maxOrthError=0.0_wp; maxProjectedError=0.0_wp
   do j=1,nroots
      if(.not.ieee_is_finite(gram(j,j)))then
         deallocate(gram); return
      endif
      maxNormError=max(maxNormError,abs(gram(j,j)-1.0_wp))
      hscale=max(1.0_wp,abs(eval(j)))
      maxProjectedError=max(maxProjectedError, &
         & abs(dot_product(C(:,j),workHC(:,j))-eval(j))/hscale)
      do i=1,nroots
         if(.not.ieee_is_finite(gram(i,j)))then
            deallocate(gram); return
         endif
         if(i/=j) maxOrthError=max(maxOrthError,abs(gram(i,j)))
      enddo
   enddo

   success=maxRelResidual<=residualTol .and. maxNormError<=normTol .and. &
      & maxOrthError<=orthTol .and. maxProjectedError<=projectedTol
   deallocate(gram)
end subroutine validateGFN1DavidsonEigensystem

!> Sparse-pair generalized block Davidson for GFN1.  Both H*V and S*V use
!> the screened AO-pair representation, so the Davidson loop no longer reads
!> dense H, dense S, or packed H0.
subroutine gfn1BlockDavidsonMF(pairH,pairS,S,nmat,matlist,seed,nroots,eval,C,success,niter,maxResidual, &
      & carry,ncarry,nrestart)
   real(wp), intent(in) :: pairH(nmat),pairS(nmat),S(:,:),seed(:,:)
   integer, intent(in) :: nmat,matlist(2,nmat)
   integer, intent(in) :: nroots
   real(wp), intent(out) :: eval(:), C(:,:)
   logical, intent(out) :: success
   integer, intent(out) :: niter
   real(wp), intent(out) :: maxResidual
   real(wp), intent(out), optional :: carry(:,:)
   integer, intent(out), optional :: ncarry,nrestart

   integer :: n, nseed, maxSub, restartTrigger, msub, iter, i, j, nadd, nkeep
   integer :: iconv, idx, best, info, m, ia, ja, restartGuard, nretain
   integer :: carryCount, restartCount
   integer, allocatable :: order(:)
   logical, allocatable :: chosen(:)
   real(wp) :: denom, floorv, rel, hc2, sc2, r2
   real(wp), allocatable :: V(:,:), AV(:,:), BV(:,:)
   real(wp), allocatable :: proj(:,:), y(:,:), theta(:)
   real(wp), allocatable :: ritz(:,:), hritz(:,:), britz(:,:), resid(:,:)
   real(wp), allocatable :: corr(:,:), restartV(:,:)
   real(wp), allocatable :: relres(:), diagH(:), diagS(:)

   success = .false.
   niter = 0
   maxResidual = huge(1.0_wp)
   restartCount = 0
   carryCount = 0
   if (present(ncarry)) ncarry = 0
   if (present(nrestart)) nrestart = 0
   n = size(S,1)
   if (size(S,2) /= n) return
   if (nroots < 1 .or. nroots >= n) return
   if (size(C,1) < n .or. size(C,2) < nroots) return
   if (size(eval) < nroots) return
   nseed = min(size(seed,2),n)
   if (size(seed,1) < n .or. nseed < nroots) return

   ! 2.1.3: keep a small low-energy guard space across SCC iterations.  It is
   ! not part of the accepted eigensystem; it only improves the next Davidson
   ! starting subspace.  The count is deliberately capped so memory remains
   ! O(N*Nocc), not O(N^2).
   restartGuard = gfn1DavidsonGuardCount(nroots,n)
   maxSub = min(n,max(nroots + 2*restartGuard + maxCorrectionBlock, &
      & min(3*nroots,nroots+256)))
   if (maxSub <= nroots) return
   restartTrigger=min(maxSub,nroots+restartGuard+max(8,maxCorrectionBlock/2))
   nseed = min(nseed,maxSub)

   allocate(V(n,maxSub),AV(n,maxSub),BV(n,maxSub))
   allocate(ritz(n,nroots),hritz(n,nroots),britz(n,nroots),resid(n,nroots))
   allocate(relres(nroots),diagH(n),diagS(n),order(nroots),chosen(nroots))
   allocate(corr(n,maxCorrectionBlock))

   diagH = 0.0_wp
   diagS = 0.0_wp
   do m=1,nmat
      ia=matlist(1,m); ja=matlist(2,m)
      if(ia==ja)then
         diagH(ia)=pairH(m)
         diagS(ia)=pairS(m)
      endif
   end do

   ! Start with all available previous-SCC low-energy vectors, including the
   ! carried guard Ritz space.  S-orthogonalization uses the same screened pair
   ! metric as every Davidson S*V, eliminating the dense-S O(N^2*Nseed) path.
   msub = nseed
   V(:,1:msub) = seed(:,1:msub)
   call metricOrthonormalizeBlockPair(pairS,nmat,matlist,V(:,1:msub),msub,nkeep)
   if (nkeep < nroots) then
      deallocate(V,AV,BV,ritz,hritz,britz,resid,relres,diagH,diagS,order,chosen,corr)
      return
   end if
   msub = nkeep
   call applyGFN1PairHamiltonianBlock(pairH,nmat,matlist,V(:,1:msub),AV(:,1:msub))
   call applyGFN1PairMetricBlock(pairS,nmat,matlist,V(:,1:msub),BV(:,1:msub))

   do iter = 1, maxDavidsonIter
      niter = iter
      allocate(proj(msub,msub),y(msub,msub),theta(msub))
      call mctc_gemm(V(:,1:msub),AV(:,1:msub),proj,transa='t')
      proj = 0.5_wp*(proj+transpose(proj))
      y = proj
      call diagonalizeSmall(y,theta,info)
      deallocate(proj)
      if (info /= 0) then
         deallocate(y,theta); exit
      end if
      call mctc_gemm(V(:,1:msub),y(:,1:nroots),ritz)
      call mctc_gemm(AV(:,1:msub),y(:,1:nroots),hritz)
      call mctc_gemm(BV(:,1:msub),y(:,1:nroots),britz)

      maxResidual = 0.0_wp; iconv = 0
      do j = 1, nroots
         resid(:,j) = hritz(:,j)-theta(j)*britz(:,j)
         r2=sum(resid(:,j)*resid(:,j)); hc2=sum(hritz(:,j)*hritz(:,j)); sc2=sum(britz(:,j)*britz(:,j))
         denom=sqrt(max(0.0_wp,hc2))+abs(theta(j))*sqrt(max(0.0_wp,sc2))+tiny(1.0_wp)
         relres(j)=sqrt(max(0.0_wp,r2))/denom
         maxResidual=max(maxResidual,relres(j))
         if (relres(j) <= davidsonResidualTol) iconv=iconv+1
      end do
      if (iconv == nroots) then
         eval(1:nroots)=theta(1:nroots)
         C(:,1:nroots)=ritz(:,1:nroots)
         success=.true.

         ! Carry the lowest extra Ritz directions to the next SCC iteration.
         ! These are not reported as converged roots and therefore cannot alter
         ! occupations/energy; they are only a better starting invariant space.
         if (present(carry)) then
            if (size(carry,1) >= n) then
               carryCount=min(size(carry,2),min(msub,nroots+restartGuard))
               if (carryCount > 0) then
                  call mctc_gemm(V(:,1:msub),y(:,1:carryCount),carry(:,1:carryCount))
               end if
            end if
         end if
         if (present(ncarry)) ncarry=carryCount
         if (present(nrestart)) nrestart=restartCount
         deallocate(y,theta)
         exit
      end if

      chosen=.false.; nadd=0
      do i=1,min(maxCorrectionBlock,nroots-iconv)
         best=0; rel=-1.0_wp
         do j=1,nroots
            if (.not.chosen(j) .and. relres(j)>davidsonResidualTol .and. relres(j)>rel) then
               best=j; rel=relres(j)
            end if
         end do
         if (best==0) exit
         chosen(best)=.true.; nadd=nadd+1; order(nadd)=best
      end do
      if (nadd==0) then
         deallocate(y,theta); exit
      end if
      do j=1,nadd
         idx=order(j)
         do i=1,n
            denom=theta(idx)*diagS(i)-diagH(i)
            floorv=max(preconditionerFloor,1.0e-6_wp*max(1.0_wp,abs(theta(idx))))
            if (abs(denom)<floorv) denom=merge(-floorv,floorv,denom<0.0_wp)
            corr(i,j)=resid(i,idx)/denom
         end do
      end do

      ! True thick restart: preserve the requested roots plus a bounded guard
      ! band of the lowest Ritz vectors instead of discarding every correction
      ! direction and restarting from only nroots.
      if (msub+nadd>maxSub .or. (msub>=restartTrigger .and. msub>nroots+restartGuard)) then
         nretain=min(msub,nroots+restartGuard)
         allocate(restartV(n,nretain))
         call mctc_gemm(V(:,1:msub),y(:,1:nretain),restartV)
         V(:,1:nretain)=restartV
         deallocate(restartV)
         call metricOrthonormalizeBlockPair(pairS,nmat,matlist,V(:,1:nretain),nretain,nkeep)
         if (nkeep < nroots) then
            deallocate(y,theta); exit
         end if
         msub=nkeep
         call applyGFN1PairHamiltonianBlock(pairH,nmat,matlist,V(:,1:msub),AV(:,1:msub))
         call applyGFN1PairMetricBlock(pairS,nmat,matlist,V(:,1:msub),BV(:,1:msub))
         restartCount=restartCount+1
      end if

      call orthogonalizeCorrectionsPair(pairS,nmat,matlist,V(:,1:msub),corr(:,1:nadd),nadd,nkeep)
      if (nkeep<=0) then
         do j=1,nadd; corr(:,j)=resid(:,order(j)); end do
         call orthogonalizeCorrectionsPair(pairS,nmat,matlist,V(:,1:msub),corr(:,1:nadd),nadd,nkeep)
      end if
      if (nkeep<=0) then
         deallocate(y,theta); exit
      end if
      nadd=min(nkeep,maxSub-msub)
      if (nadd<=0) then
         deallocate(y,theta); exit
      end if
      V(:,msub+1:msub+nadd)=corr(:,1:nadd)
      call applyGFN1PairHamiltonianBlock(pairH,nmat,matlist, &
         & V(:,msub+1:msub+nadd),AV(:,msub+1:msub+nadd))
      call applyGFN1PairMetricBlock(pairS,nmat,matlist,V(:,msub+1:msub+nadd),BV(:,msub+1:msub+nadd))
      msub=msub+nadd
      deallocate(y,theta)
   end do

   if (present(nrestart)) nrestart=restartCount
   deallocate(V,AV,BV,ritz,hritz,britz,resid,relres,diagH,diagS,order,chosen,corr)
contains

subroutine diagonalizeSmall(a,w,info)
   real(wp), intent(inout) :: a(:,:)
   real(wp), intent(out) :: w(:)
   integer, intent(out) :: info
   integer :: nn, lwork, liwork
   real(wp), allocatable :: lworkArray(:)
   integer, allocatable :: liworkArray(:)

   nn = size(a,1)
   lwork = max(1,1 + 6*nn + 2*nn*nn)
   liwork = max(1,3 + 5*nn)
   allocate(lworkArray(lwork),liworkArray(liwork))
   call lapack_syevd('v','u',nn,a,nn,w,lworkArray,lwork,liworkArray,liwork,info)
   deallocate(lworkArray,liworkArray)
end subroutine diagonalizeSmall

!> Stable block Lowdin orthogonalization in the screened S metric.  A second
!> pass is performed to suppress loss of S-orthogonality when the incoming
!> previous-SCC vectors contain nearly dependent guard directions.
subroutine metricOrthonormalizeBlockPair(pairMetric,npair,pairs,block,ncol,nkeep)
   real(wp), intent(in) :: pairMetric(npair)
   integer, intent(in) :: npair,pairs(2,npair),ncol
   real(wp), intent(inout) :: block(:,:)
   integer, intent(out) :: nkeep
   integer :: j, info, pass, ncur
   real(wp) :: scale, cutoff
   real(wp), allocatable :: sblock(:,:), gram(:,:), eig(:), tmp(:,:)

   allocate(sblock(size(block,1),ncol),gram(ncol,ncol),eig(ncol),tmp(size(block,1),ncol))
   ncur=ncol
   nkeep=0
   do pass=1,2
      call applyGFN1PairMetricBlock(pairMetric,npair,pairs,block(:,1:ncur),sblock(:,1:ncur))
      call mctc_gemm(block(:,1:ncur),sblock(:,1:ncur),gram(1:ncur,1:ncur),transa='t')
      gram(1:ncur,1:ncur)=0.5_wp*(gram(1:ncur,1:ncur)+transpose(gram(1:ncur,1:ncur)))
      call diagonalizeSmall(gram(1:ncur,1:ncur),eig(1:ncur),info)
      if (info /= 0) then
         nkeep=0; exit
      end if
      scale=maxval(abs(eig(1:ncur)))
      cutoff=max(tiny(1.0_wp),minMetricNorm*scale)
      nkeep=count(eig(1:ncur)>cutoff)
      if (nkeep<=0) exit
      call mctc_gemm(block(:,1:ncur),gram(1:ncur,ncur-nkeep+1:ncur),tmp(:,1:nkeep))
      do j=1,nkeep
         tmp(:,j)=tmp(:,j)/sqrt(eig(ncur-nkeep+j))
      end do
      block(:,1:nkeep)=tmp(:,1:nkeep)
      ncur=nkeep
   end do
   deallocate(sblock,gram,eig,tmp)
end subroutine metricOrthonormalizeBlockPair

subroutine orthogonalizeCorrectionsPair(pairMetric,npair,pairs,basis,corr,ncorr,nkeep)
   real(wp), intent(in) :: pairMetric(npair),basis(:,:)
   integer, intent(in) :: npair,pairs(2,npair),ncorr
   real(wp), intent(inout) :: corr(:,:)
   integer, intent(out) :: nkeep
   integer :: pass
   real(wp), allocatable :: scorr(:,:),coeff(:,:)

   allocate(scorr(size(corr,1),ncorr),coeff(size(basis,2),ncorr))
   do pass=1,2
      call applyGFN1PairMetricBlock(pairMetric,npair,pairs,corr(:,1:ncorr),scorr)
      call mctc_gemm(basis,scorr,coeff,transa='t')
      call mctc_gemm(basis,coeff,corr(:,1:ncorr),alpha=-1.0_wp,beta=1.0_wp)
   end do
   deallocate(scorr,coeff)
   call metricOrthonormalizeBlockPair(pairMetric,npair,pairs,corr(:,1:ncorr),ncorr,nkeep)
end subroutine orthogonalizeCorrectionsPair

end subroutine gfn1BlockDavidsonMF

subroutine gfn1BlockDavidson(H,S,seed,nroots,eval,C,success,niter,maxResidual)
   real(wp), intent(in) :: H(:,:), S(:,:), seed(:,:)
   integer, intent(in) :: nroots
   real(wp), intent(out) :: eval(:)
   real(wp), intent(out) :: C(:,:)
   logical, intent(out) :: success
   integer, intent(out) :: niter
   real(wp), intent(out) :: maxResidual

   integer :: n, nseed, maxSub, msub, iter, i, j, nadd, nkeep
   integer :: iconv, idx, best, info
   integer, allocatable :: order(:)
   logical, allocatable :: chosen(:)
   real(wp) :: denom, floorv, rel, hc2, sc2, r2
   real(wp), allocatable :: V(:,:), AV(:,:), BV(:,:)
   real(wp), allocatable :: proj(:,:), y(:,:), theta(:)
   real(wp), allocatable :: ritz(:,:), hritz(:,:), britz(:,:), resid(:,:)
   real(wp), allocatable :: corr(:,:)
   real(wp), allocatable :: relres(:), diagH(:), diagS(:)

   success = .false.
   niter = 0
   maxResidual = huge(1.0_wp)

   n = size(H,1)
   if (size(H,2) /= n .or. size(S,1) /= n .or. size(S,2) /= n) return
   if (nroots < 1 .or. nroots >= n) return
   if (size(C,1) < n .or. size(C,2) < nroots) return
   if (size(eval) < nroots) return
   nseed = min(size(seed,2), n)
   if (size(seed,1) < n .or. nseed < nroots) return

   ! Keep the reduced problem bounded.  SCC-to-SCC seeds are already close to
   ! the solution, so a modest correction space is more efficient than a large
   ! generic Davidson expansion.
   maxSub = min(n, max(nroots + 128, 3*nroots))
   if (maxSub <= nroots) return

   allocate(V(n,maxSub), AV(n,maxSub), BV(n,maxSub))
   allocate(ritz(n,nroots), hritz(n,nroots), britz(n,nroots), resid(n,nroots))
   allocate(relres(nroots), diagH(n), diagS(n))
   allocate(order(nroots), chosen(nroots))
   allocate(corr(n,maxCorrectionBlock))

   do i = 1, n
      diagH(i) = H(i,i)
      diagS(i) = 0.0_wp
   end do

   ! Start from the previous SCC eigenspace.  Orthonormalization in the S
   ! metric is deliberately repeated even though a converged prior eigenspace
   ! should already be orthonormal; it makes restart/checkpoint seeds robust.
   msub = min(nseed,maxSub)
   V(:,1:msub) = seed(:,1:msub)
   call metricOrthonormalizeBlock(S,V(:,1:msub),msub,nkeep)
   if (nkeep < nroots) return
   msub = nkeep
   call mctc_gemm(H,V(:,1:msub),AV(:,1:msub))
   call mctc_gemm(S,V(:,1:msub),BV(:,1:msub))

   do iter = 1, maxDavidsonIter
      niter = iter

      allocate(proj(msub,msub), y(msub,msub), theta(msub))
      call mctc_gemm(V(:,1:msub),AV(:,1:msub),proj,transa='t')
      ! Explicit symmetrization suppresses tiny BLAS roundoff asymmetry before
      ! the small Rayleigh-Ritz diagonalization.
      proj = 0.5_wp*(proj + transpose(proj))
      y = proj
      call diagonalizeSmall(y,theta,info)
      deallocate(proj)
      if (info /= 0) then
         deallocate(y,theta)
         return
      end if

      call mctc_gemm(V(:,1:msub),y(:,1:nroots),ritz)
      call mctc_gemm(AV(:,1:msub),y(:,1:nroots),hritz)
      call mctc_gemm(BV(:,1:msub),y(:,1:nroots),britz)

      maxResidual = 0.0_wp
      iconv = 0
      do j = 1, nroots
         resid(:,j) = hritz(:,j) - theta(j)*britz(:,j)
         r2 = sum(resid(:,j)*resid(:,j))
         hc2 = sum(hritz(:,j)*hritz(:,j))
         sc2 = sum(britz(:,j)*britz(:,j))
         denom = sqrt(max(0.0_wp,hc2)) + abs(theta(j))*sqrt(max(0.0_wp,sc2)) + tiny(1.0_wp)
         relres(j) = sqrt(max(0.0_wp,r2))/denom
         maxResidual = max(maxResidual,relres(j))
         if (relres(j) <= davidsonResidualTol) iconv = iconv + 1
      end do

      if (iconv == nroots) then
         eval(1:nroots) = theta(1:nroots)
         C(:,1:nroots) = ritz(:,1:nroots)
         success = .true.
         deallocate(y,theta)
         exit
      end if

      ! Select the largest unconverged residuals as a correction block.
      chosen = .false.
      nadd = 0
      do i = 1, min(maxCorrectionBlock,nroots-iconv)
         best = 0
         rel = -1.0_wp
         do j = 1, nroots
            if (.not.chosen(j) .and. relres(j) > davidsonResidualTol .and. relres(j) > rel) then
               best = j
               rel = relres(j)
            end if
         end do
         if (best == 0) exit
         chosen(best) = .true.
         nadd = nadd + 1
         order(nadd) = best
      end do
      if (nadd == 0) then
         deallocate(y,theta)
         exit
      end if

      do j = 1, nadd
         idx = order(j)
         do i = 1, n
            denom = theta(idx)*diagS(i) - diagH(i)
            floorv = max(preconditionerFloor,1.0e-6_wp*max(1.0_wp,abs(theta(idx))))
            if (abs(denom) < floorv) then
               if (denom < 0.0_wp) then
                  denom = -floorv
               else
                  denom = floorv
               end if
            end if
            corr(i,j) = resid(i,idx)/denom
         end do
      end do

      ! A thick restart keeps the Rayleigh-Ritz problem small.  Retain the
      ! current lowest Ritz vectors, then add freshly preconditioned residuals.
      if (msub + nadd > maxSub) then
         V(:,1:nroots) = ritz(:,1:nroots)
         AV(:,1:nroots) = hritz(:,1:nroots)
         BV(:,1:nroots) = britz(:,1:nroots)
         msub = nroots
      end if

      call orthogonalizeCorrections(S,V(:,1:msub),corr(:,1:nadd),nadd,nkeep)
      if (nkeep <= 0) then
         ! A diagonal preconditioner can occasionally map a nearly converged
         ! residual back into the current subspace.  The raw residual is still
         ! a valid Davidson expansion direction, so retry it before falling
         ! back to the full eigensolver.
         do j = 1, nadd
            corr(:,j) = resid(:,order(j))
         end do
         call orthogonalizeCorrections(S,V(:,1:msub),corr(:,1:nadd),nadd,nkeep)
      end if
      if (nkeep <= 0) then
         deallocate(y,theta)
         exit
      end if
      nadd = min(nkeep,maxSub-msub)
      if (nadd <= 0) then
         deallocate(y,theta)
         exit
      end if

      V(:,msub+1:msub+nadd) = corr(:,1:nadd)
      call mctc_gemm(H,V(:,msub+1:msub+nadd),AV(:,msub+1:msub+nadd))
      call mctc_gemm(S,V(:,msub+1:msub+nadd),BV(:,msub+1:msub+nadd))
      msub = msub + nadd

      deallocate(y,theta)
   end do

   deallocate(V,AV,BV,ritz,hritz,britz,resid,relres,diagH,diagS,order,chosen,corr)

contains

subroutine diagonalizeSmall(a,w,info)
   real(wp), intent(inout) :: a(:,:)
   real(wp), intent(out) :: w(:)
   integer, intent(out) :: info
   integer :: nn, lwork, liwork
   real(wp), allocatable :: lworkArray(:)
   integer, allocatable :: liworkArray(:)

   nn = size(a,1)
   lwork = max(1,1 + 6*nn + 2*nn*nn)
   liwork = max(1,3 + 5*nn)
   allocate(lworkArray(lwork),liworkArray(liwork))
   call lapack_syevd('v','u',nn,a,nn,w,lworkArray,lwork,liworkArray,liwork,info)
   deallocate(lworkArray,liworkArray)
end subroutine diagonalizeSmall

subroutine metricOrthonormalizeBlock(metric,block,ncol,nkeep)
   real(wp), intent(in) :: metric(:,:)
   real(wp), intent(inout) :: block(:,:)
   integer, intent(in) :: ncol
   integer, intent(out) :: nkeep
   integer :: j, info
   real(wp) :: scale, cutoff
   real(wp), allocatable :: sblock(:,:), gram(:,:), eig(:), tmp(:,:)

   allocate(sblock(size(block,1),ncol),gram(ncol,ncol),eig(ncol),tmp(size(block,1),ncol))
   call mctc_gemm(metric,block(:,1:ncol),sblock)
   call mctc_gemm(block(:,1:ncol),sblock,gram,transa='t')
   gram = 0.5_wp*(gram+transpose(gram))
   call diagonalizeSmall(gram,eig,info)
   if (info /= 0) then
      nkeep = 0
      deallocate(sblock,gram,eig,tmp)
      return
   end if
   scale = maxval(abs(eig))
   cutoff = max(tiny(1.0_wp),minMetricNorm*scale)
   nkeep = count(eig > cutoff)
   if (nkeep > 0) then
      ! SYEVD orders eigenvalues ascending; retain the well-conditioned metric
      ! directions at the upper end of the spectrum.
      call mctc_gemm(block(:,1:ncol),gram(:,ncol-nkeep+1:ncol),tmp(:,1:nkeep))
      do j = 1, nkeep
         tmp(:,j) = tmp(:,j)/sqrt(eig(ncol-nkeep+j))
      end do
      block(:,1:nkeep) = tmp(:,1:nkeep)
   end if
   deallocate(sblock,gram,eig,tmp)
end subroutine metricOrthonormalizeBlock

subroutine orthogonalizeCorrections(metric,basis,corr,ncorr,nkeep)
   real(wp), intent(in) :: metric(:,:), basis(:,:)
   real(wp), intent(inout) :: corr(:,:)
   integer, intent(in) :: ncorr
   integer, intent(out) :: nkeep
   integer :: pass, j, info
   real(wp) :: scale, cutoff
   real(wp), allocatable :: scorr(:,:), coeff(:,:), gram(:,:), eig(:), tmp(:,:)

   allocate(scorr(size(corr,1),ncorr),coeff(size(basis,2),ncorr))
   ! Two block modified-Gram-Schmidt passes in the S metric.  Using GEMM here
   ! is important: Davidson corrections are numerous enough that individual
   ! S*vector operations otherwise dominate the iterative solve.
   do pass = 1, 2
      call mctc_gemm(metric,corr(:,1:ncorr),scorr)
      call mctc_gemm(basis,scorr,coeff,transa='t')
      call mctc_gemm(basis,coeff,corr(:,1:ncorr),alpha=-1.0_wp,beta=1.0_wp)
   end do

   call mctc_gemm(metric,corr(:,1:ncorr),scorr)
   allocate(gram(ncorr,ncorr),eig(ncorr),tmp(size(corr,1),ncorr))
   call mctc_gemm(corr(:,1:ncorr),scorr,gram,transa='t')
   gram = 0.5_wp*(gram+transpose(gram))
   call diagonalizeSmall(gram,eig,info)
   if (info /= 0) then
      nkeep = 0
      deallocate(scorr,coeff,gram,eig,tmp)
      return
   end if
   scale = maxval(abs(eig))
   cutoff = max(tiny(1.0_wp),minMetricNorm*scale)
   nkeep = count(eig > cutoff)
   if (nkeep > 0) then
      call mctc_gemm(corr(:,1:ncorr),gram(:,ncorr-nkeep+1:ncorr),tmp(:,1:nkeep))
      do j = 1, nkeep
         tmp(:,j) = tmp(:,j)/sqrt(eig(ncorr-nkeep+j))
      end do
      corr(:,1:nkeep) = tmp(:,1:nkeep)
   end if
   deallocate(scorr,coeff,gram,eig,tmp)
end subroutine orthogonalizeCorrections

end subroutine gfn1BlockDavidson

end module xtb_gfn1_davidson
