! This file is part of xtb.
!
! GFN1-fast 4.0.0: canonical electronic and coupled shell-charge response
! for the analytic Hessian path. Retains the 2.0.0 closed-shell API.
!
! The routines in this module do not alter the production numerical Hessian.
! They provide first-order density and energy-weighted density response of the
! generalized eigenvalue problem
!
!     H C = S C eps,        C^T S C = I
!
! for a supplied pair of first derivatives dH and dS.  The response is the
! basic electronic building block required for a coupled SCC Hessian.  A
! shell Mulliken response is coupled to the isotropic potential Jacobian.
! Molecular geometry derivatives and full Hessian assembly are still needed.
module xtb_xtb_response
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : kB
   use xtb_mctc_blas_level3, only : mctc_gemm
   use xtb_mctc_lapack_trf, only : lapack_getrf
   use xtb_mctc_lapack_trs, only : lapack_getrs
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   public :: density_response_closed_shell
   public :: shell_charge_response
   public :: TGFN1ElectronicResponse, TGFN1CoupledResponse, density_response_spin

   ! GFN1-fast 4.0.0: cache the fixed-reference electronic susceptibility
   ! once, then apply it to nuclear derivatives and SCC charge perturbations.
   ! eps/dH/dW use atomic energy units (Eh), temperature is in kelvin.
   ! Shared alpha/beta orbitals, separate canonical particle constraints.
   ! Work arrays are local to apply, so the initialized cache can be shared
   ! by concurrent coordinate responses without changing the reference.
   type :: TGFN1ElectronicResponse
      private
      logical :: ready = .false.
      real(wp), allocatable :: c(:,:),eps(:),fp(:,:),numberWeights(:,:)
      real(wp), allocatable :: loewnerP(:,:),loewnerW(:,:),metricP(:,:),metricW(:,:)
   contains
      procedure :: initialize => initializeElectronicResponse
      procedure :: apply => applyElectronicResponse
   end type

   ! Cached dense solution of dq = dq_external + chi*dq. The supplied
   ! shift Jacobian is d(shellShift)/d(q_shell), in Eh per electron; it
   ! must include all enabled isotropic SCC/solvent terms at the reference.
   ! Geometry derivatives at fixed charge belong to apply's external dH.
   type :: TGFN1CoupledResponse
      private
      logical :: ready = .false.
      type(TGFN1ElectronicResponse) :: electronic
      integer, allocatable :: ao2sh(:),piv(:)
      real(wp), allocatable :: s(:,:),p(:,:),shiftJacobian(:,:),factor(:,:)
   contains
      procedure :: initialize => initializeCoupledResponse
      procedure :: apply => applyCoupledResponse
   end type

contains

subroutine initializeCoupledResponse(self,electronic,ao2sh,s,p,shiftJacobian,ok,rcond,threads)
   class(TGFN1CoupledResponse), intent(inout) :: self
   type(TGFN1ElectronicResponse), intent(in) :: electronic
   integer, intent(in) :: ao2sh(:)
   real(wp), intent(in) :: s(:,:),p(:,:),shiftJacobian(:,:)
   logical, intent(out) :: ok
   real(wp), intent(out), optional :: rcond
   integer, intent(in), optional :: threads
   real(wp), allocatable :: lhs(:,:),inverse(:,:),work(:,:),gram(:,:)
   real(wp) :: normA,normInverse,condition
   integer :: n,ns,i,j,info,threadCount
   logical, allocatable :: columnsOK(:)
   logical :: responseOK
   self%ready=.false.; ok=.false.
   if(present(rcond)) rcond=0.0_wp
   if(.not.electronic%ready) return
   n=size(electronic%c,1); ns=size(shiftJacobian,1)
   if(ns<1.or.size(shiftJacobian,2)/=ns.or.size(ao2sh)/=n) return
   if(size(s,1)/=n.or.size(s,2)/=n.or.size(p,1)/=n.or.size(p,2)/=n) return
   if(minval(ao2sh)<1.or.maxval(ao2sh)>ns) return
   if(.not.all(ieee_is_finite(s)).or..not.all(ieee_is_finite(p))) return
   if(.not.all(ieee_is_finite(shiftJacobian))) return
   allocate(lhs(ns,ns),inverse(ns,ns),work(n,n),gram(n,n))
   call mctc_gemm(s,electronic%c,work)
   call mctc_gemm(electronic%c,work,gram,transa='t')
   do i=1,n; gram(i,i)=gram(i,i)-1.0_wp; enddo
   if(maxval(abs(gram))>1.0e-9_wp) return
   if(maxval(abs(s-transpose(s)))>1.0e-12_wp.or.maxval(abs(p-transpose(p)))>1.0e-12_wp) return
   self%electronic=electronic; self%ao2sh=ao2sh
   self%s=s; self%p=p; self%shiftJacobian=shiftJacobian
   threadCount=1
   if(present(threads)) threadCount=max(1,min(threads,ns))
   allocate(columnsOK(ns));columnsOK=.false.
   !$omp parallel num_threads(threadCount) default(none) &
   !$omp shared(electronic,ao2sh,s,p,shiftJacobian,lhs,columnsOK,n,ns)
   call susceptibilityWorker()
   !$omp end parallel
   if(.not.all(columnsOK)) return
   self%factor=lhs
   if(allocated(self%piv)) deallocate(self%piv)
   allocate(self%piv(ns))
   call lapack_getrf(ns,ns,self%factor,ns,self%piv,info)
   if(info/=0) return
   ! Exact dense 1-norm reciprocal conditioning check; the shell matrix is
   ! much smaller than AO response matrices. Reject unstable SCC response
   ! instead of shifting the physical response denominators.
   inverse=0.0_wp
   do i=1,ns; inverse(i,i)=1.0_wp; enddo
   call lapack_getrs('n',ns,ns,self%factor,ns,self%piv,inverse,ns,info)
   if(info/=0.or..not.all(ieee_is_finite(inverse))) return
   normA=maxval(sum(abs(lhs),dim=1)); normInverse=maxval(sum(abs(inverse),dim=1))
   condition=(1.0_wp/normA)/normInverse
   if(present(rcond)) rcond=condition
   if(.not.ieee_is_finite(condition).or.condition<1.0e-12_wp) return
   self%ready=.true.; ok=.true.
contains
   ! Allocate simple scratch only after entering the OpenMP worker. Avoid
   ! private copies of nested allocatable Fortran descriptors on Windows ifx.
   recursive subroutine susceptibilityWorker()
      real(wp), allocatable :: localH(:,:),localS(:,:),localP(:,:),localW(:,:),localQ(:)
      integer :: column,allocationStatus
      logical :: valid
      allocate(localH(n,n),localS(n,n),localP(n,n),localW(n,n),localQ(ns),stat=allocationStatus)
      if(allocationStatus==0) localS=0.0_wp
      !$omp do schedule(static)
      do column=1,ns
         if(allocationStatus/=0) cycle
         call shellShiftPerturbation(ao2sh,s,shiftJacobian(:,column),localH)
         call electronic%apply(localH,localS,localP,localW,valid)
         if(.not.valid) cycle
         call shell_charge_response(ao2sh,s,p,localS,localP,localQ)
         lhs(:,column)=-localQ;lhs(column,column)=lhs(column,column)+1.0_wp
         columnsOK(column)=.true.
      enddo
      !$omp end do
   end subroutine susceptibilityWorker
end subroutine

subroutine applyCoupledResponse(self,dh,ds,dp,dw,dqsh,ok)
   class(TGFN1CoupledResponse), intent(in) :: self
   real(wp), intent(in) :: dh(:,:),ds(:,:)
   real(wp), intent(out) :: dp(:,:),dw(:,:),dqsh(:)
   logical, intent(out) :: ok
   real(wp), allocatable :: rhs(:,:),shift(:),perturbation(:,:),totalH(:,:),check(:)
   real(wp), allocatable, target :: metric(:,:)
   integer :: n,ns,info
   logical :: responseOK
   ok=.false.; dp=0.0_wp; dw=0.0_wp; dqsh=0.0_wp
   if(.not.self%ready) return
   n=size(self%s,1); ns=size(self%shiftJacobian,1)
   if(size(dqsh)/=ns) return
   if(any(shape(ds)/=[n,n])) return
   if(.not.all(ieee_is_finite(ds))) return
   allocate(rhs(ns,1),shift(ns),perturbation(n,n),totalH(n,n),check(ns),metric(n,n))
   ! Both electronic applications have the same moving AO metric. Reuse
   ! C^T dS C within this coordinate, with worker-local storage only. The
   ! perturbation buffer is scratch here and becomes the SCC shift below.
   call mctc_gemm(ds,self%electronic%c,perturbation)
   call mctc_gemm(self%electronic%c,perturbation,metric,transa='t')
   call self%electronic%apply(dh,ds,dp,dw,responseOK,cachedMetric=metric)
   if(.not.responseOK) return
   call shell_charge_response(self%ao2sh,self%s,self%p,ds,dp,check)
   rhs(:,1)=check
   call lapack_getrs('n',ns,1,self%factor,ns,self%piv,rhs,ns,info)
   if(info/=0.or..not.all(ieee_is_finite(rhs))) then
      dp=0.0_wp; dw=0.0_wp; return
   endif
   dqsh=rhs(:,1); shift=matmul(self%shiftJacobian,dqsh)
   call shellShiftPerturbation(self%ao2sh,self%s,shift,perturbation)
   totalH=dh+perturbation
   call self%electronic%apply(totalH,ds,dp,dw,responseOK,cachedMetric=metric)
   if(responseOK) then
      call shell_charge_response(self%ao2sh,self%s,self%p,ds,dp,check)
      responseOK=maxval(abs(check-dqsh))<1.0e-10_wp*max(1.0_wp,maxval(abs(dqsh)))
   endif
   if(.not.responseOK) then
      dp=0.0_wp; dw=0.0_wp; dqsh=0.0_wp; return
   endif
   ok=.true.
end subroutine

! GFN1 isotropic Hamiltonian: H_ij = H0_ij - S_ij*(v_i+v_j)/2.
pure subroutine shellShiftPerturbation(ao2sh,s,shift,dh)
   integer, intent(in) :: ao2sh(:)
   real(wp), intent(in) :: s(:,:),shift(:)
   real(wp), intent(out) :: dh(:,:)
   integer :: i,j
   do j=1,size(s,2)
      do i=1,size(s,1)
         dh(i,j)=-0.5_wp*s(i,j)*(shift(ao2sh(i))+shift(ao2sh(j)))
      enddo
   enddo
end subroutine

!> Independent-spin canonical first-order response, including dW and the
!> moving AO metric. Caller supplies the converged orbitals/occupations.
!> This is a building block: dH must include the coupled SCC potential
!> response before the result can be used for a complete molecular Hessian.
subroutine density_response_spin(c,eps,occupation,temperature,dh,ds,dp,dw,dmu,ok,gap_tol)
   real(wp), intent(in) :: c(:,:),eps(:),occupation(:,:),temperature,dh(:,:),ds(:,:)
   real(wp), intent(out) :: dp(:,:),dw(:,:),dmu(2)
   logical, intent(out) :: ok
   real(wp), intent(in), optional :: gap_tol
   type(TGFN1ElectronicResponse) :: response
   dp=0.0_wp; dw=0.0_wp; dmu=0.0_wp
   call response%initialize(c,eps,occupation,temperature,ok,gap_tol)
   if(.not.ok) return
   call response%apply(dh,ds,dp,dw,ok,dmu)
end subroutine

subroutine initializeElectronicResponse(self,c,eps,occupation,temperature,ok,gap_tol)
   class(TGFN1ElectronicResponse), intent(inout) :: self
   real(wp), intent(in) :: c(:,:),eps(:),occupation(:,:),temperature
   logical, intent(out) :: ok
   real(wp), intent(in), optional :: gap_tol
   real(wp), allocatable :: f(:,:),total(:)
   real(wp) :: kt,beta,tol,delta,lf,fi,fj,scale,chemical,mu,weight
   integer :: n,i,j,spin
   logical :: thermal,hasChemical
   self%ready=.false.; ok=.false.; n=size(c,1)
   if(n<1.or.size(c,2)/=n.or.size(eps)/=n) return
   if(size(occupation,1)/=n.or.size(occupation,2)/=2) return
   if(.not.all(ieee_is_finite(c)).or..not.all(ieee_is_finite(eps))) return
   if(.not.all(ieee_is_finite(occupation)).or..not.ieee_is_finite(temperature)) return
   if(temperature<0.0_wp.or.any(occupation<0.0_wp).or.any(occupation>1.0_wp)) return
   if(any(eps(2:n)<eps(1:n-1))) return
   tol=100.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(eps)))
   if(present(gap_tol)) then
      if(.not.ieee_is_finite(gap_tol).or.gap_tol<0.0_wp) return
      tol=gap_tol
   endif
   thermal=temperature>0.1_wp; beta=0.0_wp
   if(thermal) then
      kt=kB*temperature; beta=1.0_wp/kt
      if(.not.ieee_is_finite(beta)) return
   else
      ! A discontinuous zero-temperature occupation is differentiable only
      ! for integer filling separated from a differently occupied state.
      if(any(abs(occupation-anint(occupation))>1.0e-12_wp)) return
   endif
   f=occupation
   if(.not.thermal) f=anint(f)
   do spin=1,2
      if(any(f(2:n,spin)>f(1:n-1,spin)+1.0e-12_wp)) return
      if(.not.thermal) cycle
      hasChemical=.false.; chemical=0.0_wp
      do i=1,n
         if(f(i,spin)>1.0e-9_wp.and.f(i,spin)<1.0_wp-1.0e-9_wp) then
            mu=eps(i)-kt*log((1.0_wp-f(i,spin))/f(i,spin))
            if(hasChemical) then
               if(abs(mu-chemical)>1.0e-8_wp*max(kt,1.0e-2_wp)) return
            else
               chemical=mu; hasChemical=.true.
            endif
         endif
      enddo
   enddo
   self%c=c; self%eps=eps
   self%fp=-beta*f*(1.0_wp-f)
   self%numberWeights=self%fp
   do spin=1,2
      scale=maxval(abs(self%fp(:,spin)))
      self%numberWeights(:,spin)=0.0_wp
      if(scale>0.0_wp) then
         ! Normalize before summing, even when all derivatives are tiny.
         self%numberWeights(:,spin)=self%fp(:,spin)/scale
         weight=sum(self%numberWeights(:,spin))
         self%numberWeights(:,spin)=self%numberWeights(:,spin)/weight
      endif
   enddo
   total=f(:,1)+f(:,2)
   self%loewnerP=spread(total,2,n)*0.0_wp
   self%loewnerW=self%loewnerP; self%metricP=self%loewnerP; self%metricW=self%loewnerP
   do j=1,n
      do i=1,n
         delta=eps(i)-eps(j)
         if(.not.ieee_is_finite(delta)) return
         do spin=1,2
            fi=f(i,spin); fj=f(j,spin)
            if(thermal) then
               if(delta==0.0_wp) then
                  if(abs(fi-fj)>1.0e-12_wp) return
                  lf=0.5_wp*(self%fp(i,spin)+self%fp(j,spin))
               else
                  ! Exact Fermi divided difference in a form that avoids
                  ! subtracting almost equal occupations near degeneracy.
                  lf=-beta*max(fi,fj)*(1.0_wp-min(fi,fj))*fermiDifferenceRatio(beta*abs(delta))
               endif
            else
               lf=0.0_wp
               if(fi/=fj) then
                  if(abs(delta)<=tol) return
                  lf=(fi-fj)/delta
               endif
            endif
            self%loewnerP(i,j)=self%loewnerP(i,j)+lf
            self%loewnerW(i,j)=self%loewnerW(i,j)+0.5_wp*(fi+fj)+0.5_wp*(eps(i)+eps(j))*lf
         enddo
         self%metricP(i,j)=0.5_wp*(total(i)+total(j))
         self%metricW(i,j)=0.5_wp*(eps(i)*total(i)+eps(j)*total(j))
      enddo
   enddo
   if(.not.all(ieee_is_finite(self%loewnerP)).or..not.all(ieee_is_finite(self%loewnerW))) return
   self%ready=.true.; ok=.true.
end subroutine

subroutine applyElectronicResponse(self,dh,ds,dp,dw,ok,dmu,cachedMetric)
   class(TGFN1ElectronicResponse), intent(in) :: self
   real(wp), intent(in) :: dh(:,:),ds(:,:)
   real(wp), intent(out) :: dp(:,:),dw(:,:)
   logical, intent(out) :: ok
   real(wp), intent(out), optional :: dmu(2)
   real(wp), intent(in), optional, target :: cachedMetric(:,:)
   real(wp), allocatable, target :: ownedMetric(:,:)
   real(wp), pointer :: sigma(:,:)
   real(wp), allocatable :: k(:,:),work(:,:),pm(:,:),wm(:,:),diagonal(:)
   real(wp) :: chemical(2),tol
   integer :: n,i,j
   ok=.false.; dp=0.0_wp; dw=0.0_wp
   if(present(dmu)) dmu=0.0_wp
   if(.not.self%ready) return
   n=size(self%c,1)
   if(size(dh,1)/=n.or.size(dh,2)/=n.or.size(ds,1)/=n.or.size(ds,2)/=n) return
   if(size(dp,1)/=n.or.size(dp,2)/=n.or.size(dw,1)/=n.or.size(dw,2)/=n) return
   if(.not.all(ieee_is_finite(dh)).or..not.all(ieee_is_finite(ds))) return
   tol=100.0_wp*epsilon(1.0_wp)
   if(maxval(abs(dh-transpose(dh)))>tol*max(1.0_wp,maxval(abs(dh)))) return
   if(maxval(abs(ds-transpose(ds)))>tol*max(1.0_wp,maxval(abs(ds)))) return
   allocate(k(n,n),work(n,n),pm(n,n),wm(n,n),diagonal(n))
   if(present(cachedMetric)) then
      ! Internal coupled-response callers supply the transformed derivative
      ! of this same dS; the uncached public path remains available.
      if(any(shape(cachedMetric)/=[n,n])) return
      if(.not.all(ieee_is_finite(cachedMetric))) return
      sigma=>cachedMetric
   else
      allocate(ownedMetric(n,n));sigma=>ownedMetric
      sigma=0.0_wp
      if(any(ds/=0.0_wp)) then
         call mctc_gemm(ds,self%c,work)
         call mctc_gemm(self%c,work,sigma,transa='t')
      endif
   endif
   call mctc_gemm(dh,self%c,work)
   call mctc_gemm(self%c,work,k,transa='t')
   do j=1,n
      do i=1,n
         k(i,j)=k(i,j)-0.5_wp*(self%eps(i)+self%eps(j))*sigma(i,j)
      enddo
      diagonal(j)=k(j,j)
   enddo
   chemical=matmul(diagonal,self%numberWeights)
   pm=self%loewnerP*k-self%metricP*sigma
   wm=self%loewnerW*k-self%metricW*sigma
   do i=1,n
      pm(i,i)=pm(i,i)-sum(self%fp(i,:)*chemical)
      wm(i,i)=wm(i,i)-self%eps(i)*sum(self%fp(i,:)*chemical)
   enddo
   call mctc_gemm(self%c,pm,work)
   call mctc_gemm(work,self%c,dp,transb='t')
   call mctc_gemm(self%c,wm,work)
   call mctc_gemm(work,self%c,dw,transb='t')
   if(.not.all(ieee_is_finite(dp)).or..not.all(ieee_is_finite(dw))) then
      dp=0.0_wp; dw=0.0_wp; return
   endif
   dp=0.5_wp*(dp+transpose(dp)); dw=0.5_wp*(dw+transpose(dw))
   if(present(dmu)) dmu=chemical
   ok=.true.
end subroutine

! (1-exp(-z))/z evaluated without cancellation for z near zero.
pure real(wp) function fermiDifferenceRatio(z) result(value)
   real(wp), intent(in) :: z
   real(wp) :: term
   integer :: k
   if(z<0.01_wp) then
      value=1.0_wp; term=1.0_wp
      do k=1,8
         term=term*(-z)/real(k+1,wp); value=value+term
      enddo
   else
      value=(1.0_wp-exp(-z))/z
   endif
end function

!> Closed-shell first-order AO density response for a generalized symmetric
!> eigenvalue problem.  The occupied space is assumed to consist of the first
!> nocc eigenvectors and to carry occupation 2.0.  No approximation or
!> denominator regularisation is used: if an occupied--virtual gap is smaller
!> than gap_tol, ok is returned false so callers can fall back to the numerical
!> Hessian path.
subroutine density_response_closed_shell(C, eps, nocc, dH, dS, dP, ok, gap_tol)
   real(wp), intent(in) :: C(:, :)
   real(wp), intent(in) :: eps(:)
   integer, intent(in) :: nocc
   real(wp), intent(in) :: dH(:, :)
   real(wp), intent(in) :: dS(:, :)
   real(wp), intent(out) :: dP(:, :)
   logical, intent(out) :: ok
   real(wp), intent(in), optional :: gap_tol

   integer :: nao, nvir, i, a
   real(wp) :: denom, tol
   real(wp), allocatable :: dh_cocc(:, :), ds_cocc(:, :)
   real(wp), allocatable :: hco(:, :), sco(:, :), uvo(:, :)
   real(wp), allocatable :: xvo(:, :), yoo(:, :), soo(:, :)

   nao = size(C, 1)
   dP = 0.0_wp
   ok = .false.

   if (size(C, 2) /= nao) return
   if (size(eps) < nao) return
   if (size(dH, 1) /= nao .or. size(dH, 2) /= nao) return
   if (size(dS, 1) /= nao .or. size(dS, 2) /= nao) return
   if (size(dP, 1) /= nao .or. size(dP, 2) /= nao) return
   if (nocc < 0 .or. nocc > nao) return

   if (present(gap_tol)) then
      tol = max(0.0_wp, gap_tol)
   else
      tol = 100.0_wp*epsilon(1.0_wp)
   end if

   if (nocc == 0) then
      ok = .true.
      return
   end if

   nvir = nao - nocc

   allocate(dh_cocc(nao, nocc), ds_cocc(nao, nocc))
   allocate(hco(nao, nocc), sco(nao, nocc))
   allocate(soo(nocc, nocc), yoo(nao, nocc))

   ! Only C^T dH C_occ and C^T dS C_occ are required.  Avoid constructing
   ! full MO-basis derivative matrices because Hessian response is repeated
   ! for every nuclear Cartesian coordinate.
   dh_cocc = matmul(dH, C(:, 1:nocc))
   ds_cocc = matmul(dS, C(:, 1:nocc))
   hco = matmul(transpose(C), dh_cocc)
   sco = matmul(transpose(C), ds_cocc)
   soo = sco(1:nocc, 1:nocc)

   if (nvir > 0) then
      allocate(uvo(nvir, nocc), xvo(nao, nocc))
      do i = 1, nocc
         do a = 1, nvir
            denom = eps(i) - eps(nocc + a)
            if (abs(denom) <= tol) then
               dP = 0.0_wp
               return
            end if
            uvo(a, i) = (hco(nocc + a, i) - eps(i)*sco(nocc + a, i))/denom
         end do
      end do

      ! Orbital-relaxation part: 2*(C_v U_vo C_o^T + transpose).
      xvo = matmul(C(:, nocc+1:nao), uvo)
      dP = 2.0_wp*(matmul(xvo, transpose(C(:, 1:nocc))) &
         & + matmul(C(:, 1:nocc), transpose(xvo)))
   end if

   ! Metric response from d(C^T S C)=0.  Choosing the symmetric occupied-
   ! occupied gauge U_oo=-1/2 S^x_oo gives the gauge-invariant projector
   ! contribution -2*C_occ*S^x_oo*C_occ^T.
   yoo = matmul(C(:, 1:nocc), soo)
   dP = dP - 2.0_wp*matmul(yoo, transpose(C(:, 1:nocc)))

   ! Suppress only the insignificant antisymmetric roundoff generated by
   ! independent matrix products.  This is mathematically identical to the
   ! symmetric density response.
   dP = 0.5_wp*(dP + transpose(dP))
   ok = .true.

end subroutine density_response_closed_shell


!> First-order shell-resolved Mulliken *charge* response.
!>
!> GFN1 stores q_sh = z_sh - population_sh.  Therefore
!>
!>   dq_sh = - d population_sh
!>         = - Mulliken(dP*S + P*dS).
!>
!> The accumulation order mirrors mpopsh in scc_core.f90 so this routine can
!> later be inserted into the coupled SCC response without changing the
!> population convention.
subroutine shell_charge_response(ao2sh, S, P, dS, dP, dqsh)
   integer, intent(in) :: ao2sh(:)
   real(wp), intent(in) :: S(:, :)
   real(wp), intent(in) :: P(:, :)
   real(wp), intent(in) :: dS(:, :)
   real(wp), intent(in) :: dP(:, :)
   real(wp), intent(out) :: dqsh(:)

   integer :: nao, i, j, ish, jsh
   real(wp) :: dps

   nao = size(ao2sh)
   dqsh = 0.0_wp

   do i = 1, nao
      ish = ao2sh(i)
      do j = 1, i-1
         jsh = ao2sh(j)
         dps = dP(j,i)*S(j,i) + P(j,i)*dS(j,i)
         dqsh(ish) = dqsh(ish) - dps
         dqsh(jsh) = dqsh(jsh) - dps
      end do
      dps = dP(i,i)*S(i,i) + P(i,i)*dS(i,i)
      dqsh(ish) = dqsh(ish) - dps
   end do

end subroutine shell_charge_response

end module xtb_xtb_response
