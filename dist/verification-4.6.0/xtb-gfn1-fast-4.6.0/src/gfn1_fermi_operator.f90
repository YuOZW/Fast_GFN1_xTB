! This file is part of xtb.
! SPDX-License-Identifier: LGPL-3.0-or-later
!
! Dense, diagonalization-free reference implementation. All matrix
! functions use polynomial products and Cholesky/LDL linear solves.
! Recursive Fermi expansion: Niklasson, cond-mat/0311469.
module xtb_gfn1_fermi_operator
   use xtb_mctc_accuracy, only : dp
   use xtb_mctc_constants, only : kB
   use xtb_mctc_convert, only : autoev,evtoau
   use xtb_mctc_blas_level3, only : mctc_gemm,blas_trsm
   use xtb_mctc_lapack_gst, only : lapack_sygst
   use xtb_mctc_lapack_trf, only : lapack_potrf,lapack_sytrf
   use xtb_mctc_lapack_trs, only : lapack_potrs
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: TGFN1FermiOperator,fermi_matrix,purify_density
   type :: TGFN1FermiOperator
      real(dp), allocatable :: metric_factor(:,:)
   contains
      procedure :: initialize => initializeMetric
      procedure :: solve => solveDensity
   end type
contains

subroutine initializeMetric(self,s,ok)
   class(TGFN1FermiOperator), intent(inout) :: self
   real(dp), intent(in) :: s(:,:)
   logical, intent(out) :: ok
   integer :: n,info
   ok=.false.
   if(allocated(self%metric_factor)) deallocate(self%metric_factor)
   n=size(s,1)
   if(n<1.or.size(s,2)/=n.or..not.all(ieee_is_finite(s))) return
   self%metric_factor=s
   call lapack_potrf('u',n,self%metric_factor,n,info)
   if(info/=0) then
      deallocate(self%metric_factor)
      return
   endif
   ok=.true.
end subroutine

! H/W are in eV; P is electron density and entropy terms are in Eh.
! At T=0 use separate alpha/beta integer projectors. At finite T solve
! chemical potentials separately, including fractional/open-shell filling.
subroutine solveDensity(self,h,nel,nopen,temperature,mua,mub,p,w,ga,gb,ok,steps,numberError,maxSteps)
   class(TGFN1FermiOperator), intent(in) :: self
   real(dp), intent(in) :: h(:,:),temperature
   integer, intent(in) :: nel,nopen
   real(dp), intent(inout) :: mua,mub
   real(dp), intent(out) :: p(:,:),w(:,:),ga,gb,numberError
   logical, intent(out) :: ok
   integer, intent(out) :: steps
   integer, intent(in), optional :: maxSteps
   real(dp), allocatable :: a(:,:),da(:,:),db(:,:),d(:,:),aw(:,:)
   integer :: n,na,nb,info,sa,sb,limit
   logical :: spinOK
   ok=.false.; steps=0; ga=0.0_dp; gb=0.0_dp; numberError=huge(1.0_dp)
   p=0.0_dp; w=0.0_dp
   if(.not.allocated(self%metric_factor)) return
   n=size(h,1); limit=100
   if(present(maxSteps)) limit=maxSteps
   if(n<1.or.size(h,2)/=n.or.limit<1) return
   if(size(self%metric_factor,1)/=n.or.size(p,1)/=n.or.size(p,2)/=n) return
   if(size(w,1)/=n.or.size(w,2)/=n) return
   if(.not.all(ieee_is_finite(h)).or..not.ieee_is_finite(temperature).or.temperature<0.0_dp) return
   if(nel<0.or.nel>2*n.or.nopen<0.or.nopen>nel.or.mod(nel+nopen,2)/=0) return
   na=(nel+nopen)/2; nb=(nel-nopen)/2
   if(na>n) return
   allocate(a(n,n),da(n,n),db(n,n),d(n,n),aw(n,n))
   a=h
   call lapack_sygst(1,'u',n,a,n,self%metric_factor,n,info)
   if(info/=0) return
   call completeUpper(a)
   if(temperature<=0.1_dp) then
      call purify_density(a,na,da,spinOK,sa,limit)
   else
      call canonicalFermi(a,na,temperature,mua,da,ga,spinOK,sa,limit)
   endif
   steps=sa
   if(.not.spinOK) return
   if(na==nb) then
      db=da; gb=ga; mub=mua; sb=0
   else
      if(temperature<=0.1_dp) then
         call purify_density(a,nb,db,spinOK,sb,limit)
      else
         call canonicalFermi(a,nb,temperature,mub,db,gb,spinOK,sb,limit)
      endif
      steps=sa+sb
      if(.not.spinOK) return
   endif
   steps=sa+sb; d=da+db
   numberError=abs(trace(d)-real(nel,dp))
   if(numberError>2.0e-10_dp) return
   call mctc_gemm(a,d,aw)
   aw=0.5_dp*(aw+transpose(aw))
   p=d; w=aw
   call blas_trsm('l','u','n','n',n,n,1.0_dp,self%metric_factor,n,p,n)
   call blas_trsm('r','u','t','n',n,n,1.0_dp,self%metric_factor,n,p,n)
   call blas_trsm('l','u','n','n',n,n,1.0_dp,self%metric_factor,n,w,n)
   call blas_trsm('r','u','t','n',n,n,1.0_dp,self%metric_factor,n,w,n)
   p=0.5_dp*(p+transpose(p)); w=0.5_dp*(w+transpose(w))
   ok=all(ieee_is_finite(p)).and.all(ieee_is_finite(w)).and.ieee_is_finite(ga+gb)
end subroutine

! Trace-correcting second-order spectral projection, with rigorous initial
! Gershgorin bounds. No entry truncation or altered physical cutoff is used.
subroutine purify_density(a,nocc,d,ok,steps,maxSteps)
   real(dp), intent(in) :: a(:,:)
   integer, intent(in) :: nocc,maxSteps
   real(dp), intent(out) :: d(:,:)
   logical, intent(out) :: ok
   integer, intent(out) :: steps
   real(dp), allocatable :: square(:,:),commutator(:,:),tmp(:,:)
   real(dp) :: lo,hi,tr,idempotency
   integer :: n,k
   ok=.false.; steps=0; d=0.0_dp; n=size(a,1)
   if(n<1.or.size(a,2)/=n.or.size(d,1)/=n.or.size(d,2)/=n) return
   if(nocc<0.or.nocc>n.or.maxSteps<1.or..not.all(ieee_is_finite(a))) return
   if(nocc==0.or.nocc==n) then
      if(nocc==n) call addDiagonal(d,1.0_dp)
      ok=.true.; return
   endif
   call bounds(a,lo,hi)
   if(hi-lo<=epsilon(1.0_dp)*max(1.0_dp,abs(lo),abs(hi))) return
   allocate(square(n,n),commutator(n,n),tmp(n,n))
   d=-a/(hi-lo)
   call addDiagonal(d,hi/(hi-lo))
   do k=1,maxSteps
      steps=k
      call mctc_gemm(d,d,square)
      tr=trace(d); idempotency=maxval(abs(square-d))
      if(idempotency<2.0e-13_dp.and.abs(tr-real(nocc,dp))<2.0e-11_dp) then
         call mctc_gemm(a,d,commutator)
         call mctc_gemm(d,a,tmp)
         if(maxval(abs(commutator-tmp))>2.0e-11_dp*max(1.0_dp,maxval(abs(a)))) return
         ok=.true.; return
      endif
      if(tr>real(nocc,dp)) then
         d=square
      else
         d=2.0_dp*d-square
      endif
      d=0.5_dp*(d+transpose(d))
      if(.not.all(ieee_is_finite(d))) return
   enddo
end subroutine

subroutine canonicalFermi(a,nocc,temperature,mu,d,entropy,ok,steps,maxSteps)
   real(dp), intent(in) :: a(:,:),temperature
   integer, intent(in) :: nocc,maxSteps
   real(dp), intent(inout) :: mu
   real(dp), intent(out) :: d(:,:),entropy
   logical, intent(out) :: ok
   integer, intent(out) :: steps
   real(dp), allocatable :: trial(:,:),square(:,:),projector(:,:),af(:,:),holes(:,:)
   real(dp) :: beta,lo,hi,left,right,delta,slope,next,logPartition,rawEntropy,cut,tail
   integer :: n,iter,m,depth,count,lowCount,highCount,projectionSteps
   logical :: matrixOK
   ok=.false.; steps=0; entropy=0.0_dp; n=size(a,1)
   d=0.0_dp; call bounds(a,lo,hi)
   beta=1.0_dp/(kB*autoev*temperature)
   if(nocc==0.or.nocc==n) then
      if(nocc==n) call addDiagonal(d,1.0_dp)
      mu=merge(hi+50.0_dp/beta,lo-50.0_dp/beta,nocc==n)
      ok=.true.; return
   endif
   allocate(trial(n,n),square(n,n),projector(n,n),af(n,n),holes(n,n))
   left=lo-50.0_dp/beta; right=hi+50.0_dp/beta
   if(.not.ieee_is_finite(mu)) mu=0.5_dp*(lo+hi)
   mu=max(left,min(right,mu))
   do iter=1,30
      call fermi_matrix(a,mu,beta,d,logPartition,matrixOK,depth,maxSteps)
      steps=steps+depth
      if(.not.matrixOK) return
      delta=real(nocc,dp)-trace(d)
      if(abs(delta)<5.0e-12_dp) exit
      if(delta>0.0_dp) then
         left=mu
      else
         right=mu
      endif
      call mctc_gemm(d,d,square)
      slope=beta*max(0.0_dp,trace(d-square))
      next=0.5_dp*(left+right)
      if(slope>1.0e-14_dp) then
         if(mu+delta/slope>left.and.mu+delta/slope<right) next=mu+delta/slope
      endif
      mu=next
   enddo
   if(abs(delta)>=5.0e-12_dp) return
   ! Check expansion error against one additional temperature-doubling step.
   call fermi_matrix(a,mu,beta,trial,tail,matrixOK,m,maxSteps,depth+1)
   steps=steps+m
   if(.not.matrixOK.or.maxval(abs(d-trial))>2.0e-11_dp) return
   ! Match fermismear's 1e-9 entropy convention. The entropy of excluded
   ! tails is x*f+f with error O(f^2), <1e-18 per state at this threshold.
   ! Inertia counts determine the spectral projectors without eigenvalues.
   cut=log((1.0_dp-1.0e-9_dp)/1.0e-9_dp)/beta
   call countBelow(a,mu-cut,lowCount,matrixOK)
   if(.not.matrixOK) return
   call countBelow(a,mu+cut,highCount,matrixOK)
   if(.not.matrixOK) return
   if(lowCount==nocc.and.highCount==nocc) then
      ! Inertia certifies that no occupation belongs to the entropy window.
      ! fermismear returns exactly zero, avoiding cancellation of two large
      ! band-energy/log-partition traces and two redundant tail projectors.
      entropy=0.0_dp; ok=.true.; return
   endif
   af=a; call addDiagonal(af,-mu)
   rawEntropy=-(sum(af*d)+logPartition/beta)
   tail=0.0_dp
   count=highCount
   if(count<n) then
      call purify_density(a,count,projector,matrixOK,projectionSteps,maxSteps)
      steps=steps+projectionSteps
      if(.not.matrixOK) return
      projector=-projector; call addDiagonal(projector,1.0_dp)
      call mctc_gemm(af,d,square)
      tail=tail+sum(projector*(square+d/beta))
   endif
   count=lowCount
   if(count>0) then
      call purify_density(a,count,projector,matrixOK,projectionSteps,maxSteps)
      steps=steps+projectionSteps
      if(.not.matrixOK) return
      holes=-d; call addDiagonal(holes,1.0_dp)
      call mctc_gemm(af,holes,square)
      tail=tail+sum(projector*(-square+holes/beta))
   endif
   entropy=(rawEntropy+tail)*evtoau
   if(.not.ieee_is_finite(entropy).or.entropy>2.0e-11_dp) return
   ! A positive value of at most roundoff is zero physical entropy.
   entropy=min(0.0_dp,entropy)
   ok=.true.
end subroutine

! [6/6] Pade exponential initializes the recursive rational Fermi map at
! small beta. Each iteration doubles beta: F <- F^2/(F^2+(I-F)^2).
! log Z doubles with logdet(F^2+(I-F)^2), so entropy needs no matrix log.
subroutine fermi_matrix(a,mu,beta,d,logPartition,ok,depth,maxDepth,forceDepth)
   real(dp), intent(in) :: a(:,:),mu,beta
   real(dp), intent(out) :: d(:,:),logPartition
   logical, intent(out) :: ok
   integer, intent(out) :: depth
   integer, intent(in) :: maxDepth
   integer, intent(in), optional :: forceDepth
   real(dp), allocatable :: x(:,:),x2(:,:),even(:,:),odd(:,:),tmp(:,:),chol(:,:),rhs(:,:)
   real(dp) :: lo,hi,scale,coeff(0:6),evenLog,plusLog,denominatorLog
   integer :: n,i,k,info
   ok=.false.; depth=0; logPartition=0.0_dp; d=0.0_dp; n=size(a,1)
   if(n<1.or.size(a,2)/=n.or.size(d,1)/=n.or.size(d,2)/=n.or.maxDepth<1) return
   if(.not.all(ieee_is_finite(a)).or..not.ieee_is_finite(mu).or..not.ieee_is_finite(beta)) return
   if(beta<=0.0_dp) return
   call bounds(a,lo,hi)
   scale=beta*max(abs(lo-mu),abs(hi-mu))
   if(.not.ieee_is_finite(scale)) return
   depth=max(0,ceiling(log(max(1.0_dp,2.0_dp*scale))/log(2.0_dp)))
   if(present(forceDepth)) depth=max(depth,forceDepth)
   if(depth>min(25,maxDepth)) return
   allocate(x(n,n),x2(n,n),even(n,n),odd(n,n),tmp(n,n),chol(n,n),rhs(n,n))
   x=a; call addDiagonal(x,-mu)
   x=x*(beta/2.0_dp**depth)
   call mctc_gemm(x,x,x2)
   coeff(0)=1.0_dp
   do k=1,6; coeff(k)=coeff(k-1)*real(7-k,dp)/real((13-k)*k,dp); enddo
   even=coeff(6)*x2; call addDiagonal(even,coeff(4))
   call mctc_gemm(even,x2,tmp)
   even=tmp; call addDiagonal(even,coeff(2))
   call mctc_gemm(even,x2,tmp)
   even=tmp; call addDiagonal(even,coeff(0))
   odd=coeff(5)*x2; call addDiagonal(odd,coeff(3))
   call mctc_gemm(odd,x2,tmp)
   odd=tmp; call addDiagonal(odd,coeff(1))
   call mctc_gemm(x,odd,tmp); odd=tmp
   chol=even; call choleskyLog(chol,evenLog,info)
   if(info/=0) return
   rhs=odd
   call lapack_potrs('u',n,n,chol,n,rhs,n,info)
   if(info/=0) return
   d=-0.5_dp*rhs; call addDiagonal(d,0.5_dp)
   chol=even+odd; call choleskyLog(chol,plusLog,info)
   if(info/=0) return
   logPartition=real(n,dp)*log(2.0_dp)+evenLog-plusLog
   do k=1,depth
      call mctc_gemm(d,d,rhs)
      chol=2.0_dp*(rhs-d); call addDiagonal(chol,1.0_dp)
      call choleskyLog(chol,denominatorLog,info)
      if(info/=0) return
      call lapack_potrs('u',n,n,chol,n,rhs,n,info)
      if(info/=0) return
      d=0.5_dp*(rhs+transpose(rhs))
      logPartition=2.0_dp*logPartition+denominatorLog
      if(.not.all(ieee_is_finite(d))) return
   enddo
   ok=ieee_is_finite(logPartition)
end subroutine

subroutine countBelow(a,cut,count,ok)
   real(dp), intent(in) :: a(:,:),cut
   integer, intent(out) :: count
   logical, intent(out) :: ok
   real(dp), allocatable :: factor(:,:),work(:)
   integer, allocatable :: piv(:)
   real(dp) :: lo,hi,aa,bb,cc,det,scale,tol
   integer :: n,k,info
   ok=.false.; count=0; n=size(a,1)
   call bounds(a,lo,hi)
   if(cut<lo) then
      ok=.true.; return
   else if(cut>hi) then
      count=n; ok=.true.; return
   endif
   allocate(factor(n,n),work(max(1,64*n)),piv(n))
   factor=a; call addDiagonal(factor,-cut)
   call lapack_sytrf('l',n,factor,n,piv,work,size(work),info)
   if(info/=0) return
   tol=128.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lo-cut),abs(hi-cut))
   k=1
   do while(k<=n)
      if(piv(k)>0) then
         if(abs(factor(k,k))<tol) return
         if(factor(k,k)<0.0_dp) count=count+1
         k=k+1
      else
         if(k==n) return
         aa=factor(k,k); bb=factor(k+1,k); cc=factor(k+1,k+1)
         scale=max(abs(aa),abs(bb),abs(cc))
         if(scale<tol) return
         det=(aa/scale)*(cc/scale)-(bb/scale)**2
         if(abs(det)<tol/scale) return
         if(det<0.0_dp) then
            count=count+1
         else if(aa+cc<0.0_dp) then
            count=count+2
         endif
         k=k+2
      endif
   enddo
   ok=.true.
end subroutine

subroutine choleskyLog(a,logdet,info)
   real(dp), intent(inout) :: a(:,:)
   real(dp), intent(out) :: logdet
   integer, intent(out) :: info
   integer :: n,i
   n=size(a,1); logdet=0.0_dp
   call lapack_potrf('u',n,a,n,info)
   if(info/=0) return
   do i=1,n; logdet=logdet+2.0_dp*log(a(i,i)); enddo
end subroutine

subroutine bounds(a,lo,hi)
   real(dp), intent(in) :: a(:,:)
   real(dp), intent(out) :: lo,hi
   real(dp) :: radius
   integer :: i
   lo=huge(1.0_dp); hi=-huge(1.0_dp)
   do i=1,size(a,1)
      radius=sum(abs(a(i,:)))-abs(a(i,i))
      lo=min(lo,a(i,i)-radius); hi=max(hi,a(i,i)+radius)
   enddo
   radius=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lo),abs(hi))
   lo=lo-radius; hi=hi+radius
end subroutine

pure real(dp) function trace(a) result(t)
   real(dp), intent(in) :: a(:,:)
   integer :: i
   t=0.0_dp
   do i=1,size(a,1); t=t+a(i,i); enddo
end function

pure subroutine addDiagonal(a,shift)
   real(dp), intent(inout) :: a(:,:)
   real(dp), intent(in) :: shift
   integer :: i
   do i=1,size(a,1); a(i,i)=a(i,i)+shift; enddo
end subroutine

pure subroutine completeUpper(a)
   real(dp), intent(inout) :: a(:,:)
   integer :: i,j
   do j=1,size(a,1)
      do i=j+1,size(a,1); a(i,j)=a(j,i); enddo
   enddo
end subroutine
end module
