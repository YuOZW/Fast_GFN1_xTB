! This file is part of xtb.
! SCC subspace reuse with an S-orthonormal GFN1 reference spectrum.
module xtb_gfn1_subspace
   use xtb_mctc_accuracy, only : dp
   use xtb_mctc_blas_level3, only : mctc_gemm, blas_trsm
   use xtb_mctc_lapack_stdeigval, only : lapack_syevd
   use xtb_mctc_lapack_trf, only : lapack_potrf
   use xtb_mctc_lapack_gst, only : lapack_sygst
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: TGFN1Subspace

   type :: TGFN1Subspace
      real(dp), allocatable :: reference(:,:), reference_eval(:), seed(:,:)
      real(dp), allocatable :: metric_reference(:,:), shift_reference(:)
      real(dp), allocatable :: basis(:,:), abasis(:,:), ritz(:,:), aritz(:,:), residual(:,:)
      real(dp), allocatable :: projection(:,:), rotation(:,:), theta(:), correction(:,:), coeff(:,:)
      real(dp), allocatable :: work(:)
      integer, allocatable :: iwork(:)
      integer :: seed_cols = 0
   contains
      procedure :: clear => clearSubspace
      procedure :: remember => rememberSpectrum
      procedure :: solve => solveSubspace
      procedure :: remember_gfn1 => rememberGFN1
      procedure :: solve_gfn1 => solveGFN1
   end type

contains

subroutine clearSubspace(self)
   class(TGFN1Subspace), intent(inout) :: self
   type(TGFN1Subspace) :: empty
   ! Assignment finalizes every allocatable component without shared state.
   select type(self)
   type is(TGFN1Subspace)
      self = empty
   end select
end subroutine

subroutine rememberGFN1(self,vectors,eval,metric,shift)
   class(TGFN1Subspace), intent(inout) :: self
   real(dp), intent(in) :: vectors(:,:),eval(:),metric(:,:),shift(:)
   integer :: n
   n=size(vectors,1)
   call self%clear()
   if(size(vectors,2)/=n.or.size(eval)/=n.or.size(shift)/=n) return
   if(size(metric,1)/=n.or.size(metric,2)/=n) return
   if(.not.all(ieee_is_finite(metric)).or..not.all(ieee_is_finite(shift))) return
   call self%remember(vectors,eval)
   if(.not.allocated(self%reference)) return
   call matrix(self%metric_reference,size(vectors,1),size(vectors,2))
   call mctc_gemm(metric,vectors,self%metric_reference)
   self%shift_reference=shift
end subroutine

subroutine solveGFN1(self,shift,nroots,guard,maxiter,tolerance,eval,vectors,success,iterations,restarts,maxres)
   class(TGFN1Subspace), intent(inout) :: self
   real(dp), intent(in) :: shift(:),tolerance
   integer, intent(in) :: nroots,guard,maxiter
   real(dp), intent(out) :: eval(:),vectors(:,:)
   logical, intent(out) :: success
   integer, intent(out) :: iterations,restarts
   real(dp), intent(out) :: maxres
   real(dp), allocatable :: dc(:,:),projected(:,:)
   integer :: n,i,p
   success=.false.; iterations=0; restarts=0; maxres=huge(1.0_dp)
   if(.not.allocated(self%metric_reference).or..not.allocated(self%reference)) return
   n=size(shift); p=min(n-1,nroots+guard)
   if(nroots<1.or.nroots>=n.or.p<nroots.or.maxiter<1.or.tolerance<=0.0_dp) return
   if(size(self%reference,1)/=n.or.size(vectors,1)/=n) return
   if(size(vectors,2)<nroots.or.size(eval)<nroots) return
   if(.not.all(ieee_is_finite(shift)).or..not.ieee_is_finite(tolerance)) return
   allocate(dc(n,n),projected(n,n))
   do i=1,n
      dc(i,:)=(shift(i)-self%shift_reference(i))*self%reference(i,:)
   enddo
   ! GFN1 H changes exactly as -1/2 (S dD + dD S). In the old
   ! S-orthonormal MO basis its projection needs a single GEMM and no
   ! SYGST or full H*C transform. Screening is already included in S.
   call mctc_gemm(self%metric_reference,dc,projected,transa='t')
   projected=-0.5_dp*(projected+transpose(projected))
   do i=1,n; projected(i,i)=projected(i,i)+self%reference_eval(i); enddo
   call solveGraph(self,projected,nroots,p,maxiter,tolerance,eval,vectors,success,iterations,maxres, &
      & projectedInput=projected)
   if(success) then
      vectors(:,nroots+1:)=0.0_dp
      eval(nroots+1:)=1.0e6_dp
   endif
end subroutine

subroutine rememberSpectrum(self, vectors, eval)
   class(TGFN1Subspace), intent(inout) :: self
   real(dp), intent(in) :: vectors(:,:), eval(:)
   call self%clear()
   if(size(vectors,1)/=size(vectors,2).or.size(eval)/=size(vectors,1)) return
   if(.not.all(ieee_is_finite(vectors)).or..not.all(ieee_is_finite(eval))) return
   self%reference = vectors
   self%reference_eval = eval
   self%seed = vectors
   self%seed_cols = size(vectors,2)
end subroutine

! A is the full symmetric standard problem; all accepted roots are checked.
! The old virtual spectrum preconditions corrections, while Rayleigh--Ritz
! handles rotations inside the occupied/thermal space. Guard Ritz vectors
! are retained at a thick restart and across SCC steps, never occupied.
subroutine solveSubspace(self, a, nroots, guard, maxiter, tolerance, eval, vectors, &
      & success, iterations, restarts, maxres)
   class(TGFN1Subspace), intent(inout) :: self
   real(dp), intent(in) :: a(:,:), tolerance
   integer, intent(in) :: nroots, guard, maxiter
   real(dp), intent(out) :: eval(:), vectors(:,:)
   logical, intent(out) :: success
   integer, intent(out) :: iterations, restarts
   real(dp), intent(out) :: maxres
   integer :: n, cap, m, keep, nc, rank, j, i, k, info, nseed, nvirtual, startvirtual
   real(dp) :: scale, denom, residuals(nroots)
   logical :: separated
   real(dp), allocatable :: block(:,:), eig(:), gram(:,:), tmp(:,:), nextseed(:,:)

   success=.false.; iterations=0; restarts=0; maxres=huge(1.0_dp)
   n=size(a,1)
   if(nroots<1.or.nroots>=n.or.size(a,2)/=n.or.maxiter<1) return
   if(.not.ieee_is_finite(tolerance).or.tolerance<=0.0_dp) return
   if(.not.allocated(self%reference).or..not.allocated(self%seed)) return
   if(size(self%reference,1)/=n.or.size(self%reference,2)/=n) return
   if(size(vectors,1)/=n.or.size(vectors,2)<nroots.or.size(eval)<nroots) return
   if(.not.all(ieee_is_finite(a)).or..not.all(ieee_is_finite(self%seed))) return
   nseed=min(self%seed_cols,min(n-1,nroots+max(1,guard)))
   if(nseed<nroots) return
   call solveGraph(self,a,nroots,nseed,maxiter,tolerance,eval,vectors,success,iterations,maxres,certified=separated)
   if(success) return
   ! Residual convergence alone does not certify the lowest roots: a
   ! disconnected virtual intruder could be invisible to a low-space seed.
   if(.not.separated) return
   iterations=0; maxres=huge(1.0_dp)
   cap=min(n-1,nroots+2*max(1,guard)+64)
   call matrix(self%basis,n,cap); call matrix(self%abasis,n,cap)
   call matrix(self%ritz,n,nroots); call matrix(self%aritz,n,nroots)
   call matrix(self%residual,n,nroots)
   call matrix(self%projection,cap,cap); call matrix(self%rotation,cap,cap)
   call matrix(self%correction,n,32); call matrix(self%coeff,n,nroots)
   if(allocated(self%theta)) then
      if(size(self%theta)<cap) deallocate(self%theta)
   endif
   if(.not.allocated(self%theta)) allocate(self%theta(cap))
   m=nseed
   self%basis(:,1:m)=self%seed(:,1:m)
   allocate(gram(m,m))
   call mctc_gemm(self%basis(:,1:m),self%basis(:,1:m),gram,transa='t')
   do i=1,m; gram(i,i)=gram(i,i)-1.0_dp; enddo
   if(.not.all(ieee_is_finite(gram)).or.maxval(abs(gram))>1.0e-9_dp) return
   deallocate(gram)
   call mctc_gemm(a,self%basis(:,1:m),self%abasis(:,1:m))
   startvirtual=nseed+1; nvirtual=n-startvirtual+1

   do k=1,maxiter
      iterations=k
      call mctc_gemm(self%basis(:,1:m),self%abasis(:,1:m), &
         & self%projection(1:m,1:m),transa='t')
      self%rotation(1:m,1:m)=0.5_dp*(self%projection(1:m,1:m)+transpose(self%projection(1:m,1:m)))
      call diagonalize(self%rotation(1:m,1:m),self%theta(1:m),self%work,self%iwork,info)
      if(info/=0) return
      call mctc_gemm(self%basis(:,1:m),self%rotation(1:m,1:nroots),self%ritz)
      call mctc_gemm(self%abasis(:,1:m),self%rotation(1:m,1:nroots),self%aritz)
      maxres=0.0_dp
      do j=1,nroots
         self%residual(:,j)=self%aritz(:,j)-self%theta(j)*self%ritz(:,j)
         denom=sqrt(sum(self%aritz(:,j)**2))+abs(self%theta(j))*sqrt(sum(self%ritz(:,j)**2))+tiny(1.0_dp)
         residuals(j)=sqrt(sum(self%residual(:,j)**2))/denom
         if(.not.ieee_is_finite(residuals(j)).or..not.ieee_is_finite(self%theta(j))) return
         maxres=max(maxres,residuals(j))
      enddo
      if(maxres<=tolerance) then
         eval(1:nroots)=self%theta(1:nroots)
         vectors(:,1:nroots)=self%ritz
         keep=min(m,nroots+max(1,guard))
         allocate(nextseed(n,keep))
         call mctc_gemm(self%basis(:,1:m),self%rotation(1:m,1:keep),nextseed)
         self%seed=nextseed; self%seed_cols=keep
         success=.true.
         return
      endif
      if(k==maxiter) return
      ! Compress the entire preconditioned residual block, rather than
      ! choosing individual roots. A shared correction can improve many
      ! nearly converged occupied vectors at once.
      call mctc_gemm(self%reference(:,startvirtual:n),self%residual, &
         & self%coeff(1:nvirtual,1:nroots),transa='t')
      do j=1,nroots
         do i=1,nvirtual
            denom=self%reference_eval(startvirtual+i-1)-self%theta(j)
            if(abs(denom)<1.0e-4_dp) denom=sign(1.0e-4_dp,denom)
            self%coeff(i,j)=self%coeff(i,j)/denom
         enddo
      enddo
      allocate(gram(nroots,nroots),eig(nroots))
      call mctc_gemm(self%coeff(1:nvirtual,1:nroots),self%coeff(1:nvirtual,1:nroots),gram,transa='t')
      gram=0.5_dp*(gram+transpose(gram))
      call diagonalize(gram,eig,self%work,self%iwork,info)
      if(info/=0) return
      scale=maxval(abs(eig)); nc=0
      do j=nroots,max(1,nroots-31),-1
         if(eig(j)>max(1.0e-30_dp,1.0e-12_dp*scale)) nc=nc+1
      enddo
      if(nc<1) return
      allocate(tmp(nvirtual,nc))
      call mctc_gemm(self%coeff(1:nvirtual,1:nroots),gram(:,nroots-nc+1:nroots),tmp)
      do j=1,nc
         tmp(:,j)=tmp(:,j)/sqrt(eig(nroots-nc+j))
      enddo
      call mctc_gemm(self%reference(:,startvirtual:n),tmp,self%correction(:,1:nc))
      deallocate(gram,eig,tmp)
      if(m+nc>cap) then
         keep=min(m,nroots+max(1,guard))
         allocate(tmp(n,keep))
         call mctc_gemm(self%basis(:,1:m),self%rotation(1:m,1:keep),tmp)
         self%basis(:,1:keep)=tmp
         deallocate(tmp)
         m=keep; restarts=restarts+1
         call mctc_gemm(a,self%basis(:,1:m),self%abasis(:,1:m))
      endif
      allocate(block(n,nc),gram(nc,nc),eig(nc),tmp(n,nc))
      block=self%correction(:,1:nc)
      do i=1,2
         call mctc_gemm(self%basis(:,1:m),block,self%projection(1:m,1:nc),transa='t')
         call mctc_gemm(self%basis(:,1:m),self%projection(1:m,1:nc),tmp)
         block=block-tmp
      enddo
      call mctc_gemm(block,block,gram,transa='t')
      gram=0.5_dp*(gram+transpose(gram))
      call diagonalize(gram,eig,self%work,self%iwork,info)
      if(info/=0) return
      scale=maxval(abs(eig)); rank=0
      do j=1,nc
         if(eig(j)>max(1.0e-30_dp,1.0e-12_dp*scale)) then
            rank=rank+1
            self%projection(1:nc,rank)=gram(:,j)/sqrt(eig(j))
         endif
      enddo
      rank=min(rank,cap-m)
      if(rank<1) return
      call mctc_gemm(block,self%projection(1:nc,1:rank),self%basis(:,m+1:m+rank))
      call mctc_gemm(a,self%basis(:,m+1:m+rank),self%abasis(:,m+1:m+rank))
      m=m+rank
      deallocate(block,gram,eig,tmp)
   enddo
end subroutine

! A graph of the low-energy invariant subspace gives a much smaller
! correction problem when SCC Hamiltonian changes are small. Solve the
! Riccati equation using diagonal spectral denominators and retain the
! off-diagonal update in its residual. The residual gate below is mandatory;
! unsafe gap/coupling or failed iteration returns to block Davidson/full.
subroutine solveGraph(self,a,nroots,p,maxiter,tolerance,eval,vectors,success,iterations,maxres,projectedInput,certified)
   class(TGFN1Subspace), intent(inout) :: self
   real(dp), intent(in) :: a(:,:),tolerance
   integer, intent(in) :: nroots,p,maxiter
   real(dp), intent(out) :: eval(:),vectors(:,:)
   logical, intent(out) :: success
   integer, intent(out) :: iterations
   real(dp), intent(out) :: maxres
   real(dp), intent(in), optional :: projectedInput(:,:)
   logical, intent(out), optional :: certified
   real(dp), allocatable :: tmp(:,:),projected(:,:),low(:,:),high(:,:),b(:,:),z(:,:),next(:,:)
   real(dp), allocatable :: le(:),he(:),g(:,:),effective(:,:),w(:),lowvec(:,:)
   real(dp), allocatable :: c(:,:),ac(:,:),res(:,:),small(:,:)
   real(dp) :: gap,denom,error
   integer :: n,q,i,j,k,info
   success=.false.; iterations=0; maxres=huge(1.0_dp)
   if(present(certified)) certified=.false.
   n=size(a,1); q=n-p
   if(p<1.or.q<1.or.nroots>p.or.maxiter<1) return
   if(.not.all(ieee_is_finite(a))) return
   allocate(projected(n,n),low(p,p),high(q,q),b(q,p),le(p),he(q))
   if(present(projectedInput)) then
      projected=projectedInput
   else
      allocate(tmp(n,n))
      call mctc_gemm(a,self%reference,tmp)
      call mctc_gemm(self%reference,tmp,projected,transa='t')
      deallocate(tmp)
   endif
   low=projected(1:p,1:p); high=projected(p+1:n,p+1:n)
   low=0.5_dp*(low+transpose(low)); high=0.5_dp*(high+transpose(high))
   ! The previous full MO spectrum is already diagonal. Retain the small
   ! off-diagonal SCC update in the fixed-point residual instead of solving
   ! both diagonal blocks again. Unsafe updates fail the residual gate.
   do i=1,p; le(i)=low(i,i); low(i,i)=0.0_dp; enddo
   do i=1,q; he(i)=high(i,i); high(i,i)=0.0_dp; enddo
   ! Gershgorin bounds certify separation of the entire diagonal blocks,
   ! including off-diagonal updates. Diagonal values alone cannot exclude
   ! an uncomputed virtual eigenvalue moving into the occupied spectrum.
   gap=minval(he-sum(abs(high),dim=2))-maxval(le+sum(abs(low),dim=2))
   if(gap<=1.0e-6_dp) return
   b=projected(p+1:n,1:p)
   if(sqrt(sum(b*b))>0.25_dp*gap) return
   if(present(certified)) certified=.true.
   allocate(z(q,p),next(q,p),tmp(q,p),small(p,p),g(p,p),effective(p,p),w(p))
   z=0.0_dp
   do k=1,maxiter
      iterations=k
      call mctc_gemm(b,z,small,transa='t')
      call mctc_gemm(z,small,next)
      next=next-b
      call mctc_gemm(high,z,tmp)
      next=next-tmp
      call mctc_gemm(z,low,tmp)
      next=next+tmp
      if(.not.all(ieee_is_finite(next))) return
      do j=1,p
         do i=1,q
            next(i,j)=next(i,j)/(he(i)-le(j))
         enddo
      enddo
      error=maxval(abs(next-z))
      z=next
      if(error<=0.01_dp*tolerance) exit
   enddo
   if(error>0.01_dp*tolerance) return
   call mctc_gemm(z,z,g,transa='t')
   do i=1,p; g(i,i)=g(i,i)+1.0_dp; enddo
   call lapack_potrf('u',p,g,p,info)
   if(info/=0) return
   call mctc_gemm(b,z,small,transa='t')
   effective=small+transpose(small)
   effective=effective+projected(1:p,1:p)
   ! D*z = D_off*z + diag(D)*z.
   call mctc_gemm(high,z,next)
   do i=1,q; next(i,:)=next(i,:)+he(i)*z(i,:); enddo
   call mctc_gemm(z,next,small,transa='t')
   effective=effective+small
   call lapack_sygst(1,'u',p,effective,p,g,p,info)
   if(info/=0) return
   call diagonalize(effective,w,self%work,self%iwork,info)
   if(info/=0) return
   allocate(lowvec(n,p),c(n,p),ac(n,nroots),res(n,nroots))
   call blas_trsm('l','u','n','n',p,p,1.0_dp,g,p,effective,p)
   c(1:p,:)=effective
   call mctc_gemm(z,effective,c(p+1:n,:))
   call mctc_gemm(projected,c(:,1:nroots),ac)
   maxres=0.0_dp
   do j=1,nroots
      res(:,j)=ac(:,j)-w(j)*c(:,j)
      denom=sqrt(sum(ac(:,j)**2))+abs(w(j))+tiny(1.0_dp)
      error=sqrt(sum(res(:,j)**2))/denom
      if(.not.ieee_is_finite(error)) return
      maxres=max(maxres,error)
   enddo
   if(.not.ieee_is_finite(maxres).or.maxres>tolerance) return
   call mctc_gemm(c,c,g,transa='t')
   do i=1,p; g(i,i)=g(i,i)-1.0_dp; enddo
   if(.not.all(ieee_is_finite(g)).or.maxval(abs(g))>1.0e-10_dp) return
   eval(1:nroots)=w(1:nroots)
   call mctc_gemm(self%reference,c,lowvec)
   vectors(:,1:nroots)=lowvec(:,1:nroots)
   self%seed=lowvec; self%seed_cols=p
   success=.true.
end subroutine

subroutine matrix(a,n,m)
   real(dp), allocatable, intent(inout) :: a(:,:)
   integer, intent(in) :: n,m
   if(allocated(a)) then
      if(size(a,1)/=n.or.size(a,2)/=m) deallocate(a)
   endif
   if(.not.allocated(a)) allocate(a(n,m))
end subroutine

subroutine diagonalize(a,w,work,iwork,info)
   real(dp), intent(inout) :: a(:,:)
   real(dp), intent(out) :: w(:)
   real(dp), allocatable, intent(inout) :: work(:)
   integer, allocatable, intent(inout) :: iwork(:)
   integer, intent(out) :: info
   integer :: n,lw,liw
   real(dp), allocatable :: dense(:,:)
   n=size(a,1); lw=1+6*n+2*n*n; liw=3+5*n
   if(allocated(work)) then
      if(size(work)<lw) deallocate(work)
   endif
   if(allocated(iwork)) then
      if(size(iwork)<liw) deallocate(iwork)
   endif
   if(.not.allocated(work)) allocate(work(lw))
   if(.not.allocated(iwork)) allocate(iwork(liw))
   if(is_contiguous(a)) then
      call lapack_syevd('v','u',n,a,n,w,work,size(work),iwork,size(iwork),info)
   else
      ! Active Ritz blocks can have a larger parent leading dimension.
      ! Make that packing explicit instead of a hidden LAPACK temporary.
      allocate(dense(n,n)); dense=a
      call lapack_syevd('v','u',n,dense,n,w,work,size(work),iwork,size(iwork),info)
      a=dense
   endif
end subroutine
end module
