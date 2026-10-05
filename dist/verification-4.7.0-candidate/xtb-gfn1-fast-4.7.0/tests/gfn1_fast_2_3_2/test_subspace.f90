program test_subspace
   use xtb_mctc_accuracy, only : dp
   use xtb_gfn1_subspace, only : TGFN1Subspace
   use xtb_mctc_lapack_stdeigval, only : lapack_syevd
   use xtb_mctc_lapack_trf, only : lapack_potrf
   use xtb_mctc_lapack_gst, only : lapack_sygst
   use xtb_mctc_blas_level3, only : blas_trsm
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   integer, parameter :: n=64,nroots=12
   type(TGFN1Subspace) :: solver
   real(dp) :: identity(n,n),old(n),a(n,n),exact(n,n),w(n),c(n,n),eval(n),res
   real(dp) :: work(1+6*n+2*n*n),error,metric(n,n),l(n,n),q(n,n),reference(n,n)
   real(dp) :: shift(n),zero(n),sc(n,n),h(n,n),cholesky(n,n),density(n,n),r(n,nroots)
   integer :: iwork(3+5*n),i,j,info,iterations,restarts
   logical :: ok
   identity=0.0_dp; zero=0.0_dp
   do i=1,n
      identity(i,i)=1.0_dp; old(i)=real(i,dp)/8.0_dp
   enddo
   do j=1,n
      do i=1,n
         a(i,j)=1.0e-4_dp*sin(real(i+j,dp))
      enddo
      a(j,j)=a(j,j)+old(j)
   enddo
   exact=a
   call diagonalize(exact,w)
   call solver%remember(identity,old)
   call solver%solve(a,nroots,4,30,1.0e-12_dp,eval,c,ok,iterations,restarts,res)
   call require(ok,'standard symmetric problem convergence')
   call require(maxval(abs(eval(1:nroots)-w(1:nroots)))<1.0e-10_dp,'standard eigenvalue parity')
   ! A deliberately unattainable tolerance forces correction expansion and
   ! thick restart before the bounded solver rejects the candidate.
   do j=1,n
      do i=1,n
         a(i,j)=0.02_dp*sin(real(i*j,dp))
      enddo
      if(j>nroots+4) old(j)=old(j)+10.0_dp
      a(j,j)=a(j,j)+old(j)
   enddo
   call solver%remember(identity,old)
   call solver%solve(a,nroots,4,8,1.0e-30_dp,eval,c,ok,iterations,restarts,res)
   call require(.not.ok.and.restarts>0,'bounded block corrections and thick restart')
   do i=1,n; old(i)=real(i,dp)/8.0_dp; enddo
   a=0.0_dp
   do i=1,n; a(i,i)=old(i); enddo
   a(17,18)=3.0_dp; a(18,17)=3.0_dp
   call solver%remember(identity,old)
   call solver%solve(a,nroots,4,8,1.0e-12_dp,eval,c,ok,iterations,restarts,res)
   call require(.not.ok,'disconnected virtual intruder rejected')
   call solver%clear()
   call solver%solve(a,nroots,4,30,1.0e-12_dp,eval,c,ok,iterations,restarts,res)
   call require(.not.ok,'missing reference rejected')

   ! Non-diagonal SPD metric and rotated MOs independently verify the
   ! exact GFN1 shift projection and AO generalized residual/density.
   l=0.0_dp; q=identity
   do i=1,n
      l(i,i)=1.2_dp+0.01_dp*i
      do j=1,i-1
         l(i,j)=0.01_dp*sin(real(i+2*j,dp))
      enddo
   enddo
   do i=1,n/2
      j=i+n/2
      q(i,i)=cos(0.3_dp); q(j,j)=cos(0.3_dp)
      q(j,i)=sin(0.3_dp); q(i,j)=-sin(0.3_dp)
   enddo
   metric=matmul(l,transpose(l)); reference=q
   call blas_trsm('l','l','t','n',n,n,1.0_dp,l,n,reference,n)
   sc=matmul(metric,reference)
   h=sc
   do j=1,n; h(:,j)=h(:,j)*old(j); enddo
   h=matmul(h,transpose(sc))
   do i=1,n; shift(i)=1.0e-4_dp*sin(real(i,dp)); enddo
   do j=1,n
      do i=1,n
         h(i,j)=h(i,j)-0.5_dp*metric(i,j)*(shift(i)+shift(j))
      enddo
   enddo
   exact=h; cholesky=metric
   call lapack_potrf('u',n,cholesky,n,info)
   call require(info==0,'reference Cholesky')
   call lapack_sygst(1,'u',n,exact,n,cholesky,n,info)
   call require(info==0,'reference generalized reduction')
   call diagonalize(exact,w)
   call blas_trsm('l','u','n','n',n,n,1.0_dp,cholesky,n,exact,n)
   call solver%remember_gfn1(reference,old,metric,zero)
   call solver%solve_gfn1(shift,nroots,4,30,1.0e-12_dp,eval,c,ok,iterations,restarts,res)
   call require(ok,'GFN1 shifted generalized convergence')
   call require(maxval(abs(eval(1:nroots)-w(1:nroots)))<1.0e-10_dp,'GFN1 eigenvalue parity')
   r=matmul(h,c(:,1:nroots))
   sc(:,1:nroots)=matmul(metric,c(:,1:nroots))
   do j=1,nroots; r(:,j)=r(:,j)-eval(j)*sc(:,j); enddo
   call require(maxval(abs(r))<1.0e-11_dp,'independent AO generalized residual')
   density=matmul(c(:,1:nroots),transpose(c(:,1:nroots))) &
      & -matmul(exact(:,1:nroots),transpose(exact(:,1:nroots)))
   call require(maxval(abs(density))<1.0e-10_dp,'occupied density parity')
   call require(all(c(:,nroots+1:)==0.0_dp),'unused vectors cleared')
   call solver%solve_gfn1(shift,nroots,4,0,1.0e-12_dp,eval,c,ok,iterations,restarts,res)
   call require(.not.ok,'iteration limit rejected')
   shift(1)=ieee_value(0.0_dp,ieee_quiet_nan)
   call solver%solve_gfn1(shift,nroots,4,30,1.0e-12_dp,eval,c,ok,iterations,restarts,res)
   call require(.not.ok,'nonfinite shift rejected')
   call solver%remember_gfn1(reference,old,metric,zero(1:n-1))
   call solver%solve_gfn1(zero,nroots,4,30,1.0e-12_dp,eval,c,ok,iterations,restarts,res)
   call require(.not.ok,'mismatched reference shape rejected')
contains
   subroutine require(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if(.not.condition) then
         print *, 'FAIL ',label
         error stop 1
      endif
      print *, 'PASS ',label
   end subroutine
   subroutine diagonalize(matrix,eigenvalues)
      real(dp), intent(inout) :: matrix(n,n)
      real(dp), intent(out) :: eigenvalues(n)
      call lapack_syevd('v','u',n,matrix,n,eigenvalues,work,size(work),iwork,size(iwork),info)
      call require(info==0,'reference SYEVD')
   end subroutine
end program
