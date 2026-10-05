program test_fermi_operator
   use xtb_mctc_accuracy, only : dp
   use xtb_mctc_constants, only : kB
   use xtb_mctc_convert, only : autoev
   use xtb_mctc_blas_level3, only : blas_trsm
   use xtb_gfn1_fermi_operator, only : TGFN1FermiOperator,purify_density,fermi_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
   implicit none
   integer, parameter :: n=12
   type(TGFN1FermiOperator) :: solver
   real(dp) :: s(n,n),l(n,n),q(n,n),c(n,n),sc(n,n),h(n,n),e(n),p(n,n),w(n,n)
   real(dp) :: reference(n,n),rw(n,n),f(n),fa(n),fb(n),mu,mub,ga,gb,countError,temp,sa,sb
   real(dp) :: maximumP,maximumW,part,check(n,n),identity(n,n)
   integer :: i,j,k,steps,spin,na,nb
   logical :: ok
   l=0.0_dp; q=0.0_dp; identity=0.0_dp
   do i=1,n
      l(i,i)=1.2_dp+0.01_dp*i; q(i,i)=1.0_dp; identity(i,i)=1.0_dp
      e(i)=0.08_dp*(real(i,dp)-6.5_dp)
      do j=1,i-1; l(i,j)=0.02_dp*sin(real(i+j,dp)); enddo
   enddo
   do i=1,n/2
      j=i+n/2
      q(i,i)=cos(0.3_dp); q(j,j)=cos(0.3_dp)
      q(j,i)=sin(0.3_dp); q(i,j)=-sin(0.3_dp)
   enddo
   s=matmul(l,transpose(l)); c=q
   call blas_trsm('l','l','t','n',n,n,1.0_dp,l,n,c,n)
   sc=matmul(s,c); h=sc
   do i=1,n; h(:,i)=h(:,i)*e(i); enddo
   h=matmul(h,transpose(sc))
   call solver%initialize(s,ok)
   call require(ok,'SPD metric initialized')
   do spin=0,2,2
      na=(n+spin)/2; nb=(n-spin)/2
      do k=1,4
         select case(k)
         case(1); temp=0.0_dp
         case(2); temp=50.0_dp
         case(3); temp=300.0_dp
         case(4); temp=1000.0_dp
         end select
         mu=0.5_dp*(e(na)+e(na+1)); mub=0.5_dp*(e(nb)+e(nb+1))
         call exactOccupations(na,temp,mu,fa,sa)
         call exactOccupations(nb,temp,mub,fb,sb)
         f=fa+fb
         reference=c; rw=c
         do i=1,n
            reference(:,i)=reference(:,i)*f(i)
            rw(:,i)=rw(:,i)*(f(i)*e(i))
         enddo
         reference=matmul(reference,transpose(c)); rw=matmul(rw,transpose(c))
         call solver%solve(h,n,spin,temp,mu,mub,p,w,ga,gb,ok,steps,countError)
         print *, 'T=',temp,' spin=',spin,' success=',ok,' expansion steps=',steps
         call require(ok,'density solve convergence')
         maximumP=maxval(abs(p-reference)); maximumW=maxval(abs(w-rw))
         print *, 'max dP=',maximumP,' max dW=',maximumW,' dS=',abs(ga+gb-sa-sb)
         call require(maximumP<1.0e-10_dp,'generalized AO density parity')
         call require(maximumW<1.0e-10_dp,'energy-weighted density parity')
         call require(abs(ga+gb-sa-sb)<1.0e-11_dp,'fermismear clipped entropy parity')
         check=matmul(h,p)-matmul(s,w)
         call require(maxval(abs(check))<1.0e-10_dp,'HP=SW identity')
         call require(abs(sum(s*p)-real(n,dp))<1.0e-10_dp,'AO electron count')
      enddo
   enddo
   call solver%solve(h,n,0,300.0_dp,mu,mub,p,w,ga,gb,ok,steps,countError,1)
   call require(.not.ok,'bounded expansion rejection')
   check=0.0_dp
   call purify_density(check,n/2,p,ok,steps,100)
   call require(.not.ok,'zero-temperature degenerate boundary rejected')
   call solver%initialize(identity,ok)
   mu=0.0_dp; mub=0.0_dp
   call solver%solve(check,n,0,300.0_dp,mu,mub,p,w,ga,gb,ok,steps,countError)
   call require(ok.and.maxval(abs(p-identity))<1.0e-12_dp,'fractional degenerate finite-temperature density')
   call require(abs(ga+gb+2.0_dp*n*log(2.0_dp)*kB*300.0_dp)<1.0e-11_dp,'degenerate entropy')
   s(1,1)=-1.0_dp
   call solver%initialize(s,ok)
   call require(.not.ok,'indefinite metric rejected')
   call solver%initialize(identity,ok)
   h(1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
   call solver%solve(h,n,0,300.0_dp,mu,mub,p,w,ga,gb,ok,steps,countError)
   call require(.not.ok,'nonfinite Hamiltonian rejected')
contains
   subroutine exactOccupations(nelec,t,chemical,focc,entropy)
      integer, intent(in) :: nelec
      real(dp), intent(in) :: t
      real(dp), intent(inout) :: chemical
      real(dp), intent(out) :: focc(n),entropy
      real(dp) :: beta,total,slope,change,x
      integer :: iter,i
      focc=0.0_dp; entropy=0.0_dp
      if(t<=0.1_dp) then
         focc(1:nelec)=1.0_dp; return
      endif
      beta=1.0_dp/(kB*autoev*t)
      do iter=1,100
         do i=1,n
            x=(e(i)-chemical)*beta
            if(x<50.0_dp) then
               focc(i)=1.0_dp/(1.0_dp+exp(x))
            else
               focc(i)=0.0_dp
            endif
         enddo
         total=sum(focc); slope=beta*sum(focc*(1.0_dp-focc))
         if(abs(total-real(nelec,dp))<1.0e-13_dp) exit
         chemical=chemical+(real(nelec,dp)-total)/slope
      enddo
      do i=1,n
         if(focc(i)>1.0e-9_dp.and.1.0_dp-focc(i)>1.0e-9_dp) then
            entropy=entropy+focc(i)*log(focc(i))+(1.0_dp-focc(i))*log(1.0_dp-focc(i))
         endif
      enddo
      entropy=entropy*kB*t
   end subroutine
   subroutine require(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if(.not.condition) then
         print *, 'FAIL ',label
         error stop 1
      endif
      print *, 'PASS ',label
   end subroutine
end program
