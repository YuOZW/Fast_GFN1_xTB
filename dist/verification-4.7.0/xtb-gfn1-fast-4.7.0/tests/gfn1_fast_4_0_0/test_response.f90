program test_response
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : kB
   use xtb_mctc_blas_level3, only : blas_trsm
   use xtb_mctc_lapack_geneigval, only : lapack_sygv
   use xtb_xtb_response, only : TGFN1ElectronicResponse,TGFN1CoupledResponse,density_response_spin, &
      & density_response_closed_shell,shell_charge_response
   use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
   implicit none
   integer, parameter :: n=6
   type(TGFN1ElectronicResponse) :: response,uninitialized
   real(wp) :: c(n,n),e(n),s(n,n),h(n,n),l(n,n),q(n,n),sc(n,n),dh(n,n),ds(n,n)
   real(wp) :: f(n,2),p(n,n),w(n,n),dp(n,n),dw(n,n),pm(n,n),pp(n,n),wm(n,n),ww(n,n)
   real(wp) :: dmu(2),old(n,n),temp,step,ep,ew,previous,identity(n,n),dq(3),fdq(3),counts(2)
   real(wp) :: shellP(3),shellM(3),zero(n,n),reference(n,n),r(n,n),bad(n,n)
   integer :: i,j,k,spin,number(2),ao2sh(n)
   logical :: ok
   l=0.0_wp; q=0.0_wp; identity=0.0_wp; zero=0.0_wp
   do i=1,n
      l(i,i)=1.1_wp+0.02_wp*i; q(i,i)=1.0_wp; identity(i,i)=1.0_wp
      e(i)=0.002_wp*(real(i,wp)-3.4_wp)
      do j=1,i-1; l(i,j)=0.025_wp*sin(real(i+j,wp)); enddo
      do j=1,n
         dh(i,j)=0.001_wp*cos(real(i+j,wp))
         ds(i,j)=0.02_wp*sin(real(i+j,wp))
      enddo
   enddo
   do i=1,n/2
      j=i+n/2
      q(i,i)=cos(0.4_wp); q(j,j)=cos(0.4_wp)
      q(j,i)=sin(0.4_wp); q(i,j)=-sin(0.4_wp)
   enddo
   s=matmul(l,transpose(l)); c=q
   call blas_trsm('l','l','t','n',n,n,1.0_wp,l,n,c,n)
   sc=matmul(s,c); h=sc
   do i=1,n; h(:,i)=h(:,i)*e(i); enddo
   h=matmul(h,transpose(sc)); ao2sh=[1,1,2,2,3,3]

   do spin=0,1
      number=[3+spin,3-spin]
      do k=1,4
         select case(k)
         case(1); temp=0.0_wp
         case(2); temp=300.0_wp
         case(3); temp=3000.0_wp
         case(4); temp=30000.0_wp
         end select
         call occupations(e,number,temp,f)
         call makeDensity(c,e,f,p,w)
         call response%initialize(c,e,f,temp,ok)
         call require(ok,'response cache initialized')
         call response%apply(dh,ds,dp,dw,ok,dmu)
         call require(ok,'P/W response evaluated')
         reference=matmul(transpose(c),matmul(ds,c))
         call response%apply(dh,ds,pp,ww,ok,counts,cachedMetric=reference)
         call require(ok.and.maxval(abs(pp-dp))<1.0e-12_wp.and.maxval(abs(ww-dw))<1.0e-12_wp, &
            & 'cached overlap metric preserves P/W response')
         call require(maxval(abs(counts-dmu))<1.0e-12_wp,'cached canonical chemical-potential response')
         call require(abs(sum(s*dp)+sum(ds*p))<1.0e-11_wp,'canonical AO electron response')
         reference=matmul(dh,p)+matmul(h,dp)-matmul(ds,w)-matmul(s,dw)
         call require(maxval(abs(reference))<1.0e-11_wp,'d(HP=SW) identity')
         if(spin==0.and.temp==0.0_wp) then
            call density_response_closed_shell(c,e,3,dh,ds,old,ok)
            call require(ok.and.maxval(abs(old-dp))<1.0e-12_wp,'closed-shell kernel compatibility')
         endif
         previous=huge(1.0_wp)
         do j=1,2
            step=2.0e-4_wp/real(2**(j-1),wp)
            call displaced(h,s,number,temp,pm,wm,-step)
            call displaced(h,s,number,temp,pp,ww,step)
            ep=maxval(abs((pp-pm)/(2*step)-dp))
            ew=maxval(abs((ww-wm)/(2*step)-dw))
            print *, 'spin=',spin,' T=',temp,' step=',step,' dP error=',ep,' dW error=',ew
            call require(ep<2.0e-8_wp.and.ew<2.0e-10_wp,'generalized eigensolve finite difference')
            if(j==2.and.previous>1.0e-9_wp) then
               call require(ep<0.4_wp*previous,'second-order finite-difference convergence')
            endif
            previous=ep
         enddo
         call shell_charge_response(ao2sh,s,p,ds,dp,dq)
         call shell_charge_response(ao2sh,s+step*ds,pp,zero,pp,shellP)
         call shell_charge_response(ao2sh,s-step*ds,pm,zero,pm,shellM)
         fdq=(shellP-shellM)/(2*step)
         call require(maxval(abs(fdq-dq))<3.0e-7_wp,'shell Mulliken response finite difference')
         call require(abs(sum(dq))<1.0e-11_wp,'charge response conservation')
         if(k==1.or.k==3) call testCoupled()
      enddo
   enddo

   ! Exactly degenerate fractional states require no denominator shift and
   ! are invariant under rotations of their reference eigenvectors.
   e=0.0_wp; c=identity; s=identity; h=zero; temp=300.0_wp
   f=0.5_wp; number=[3,3]
   call density_response_spin(c,e,f,temp,dh,ds,dp,dw,dmu,ok)
   call require(ok,'finite-temperature exact degeneracy accepted')
   call density_response_spin(q,e,f,temp,dh,ds,pp,ww,counts,ok)
   call require(ok.and.maxval(abs(dp-pp))<1.0e-11_wp.and.maxval(abs(dw-ww))<1.0e-12_wp, &
      & 'degenerate gauge invariance')
   step=1.0e-4_wp
   call displaced(h,s,number,temp,pm,wm,-step)
   call displaced(h,s,number,temp,pp,ww,step)
   call require(maxval(abs((pp-pm)/(2*step)-dp))<1.0e-7_wp,'degenerate response finite difference')

   ! Near degeneracy exercises the stable divided difference rather than
   ! a subtraction of two almost equal occupations.
   e=[-0.01_wp,-0.003_wp,-1.0e-13_wp,1.0e-13_wp,0.003_wp,0.01_wp]
   call occupations(e,number,temp,f)
   h=zero
   do i=1,n; h(i,i)=e(i); enddo
   call density_response_spin(identity,e,f,temp,dh,ds,dp,dw,dmu,ok)
   call require(ok,'near-degenerate finite-temperature response')
   step=1.0e-4_wp
   call displaced(h,s,number,temp,pm,wm,-step)
   call displaced(h,s,number,temp,pp,ww,step)
   call require(maxval(abs((pp-pm)/(2*step)-dp))<1.0e-7_wp,'near-degenerate finite difference')

   e=0.0_wp; f=0.0_wp; f(1:3,:)=1.0_wp
   call response%initialize(identity,e,f,0.0_wp,ok)
   call require(.not.ok,'zero-temperature occupied/virtual degeneracy rejected')
   call response%apply(dh,ds,dp,dw,ok)
   call require(.not.ok.and.maxval(abs(dp))+maxval(abs(dw))==0.0_wp,'failed reinitialization invalidates cache')
   call uninitialized%apply(dh,ds,dp,dw,ok)
   call require(.not.ok,'missing reference rejected')
   f=0.5_wp
   call response%initialize(identity,e,f,0.0_wp,ok)
   call require(.not.ok,'zero-temperature fractional filling rejected')
   f(2,1)=0.8_wp
   call response%initialize(identity,e,f,300.0_wp,ok)
   call require(.not.ok,'inconsistent occupations rejected')
   f=0.5_wp
   call response%initialize(identity,e,f,300.0_wp,ok)
   call require(ok,'valid reference restored')
   bad=dh; bad(1,1)=ieee_value(0.0_wp,ieee_quiet_nan)
   call response%apply(bad,ds,dp,dw,ok)
   call require(.not.ok,'nonfinite perturbation rejected')
   bad=dh; bad(1,2)=bad(1,2)+0.01_wp
   call response%apply(bad,ds,dp,dw,ok)
   call require(.not.ok,'nonsymmetric perturbation rejected')
   call response%apply(dh(1:n-1,1:n-1),ds,dp,dw,ok)
   call require(.not.ok,'shape mismatch rejected')
   call response%apply(dh,ds,dp,dw,ok,cachedMetric=bad(1:n-1,1:n-1))
   call require(.not.ok,'cached metric shape mismatch rejected')
   bad=zero;bad(1,1)=ieee_value(0.0_wp,ieee_quiet_nan)
   call response%apply(dh,ds,dp,dw,ok,cachedMetric=bad)
   call require(.not.ok,'nonfinite cached metric rejected')
contains
   subroutine testCoupled()
      type(TGFN1CoupledResponse) :: coupled
      real(wp) :: jacobian(3,3),rcond,cdp(n,n),cdw(n,n),cq(3),qplus(3),qminus(3)
      real(wp) :: pplus(n,n),pminus(n,n),wplus(n,n),wminus(n,n),xstep,err
      real(wp) :: energyPlus,energyMinus,hessianAA,hessianAB,hessianBA,shift(3),dH2(n,n),dS2(n,n)
      real(wp) :: dp2(n,n),dw2(n,n),dq2(3),gradientPlus,gradientMinus
      real(wp) :: susceptibility(3,3),basisShift(3),perturbation(n,n),smallIdentity(3,3),eval(3),work(128)
      integer :: a,b,info,i
      logical :: good
      jacobian=0.0002_wp
      do a=1,3; jacobian(a,a)=0.0008_wp; enddo
      call coupled%initialize(response,ao2sh,s,p,jacobian,good,rcond)
      call require(good.and.rcond>0.1_wp,'coupled shell-response factorization')
      call coupled%apply(dh,ds,cdp,cdw,cq,good)
      call require(good,'coupled P/W/shell-charge response')
      call require(abs(sum(cq))<1.0e-11_wp,'coupled charge conservation')
      xstep=1.0e-4_wp
      call selfConsistent(-xstep,jacobian,pminus,wminus,qminus,energy=energyMinus)
      call selfConsistent(xstep,jacobian,pplus,wplus,qplus,energy=energyPlus)
      err=maxval(abs((pplus-pminus)/(2*xstep)-cdp))
      print *, 'coupled spin=',spin,' T=',temp,' dP error=',err, &
         & ' dW error=',maxval(abs((wplus-wminus)/(2*xstep)-cdw))
      call require(err<2.0e-8_wp,'self-consistent P finite difference')
      call require(maxval(abs((wplus-wminus)/(2*xstep)-cdw))<2.0e-10_wp,'self-consistent W finite difference')
      call require(maxval(abs((qplus-qminus)/(2*xstep)-cq))<3.0e-8_wp,'self-consistent charge finite difference')
      call require(abs((energyPlus-energyMinus)/(2*xstep)-(sum(p*dh)-sum(w*ds)))<2.0e-8_wp, &
         & 'self-consistent free-energy gradient')
      hessianAA=electronicHessian(dh,ds,cdp,cdw,cq,jacobian)
      gradientPlus=electronicGradient(dh,ds,pplus,wplus,qplus,jacobian)
      gradientMinus=electronicGradient(dh,ds,pminus,wminus,qminus,jacobian)
      call require(abs((gradientPlus-gradientMinus)/(2*xstep)-hessianAA)<2.0e-8_wp, &
         & 'electronic Hessian from gradient differentiation')
      do b=1,n
         do a=1,n
            dH2(a,b)=0.001_wp*sin(real(a+b,wp))
            dS2(a,b)=0.015_wp*cos(real(a+b,wp))
         enddo
      enddo
      call coupled%apply(dH2,dS2,dp2,dw2,dq2,good)
      call require(good,'second coordinate response reuses factorization')
      hessianAB=electronicHessian(dh,ds,dp2,dw2,dq2,jacobian)
      hessianBA=electronicHessian(dH2,dS2,cdp,cdw,cq,jacobian)
      print *, 'electronic Hessian mixed symmetry error=',abs(hessianAB-hessianBA)
      call require(abs(hessianAB-hessianBA)<1.0e-11_wp,'unsymmetrized electronic Hessian reciprocity')
      call selfConsistent(-xstep,jacobian,pminus,wminus,qminus,dH2,dS2)
      call selfConsistent(xstep,jacobian,pplus,wplus,qplus,dH2,dS2)
      gradientPlus=electronicGradient(dh,ds,pplus,wplus,qplus,jacobian)
      gradientMinus=electronicGradient(dh,ds,pminus,wminus,qminus,jacobian)
      call require(abs((gradientPlus-gradientMinus)/(2*xstep)-hessianAB)<2.0e-8_wp, &
         & 'mixed electronic Hessian finite difference')
      ! Build a deliberately critical attractive kernel by finding the
      ! nonzero susceptibility eigenvalue independently. The response must
      ! reject this singular linearization without regularizing it.
      if(spin==0.and.temp==0.0_wp) then
         do a=1,3
            basisShift=0.0_wp; basisShift(a)=1.0_wp
            do b=1,n
               do i=1,n
                  perturbation(i,b)=-0.5_wp*s(i,b)*(basisShift(ao2sh(i))+basisShift(ao2sh(b)))
               enddo
            enddo
            call response%apply(perturbation,zero,cdp,cdw,good)
            call require(good,'independent shell susceptibility')
            call shell_charge_response(ao2sh,s,p,zero,cdp,susceptibility(:,a))
         enddo
         smallIdentity=0.0_wp
         do a=1,3; smallIdentity(a,a)=1.0_wp; enddo
         call lapack_sygv(1,'n','u',3,susceptibility,3,smallIdentity,3,eval,work,size(work),info)
         call require(info==0.and.eval(1)<-1.0e-6_wp,'critical-kernel construction')
         jacobian=0.0_wp
         do a=1,3; jacobian(a,a)=1.0_wp/eval(1); enddo
         call coupled%initialize(response,ao2sh,s,p,jacobian,good,rcond)
         call require(.not.good,'singular coupled SCC response rejected')
         call coupled%apply(dh,ds,cdp,cdw,cq,good)
         call require(.not.good,'failed coupled cache invalidated')
      endif
   end subroutine
   real(wp) function electronicGradient(dha,dsa,pp,ww,charge,jacobian) result(value)
      real(wp), intent(in) :: dha(n,n),dsa(n,n),pp(n,n),ww(n,n),charge(3),jacobian(3,3)
      real(wp) :: shift(3)
      integer :: a,b
      value=sum(pp*dha)-sum(ww*dsa); shift=matmul(jacobian,charge)
      do b=1,n
         do a=1,n
            value=value-0.5_wp*pp(a,b)*dsa(a,b)*(shift(ao2sh(a))+shift(ao2sh(b)))
         enddo
      enddo
   end function
   real(wp) function electronicHessian(dha,dsa,pb,wb,qb,jacobian) result(value)
      real(wp), intent(in) :: dha(n,n),dsa(n,n),pb(n,n),wb(n,n),qb(3),jacobian(3,3)
      real(wp) :: shift(3)
      integer :: a,b
      value=sum(pb*dha)-sum(wb*dsa); shift=matmul(jacobian,qb)
      do b=1,n
         do a=1,n
            value=value-0.5_wp*p(a,b)*dsa(a,b)*(shift(ao2sh(a))+shift(ao2sh(b)))
         enddo
      enddo
   end function
   subroutine selfConsistent(displacement,jacobian,pp,ww,charge,directionH,directionS,energy)
      real(wp), intent(in) :: displacement,jacobian(3,3)
      real(wp), intent(out) :: pp(n,n),ww(n,n),charge(3)
      real(wp), intent(in), optional :: directionH(n,n),directionS(n,n)
      real(wp), intent(out), optional :: energy
      real(wp) :: ss(n,n),hh(n,n),eval(n),cc(n,n),b(n,n),work(256),occ(n,2)
      real(wp) :: referenceCharge(3),newCharge(3),shift(3),h0(n,n),entropy,occupation
      integer :: a,bindex,iteration,info
      call shell_charge_response(ao2sh,s,p,zero,p,referenceCharge)
      ss=s+displacement*ds; h0=h+displacement*dh
      if(present(directionH)) h0=h+displacement*directionH
      if(present(directionS)) ss=s+displacement*directionS
      charge=0.0_wp
      do iteration=1,200
         shift=matmul(jacobian,charge)
         hh=h0
         do bindex=1,n
            do a=1,n
               hh(a,bindex)=hh(a,bindex)-0.5_wp*ss(a,bindex)*(shift(ao2sh(a))+shift(ao2sh(bindex)))
            enddo
         enddo
         cc=hh; b=ss
         call lapack_sygv(1,'v','u',n,cc,n,b,n,eval,work,size(work),info)
         call require(info==0,'self-consistent generalized eigensolve')
         call occupations(eval,number,temp,occ)
         call makeDensity(cc,eval,occ,pp,ww)
         call shell_charge_response(ao2sh,ss,pp,zero,pp,newCharge)
         newCharge=newCharge-referenceCharge
         if(maxval(abs(newCharge-charge))<1.0e-14_wp) exit
         charge=0.5_wp*(newCharge+charge)
      enddo
      call require(iteration<=200,'independent nonlinear SCC converged')
      charge=newCharge
      if(present(energy)) then
         entropy=0.0_wp
         do bindex=1,2
            do a=1,n
               occupation=occ(a,bindex)
               if(occupation>0.0_wp.and.occupation<1.0_wp) entropy=entropy &
                  & +occupation*log(occupation)+(1.0_wp-occupation)*log(1.0_wp-occupation)
            enddo
         enddo
         energy=sum(pp*h0)+0.5_wp*dot_product(charge,matmul(jacobian,charge))+kB*temp*entropy
      endif
   end subroutine
   subroutine displaced(hh,ss,nel,t,pp,ww,perturbation)
      real(wp), intent(in) :: hh(n,n),ss(n,n),t,perturbation
      integer, intent(in) :: nel(2)
      real(wp), intent(out) :: pp(n,n),ww(n,n)
      real(wp) :: cc(n,n),b(n,n),eval(n),work(256),occ(n,2)
      integer :: info
      cc=hh+perturbation*dh; b=ss+perturbation*ds
      call lapack_sygv(1,'v','u',n,cc,n,b,n,eval,work,size(work),info)
      call require(info==0,'independent generalized eigensolve')
      call occupations(eval,nel,t,occ)
      call makeDensity(cc,eval,occ,pp,ww)
   end subroutine
   subroutine makeDensity(cc,eval,occ,pp,ww)
      real(wp), intent(in) :: cc(n,n),eval(n),occ(n,2)
      real(wp), intent(out) :: pp(n,n),ww(n,n)
      real(wp) :: tmp(n,n),total(n)
      integer :: a
      total=occ(:,1)+occ(:,2); tmp=cc
      do a=1,n; tmp(:,a)=tmp(:,a)*total(a); enddo
      pp=matmul(tmp,transpose(cc))
      do a=1,n; tmp(:,a)=tmp(:,a)*eval(a); enddo
      ww=matmul(tmp,transpose(cc))
   end subroutine
   subroutine occupations(eval,nel,t,occ)
      real(wp), intent(in) :: eval(n),t
      integer, intent(in) :: nel(2)
      real(wp), intent(out) :: occ(n,2)
      real(wp) :: left,right,mu,x
      integer :: sp,a,iter
      occ=0.0_wp
      do sp=1,2
         if(t<=0.1_wp) then
            occ(1:nel(sp),sp)=1.0_wp; cycle
         endif
         left=eval(1)-100.0_wp*kB*t; right=eval(n)+100.0_wp*kB*t
         do iter=1,200
            mu=0.5_wp*(left+right)
            do a=1,n
               x=(eval(a)-mu)/(kB*t)
               if(x>0.0_wp) then
                  occ(a,sp)=exp(-x)/(1.0_wp+exp(-x))
               else
                  occ(a,sp)=1.0_wp/(1.0_wp+exp(x))
               endif
            enddo
            if(abs(sum(occ(:,sp))-real(nel(sp),wp))<1.0e-14_wp) exit
            if(sum(occ(:,sp))>real(nel(sp),wp)) then
               right=mu
            else
               left=mu
            endif
         enddo
         call require(abs(sum(occ(:,sp))-real(nel(sp),wp))<1.0e-13_wp,'independent canonical occupations')
      enddo
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
