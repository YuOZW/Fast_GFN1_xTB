! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> TODO
module xtb_xtb_halogen
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_lin, only : lin
   use xtb_xtb_data
   use xtb_type_molecule, only : TMolecule
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   public :: xbpot, xbMolecularHessian


contains


! Differentiate the same B-X...A correction and nearest-neighbour selection
! used by GFN1 SCF. Branch boundaries have no unique Cartesian Hessian and
! must be handled by the caller's numerical fallback. Contributions are fresh.
subroutine xbMolecularHessian(mol,halData,a,energy,gradient,hessian,ok,nterms,guardRadius)
   type(TMolecule), intent(in) :: mol
   type(THalogenData), intent(in) :: halData
   real(wp), intent(in) :: a
   real(wp), intent(out) :: energy,gradient(:,:),hessian(:,:)
   logical, intent(out) :: ok
   integer, intent(out), optional :: nterms
   real(wp), intent(in), optional :: guardRadius
   real(wp) :: u(3),v(3),xx,yy,zz,distance,nearest,nextNearest,r0,cc,branchTolerance
   real(wp) :: value,df(3),d2f(3,3),dvars(3,9),gg(9),hh(9,9),window
   integer :: n,x,acceptor,neighbour,m,i,j,k,atoms(3),coordinates(9),terms
   energy=0.0_wp; gradient=0.0_wp; hessian=0.0_wp; ok=.false.; terms=0
   if(present(nterms)) nterms=0
   n=mol%n
   if(n<1.or.mol%npbc/=0) return
   if(any(shape(gradient)/=[3,n]).or.any(shape(hessian)/=[3*n,3*n])) return
   if(.not.all(ieee_is_finite(mol%xyz))) return
   if(.not.allocated(halData%bondStrength).or..not.allocated(halData%atomicRad)) return
   if(minval(mol%at)<1.or.maxval(mol%at)>size(halData%bondStrength)) return
   if(maxval(mol%at)>size(halData%atomicRad)) return
   if(.not.all(ieee_is_finite(halData%bondStrength(mol%at)))) return
   if(.not.all(ieee_is_finite(halData%atomicRad(mol%at)))) return
   if(.not.ieee_is_finite(halData%radScale).or..not.ieee_is_finite(halData%dampingPar)) return
   if(.not.ieee_is_finite(a).or.a<=0.0_wp.or.halData%radScale<=0.0_wp) return
   window=0.0_wp
   if(present(guardRadius)) window=guardRadius
   if(.not.ieee_is_finite(window).or.window<0.0_wp) return
   do x=1,n
      if(.not.isHalogen(mol%at(x))) cycle
      cc=halData%bondStrength(mol%at(x))
      if(cc==0.0_wp) cycle
      nearest=huge(1.0_wp); nextNearest=huge(1.0_wp); neighbour=0
      do m=1,n
         if(m==x) cycle
         distance=sum((mol%xyz(:,m)-mol%xyz(:,x))**2)
         if(distance<nearest) then
            nextNearest=nearest; nearest=distance; neighbour=m
         else
            nextNearest=min(nextNearest,distance)
         endif
      enddo
      do acceptor=1,n
         if(.not.isAcceptor(mol%at(acceptor))) cycle
         u=mol%xyz(:,acceptor)-mol%xyz(:,x); xx=sum(u*u)
         branchTolerance=1.0e-10_wp*max(1.0_wp,xx)
         if(abs(xx-400.0_wp)<=branchTolerance) return
         if(abs(sqrt(xx)-20.0_wp)<=window+1.0e-10_wp*20.0_wp) return
         if(xx>=400.0_wp) cycle
         if(neighbour==0.or.xx<=1.0e-16_wp.or.nearest<=1.0e-16_wp) return
         if(nextNearest-nearest<=1.0e-10_wp*max(1.0_wp,nearest)) return
         ! Moving the halogen can change each neighbour distance by window.
         if(sqrt(nextNearest)-sqrt(nearest)<=2.0_wp*window) return
         ! If the acceptor itself is the nearest atom, the angle factor is
         ! identically zero in a neighbourhood of this reference geometry.
         if(neighbour==acceptor) cycle
         v=mol%xyz(:,neighbour)-mol%xyz(:,x); yy=sum(v*v); zz=dot_product(u,v)
         r0=halData%radScale*(halData%atomicRad(mol%at(x))+halData%atomicRad(mol%at(acceptor)))
         if(r0<=0.0_wp.or..not.ieee_is_finite(r0)) return
         if(r0>1.0e8_wp*sqrt(xx)) return
         call halogenScalarDerivatives(xx,yy,zz,r0,cc,halData%dampingPar,a,value,df,d2f)
         if(.not.ieee_is_finite(value).or..not.all(ieee_is_finite(df)).or..not.all(ieee_is_finite(d2f))) return
         ! Scalar variables x=|A-X|^2, y=|B-X|^2, z=(A-X).(B-X).
         dvars=0.0_wp
         dvars(1,1:3)=-2.0_wp*u; dvars(1,4:6)=2.0_wp*u
         dvars(2,1:3)=-2.0_wp*v; dvars(2,7:9)=2.0_wp*v
         dvars(3,1:3)=-u-v; dvars(3,4:6)=v; dvars(3,7:9)=u
         gg=matmul(transpose(dvars),df)
         hh=matmul(transpose(dvars),matmul(d2f,dvars))
         do k=1,3
            hh(k,k)=hh(k,k)+2.0_wp*(df(1)+df(2)+df(3))
            hh(k,k+3)=hh(k,k+3)-2.0_wp*df(1)-df(3)
            hh(k+3,k)=hh(k+3,k)-2.0_wp*df(1)-df(3)
            hh(k,k+6)=hh(k,k+6)-2.0_wp*df(2)-df(3)
            hh(k+6,k)=hh(k+6,k)-2.0_wp*df(2)-df(3)
            hh(k+3,k+3)=hh(k+3,k+3)+2.0_wp*df(1)
            hh(k+6,k+6)=hh(k+6,k+6)+2.0_wp*df(2)
            hh(k+3,k+6)=hh(k+3,k+6)+df(3)
            hh(k+6,k+3)=hh(k+6,k+3)+df(3)
         enddo
         atoms=[x,acceptor,neighbour]
         do i=1,3
            gradient(:,atoms(i))=gradient(:,atoms(i))+gg(3*i-2:3*i)
            coordinates(3*i-2:3*i)=[3*atoms(i)-2,3*atoms(i)-1,3*atoms(i)]
         enddo
         do j=1,9
            do i=1,9
               hessian(coordinates(i),coordinates(j))=hessian(coordinates(i),coordinates(j))+hh(i,j)
            enddo
         enddo
         energy=energy+value; terms=terms+1
      enddo
   enddo
   if(present(nterms)) nterms=terms
   ok=all(ieee_is_finite(gradient)).and.all(ieee_is_finite(hessian)).and.ieee_is_finite(energy)
end subroutine xbMolecularHessian

! f(x,y,z)=c*((r0/sqrt(x))^a-d*(r0/sqrt(x))^(a/2)) /
! (1+(r0/sqrt(x))^a) * (1/2-z/(2*sqrt(x*y)))^6.
! Analytic chain rule in three scalar variables; no coordinate differences.
pure subroutine halogenScalarDerivatives(x,y,z,r0,c,d,a,value,df,d2f)
   real(wp), intent(in) :: x,y,z,r0,c,d,a
   real(wp), intent(out) :: value,df(3),d2f(3,3)
   real(wp) :: angular,angleD(3),angleH(3,3),invNorm,angularD(3),angularH(3,3)
   real(wp) :: t,denominator,numerator,radial,radialT,radialTT,tx,txx,radialX,radialXX,alpha
   integer :: i,j
   invNorm=1.0_wp/(sqrt(x)*sqrt(y))
   angular=0.5_wp-0.5_wp*z*invNorm
   angleD=[0.25_wp*z*invNorm/x,0.25_wp*z*invNorm/y,-0.5_wp*invNorm]
   angleH=0.0_wp
   angleH(1,1)=-0.375_wp*z*invNorm/(x*x)
   angleH(2,2)=-0.375_wp*z*invNorm/(y*y)
   angleH(1,2)=-0.125_wp*z*invNorm/(x*y); angleH(2,1)=angleH(1,2)
   angleH(1,3)=0.25_wp*invNorm/x; angleH(3,1)=angleH(1,3)
   angleH(2,3)=0.25_wp*invNorm/y; angleH(3,2)=angleH(2,3)
   angularD=6.0_wp*angular**5*angleD
   do j=1,3
      do i=1,3
         angularH(i,j)=30.0_wp*angular**4*angleD(i)*angleD(j)+6.0_wp*angular**5*angleH(i,j)
      enddo
   enddo
   t=(r0/sqrt(x))**(0.5_wp*a); denominator=1.0_wp+t*t
   numerator=2.0_wp*t-d+d*t*t
   radial=c*(t*t-d*t)/denominator
   radialT=c*numerator/(denominator*denominator)
   radialTT=c*((2.0_wp+2.0_wp*d*t)/(denominator*denominator) &
      & -4.0_wp*t*numerator/(denominator*denominator*denominator))
   alpha=0.25_wp*a; tx=-alpha*t/x; txx=alpha*(alpha+1.0_wp)*t/(x*x)
   radialX=radialT*tx; radialXX=radialTT*tx*tx+radialT*txx
   value=radial*angular**6
   df=radial*angularD; df(1)=df(1)+radialX*angular**6
   d2f=radial*angularH
   d2f(1,:)=d2f(1,:)+radialX*angularD
   d2f(:,1)=d2f(:,1)+radialX*angularD
   d2f(1,1)=d2f(1,1)+radialXX*angular**6
end subroutine halogenScalarDerivatives

pure elemental logical function isHalogen(at) result(found)
   integer, intent(in) :: at
   found=at==17.or.at==35.or.at==53.or.at==85
end function

pure elemental logical function isAcceptor(at) result(found)
   integer, intent(in) :: at
   found=at==7.or.at==8.or.at==15.or.at==16
end function


subroutine xbpot(halData,n,at,xyz,xblist,nxb,a,exb,g)
   type(THalogenData), intent(in) :: halData
   integer, intent(in) :: n
   integer, intent(in) :: at(:)
   integer, intent(in) :: nxb
   integer, intent(in) :: xblist(:,:)
   real(wp), intent(in) :: xyz(:,:)
   real(wp), intent(inout) :: g(:,:)
   real(wp), intent(inout) :: exb
   real(wp), intent(in) :: a

   integer :: m,k,AA,B,X,ati,atj
   real(wp) :: cc,r0ax,t13,t14,t16
   real(wp) :: d2ax,rax,term,aterm,xy,d2bx,d2ab,alp,lj2
   real(wp) :: er,el,step,dxa(3),dxb(3),dba(3),dcosterm
   real(wp) :: dtermlj,termlj,prefactor,numerator,denominator,rbx
   alp=6.0_wp
   lj2=0.50_wp*a

   exb = 0.0_wp
   if(nxb.lt.1) return

   ! B-X...A
   !$omp parallel do schedule(runtime) default(none) reduction(+:exb) &
   !$omp shared(nxb, xblist, at, halData, xyz, lj2, alp, a) &
   !$omp private(X, AA, B, ati, atj, cc, r0ax, dxa, dxb, dba, d2ax, &
   !$omp& d2bx, d2ab, rax, XY, aterm, t13, t14, TERM)
   do k=1,nxb
      X =xblist(1,k)
      AA=xblist(2,k)
      B=xblist(3,k)
      ati=at(X)
      atj=at(AA)
      cc=halData%bondStrength(ati)
      r0ax=halData%radScale*(halData%atomicRad(ati)+halData%atomicRad(atj))
      dxa=xyz(:,AA)-xyz(:,X)   ! acceptor - halogen
      dxb=xyz(:, B)-xyz(:,X)   ! neighbor - halogen
      dba=xyz(:,AA)-xyz(:, B)  ! acceptor - neighbor
      d2ax=sum(dxa*dxa)
      d2bx=sum(dxb*dxb)
      d2ab=sum(dba*dba)
      rax=sqrt(d2ax)
      ! angle part. term = cos angle B-X-A
      XY = SQRT(D2BX*D2AX)
      TERM = (D2BX+D2AX-D2AB) / XY
      aterm = (0.5_wp-0.25_wp*term)**alp
      t13 = r0ax/rax
      t14 = t13**a
      exb = exb +  aterm*cc*(t14-halData%dampingPar*t13**lj2) / (1.0_wp+t14)
   enddo

   ! analytic gradient
   !$omp parallel do schedule(runtime) default(none) reduction(+:g) &
   !$omp shared(nxb, xblist, at, halData, xyz, lj2, alp) &
   !$omp private(X, AA, B, ati, atj, cc, r0ax, dxa, dxb, dba, d2ax, &
   !$omp& d2bx, d2ab, rax, XY, aterm, numerator, denominator, termLJ, &
   !$omp& dtermlj, prefactor, dcosterm, t13, t14, rbx, TERM)
   do k=1,nxb
      X =xblist(1,k)
      AA=xblist(2,k)
      B=xblist(3,k)
      ati=at(X)
      atj=at(AA)
      cc=halData%bondStrength(ati)
      r0ax=halData%radScale*(halData%atomicRad(ati)+halData%atomicRad(atj))

      dxa=xyz(:,AA)-xyz(:,X)   ! acceptor - halogen
      dxb=xyz(:, B)-xyz(:,X)   ! neighbor - halogen
      dba=xyz(:,AA)-xyz(:, B)  ! acceptor - neighbor

      d2ax=sum(dxa*dxa)
      d2bx=sum(dxb*dxb)
      d2ab=sum(dba*dba)
      rax=sqrt(d2ax)+1.0e-18_wp
      rbx=sqrt(d2bx)+1.0e-18_wp

      XY = SQRT(D2BX*D2AX)
      TERM = (D2BX+D2AX-D2AB) / XY
      ! now compute angular damping function
      aterm = (0.5_wp-0.25_wp*term)**alp

      ! set up weighted inverted distance and compute the modified Lennard-Jones potential
      t14 = (r0ax/rax)**lj2 ! (rov/r)^lj2 ; lj2 = 6 in GFN1
      numerator = (t14*t14 - halData%dampingPar*t14)
      denominator = (1.0_wp + t14*t14)
      termLJ= numerator/denominator

      ! ----
      ! LJ derivative
      ! denominator part
      dtermlj=2.0_wp*lj2*numerator*t14*t14/(rax*denominator*denominator)
      ! numerator part
      dtermlj=dtermlj+lj2*t14*(halData%dampingPar - 2.0_wp*t14)/(rax*denominator)
      ! scale w/ angular damping term
      dtermlj=dtermlj*aterm*cc/rax
      ! gradient for the acceptor
      g(:,AA)=g(:,AA)+dtermlj*dxa(:)
      ! halogen gradient
      g(:,X)=g(:,X)-dtermlj*dxa(:)
      ! ----
      ! cosine term derivative
      prefactor=-0.250_wp*alp*(0.5_wp-0.25_wp*term)**(alp-1.0_wp)
      prefactor=prefactor*cc*termlj
      ! AX part
      dcosterm=2.0_wp/rbx - term/rax
      dcosterm=dcosterm*prefactor/rax
      ! gradient for the acceptor
      g(:,AA)=g(:,AA)+dcosterm*dxa(:)
      ! halogen gradient
      g(:,X)=g(:,X)-dcosterm*dxa(:)
      ! BX part
      dcosterm=2.0_wp/rax - term/rbx
      dcosterm=dcosterm*prefactor/rbx
      ! gradient for the acceptor
      g(:,B)=g(:,B)+dcosterm*dxb(:)
      ! halogen gradient
      g(:,X)=g(:,X)-dcosterm*dxb(:)
      ! AB part
      t13=2.0_wp*prefactor/xy
      ! acceptor
      g(:,AA)=g(:,AA)-t13*dba(:)
      ! neighbor
      g(:,B)=g(:,B)+t13*dba(:)

   enddo
end subroutine xbpot


end module xtb_xtb_halogen
