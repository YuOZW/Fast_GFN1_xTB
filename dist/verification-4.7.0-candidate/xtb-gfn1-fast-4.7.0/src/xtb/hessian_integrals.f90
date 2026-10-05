! This file is part of xtb.
! SPDX-License-Identifier: LGPL-3.0-or-later
module xtb_xtb_hessian_integrals
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : evtoau
   use xtb_type_molecule, only : TMolecule
   use xtb_type_basisset, only : TBasisset
   use xtb_xtb_data, only : THamiltonianData
   use xtb_intgrad, only : get_hess_overlap,dtrf2
   use xtb_scc_core, only : h0scal
   !$ use omp_lib, only : omp_get_thread_num
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: buildGFN1IntegralResponse,shellPolynomialHessian
contains

! Molecular GFN1 H0/S coordinate derivatives, with linear CN self-energies.
! selfEnergy/dSEdcn are in the native eV convention; output dH and geometric
! Hessian use Eh. p is density; pew is energy-weighted density in Eh.
! geometricHessian = Tr(P d2H0)-Tr(W d2S), with P/W held fixed.
! Stream the second derivatives into a 3N x 3N matrix, avoiding the
! prohibitively large NAO x NAO x 3N x 3N integral tensor.
! SCC potential/Coulomb and other energy terms must be assembled separately.
subroutine buildGFN1IntegralResponse(mol,basis,hData,nShell,intcut,selfEnergy,dSEdcn,dcndr,hcn, &
      & p,pew,dS,dH,geometricHessian,ok,threads)
   type(TMolecule), intent(in) :: mol
   type(TBasisset), intent(in) :: basis
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   real(wp), intent(in) :: intcut,selfEnergy(:,:),dSEdcn(:,:),dcndr(:,:,:),hcn(:,:,:,:,:)
   real(wp), intent(in) :: p(:,:),pew(:,:)
   real(wp), intent(out) :: dS(:,:,:),dH(:,:,:),geometricHessian(:,:)
   logical, intent(out) :: ok
   integer, intent(in), optional :: threads
   integer, parameter :: nc(0:2)=[1,3,6],ns(0:2)=[1,3,5]
   real(wp) :: ss(6,6),sg(3,6,6),sh(3,3,6,6),tmp(6,6),poly,gpoly(3),hpoly(3,3)
   real(wp) :: rij(3),point(3),r2,zi,zj,zeta,km,factor,hav,value,first(3),second(3,3),block(3,3),pn,wn
   real(wp), allocatable :: dcn(:,:),cnWeight(:)
   real(wp), allocatable, target :: geometryWork(:,:,:),weightWork(:,:)
   logical, allocatable :: threadReady(:)
   integer :: threadCount,threadIndex,allocationStatus
   integer :: n,nao,n3,iat,jat,ish,jsh,iz,jz,li,lj,iao,jao,ii,jj,a,b,ci,cj,di,dj
   dS=0.0_wp; dH=0.0_wp; geometricHessian=0.0_wp; ok=.false.
   n=mol%n; nao=basis%nao; n3=3*n
   if(n<1.or.nao<1.or.mol%npbc/=0.or.basis%n/=n) return
   if(any(shape(p)/=[nao,nao]).or.any(shape(pew)/=[nao,nao])) return
   if(any(shape(dS)/=[nao,nao,n3]).or.any(shape(dH)/=[nao,nao,n3])) return
   if(any(shape(geometricHessian)/=[n3,n3])) return
   if(any(shape(dcndr)/=[3,n,n]).or.any(shape(hcn)/=[3,n,3,n,n])) return
   if(size(selfEnergy,2)/=n.or.size(dSEdcn,2)/=n) return
   if(any(mol%at<1).or.any(mol%at>size(nShell))) return
   if(.not.ieee_is_finite(intcut).or.intcut<=0.0_wp) return
   if(.not.all(ieee_is_finite(p)).or..not.all(ieee_is_finite(pew))) return
   if(.not.all(ieee_is_finite(mol%xyz))) return
   if(.not.all(ieee_is_finite(selfEnergy)).or..not.all(ieee_is_finite(dSEdcn))) return
   if(maxval(abs(p-transpose(p)))>1.0e-12_wp*max(1.0_wp,maxval(abs(p)))) return
   if(maxval(abs(pew-transpose(pew)))>1.0e-12_wp*max(1.0_wp,maxval(abs(pew)))) return
   if(.not.all(ieee_is_finite(dcndr)).or..not.all(ieee_is_finite(hcn))) return
   do iat=1,n
      iz=mol%at(iat)
      if(nShell(iz)>min(size(selfEnergy,1),size(dSEdcn,1))) return
      do ish=1,nShell(iz)
         li=hData%angShell(ish,iz)
         if(li<0.or.li>2) return
      enddo
   enddo
   threadCount=1
   if(present(threads)) threadCount=max(1,min(threads,n))
   do iat=1,n
      do jat=1,iat-1
         if(sum((mol%xyz(:,iat)-mol%xyz(:,jat))**2)<1.0e-12_wp) return
      enddo
   enddo
   allocate(dcn(n3,n),cnWeight(n),geometryWork(n3,n3,threadCount),weightWork(n,threadCount), &
      & threadReady(threadCount),stat=allocationStatus)
   if(allocationStatus/=0) return
   dcn=reshape(dcndr,[n3,n]);cnWeight=0.0_wp
   geometryWork=0.0_wp;weightWork=0.0_wp;threadReady=.true.
   !$omp parallel num_threads(threadCount) default(none) &
   !$omp shared(mol,basis,hData,nShell,intcut,selfEnergy,dSEdcn,dcn,p,pew,dS,dH,n,nao,n3, &
   !$omp& geometryWork,weightWork,threadReady) private(threadIndex)
   threadIndex=1
   !$ threadIndex=omp_get_thread_num()+1
   call integralWorker(threadIndex)
   !$omp end parallel
   if(.not.all(threadReady)) return
   do threadIndex=1,threadCount
      geometricHessian=geometricHessian+geometryWork(:,:,threadIndex)
      cnWeight=cnWeight+weightWork(:,threadIndex)
   enddo
   do iat=1,n
      geometricHessian=geometricHessian+cnWeight(iat)*reshape(hcn(:,:,:,:,iat),[n3,n3])
   enddo
   ok=all(ieee_is_finite(dS)).and.all(ieee_is_finite(dH)).and.all(ieee_is_finite(geometricHessian))
contains
   recursive subroutine integralWorker(index)
      integer, intent(in) :: index
      real(wp), pointer :: localGeometry(:,:),localWeight(:)
      real(wp), allocatable :: havFirst(:),hfirst(:)
      real(wp) :: ss(6,6),sg(3,6,6),sh(3,3,6,6),tmp(6,6),poly,gpoly(3),hpoly(3,3)
      real(wp) :: rij(3),point(3),r2,zi,zj,zeta,km,factor,hav,value,first(3),second(3,3),block(3,3),pn,wn
      integer :: iat,jat,ish,jsh,iz,jz,li,lj,iao,jao,ii,jj,a,b,ci,cj,di,dj,allocationStatus
      allocate(havFirst(n3),hfirst(n3),stat=allocationStatus)
      threadReady(index)=allocationStatus==0
      localGeometry=>geometryWork(:,:,index);localWeight=>weightWork(:,index)
      point=0.0_wp
      !$omp do schedule(dynamic,1)
   do iat=1,n
      if(allocationStatus/=0) cycle
      iz=mol%at(iat)
      do jat=1,iat-1
         jz=mol%at(jat); rij=mol%xyz(:,iat)-mol%xyz(:,jat); r2=sum(rij**2)
         if(r2>2000.0_wp) cycle
         do ish=1,nShell(iz)
            li=hData%angShell(ish,iz)
            do jsh=1,nShell(jz)
               lj=hData%angShell(jsh,jz)
               zi=hData%slaterExponent(ish,iz); zj=hData%slaterExponent(jsh,jz)
               zeta=(2.0_wp*sqrt(zi*zj)/(zi+zj))**hData%wExp
               call h0scal(hData,li+1,lj+1,iz,jz,hData%valenceShell(ish,iz)/=0, &
                  & hData%valenceShell(jsh,jz)/=0,km)
               factor=0.5_wp*km*zeta*evtoau
               hav=factor*(selfEnergy(ish,iat)+selfEnergy(jsh,jat))
               havFirst=factor*(dSEdcn(ish,iat)*dcn(:,iat)+dSEdcn(jsh,jat)*dcn(:,jat))
               call shellPolynomialHessian(hData%shellPoly(li+1,iz),hData%shellPoly(lj+1,jz), &
                  & hData%atomicRad(iz),hData%atomicRad(jz),rij,poly,gpoly,hpoly)
               call get_hess_overlap(basis%caoshell(ish,iat),basis%caoshell(jsh,jat),nc(li),nc(lj), &
                  & li,lj,mol%xyz(:,iat),mol%xyz(:,jat),point,intcut, &
                  & basis%nprim,basis%primcount,basis%alp,basis%cont,ss,sg,sh)
               call dtrf2(ss,li,lj)
               do a=1,3
                  tmp=sg(a,:,:); call dtrf2(tmp,li,lj); sg(a,:,:)=tmp
                  do b=1,3
                     tmp=sh(a,b,:,:); call dtrf2(tmp,li,lj); sh(a,b,:,:)=tmp
                  enddo
               enddo
               do ii=1,ns(li)
                  iao=ii+basis%saoshell(ish,iat)
                  do jj=1,ns(lj)
                     jao=jj+basis%saoshell(jsh,jat)
                     value=poly*ss(jj,ii)
                     first=poly*sg(:,jj,ii)+gpoly*ss(jj,ii)
                     second=poly*sh(:,:,jj,ii)+hpoly*ss(jj,ii) &
                        & +spread(gpoly,2,3)*spread(sg(:,jj,ii),1,3) &
                        & +spread(sg(:,jj,ii),2,3)*spread(gpoly,1,3)
                     hfirst=havFirst*value
                     pn=2.0_wp*p(jao,iao); wn=2.0_wp*pew(jao,iao)
                     localWeight(iat)=localWeight(iat)+pn*value*factor*dSEdcn(ish,iat)
                     localWeight(jat)=localWeight(jat)+pn*value*factor*dSEdcn(jsh,jat)
                     do a=1,3
                        ci=3*(iat-1)+a; cj=3*(jat-1)+a
                        dS(jao,iao,ci)=sg(a,jj,ii); dS(jao,iao,cj)=-sg(a,jj,ii)
                        hfirst(ci)=hfirst(ci)+hav*first(a); hfirst(cj)=hfirst(cj)-hav*first(a)
                        localGeometry(:,ci)=localGeometry(:,ci)+pn*first(a)*havFirst
                        localGeometry(ci,:)=localGeometry(ci,:)+pn*first(a)*havFirst
                        localGeometry(:,cj)=localGeometry(:,cj)-pn*first(a)*havFirst
                        localGeometry(cj,:)=localGeometry(cj,:)-pn*first(a)*havFirst
                        do b=1,3
                           di=3*(iat-1)+b; dj=3*(jat-1)+b
                           block(a,b)=pn*hav*second(a,b)-wn*sh(a,b,jj,ii)
                           localGeometry(ci,di)=localGeometry(ci,di)+block(a,b)
                           localGeometry(cj,dj)=localGeometry(cj,dj)+block(a,b)
                           localGeometry(ci,dj)=localGeometry(ci,dj)-block(a,b)
                           localGeometry(cj,di)=localGeometry(cj,di)-block(a,b)
                        enddo
                     enddo
                     dS(iao,jao,:)=dS(jao,iao,:)
                     dH(jao,iao,:)=hfirst; dH(iao,jao,:)=hfirst
                  enddo
               enddo
            enddo
         enddo
      enddo
      do ish=1,nShell(iz)
         li=hData%angShell(ish,iz)
         do ii=1,ns(li)
            iao=ii+basis%saoshell(ish,iat)
            dH(iao,iao,:)=evtoau*dSEdcn(ish,iat)*dcn(:,iat)
            localWeight(iat)=localWeight(iat)+p(iao,iao)*evtoau*dSEdcn(ish,iat)
         enddo
      enddo
   enddo
      !$omp end do
   end subroutine integralWorker

end subroutine

! GFN1's radius convention and R^0.5 polynomial are preserved exactly.
pure subroutine shellPolynomialHessian(iPoly,jPoly,iRad,jRad,rij,value,gradient,hessian)
   real(wp), intent(in) :: iPoly,jPoly,iRad,jRad,rij(3)
   real(wp), intent(out) :: value,gradient(3),hessian(3,3)
   real(wp) :: r,rad,k1,k2,root,first,second,u(3)
   integer :: a
   r=sqrt(sum(rij**2)); rad=iRad+jRad; k1=0.01_wp*iPoly; k2=0.01_wp*jPoly
   root=sqrt(r/rad); value=(1.0_wp+k1*root)*(1.0_wp+k2*root)
   first=0.5_wp*(k1+k2)*root/r+k1*k2/rad
   second=-0.25_wp*(k1+k2)*root/r**2
   u=rij/r; gradient=first*u
   hessian=(second-first/r)*spread(u,2,3)*spread(u,1,3)
   do a=1,3; hessian(a,a)=hessian(a,a)+first/r; enddo
end subroutine
end module
