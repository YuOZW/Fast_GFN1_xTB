program test_geometry
   use xtb_mctc_accuracy, only : wp
   use xtb_intgrad, only : get_hess_overlap,get_grad_overlap,dtrf2
   use xtb_disp_coordinationnumber, only : cnType,getCoordinationNumber,getCoordinationNumberHessian
   use xtb_type_molecule, only : TMolecule,init
   implicit none
   integer, parameter :: nc(0:2)=[1,3,6]
   integer :: nprim(12),primcount(12),li,lj,i,j,a,b,geometry,cf,atom
   real(wp) :: alp(24),cont(24),ri(3),rj(3),point(3),s(6,6),g(3,6,6),hh(3,3,6,6)
   real(wp) :: oldS(6,6),oldG(3,6,6),gp(3,6,6),gm(3,6,6),tmp(6,6),displaced(3),error,step
   real(wp) :: xyz(3,4),cn(4),dcn(3,4,4),hcn(3,4,3,4,4),oldCN(4),oldDCN(3,4,4)
   real(wp) :: plus(3,4,4),minus(3,4,4),strain(3,3,4),trans(3,1),fd(3,4,4)
   type(TMolecule) :: mol,moved
   logical :: ok
   point=0.0_wp; nprim=2
   do i=1,12
      primcount(i)=2*(i-1); alp(2*i-1)=0.7_wp; alp(2*i)=1.3_wp
      cont(2*i-1)=0.5_wp*sin(real(i,wp)); cont(2*i)=0.4_wp*cos(real(i,wp))
   enddo
   step=2.0e-5_wp
   do geometry=1,3
      ri=[0.1_wp,0.2_wp,-0.3_wp]; rj=[1.3_wp,-0.5_wp,0.7_wp]
      if(geometry==2) rj=ri
      if(geometry==3) rj=[5.0_wp,-2.0_wp,1.0_wp]
      do li=0,2
         do lj=0,2
            call get_hess_overlap(0,6,nc(li),nc(lj),li,lj,ri,rj,point,25.0_wp, &
               & nprim,primcount,alp,cont,s,g,hh)
            call get_grad_overlap(0,6,nc(li),nc(lj),li,lj,ri,rj,point,25.0_wp, &
               & nprim,primcount,alp,cont,oldS,oldG)
            call require(maxval(abs(s-oldS))+maxval(abs(g-oldG))<2.0e-13_wp,'overlap value/gradient compatibility')
            error=0.0_wp
            do a=1,3
               displaced=ri; displaced(a)=displaced(a)+step
               call get_grad_overlap(0,6,nc(li),nc(lj),li,lj,displaced,rj,point,25.0_wp, &
                  & nprim,primcount,alp,cont,oldS,gp)
               displaced=ri; displaced(a)=displaced(a)-step
               call get_grad_overlap(0,6,nc(li),nc(lj),li,lj,displaced,rj,point,25.0_wp, &
                  & nprim,primcount,alp,cont,oldS,gm)
               do b=1,3
                  tmp=(gp(b,:,:)-gm(b,:,:))/(2.0_wp*step)-hh(b,a,:,:)
                  call dtrf2(tmp,li,lj)
                  error=max(error,maxval(abs(tmp)))
               enddo
               ! Mixed-center derivative is -ii, from translational invariance.
               displaced=rj; displaced(a)=displaced(a)+step
               call get_grad_overlap(0,6,nc(li),nc(lj),li,lj,ri,displaced,point,25.0_wp, &
                  & nprim,primcount,alp,cont,oldS,gp)
               displaced=rj; displaced(a)=displaced(a)-step
               call get_grad_overlap(0,6,nc(li),nc(lj),li,lj,ri,displaced,point,25.0_wp, &
                  & nprim,primcount,alp,cont,oldS,gm)
               do b=1,3
                  error=max(error,maxval(abs((gp(b,:,:)-gm(b,:,:))/(2.0_wp*step)+hh(b,a,:,:))))
               enddo
            enddo
            print *, 'geometry=',geometry,' shells=',li,lj,' Hessian error=',error
            call require(error<1.0e-8_wp,'Cartesian/spherical overlap Hessian finite difference')
         enddo
      enddo
   enddo
   rj=ri+[50.0_wp,0.0_wp,0.0_wp]
   call get_hess_overlap(0,6,6,6,2,2,ri,rj,point,25.0_wp,nprim,primcount,alp,cont,s,g,hh)
   call require(maxval(abs(s))+maxval(abs(g))+maxval(abs(hh))==0.0_wp,'overlap distance screening')

   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]; xyz(:,2)=[1.43_wp,0.2_wp,0.95_wp]
   xyz(:,3)=[-1.2_wp,0.1_wp,1.15_wp]; xyz(:,4)=[4.0_wp,-0.3_wp,0.1_wp]
   call init(mol,[8,1,1,14],xyz); trans=0.0_wp
   do cf=1,4
      call getCoordinationNumberHessian(mol,cf,cn,dcn,hcn,ok)
      call require(ok,'molecular CN Hessian')
      call getCoordinationNumber(mol,trans,40.0_wp,cf,oldCN,oldDCN,strain)
      call require(maxval(abs(cn-oldCN))+maxval(abs(dcn-oldDCN))<1.0e-13_wp,'production CN/gradient compatibility')
      error=0.0_wp
      do atom=1,4
         do a=1,3
            call moved%copy(mol); moved%xyz(a,atom)=mol%xyz(a,atom)+step
            call getCoordinationNumber(moved,trans,40.0_wp,cf,oldCN,plus,strain)
            moved%xyz(a,atom)=mol%xyz(a,atom)-step
            call getCoordinationNumber(moved,trans,40.0_wp,cf,oldCN,minus,strain)
            fd=(plus-minus)/(2.0_wp*step)-hcn(:,:,a,atom,:)
            error=max(error,maxval(abs(fd)))
         enddo
      enddo
      print *, 'CN type=',cf,' Hessian error=',error
      call require(error<1.0e-7_wp,'CN Hessian vs production gradient finite difference')
      call require(maxval(abs(sum(hcn,dim=2)))<1.0e-12_wp,'CN translational identity')
      do i=1,4
         do j=1,4
            do a=1,3
               do b=1,3
                  call require(maxval(abs(hcn(a,i,b,j,:)-hcn(b,j,a,i,:)))<1.0e-12_wp,'CN Hessian symmetry')
               enddo
            enddo
         enddo
      enddo
   enddo
   call getCoordinationNumberHessian(mol,0,cn,dcn,hcn,ok)
   call require(.not.ok,'unsupported CN type rejected')
   moved%npbc=3
   call getCoordinationNumberHessian(moved,cnType%exp,cn,dcn,hcn,ok)
   call require(.not.ok,'periodic CN Hessian defers to fallback')
contains
   subroutine require(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if(.not.condition) then
         print *, 'FAIL ',label
         error stop 1
      endif
   end subroutine
end program
