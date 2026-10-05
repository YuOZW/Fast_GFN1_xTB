program test_contraction
   use xtb_mctc_accuracy, only : wp
   use xtb_intgrad, only : get_grad_overlap, get_overlap, dtrf2, dtrf2_weights
   implicit none
   integer, parameter :: nc(0:2)=[1,3,6],ns(0:2)=[1,3,5]
   integer :: nprim(12),primcount(12),li,lj,i,j,k
   real(wp) :: alp(24),cont(24),s(6,6),sg(3,6,6),weights(6,6),cart(6,6),tmp(6,6)
   real(wp) :: ri(3),rj(3),point(3),expected(3),actual(3),plus(6,6),minus(6,6),displaced(3),fd
   real(wp), parameter :: step=1.0e-5_wp
   nprim=2; point=0.0_wp; ri=[0.1_wp,0.2_wp,-0.3_wp]; rj=[1.3_wp,-0.5_wp,0.7_wp]
   do i=1,12
      primcount(i)=2*(i-1)
      alp(2*i-1)=0.7_wp; alp(2*i)=1.3_wp
      cont(2*i-1)=0.5_wp*sin(real(i,wp)); cont(2*i)=0.4_wp*cos(real(i,wp))
   enddo
   do li=0,2
      do lj=0,2
         weights=0.0_wp
         do i=1,ns(li)
            do j=1,ns(lj)
               weights(j,i)=sin(real(i+3*j,wp))
            enddo
         enddo
         call get_grad_overlap(0,6,nc(li),nc(lj),li,lj,ri,rj,point,25.0_wp, &
            & nprim,primcount,alp,cont,s,sg)
         do k=1,3
            tmp=sg(k,:,:)
            call dtrf2(tmp,li,lj)
            expected(k)=sum(tmp(1:ns(lj),1:ns(li))*weights(1:ns(lj),1:ns(li)))
         enddo
         cart=weights
         call dtrf2_weights(cart,li,lj)
         call get_grad_overlap(0,6,nc(li),nc(lj),li,lj,ri,rj,point,25.0_wp, &
            & nprim,primcount,alp,cont,s,sg,cart,actual)
         if(maxval(abs(expected-actual))>1.0e-12_wp) error stop 'spherical/Cartesian contraction mismatch'
         do k=1,3
            displaced=ri; displaced(k)=displaced(k)+step; plus=0.0_wp
            call get_overlap(0,6,nc(li),nc(lj),li,lj,displaced,rj,point,25.0_wp, &
               & nprim,primcount,alp,cont,plus)
            call dtrf2(plus,li,lj)
            displaced=ri; displaced(k)=displaced(k)-step; minus=0.0_wp
            call get_overlap(0,6,nc(li),nc(lj),li,lj,displaced,rj,point,25.0_wp, &
               & nprim,primcount,alp,cont,minus)
            call dtrf2(minus,li,lj)
            fd=sum((plus(1:ns(lj),1:ns(li))-minus(1:ns(lj),1:ns(li))) &
               & *weights(1:ns(lj),1:ns(li)))/(2.0_wp*step)
            if(abs(fd-actual(k))>1.0e-8_wp) error stop 'contracted derivative finite difference mismatch'
         enddo
         print *, 'PASS primitive contraction and finite differences, shell l=',li,lj
      enddo
   enddo
end program
