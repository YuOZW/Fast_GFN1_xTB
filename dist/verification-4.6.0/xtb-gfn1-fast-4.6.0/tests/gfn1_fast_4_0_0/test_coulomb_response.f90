program test_coulomb_response
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_coulomb, only : TxTBCoulomb
   use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
   implicit none
   type(TxTBCoulomb) :: coulomb,empty
   integer, parameter :: nsh=3,nat=2
   integer :: sh2at(nsh),i,j,mode
   real(wp) :: qat(nat),qsh(nsh),kernel(nsh,nsh),qa(nat),qs(nsh),atomic(nat),shell(nsh)
   real(wp) :: plus(nsh),minus(nsh),step,err
   logical :: ok
   sh2at=[1,1,2]; qat=[0.2_wp,-0.2_wp]; qsh=[0.1_wp,0.1_wp,-0.2_wp]
   allocate(coulomb%jmat(nsh,nsh))
   coulomb%jmat=ieee_value(0.0_wp,ieee_quiet_nan)
   do j=1,nsh
      do i=j,nsh
         coulomb%jmat(i,j)=0.2_wp/real(i+j,wp)
      enddo
   enddo
   do mode=1,3
      if(mode==2) coulomb%thirdOrder%atomicGam=[0.04_wp,0.07_wp]
      if(mode==3) coulomb%thirdOrder%shellGam=[0.03_wp,-0.05_wp,0.02_wp]
      call coulomb%getResponseKernel(qat,qsh,sh2at,kernel,ok)
      call require(ok,'Coulomb response kernel uses lower triangle')
      step=1.0e-5_wp; err=0.0_wp
      do j=1,nsh
         qs=qsh; qa=qat; qs(j)=qs(j)+step; qa(sh2at(j))=qa(sh2at(j))+step
         atomic=0.0_wp; shell=0.0_wp
         call coulomb%addShift(qa,qs,atomic,shell)
         plus=shell+atomic(sh2at)
         qs=qsh; qa=qat; qs(j)=qs(j)-step; qa(sh2at(j))=qa(sh2at(j))-step
         atomic=0.0_wp; shell=0.0_wp
         call coulomb%addShift(qa,qs,atomic,shell)
         minus=shell+atomic(sh2at)
         err=max(err,maxval(abs((plus-minus)/(2.0_wp*step)-kernel(:,j))))
      enddo
      print *, 'third-order mode=',mode,' kernel derivative error=',err
      call require(err<1.0e-10_wp,'actual addShift finite difference')
   enddo
   call empty%getResponseKernel(qat,qsh,sh2at,kernel,ok)
   call require(.not.ok,'missing Coulomb matrix rejected')
   sh2at(1)=0
   call coulomb%getResponseKernel(qat,qsh,sh2at,kernel,ok)
   call require(.not.ok,'invalid shell mapping rejected')
   sh2at(1)=1; coulomb%jmat(2,1)=ieee_value(0.0_wp,ieee_quiet_nan)
   call coulomb%getResponseKernel(qat,qsh,sh2at,kernel,ok)
   call require(.not.ok,'nonfinite authoritative Coulomb entry rejected')
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
end program
