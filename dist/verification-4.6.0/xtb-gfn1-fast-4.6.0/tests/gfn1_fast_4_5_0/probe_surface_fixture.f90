! Select a deterministic smooth angular-quadrature fixture by original
! surface branch margins alone, before any molecular Hessian comparison.
program probe_surface_fixture
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : TSolvModel,init,newBornModel
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_response, only : getSurfaceResponse
   implicit none
   type(TEnvironment) :: env
   type(TSolvModel) :: model
   type(TBorn) :: born(2)
   integer :: atoms(8),trial,i,a,m
   real(wp) :: base(3,8),xyz(3,8),surface(8),g(24,8),h(24,24,8)
   logical :: ok,good
   call mctc_init('probe_surface_fixture',10,.true.);call init(env)
   atoms=[35,6,1,1,1,8,1,1]
   base(:,1)=0.0_wp;base(:,2)=[-1.94_wp,0.04_wp,-0.05_wp]
   base(:,3)=base(:,2)+[-0.36_wp,1.026_wp,0.03_wp]
   base(:,4)=base(:,2)+[-0.37_wp,-0.51_wp,0.89_wp]
   base(:,5)=base(:,2)+[-0.35_wp,-0.52_wp,-0.88_wp]
   base(:,6)=[3.2_wp,0.17_wp,-0.13_wp]
   base(:,7)=base(:,6)+[0.57_wp,0.74_wp,0.12_wp]
   base(:,8)=base(:,6)+[0.58_wp,-0.75_wp,-0.09_wp]
   do m=1,2
      call init(model,env,TSolvInput(solvent='water',alpb=m==2,kernel=1),1)
      call newBornModel(model,env,born(m),atoms)
   enddo
   do trial=1,100
      do i=1,8
         do a=1,3
            xyz(a,i)=aatoau*(base(a,i)+0.1_wp*sin(real(19*trial+7*i+3*a,wp)))
         enddo
      enddo
      good=.true.
      do m=1,2
         call born(m)%update(env,atoms,xyz)
         call getSurfaceResponse(born(m),xyz,surface,g,h,ok,guardRadius=2.0e-4_wp)
         good=good.and.ok
      enddo
      if(good) then
         print *, 'smooth original surface margin trial:',trial
         do i=1,8
            write(*,'(i3,3f23.16)') atoms(i),xyz(:,i)/aatoau
         enddo
         stop
      endif
   enddo
   error stop 'no smooth fixture within bounded search'
end program
