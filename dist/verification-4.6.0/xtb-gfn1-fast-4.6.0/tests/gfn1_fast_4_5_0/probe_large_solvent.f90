! Inspect actual large-molecule solvent branch acceptance before timing a
! complete Hessian. The Born model and coordinates use the production reader.
program probe_large_solvent
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule
   use xtb_io_reader, only : readMolecule
   use xtb_mctc_filetypes, only : fileType
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : TSolvModel,init,newBornModel
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_response, only : getSurfaceResponse
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TSolvModel) :: model
   type(TBorn) :: born
   real(wp), allocatable :: surface(:),g(:,:),h(:,:,:)
   real(wp) :: windows(5),t0,t1
   integer :: m,k,unit,n
   logical :: ok,failed
   character(len=512) :: path
   call mctc_init('probe_large_solvent',10,.true.);call init(env)
   call get_command_argument(1,path)
   open(newunit=unit,file=trim(path),status='old')
   call readMolecule(env,mol,unit,fileType%xyz);close(unit)
   call env%check(failed);if(failed) error stop 'large fixture read'
   n=mol%n;allocate(surface(n),g(3*n,n),h(3*n,3*n,n))
   windows=[5.0e-3_wp,5.0e-4_wp,1.0e-4_wp,1.0e-5_wp,0.0_wp]
   do m=1,2
      call init(model,env,TSolvInput(solvent='water',alpb=m==2,kernel=gbKernel%still),1)
      call newBornModel(model,env,born,mol%at);call born%update(env,mol%at,mol%xyz)
      do k=1,5
         call cpu_time(t0)
         call getSurfaceResponse(born,mol%xyz,surface,g,h,ok,guardRadius=windows(k))
         call cpu_time(t1)
         print *, 'model/window/accepted/CPU:',m,windows(k),ok,t1-t0
         flush(6)
         if(k==5.and..not.ok) error stop 'local surface response rejected'
      enddo
   enddo
end program
