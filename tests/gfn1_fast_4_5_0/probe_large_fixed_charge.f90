! Full fixed-bare-charge solvent Hessian vs the original functional at the
! actual converged large-molecule charge distribution. No SCC in displacements.
program probe_large_fixed_charge
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule
   use xtb_io_reader, only : readMolecule
   use xtb_mctc_filetypes, only : fileType
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator,newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : newBornModel
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_cm5, only : calc_cm5
   use xtb_solv_response, only : buildGFN1SolventResponse
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TxTBCalculator) :: calc
   type(TRestart) :: center
   type(scc_results) :: results
   type(TBorn) :: born
   real(wp), allocatable :: q(:),moved(:,:),g(:,:),h(:,:),v(:),dv(:,:),amat(:,:)
   real(wp), allocatable :: gp(:,:),gm(:,:),fd(:,:),cm(:),dc(:,:,:),charges(:)
   real(wp) :: e,sigma(3,3),gap,step,t0,t1
   integer :: m,k,x,a,iat,n,unit,coord
   logical :: ok,failed
   character(len=512) :: path
   call mctc_init('probe_large_fixed_charge',10,.true.);call init(env)
   call get_command_argument(1,path);open(newunit=unit,file=trim(path),status='old')
   call readMolecule(env,mol,unit,fileType%xyz);close(unit)
   n=mol%n
   allocate(q(n),moved(3,n),g(3,n),h(3*n,3*n),v(n),dv(n,3*n),amat(n,n), &
      & gp(3,n),gm(3,n),fd(3*n,3*n),cm(n),dc(3,n,n),charges(n))
   do m=1,2
      call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
      call addSolvationModel(env,calc,TSolvInput(solvent='water',alpb=m==2,kernel=m))
      call newWavefunction(env,mol,calc,center)
      call calc%singlepoint(env,mol,center,0,.false.,e,g,sigma,gap,results)
      call env%check(failed);if(failed.or..not.results%converged) error stop 'large SCC convergence'
      q=center%wfn%q
      call newBornModel(calc%solvation,env,born,mol%at);call born%update(env,mol%at,mol%xyz)
      call cpu_time(t0)
      call buildGFN1SolventResponse(born,mol%at,mol%xyz,q,e,g,h,v,dv,amat,ok)
      call cpu_time(t1)
      if(.not.ok) error stop 'large local solvent response rejected'
      call parent(mol%xyz,gp)
      print *, 'large solvent model/local CPU/Gradient parity:',m,t1-t0,maxval(abs(g-gp))
      if(maxval(abs(g-gp))>1.0e-11_wp) error stop 'original large Gradient parity'
      do k=1,3
         step=5.0e-4_wp/real(5**(k-1),wp)
         call cpu_time(t0)
         do x=1,3*n
            a=mod(x-1,3)+1;iat=(x-1)/3+1
            moved=mol%xyz;moved(a,iat)=moved(a,iat)+step;call parent(moved,gp)
            moved(a,iat)=mol%xyz(a,iat)-step;call parent(moved,gm)
            fd(:,x)=reshape(gp-gm,[3*n])/(2.0_wp*step)
         enddo
         call cpu_time(t1)
         coord=maxloc(maxval(abs(fd-h),dim=1),dim=1)
         print *, 'model/step/full fixed-q Hessian error/column/CPU:',m,step,maxval(abs(fd-h)),coord,t1-t0
         flush(6)
      enddo
   enddo
contains
   subroutine parent(coordinates,gradient)
      real(wp), intent(in) :: coordinates(:,:)
      real(wp), intent(out) :: gradient(:,:)
      integer :: i
      call born%update(env,mol%at,coordinates)
      call calc_cm5(n,mol%at,coordinates,cm,dc);charges=q+cm
      v=matmul(born%bornMat,charges);gradient=0.0_wp
      call born%addGradient(env,mol%at,coordinates,charges,center%wfn%qsh,gradient)
      do i=1,n
         gradient=gradient+dc(:,:,i)*v(i)
      enddo
   end subroutine
end program
