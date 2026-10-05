! All 339 Cartesian columns of the raw taxol solvent Hessian, independently
! differentiated through converged nonlinear SCC. No projection/symmetrizing.
program test_large_molecular_solvent
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
   use xtb_xtb_hessian_response, only : buildGFN1SolvatedHessian
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol,moved
   type(TxTBCalculator) :: calc
   type(TRestart) :: center,plus,minus
   type(scc_results) :: results
   type(TBorn) :: born
   real(wp), allocatable :: g(:,:),analyticG(:,:),h(:,:),gp(:,:),gm(:,:),fd(:,:),translation(:)
   real(wp) :: e,ep,em,sigma(3,3),gap,step,t0,t1,err,eg
   integer :: m,x,a,iat,n,unit
   logical :: ok,failed
   character(len=512) :: path
   call mctc_init('test_large_molecular_solvent',10,.true.);call init(env)
   call get_command_argument(1,path);open(newunit=unit,file=trim(path),status='old')
   call readMolecule(env,mol,unit,fileType%xyz);close(unit)
   n=mol%n;step=1.0e-5_wp
   allocate(g(3,n),analyticG(3,n),h(3*n,3*n),gp(3,n),gm(3,n),fd(3*n,3*n),translation(3*n))
   do m=1,2
      call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
      call addSolvationModel(env,calc,TSolvInput(solvent='water',alpb=m==2,kernel=m))
      calc%maxiter=500
      call newWavefunction(env,mol,calc,center)
      call calc%singlepoint(env,mol,center,0,.false.,e,g,sigma,gap,results)
      call env%check(failed);if(failed.or..not.results%converged) error stop 'large SCC convergence'
      call newBornModel(calc%solvation,env,born,mol%at);call born%update(env,mol%at,mol%xyz)
      call cpu_time(t0)
      call buildGFN1SolvatedHessian(env,mol,calc%basis,calc%xtbData,center%wfn,calc%etemp, &
         & calc%accuracy,born,analyticG,h,ok,guardRadius=step)
      call cpu_time(t1)
      if(.not.ok) error stop 'large guarded solvent response rejected'
      print *, 'large model/E/G/analytic CPU/Gradient error:',m,e,norm2(g),t1-t0,maxval(abs(g-analyticG))
      flush(6)
      if(maxval(abs(g-analyticG))>1.0e-9_wp) error stop 'large original Gradient parity'
      if(maxval(abs(h-transpose(h)))>1.0e-7_wp) error stop 'large raw reciprocity'
      do a=1,3
         translation=0.0_wp;translation(a:3*n:3)=1.0_wp
         if(maxval(abs(matmul(h,translation)))>1.0e-7_wp) error stop 'large raw translation'
      enddo
      call cpu_time(t0);eg=0.0_wp
      do x=1,3*n
         a=mod(x-1,3)+1;iat=(x-1)/3+1
         call moved%copy(mol);moved%xyz(a,iat)=mol%xyz(a,iat)+step
         call plus%copy(center)
         call calc%singlepoint(env,moved,plus,0,.true.,ep,gp,sigma,gap,results)
         call env%check(failed);if(failed.or..not.results%converged) error stop 'large positive SCC convergence'
         moved%xyz(a,iat)=mol%xyz(a,iat)-step
         call minus%copy(center)
         call calc%singlepoint(env,moved,minus,0,.true.,em,gm,sigma,gap,results)
         call env%check(failed);if(failed.or..not.results%converged) error stop 'large negative SCC convergence'
         fd(:,x)=reshape(gp-gm,[3*n])/(2.0_wp*step)
         eg=max(eg,abs((ep-em)/(2.0_wp*step)-g(a,iat)))
         if(mod(x,50)==0) then
            print *, 'large model/completed columns/H difference:',m,x,maxval(abs(fd(:,:x)-h(:,:x)))
            flush(6)
         endif
      enddo
      call cpu_time(t1);err=maxval(abs(fd-h))
      print *, 'large full SCC model/columns/step/H error/E-G error/FD CPU:',m,3*n,step,err,eg,t1-t0
      flush(6)
      if(err>2.0e-6_wp) error stop 'large full nonlinear SCC Hessian differences'
   enddo
   print *, 'PASS two large solvated SCC Hessians, all 339 raw columns'
end program
