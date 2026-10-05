! Only established calculator APIs: the same diagnostic links against stock
! and current libraries to identify pre-existing Energy/occupation effects.
program probe_thermal_energy
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator,newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol,moved
   type(TxTBCalculator) :: calc
   type(TRestart) :: center,plus,minus
   type(scc_results) :: results
   real(wp) :: xyz(3,3),e,g(3,3),sigma(3,3),gap,gp(3,3),gm(3,3),ep,em,step,fd(9),dn(9)
   integer :: atoms(3),level,x,a,iat
   logical :: failed
   call mctc_init('probe_thermal_energy',10,.true.);call init(env)
   atoms=[8,1,1];xyz(:,1)=0.0_wp
   xyz(:,2)=aatoau*[0.758602_wp,0.0_wp,0.504284_wp]
   xyz(:,3)=aatoau*[-0.3_wp,0.67_wp,0.62_wp]
   call init(mol,atoms,xyz,chrg=1.0_wp,uhf=1)
   call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
   call addSolvationModel(env,calc,TSolvInput(solvent='water',alpb=.false.,kernel=gbKernel%still))
   calc%maxiter=500;calc%etemp=1000.0_wp
   call newWavefunction(env,mol,calc,center)
   call calc%singlepoint(env,mol,center,0,.false.,e,g,sigma,gap,results)
   call env%check(failed);if(failed.or..not.results%converged) error stop 'reference convergence'
   print *, 'reference energy:',e
   print *, 'occupations:',center%wfn%focc
   do level=1,4
      step=2.0e-4_wp/real(2**(level-1),wp)
      do x=1,9
         a=mod(x-1,3)+1;iat=(x-1)/3+1
         call moved%copy(mol);moved%xyz(a,iat)=mol%xyz(a,iat)+step
         call plus%copy(center);call calc%singlepoint(env,moved,plus,0,.true.,ep,gp,sigma,gap,results)
         call env%check(failed);if(failed.or..not.results%converged) error stop 'positive convergence'
         moved%xyz(a,iat)=mol%xyz(a,iat)-step
         call minus%copy(center);call calc%singlepoint(env,moved,minus,0,.true.,em,gm,sigma,gap,results)
         call env%check(failed);if(failed.or..not.results%converged) error stop 'negative convergence'
         fd(x)=(ep-em)/(2.0_wp*step)
         dn(x)=sum(plus%wfn%focc-minus%wfn%focc)/(2.0_wp*step)
      enddo
      print *, 'step:',step
      do x=1,9
         a=mod(x-1,3)+1;iat=(x-1)/3+1
         write(*,'(i3,3es25.16)') x,fd(x),fd(x)-g(a,iat),dn(x)
      enddo
   enddo
end program
