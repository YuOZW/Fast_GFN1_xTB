program probe_solvent_gradient
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator,newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel,addGradientStill
   use xtb_solv_model, only : newBornModel
   use xtb_solv_gbsa, only : TBorn,addGradientHBond
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol,moved
   type(TxTBCalculator) :: calc
   type(TRestart) :: center,plus,minus
   type(scc_results) :: results
   type(TBorn) :: born
   real(wp) :: xyz(3,3),g(3,3),gp(3,3),gm(3,3),fd(3,3),fixedG(3,3),fixedFD(3,3)
   real(wp) :: energy,gap,sigma(3,3),ep,em,step,qat(3),sp,sm
   real(wp) :: partsP(4),partsM(4),componentG(3,3,3),partsFD(3),dum
   integer :: coordinate,a,atom,level
   logical :: failed,stable
   character(len=32) :: mode
   call mctc_init('probe_solvent_gradient',10,.true.)
   call get_command_argument(1,mode); stable=trim(mode)=='stable'
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]
   xyz(:,2)=[0.758602_wp,0.0_wp,0.504284_wp]
   xyz(:,3)=[-0.3_wp,0.65_wp,0.62_wp]
   if(stable) xyz(2,3)=0.67_wp
   call init(env); call init(mol,[8,1,1],aatoau*xyz)
   call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
   call addSolvationModel(env,calc,TSolvInput(solvent='water',alpb=.false.,kernel=gbKernel%still))
   calc%maxiter=500
   call newWavefunction(env,mol,calc,center)
   call calc%singlepoint(env,mol,center,0,.false.,energy,g,sigma,gap,results)
   call env%check(failed); if(failed.or..not.results%converged) error stop 'reference convergence'
   call newBornModel(calc%solvation,env,born,mol%at)
   call born%update(env,mol%at,mol%xyz)
   call born%updateCM5Cache(mol%at,mol%xyz)
   qat=center%wfn%q+born%cm5aCache
   fixedG=0.0_wp
   call born%addGradient(env,mol%at,mol%xyz,qat,center%wfn%qsh,fixedG)
   componentG=0.0_wp
   call addGradientStill(born%nat,born%ntpair,born%ppind,born%ddpair,qat,born%keps, &
      & born%brad,born%brdr,dum,componentG(:,:,1))
   call addGradientHBond(born%nat,born%at,qat,born%hbw,born%dhbdw,born%dsdrt,dum,componentG(:,:,2))
   componentG(:,:,3)=born%dsdr
   if(maxval(abs(sum(componentG,dim=3)-fixedG))>1.0e-12_wp) error stop 'solvent component assembly'
   do level=1,4
      step=4.0e-4_wp/real(2**(level-1),wp)
      do coordinate=1,9
         a=mod(coordinate-1,3)+1; atom=(coordinate-1)/3+1
         call moved%copy(mol); moved%xyz(a,atom)=mol%xyz(a,atom)+step
         call plus%copy(center)
         call calc%singlepoint(env,moved,plus,0,.true.,ep,gp,sigma,gap,results)
         call env%check(failed); if(failed.or..not.results%converged) error stop 'positive convergence'
         call born%update(env,moved%at,moved%xyz)
         call born%getEnergy(env,qat,center%wfn%qsh,sp)
         call born%getEnergyParts(env,qat,center%wfn%qsh,partsP(1),partsP(2),partsP(3),partsP(4))
         moved%xyz(a,atom)=mol%xyz(a,atom)-step
         call minus%copy(center)
         call calc%singlepoint(env,moved,minus,0,.true.,em,gm,sigma,gap,results)
         call env%check(failed); if(failed.or..not.results%converged) error stop 'negative convergence'
         call born%update(env,moved%at,moved%xyz)
         call born%getEnergy(env,qat,center%wfn%qsh,sm)
         call born%getEnergyParts(env,qat,center%wfn%qsh,partsM(1),partsM(2),partsM(3),partsM(4))
         fd(a,atom)=(ep-em)/(2.0_wp*step)
         fixedFD(a,atom)=(sp-sm)/(2.0_wp*step)
         partsFD=(partsP(1:3)-partsM(1:3))/(2.0_wp*step)
         print *, 'coordinate/step/full/fixed error:',coordinate,step,fd(a,atom)-g(a,atom), &
            & fixedFD(a,atom)-fixedG(a,atom)
         print *, 'Born/HB/SASA errors:',partsFD-componentG(a,atom,:)
         if(abs(partsFD(1)-componentG(a,atom,1))>1.0e-9_wp) error stop 'fixed-charge Born derivative'
      enddo
      print *, 'strict SCC step / max error:',step,maxval(abs(g-fd))
      print *, 'fixed charge step / max error:',step,maxval(abs(fixedG-fixedFD))
      if(stable.and.maxval(abs(g-fd))>1.0e-7_wp) error stop 'smooth nodal GBSA Energy derivative'
      if(stable.and.level==2.and.maxval(abs(g-fd))>2.0e-8_wp) error stop 'smooth GBSA refined derivative'
      if(.not.stable.and.level==2) then
         if(maxval(abs(fixedG-fixedFD))<1.0e-7_wp) error stop 'screened surface boundary not exercised'
         if(maxval(abs((g-fd)-(fixedG-fixedFD)))>1.0e-8_wp) error stop 'boundary error not solvent-local'
      endif
   enddo
   if(stable) then
      print *, 'PASS strict all-coordinate nodal GBSA Energy derivatives'
   else
      print *, 'PASS independent fixed-charge screened surface discontinuity reproduction'
   endif
end program
