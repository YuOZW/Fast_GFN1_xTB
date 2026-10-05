program probe_nodal_gradient
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator,newWavefunction
   use xtb_xtb_hessian_response, only : buildGFN1GasHessian
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol,moved
   type(TxTBCalculator) :: calc
   type(TRestart) :: center,plus,minus
   type(TRestart) :: rotatedReference
   type(scc_results) :: results
   real(wp) :: xyz(3,3),g(3,3),ag(3,3),gp(3,3),gm(3,3),fd(3,3),h(9,9),fdH(9,9)
   real(wp) :: rotation(3,3),rotG(3,3),rotH(9,9),transform(9,9)
   real(wp) :: energy,gap,sigma(3,3),ep,em,step,centerEnergy,previous,errorH
   integer :: coordinate,a,atom,level
   logical :: ok,failed,legacy
   character(len=32) :: mode
   call mctc_init('probe_nodal_gradient',10,.true.)
   call get_command_argument(1,mode); legacy=trim(mode)=='legacy'
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]
   xyz(:,2)=[0.758602_wp,0.0_wp,0.504284_wp]
   xyz(:,3)=[-0.3_wp,0.65_wp,0.62_wp]
   call init(env); call init(mol,[8,1,1],aatoau*xyz)
   call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
   call env%check(failed); if(failed) error stop 'calculator construction'
   calc%maxiter=500; calc%etemp=300.0_wp
   call newWavefunction(env,mol,calc,center)
   call calc%singlepoint(env,mol,center,0,.false.,energy,g,sigma,gap,results)
   call env%check(failed); if(failed.or..not.results%converged) error stop 'reference convergence'
   centerEnergy=energy
   call buildGFN1GasHessian(env,mol,calc%basis,calc%xtbData,center%wfn,calc%etemp,calc%accuracy,ag,h,ok)
   if(.not.ok) error stop 'gas Hessian reference'
   previous=huge(1.0_wp)
   do level=1,2
   step=4.0e-4_wp/real(2**(level-1),wp)
   do coordinate=1,9
      a=mod(coordinate-1,3)+1; atom=(coordinate-1)/3+1
      call moved%copy(mol); moved%xyz(a,atom)=mol%xyz(a,atom)+step
      call plus%copy(center)
      call calc%singlepoint(env,moved,plus,0,.true.,ep,gp,sigma,gap,results)
      call env%check(failed); if(failed.or..not.results%converged) error stop 'positive convergence'
      moved%xyz(a,atom)=mol%xyz(a,atom)-step
      call minus%copy(center)
      call calc%singlepoint(env,moved,minus,0,.true.,em,gm,sigma,gap,results)
      call env%check(failed); if(failed.or..not.results%converged) error stop 'negative convergence'
      fd(a,atom)=(ep-em)/(2.0_wp*step)
      fdH(:,coordinate)=reshape((gp-gm)/(2.0_wp*step),[9])
   enddo
   if(.not.legacy) then
      errorH=maxval(abs(h-fdH))
      print *, 'step=',step,' raw Hessian vs production Gradient FD=',errorH
      if(maxval(abs(g-ag))>1.0e-10_wp) error stop 'nodal production Gradient vs analytic'
      if(maxval(abs(g-fd))>5.0e-8_wp) error stop 'nodal Gradient vs all-coordinate Energy FD'
      if(errorH>5.0e-7_wp) error stop 'nodal Hessian vs all-coordinate Gradient FD'
      if(level==2.and.errorH>0.4_wp*previous+1.0e-8_wp) error stop 'nodal Hessian quadratic convergence'
      previous=errorH
   endif
   enddo
   print *, 'production vs analytic Gradient = ',maxval(abs(g-ag))
   print *, 'production vs energy FD = ',maxval(abs(g-fd))
   print *, 'analytic vs energy FD = ',maxval(abs(ag-fd))
   print *, 'production Gradient = ',g
   print *, 'analytic Gradient = ',ag
   print *, 'energy FD = ',fd
   if(legacy) then
      if(maxval(abs(g-ag))<1.0e-2_wp) error stop 'legacy nodal defect was not exercised'
      if(maxval(abs(ag-fd))>5.0e-8_wp) error stop 'independent analytic reference differs from Energy FD'
      print *, 'PASS legacy defect reproduction against independent Energy derivatives'
   else
      rotation=reshape([0.36_wp,0.8_wp,-0.48_wp,-0.48_wp,0.6_wp,0.64_wp,0.8_wp,0.0_wp,0.6_wp],[3,3])
      call moved%copy(mol); moved%xyz=matmul(rotation,mol%xyz)
      call newWavefunction(env,moved,calc,rotatedReference)
      call calc%singlepoint(env,moved,rotatedReference,0,.false.,energy,rotG,sigma,gap,results)
      call env%check(failed); if(failed.or..not.results%converged) error stop 'rotated SCC convergence'
      if(abs(energy-centerEnergy)>1.0e-10_wp) error stop 'rotational Energy invariance'
      if(maxval(abs(rotG-matmul(rotation,g)))>1.0e-9_wp) error stop 'rotational Gradient covariance'
      call buildGFN1GasHessian(env,moved,calc%basis,calc%xtbData,rotatedReference%wfn, &
         & calc%etemp,calc%accuracy,gp,rotH,ok)
      if(.not.ok) error stop 'rotated Hessian reference'
      transform=0.0_wp
      do atom=1,3
         transform(3*atom-2:3*atom,3*atom-2:3*atom)=rotation
      enddo
      if(maxval(abs(rotH-matmul(transform,matmul(h,transpose(transform)))))>1.0e-7_wp) &
         & error stop 'rotational Hessian covariance'
      print *, 'PASS nodal Gradient, full Hessian, Energy derivatives and rotational covariance'
   endif
end program
