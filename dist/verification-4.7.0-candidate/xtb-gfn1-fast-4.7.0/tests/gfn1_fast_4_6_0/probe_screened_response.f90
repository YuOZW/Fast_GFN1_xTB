! Compare the response Gradient with the actual converged production Gradient
! at an exact overlap node and at small nonzero, SCF-screened offsets.
! Diagnose mode records the old behavior without claiming it is correct.
program probe_screened_response
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator,newWavefunction
   use xtb_xtb_hessian_response, only : buildGFN1GasHessian
   use xtb_intgrad, only : get_grad_overlap
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TxTBCalculator) :: calc
   type(TRestart) :: center
   type(scc_results) :: results
   real(wp), parameter :: accuracy=1.0e-4_wp
   real(wp), parameter :: offsets(4)=[0.0_wp,1.0e-13_wp,-1.0e-13_wp,1.0e-5_wp]
   real(wp) :: xyz(3,3),g(3,3),ag(3,3),h(9,9),energy,gap,sigma(3,3)
   real(wp) :: point(3),ss(6,6),sg(3,6,6),raw,error
   integer :: k,accepted,rejected
   logical :: ok,failed,diagnose
   character(len=32) :: mode
   call mctc_init('probe_screened_response',10,.true.)
   call init(env);call get_command_argument(1,mode);diagnose=trim(mode)=='diagnose'
   accepted=0;rejected=0
   do k=1,size(offsets)
      xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]
      xyz(:,2)=[0.758602_wp,offsets(k),0.504284_wp]
      xyz(:,3)=[-0.3_wp,0.65_wp,0.62_wp]
      call init(mol,[8,1,1],aatoau*xyz)
      call newXTBCalculator(env,mol,calc,method=1,accuracy=accuracy)
      call env%check(failed);if(failed) error stop 'screened calculator construction'
      calc%maxiter=500
      call newWavefunction(env,mol,calc,center)
      call calc%singlepoint(env,mol,center,0,.false.,energy,g,sigma,gap,results)
      call env%check(failed)
      if(failed.or..not.results%converged) error stop 'screened SCC convergence'
      point=0.0_wp
      call get_grad_overlap(calc%basis%caoshell(2,1),calc%basis%caoshell(1,2), &
         & 3,1,1,0,mol%xyz(:,1),mol%xyz(:,2),point, &
         & max(20.0_wp,25.0_wp-10.0_wp*log10(accuracy)), &
         & calc%basis%nprim,calc%basis%primcount,calc%basis%alp,calc%basis%cont,ss,sg)
      raw=ss(1,2)
      if(k==1.and.raw/=0.0_wp) error stop 'exact overlap node missing'
      if(k==2.or.k==3) then
         if(raw==0.0_wp.or.abs(raw)>=1.0e-8_wp*accuracy) error stop 'nonzero screened overlap missing'
      endif
      if(k==4.and.abs(raw)<=1.0e-8_wp*accuracy) error stop 'ordinary unscreened overlap missing'
      call buildGFN1GasHessian(env,mol,calc%basis,calc%xtbData,center%wfn, &
         & calc%etemp,calc%accuracy,ag,h,ok)
      error=maxval(abs(ag-g))
      print *, 'case/offset/raw overlap/response accepted/Gradient difference:',k,offsets(k),raw,ok,error
      flush(6)
      if(ok) then
         accepted=accepted+1
         if(.not.diagnose.and.error>1.0e-10_wp) error stop 'screened production/response Gradient mismatch'
      else
         rejected=rejected+1
      endif
      if(.not.diagnose) then
         if((k==1.or.k==4).and..not.ok) error stop 'smooth nodal/ordinary response rejected'
         if((k==2.or.k==3).and.ok) error stop 'nonvariational screened node must fall back'
      endif
   enddo
   print *, 'screened response accepted/rejected counts:',accepted,rejected
   if(.not.diagnose) print *, 'PASS exact/ordinary response and screened-node fallback'
end program
