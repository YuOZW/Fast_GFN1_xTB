program test_calculator_hessian
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator,newWavefunction
   use xtb_xtb_hessian_response, only : buildGFN1GasHessian
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy,getGFN1FastPolicy
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TxTBCalculator) :: calc
   type(TRestart) :: chk
   type(scc_results) :: results
   type(TGFN1FastPolicy) :: policy
   real(wp) :: xyz(3,3),energy,gap,sigma(3,3),gradient(3,3),expected(9,9),h(9,9),dip(3,9),pol(6,9)
   real(wp) :: referenceP(8,8),referenceQ(6)
   logical :: ok,failed
   integer :: col
   call mctc_init('test_calculator_hessian',10,.true.)
   call init(env)
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]
   xyz(:,2)=[0.758602_wp,0.13_wp,0.504284_wp]
   xyz(:,3)=[-0.758602_wp,-0.18_wp,0.62_wp]
   call init(mol,[8,1,1],xyz*aatoau)
   call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
   call newWavefunction(env,mol,calc,chk)
   call calc%singlepoint(env,mol,chk,0,.false.,energy,gradient,sigma,gap,results)
   call env%check(failed)
   call require(.not.failed.and.results%converged,'calculator reference SCC')
   referenceP=chk%wfn%p; referenceQ=chk%wfn%qsh
   call getGFN1FastPolicy(policy)
   call require(policy%analyticHessian.and.policy%profile,'runner enables analytic Hessian/profile')
   call buildGFN1GasHessian(env,mol,calc%basis,calc%xtbData,chk%wfn,calc%etemp,calc%accuracy, &
      & gradient,expected,ok)
   call require(ok,'complete molecular reference')
   !$omp parallel num_threads(2) default(none) shared(env,mol,calc,chk,gradient,h,ok)
   !$omp single
   call buildGFN1GasHessian(env,mol,calc%basis,calc%xtbData,chk%wfn,calc%etemp,calc%accuracy,gradient,h,ok)
   !$omp end single
   !$omp end parallel
   call require(ok.and.maxval(abs(h-expected))<1.0e-12_wp,'nested analytic API preserves full Hessian')
   h=1.0_wp; dip=2.0_wp
   call calc%hessian(env,mol,chk,[2],5.0e-4_wp,h,dip)
   do col=1,9
      if(col>=4.and.col<=6) then
         call require(maxval(abs(h(:,col)-1.0_wp-expected(:,col)))<1.0e-12_wp,'selected additive analytic columns')
         call require(maxval(abs(dip(:,col)))==0.0_wp,'selected dipole output')
      else
         call require(maxval(abs(h(:,col)-1.0_wp))==0.0_wp,'unselected Hessian columns preserved')
         call require(maxval(abs(dip(:,col)-2.0_wp))==0.0_wp,'unselected dipole columns preserved')
      endif
   enddo
   call require(maxval(abs(chk%wfn%p-referenceP))==0.0_wp,'reference density preserved')
   call require(maxval(abs(chk%wfn%qsh-referenceQ))==0.0_wp,'reference charges preserved')
   print *, 'PASS partial additive analytic calculator Hessian'
   h=1.0_wp; dip=2.0_wp; pol=3.0_wp
   call calc%hessian(env,mol,chk,[2],5.0e-4_wp,h,dip,pol)
   do col=1,9
      if(col>=4.and.col<=6) then
         call require(maxval(abs(h(:,col)-1.0_wp-expected(:,col)))<2.0e-6_wp,'polarizability numerical fallback')
         call require(maxval(abs(pol(:,col)))==0.0_wp,'selected polarizability output')
      else
         call require(maxval(abs(h(:,col)-1.0_wp))==0.0_wp,'fallback unselected Hessian preserved')
         call require(maxval(abs(pol(:,col)-3.0_wp))==0.0_wp,'fallback unselected polarizability preserved')
      endif
   enddo
   print *, 'PASS requested polarizability numerical fallback'
contains
   subroutine require(condition,label)
      logical,intent(in) :: condition
      character(len=*),intent(in) :: label
      if(.not.condition) then
         print *, 'FAIL ',label
         error stop 1
      endif
   end subroutine
end program
