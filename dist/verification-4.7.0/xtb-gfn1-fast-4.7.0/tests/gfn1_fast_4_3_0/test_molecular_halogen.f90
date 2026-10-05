program test_molecular_halogen
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator,newWavefunction
   use xtb_xtb_hessian_response, only : buildGFN1GasHessian
   implicit none
   integer :: donors(4),i
   real(wp) :: lengths(4)
   call mctc_init('test_molecular_halogen',10,.true.)
   donors=[17,35,53,85]; lengths=[1.78_wp,1.94_wp,2.15_wp,2.25_wp]
   do i=1,4
      call check_molecule(donors(i),lengths(i),300.0_wp,0.0_wp,0)
   enddo
   call check_molecule(35,1.94_wp,0.0_wp,0.0_wp,0)
   call check_molecule(35,1.94_wp,30000.0_wp,0.0_wp,0)
   call check_molecule(53,2.15_wp,3000.0_wp,1.0_wp,1)
contains
   subroutine check_molecule(donor,bond,temperature,charge,spin)
      integer, intent(in) :: donor,spin
      real(wp), intent(in) :: bond,temperature,charge
      type(TEnvironment) :: env
      type(TMolecule) :: mol,moved
      type(TxTBCalculator) :: calc
      type(TRestart) :: center,plus,minus
      type(scc_results) :: results
      real(wp) :: xyz(3,8),gradient(3,8),analyticG(3,8),gp(3,8),gm(3,8)
      real(wp) :: hessian(24,24),fd(24,24),translation(24),sigma(3,3)
      real(wp) :: energy,gap,step,error,previous,exb,gradientError
      integer :: coord,a,atom,level,terms
      logical :: ok,failed
      xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]
      xyz(:,2)=[-bond,0.04_wp,-0.05_wp]
      xyz(:,3)=xyz(:,2)+[-0.36_wp,1.026_wp,0.03_wp]
      xyz(:,4)=xyz(:,2)+[-0.37_wp,-0.51_wp,0.89_wp]
      xyz(:,5)=xyz(:,2)+[-0.35_wp,-0.52_wp,-0.88_wp]
      xyz(:,6)=[3.2_wp,0.17_wp,-0.13_wp]
      xyz(:,7)=xyz(:,6)+[0.57_wp,0.74_wp,0.12_wp]
      xyz(:,8)=xyz(:,6)+[0.58_wp,-0.75_wp,-0.09_wp]
      call init(env)
      call init(mol,[donor,6,1,1,1,8,1,1],aatoau*xyz,chrg=charge,uhf=spin)
      call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
      call env%check(failed); call require(.not.failed,'halogen calculator construction')
      calc%etemp=temperature; calc%maxiter=500
      call newWavefunction(env,mol,calc,center)
      call env%check(failed); call require(.not.failed,'halogen wavefunction initialization')
      call calc%singlepoint(env,mol,center,0,.false.,energy,gradient,sigma,gap,results)
      call env%check(failed); call require(.not.failed.and.results%converged,'reference SCC convergence')
      exb=results%e_xb
      call buildGFN1GasHessian(env,mol,calc%basis,calc%xtbData,center%wfn,temperature,calc%accuracy, &
         & analyticG,hessian,ok,terms)
      call require(ok,'gas Hessian including halogen correction')
      if(donor==35.or.donor==53.or.donor==85) then
         call require(abs(exb)>1.0e-8_wp.and.terms==1,'active production halogen correction')
      else
         call require(abs(exb)<1.0e-14_wp.and.terms==0,'default zero Cl correction')
      endif
      gradientError=maxval(abs(gradient-analyticG))
      call require(gradientError<2.0e-9_wp,'complete Gradient vs production singlepoint')
      call require(maxval(abs(hessian-transpose(hessian)))<1.0e-7_wp,'raw Hessian reciprocity')
      do a=1,3
         translation=0.0_wp; translation(a:24:3)=1.0_wp
         call require(maxval(abs(matmul(hessian,translation)))<1.0e-9_wp,'raw Hessian translation')
      enddo
      previous=huge(1.0_wp)
      do level=1,2
         step=4.0e-4_wp/real(2**(level-1),wp)
         do coord=1,24
            a=mod(coord-1,3)+1; atom=(coord-1)/3+1
            call moved%copy(mol); moved%xyz(a,atom)=mol%xyz(a,atom)+step
            call plus%copy(center)
            call calc%singlepoint(env,moved,plus,0,.true.,energy,gp,sigma,gap,results)
            call env%check(failed); call require(.not.failed.and.results%converged,'positive displacement SCC')
            moved%xyz(a,atom)=mol%xyz(a,atom)-step
            call minus%copy(center)
            call calc%singlepoint(env,moved,minus,0,.true.,energy,gm,sigma,gap,results)
            call env%check(failed); call require(.not.failed.and.results%converged,'negative displacement SCC')
            fd(:,coord)=reshape((gp-gm)/(2.0_wp*step),[24])
         enddo
         error=maxval(abs(hessian-fd))
         print *, 'halogen=',donor,' temperature=',temperature,' charge=',charge,' step=',step,' H error=',error
         call require(error<1.0e-6_wp,'all-coordinate production Gradient differences')
         if(level==2) call require(error<0.4_wp*previous+1.0e-8_wp,'quadratic finite-difference convergence')
         previous=error
      enddo
      call require(error<3.0e-7_wp,'smaller-step molecular Hessian tolerance')
      print *, 'PASS molecular halogen=',donor,' temperature=',temperature,' charge=',charge, &
         & ' XB energy=',exb,' terms=',terms,' Gradient error=',gradientError
   end subroutine
   subroutine require(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if(.not.condition) then
         print *, 'FAIL ',label
         error stop 1
      endif
   end subroutine
end program
