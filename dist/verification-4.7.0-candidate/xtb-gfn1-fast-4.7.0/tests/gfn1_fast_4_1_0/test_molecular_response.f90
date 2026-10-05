program test_molecular_response
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator,newWavefunction
   use xtb_xtb_hessian_response, only : buildGFN1ElectrostaticHessian,buildGFN1GasHessian
   use xtb_xtb_repulsion, only : repulsionEnGrad
   use xtb_disp_dftd3, only : d3_gradient
   use xtb_disp_coordinationnumber, only : getCoordinationNumber,cnType
   implicit none
   real(wp) :: xyz(3,8)
   call mctc_init('test_molecular_response',10,.true.)
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]; xyz(:,2)=[0.758602_wp,0.13_wp,0.504284_wp]
   xyz(:,3)=[-0.758602_wp,-0.18_wp,0.62_wp]
   call test_molecule('water 0K',[8,1,1],aatoau*xyz(:,1:3),0.0_wp,0.0_wp,0)
   call test_molecule('water 30000K',[8,1,1],aatoau*xyz(:,1:3),30000.0_wp,0.0_wp,0)
   call test_molecule('water open shell',[8,1,1],aatoau*xyz(:,1:3),3000.0_wp,1.0_wp,1)
   xyz(:,1)=[-1.15_wp,0.0_wp,0.0_wp]; xyz(:,2)=[1.15_wp,0.1_wp,0.05_wp]
   xyz(:,3)=[-1.65_wp,1.15_wp,0.1_wp]; xyz(:,4)=[-1.65_wp,-0.6_wp,1.0_wp]
   xyz(:,5)=[-1.7_wp,-0.5_wp,-1.05_wp]; xyz(:,6)=[1.7_wp,1.2_wp,0.15_wp]
   xyz(:,7)=[1.65_wp,-0.5_wp,1.1_wp]; xyz(:,8)=[1.7_wp,-0.4_wp,-1.0_wp]
   call test_molecule('disilane 300K',[14,14,1,1,1,1,1,1],aatoau*xyz,300.0_wp,0.0_wp,0)
contains
   subroutine test_molecule(label,at,xyz,temperature,charge,spin)
      character(len=*), intent(in) :: label
      integer, intent(in) :: at(:),spin
      real(wp), intent(in) :: xyz(:,:),temperature,charge
      type(TEnvironment) :: env
      type(TMolecule) :: mol,moved
      type(TxTBCalculator) :: calc
      type(TRestart) :: center,plus,minus
      type(scc_results) :: results
      real(wp), allocatable :: gradient(:,:),analyticGradient(:,:),hessian(:,:),dq(:,:),dp(:,:,:)
      real(wp), allocatable :: gp(:,:),gm(:,:),fd(:,:),translation(:)
      real(wp), allocatable :: fullGradient(:,:),fullHessian(:,:),fullGp(:,:),fullGm(:,:),fullFD(:,:)
      real(wp) :: energy,gap,sigma(3,3),step,errorH,errorQ,errorP,reciprocity,condition,errorFullH
      integer :: n,nao,ns,n3,coord,a,atom,level
      logical :: ok,failed
      call init(env); call init(mol,at,xyz,chrg=charge,uhf=spin)
      call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
      call env%check(failed); call require(.not.failed,'GFN1 molecular calculator construction')
      calc%etemp=temperature; calc%maxiter=500
      call newWavefunction(env,mol,calc,center)
      call env%check(failed); call require(.not.failed,'GFN1 molecular wavefunction initialization')
      n=mol%n; n3=3*n; nao=calc%basis%nao; ns=calc%basis%nshell
      allocate(gradient(3,n),analyticGradient(3,n),hessian(n3,n3),dq(ns,n3),dp(nao,nao,n3))
      allocate(gp(3,n),gm(3,n),fd(n3,n3),translation(n3))
      allocate(fullGradient(3,n),fullHessian(n3,n3),fullGp(3,n),fullGm(3,n),fullFD(n3,n3))
      call calc%singlepoint(env,mol,center,0,.false.,energy,gradient,sigma,gap,results)
      call env%check(failed); call require(.not.failed.and.results%converged,'reference SCC convergence')
      call buildGFN1GasHessian(env,mol,calc%basis,calc%xtbData,center%wfn,temperature,calc%accuracy, &
         & fullGradient,fullHessian,ok)
      call require(ok,'complete gas-phase analytic GFN1 Hessian')
      print *, label,' complete Gradient error=',maxval(abs(gradient-fullGradient))
      call require(maxval(abs(gradient-fullGradient))<1.0e-9_wp,'complete production molecular Gradient')
      call require(maxval(abs(fullHessian-transpose(fullHessian)))<1.0e-7_wp,'complete Hessian reciprocity')
      call subtract_classical(mol,calc,gradient)
      call buildGFN1ElectrostaticHessian(env,mol,calc%basis,calc%xtbData,center%wfn, &
         & temperature,calc%accuracy,analyticGradient,hessian,ok,dq,dp,condition)
      call require(ok,'actual molecular coupled electrostatic Hessian')
      print *, label,' Gradient error=',maxval(abs(gradient-analyticGradient)),' rcond=',condition
      call require(maxval(abs(gradient-analyticGradient))<1.0e-9_wp,'actual production electronic/ES Gradient')
      reciprocity=maxval(abs(hessian-transpose(hessian)))
      print *, label,' unsymmetrized Hessian reciprocity=',reciprocity
      call require(reciprocity<1.0e-7_wp,'self-consistent molecular Hessian reciprocity')
      do a=1,3
         translation=0.0_wp; translation(a:n3:3)=1.0_wp
         call require(maxval(abs(matmul(hessian,translation)))<1.0e-9_wp,'molecular Hessian translation')
         call require(maxval(abs(matmul(dq,translation)))<1.0e-9_wp,'shell-charge response translation')
         call require(maxval(abs(matmul(fullHessian,translation)))<1.0e-9_wp,'complete molecular Hessian translation')
      enddo
      do level=1,2
         step=5.0e-4_wp/real(2**(level-1),wp); errorQ=0.0_wp; errorP=0.0_wp
         do coord=1,n3
            a=mod(coord-1,3)+1; atom=(coord-1)/3+1
            call moved%copy(mol); moved%xyz(a,atom)=mol%xyz(a,atom)+step
            call plus%copy(center)
            call calc%singlepoint(env,moved,plus,0,.true.,energy,gp,sigma,gap,results)
            call env%check(failed); call require(.not.failed.and.results%converged,'positive displacement SCC')
            fullGp=gp
            call subtract_classical(moved,calc,gp)
            moved%xyz(a,atom)=mol%xyz(a,atom)-step
            call minus%copy(center)
            call calc%singlepoint(env,moved,minus,0,.true.,energy,gm,sigma,gap,results)
            call env%check(failed); call require(.not.failed.and.results%converged,'negative displacement SCC')
            fullGm=gm
            call subtract_classical(moved,calc,gm)
            fd(:,coord)=reshape((gp-gm)/(2.0_wp*step),[n3])
            fullFD(:,coord)=reshape((fullGp-fullGm)/(2.0_wp*step),[n3])
            errorQ=max(errorQ,maxval(abs((plus%wfn%qsh-minus%wfn%qsh)/(2.0_wp*step)-dq(:,coord))))
            errorP=max(errorP,maxval(abs((plus%wfn%p-minus%wfn%p)/(2.0_wp*step)-dp(:,:,coord))))
         enddo
         errorH=maxval(abs(fd-hessian))
         errorFullH=maxval(abs(fullFD-fullHessian))
         print *, label,' step=',step,' Hessian error=',errorH,' dq error=',errorQ,' dP error=',errorP
         call require(errorH<2.0e-6_wp,'all-coordinate self-consistent electronic/ES Hessian differences')
         call require(errorQ<2.0e-6_wp.and.errorP<2.0e-6_wp,'actual SCC shell-charge and density response differences')
         print *, label,' complete Hessian error=',errorFullH
         call require(errorFullH<2.0e-6_wp,'complete gas-phase Hessian vs production Gradient differences')
      enddo
      call require(errorH<5.0e-7_wp,'smaller-step self-consistent electrostatic Hessian tolerance')
      call require(errorFullH<5.0e-7_wp,'smaller-step complete gas-phase Hessian tolerance')
      call moved%copy(mol); moved%npbc=3
      call buildGFN1ElectrostaticHessian(env,moved,calc%basis,calc%xtbData,center%wfn, &
         & temperature,calc%accuracy,analyticGradient,fd,ok)
      call require(.not.ok,'periodic molecular response defers to numerical fallback')
      print *, 'PASS ',label,' actual molecular coupled response'
   end subroutine

   ! Independent reference: remove only geometry-dependent classical terms
   ! from the production singlepoint Gradient, retaining the full SCC response.
   ! These fixtures contain no halogens, solvent, external fields or constraints.
   subroutine subtract_classical(mol,calc,gradient)
      type(TMolecule), intent(in) :: mol
      type(TxTBCalculator), intent(in) :: calc
      real(wp), intent(inout) :: gradient(:,:)
      real(wp) :: other(3,mol%n),cn(mol%n),dcn(3,mol%n,mol%n),strain(3,3,mol%n)
      real(wp) :: trans(3,1),sigma(3,3),energy
      trans=0.0_wp; sigma=0.0_wp; other=0.0_wp; energy=0.0_wp
      call repulsionEnGrad(mol,calc%xtbData%repulsion,trans,40.0_wp,energy,other,sigma)
      call getCoordinationNumber(mol,trans,40.0_wp,cnType%exp,cn,dcn,strain)
      call d3_gradient(mol,trans,calc%xtbData%dispersion%dpar,4.0_wp,60.0_wp,cn,dcn,strain,energy,other,sigma)
      gradient=gradient-other
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
