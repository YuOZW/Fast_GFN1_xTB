! Original D3 Gradient on both sides of its real 60-bohr hard cutoff.
! No SCC or surrogate pair formula is used to differentiate the reference.
program probe_cutoff_response
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator
   use xtb_disp_coordinationnumber, only : cnType,getCoordinationNumberHessian
   use xtb_disp_dftd3, only : d3_gradient,d3MolecularHessian
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: mol,moved
   type(TxTBCalculator) :: calc
   real(wp), parameter :: distances(5)=[59.9_wp,59.9998_wp,60.0_wp,60.0002_wp,60.1_wp],step=5.0e-4_wp
   real(wp) :: xyz(3,2),g(3,2),ag(3,2),gp(3,2),gm(3,2),h(6,6),fd(6,6)
   real(wp) :: cn(2),dcn(3,2,2),hcn(3,2,3,2,2),energy,error
   integer :: k,x,a,atom
   logical :: ok,good,failed,diagnose
   character(len=32) :: mode
   call mctc_init('probe_cutoff_response',10,.true.);call init(env)
   call get_command_argument(1,mode);diagnose=trim(mode)=='diagnose'
   do k=1,size(distances)
      xyz=0.0_wp;xyz(1,2)=distances(k)
      call init(mol,[47,47],xyz)
      call newXTBCalculator(env,mol,calc,method=1,accuracy=1.0e-7_wp)
      call env%check(failed);if(failed) error stop 'cutoff parameter construction'
      call getCoordinationNumberHessian(mol,cnType%exp,cn,dcn,hcn,good)
      if(.not.good) error stop 'smooth CN reference'
      call referenceGradient(mol,g)
      call d3MolecularHessian(mol,calc%xtbData%dispersion%dpar,4.0_wp,60.0_wp, &
         & cn,dcn,hcn,energy,ag,h,ok,guardRadius=step)
      do x=1,6
         a=mod(x-1,3)+1;atom=(x-1)/3+1
         call moved%copy(mol);moved%xyz(a,atom)=mol%xyz(a,atom)+step
         call referenceGradient(moved,gp)
         moved%xyz(a,atom)=mol%xyz(a,atom)-step
         call referenceGradient(moved,gm)
         fd(:,x)=reshape(gp-gm,[6])/(2.0_wp*step)
      enddo
      error=maxval(abs(h-fd))
      print *, 'D3 distance/accepted/original G norm/G difference/H FD difference:', &
         & distances(k),ok,norm2(g),maxval(abs(g-ag)),error
      flush(6)
      if(k>=2.and.k<=4) then
         if(maxval(abs(fd))<=1.0e-10_wp) error stop 'actual cutoff Gradient jump missing'
         if(.not.diagnose.and.ok) error stop 'D3 cutoff point must reject analytic response'
      else
         if(.not.ok) error stop 'smooth D3 response rejected'
         if(maxval(abs(g-ag))>1.0e-12_wp) error stop 'smooth original D3 Gradient parity'
         if(error>1.0e-9_wp) error stop 'smooth original D3 Gradient FD parity'
      endif
   enddo
   if(.not.diagnose) print *, 'PASS smooth D3 responses and actual hard-cutoff rejection'
contains
   subroutine referenceGradient(structure,gradient)
      type(TMolecule), intent(in) :: structure
      real(wp), intent(out) :: gradient(3,2)
      real(wp) :: c(2),dc(3,2,2),hc(3,2,3,2,2),trans(3,1),strain(3,3,2),sigma(3,3),e
      logical :: valid
      call getCoordinationNumberHessian(structure,cnType%exp,c,dc,hc,valid)
      if(.not.valid) error stop 'displaced CN reference'
      trans=0.0_wp;strain=0.0_wp;sigma=0.0_wp;e=0.0_wp;gradient=0.0_wp
      call d3_gradient(structure,trans,calc%xtbData%dispersion%dpar,4.0_wp,60.0_wp, &
         & c,dc,strain,e,gradient,sigma)
   end subroutine
end program
