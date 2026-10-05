program test_coulomb_geometry
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_molecule, only : TMolecule,init
   use xtb_type_environment, only : TEnvironment,init
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator
   use xtb_coulomb_klopmanohno, only : TKlopmanOhno,init,gamAverage
   use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
   implicit none
   real(wp) :: xyz(3,4)
   call mctc_init('test_coulomb_geometry',10,.true.)
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]; xyz(:,2)=[0.758602_wp,0.13_wp,0.504284_wp]
   xyz(:,3)=[-0.758602_wp,-0.18_wp,0.62_wp]
   call test_molecule('water',[8,1,1],aatoau*xyz(:,1:3))
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]; xyz(:,2)=[2.3_wp,0.1_wp,0.05_wp]
   xyz(:,3)=[-0.5_wp,1.15_wp,0.1_wp]; xyz(:,4)=[2.85_wp,-0.5_wp,1.1_wp]
   call test_molecule('Si2H2',[14,14,1,1],aatoau*xyz)
contains
   subroutine test_molecule(label,at,xyz)
      character(len=*), intent(in) :: label
      integer, intent(in) :: at(:)
      real(wp), intent(in) :: xyz(:,:)
      type(TEnvironment) :: env
      type(TMolecule) :: mol,moved
      type(TxTBCalculator) :: calc
      type(TKlopmanOhno) :: coulomb,missing
      integer, allocatable :: idnum(:)
      real(wp), allocatable :: q(:),ds(:,:),hh(:,:),g(:,:),oldg(:,:),gp(:,:),gm(:,:),jp(:,:),jm(:,:),fd(:,:)
      real(wp), allocatable :: djdr(:,:,:),trace(:,:),strain(:,:,:),translation(:)
      real(wp) :: exponents(3),step,err_shift,err_hh,previous_error,err_g
      integer :: n,ns,n3,i,a,atom,coord,level,model
      logical :: ok,failed
      call init(env); call init(mol,at,xyz); call newXTBCalculator(env,mol,calc,method=1)
      call env%check(failed); call require(.not.failed,'GFN1 parameter construction')
      n=mol%n; ns=calc%basis%nshell; n3=3*n
      allocate(idnum(maxval(mol%id)))
      do i=1,n; idnum(mol%id(i))=mol%at(i); enddo
      allocate(q(ns),ds(ns,n3),hh(n3,n3),g(3,n),oldg(3,n),gp(3,n),gm(3,n),jp(ns,ns),jm(ns,ns),fd(n3,n3))
      allocate(djdr(3,n,ns),trace(3,ns),strain(3,3,ns),translation(n3))
      do i=1,ns; q(i)=0.2_wp*sin(real(i,wp)); enddo
      q=q-sum(q)/real(ns,wp)
      exponents=[calc%xtbData%coulomb%gExp,1.3_wp,3.0_wp]
      do model=1,3
         call init(coulomb,env,mol,gamAverage%harmonic,calc%xtbData%coulomb%shellHardness, &
            & exponents(model),num=idnum,nshell=calc%xtbData%nShell)
         call env%check(failed); call require(.not.failed,'molecular Klopman-Ohno construction')
         call coulomb%getMolecularHessianResponse(mol,q,ds,hh,ok,g)
         call require(ok,'molecular Coulomb Hessian response')
         call coulomb%getCoulombDerivs(mol,q,djdr,trace,strain)
         oldg=reshape(matmul(reshape(djdr,[n3,ns]),q),[3,n])
         err_g=maxval(abs(oldg-g))
         print *, label,' NSH=',ns,' gExp=',exponents(model),' Gradient error=',err_g
         call require(err_g<1.0e-13_wp,'Coulomb Gradient vs existing production derivative')
         call require(maxval(abs(reshape(matmul(transpose(ds),q),[3,n])-2.0_wp*g))<1.0e-13_wp, &
            & 'fixed-charge scalar energy gradient contraction')
         call require(maxval(abs(hh-transpose(hh)))<1.0e-13_wp,'Coulomb Hessian reciprocity')
         do a=1,3
            translation=0.0_wp; translation(a:n3:3)=1.0_wp
            call require(maxval(abs(matmul(hh,translation)))<1.0e-13_wp,'Coulomb Hessian translation')
            call require(maxval(abs(matmul(ds,translation)))<1.0e-13_wp,'Coulomb potential response translation')
         enddo
         previous_error=huge(1.0_wp)
         do level=1,2
            step=4.0e-4_wp/real(2**(level-1),wp); err_shift=0.0_wp
            do coord=1,n3
               a=mod(coord-1,3)+1; atom=(coord-1)/3+1
               call moved%copy(mol); moved%xyz(a,atom)=mol%xyz(a,atom)+step
               call coulomb%getCoulombMatrix(moved,jp)
               call coulomb%getCoulombDerivs(moved,q,djdr,trace,strain)
               gp=reshape(matmul(reshape(djdr,[n3,ns]),q),[3,n])
               moved%xyz(a,atom)=mol%xyz(a,atom)-step
               call coulomb%getCoulombMatrix(moved,jm)
               call coulomb%getCoulombDerivs(moved,q,djdr,trace,strain)
               gm=reshape(matmul(reshape(djdr,[n3,ns]),q),[3,n])
               err_shift=max(err_shift,maxval(abs(matmul((jp-jm)/(2.0_wp*step),q)-ds(:,coord))))
               fd(:,coord)=reshape((gp-gm)/(2.0_wp*step),[n3])
            enddo
            err_hh=maxval(abs(fd-hh))
            print *, label,' step=',step,' d(Jq)=',err_shift,' Coulomb Hessian=',err_hh
            call require(err_shift<1.0e-9_wp.and.err_hh<1.0e-9_wp,'all-coordinate Coulomb finite differences')
            if(level==2) then
               call require(err_hh<0.35_wp*previous_error,'Coulomb second-order FD convergence')
               call require(err_hh<1.0e-10_wp,'smaller-step Coulomb Hessian tolerance')
            endif
            previous_error=err_hh
         enddo
      enddo
      ! Inactive/not constructed objects and unsupported geometries use fallback.
      call missing%getMolecularHessianResponse(mol,q,ds,hh,ok)
      call require(.not.ok,'uninitialized Coulomb response rejected')
      call moved%copy(mol); moved%npbc=3
      call coulomb%getMolecularHessianResponse(moved,q,ds,hh,ok)
      call require(.not.ok,'periodic Coulomb response defers to fallback')
      call moved%copy(mol); moved%xyz(:,2)=mol%xyz(:,1)
      call coulomb%getMolecularHessianResponse(moved,q,ds,hh,ok)
      call require(.not.ok,'coincident centers reject singular radial derivative')
      q(1)=ieee_value(0.0_wp,ieee_quiet_nan)
      call coulomb%getMolecularHessianResponse(mol,q,ds,hh,ok)
      call require(.not.ok,'nonfinite charges rejected without floating-point trap')
      print *, 'PASS ',label,' Coulomb geometry response'
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
