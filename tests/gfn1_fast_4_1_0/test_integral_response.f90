program test_integral_response
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : evtoau, aatoau
   use xtb_type_molecule, only : TMolecule, init
   use xtb_type_environment, only : TEnvironment, init
   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator
   use xtb_xtb_hamiltonian, only : getSelfEnergy, build_SH0_GFN1, build_dSH0_GFN1_noreset
   use xtb_xtb_hessian_integrals, only : buildGFN1IntegralResponse
   use xtb_disp_coordinationnumber, only : cnType, getCoordinationNumber, getCoordinationNumberHessian
   use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
   implicit none
   real(wp) :: xyz(3,8),rotation(3,3),rz(3,3),ry(3,3)
   call mctc_init('test_integral_response',10,.true.)
   rz=0.0_wp; ry=0.0_wp
   rz(1,1)=cos(0.37_wp); rz(2,2)=rz(1,1); rz(1,2)=-sin(0.37_wp); rz(2,1)=-rz(1,2); rz(3,3)=1.0_wp
   ry(1,1)=cos(0.61_wp); ry(3,3)=ry(1,1); ry(1,3)=sin(0.61_wp); ry(3,1)=-ry(1,3); ry(2,2)=1.0_wp
   rotation=matmul(ry,rz)
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]
   xyz(:,2)=[0.758602_wp,0.0_wp,0.504284_wp]
   xyz(:,3)=[-0.758602_wp,0.0_wp,0.504284_wp]
   call test_molecule('water',[8,1,1],aatoau*matmul(rotation,xyz(:,1:3)))
   call test_molecule('water nodal overlap',[8,1,1],aatoau*xyz(:,1:3),nodal=.true.)
   xyz(:,1)=[-1.15_wp,0.0_wp,0.0_wp]; xyz(:,2)=[1.15_wp,0.1_wp,0.05_wp]
   xyz(:,3)=[-1.65_wp,1.15_wp,0.1_wp]; xyz(:,4)=[-1.65_wp,-0.6_wp,1.0_wp]
   xyz(:,5)=[-1.7_wp,-0.5_wp,-1.05_wp]; xyz(:,6)=[1.7_wp,1.2_wp,0.15_wp]
   xyz(:,7)=[1.65_wp,-0.5_wp,1.1_wp]; xyz(:,8)=[1.7_wp,-0.4_wp,-1.0_wp]
   call test_molecule('disilane',[14,14,1,1,1,1,1,1],aatoau*matmul(rotation,xyz))
contains
   subroutine test_molecule(label,at,xyz,nodal)
      character(len=*), intent(in) :: label
      integer, intent(in) :: at(:)
      real(wp), intent(in) :: xyz(:,:)
      logical, intent(in), optional :: nodal
      type(TEnvironment) :: env
      type(TMolecule) :: mol,moved
      type(TxTBCalculator) :: calc
      real(wp), allocatable :: cn(:),dcn(:,:,:),hcn(:,:,:,:,:),se(:,:),dse(:,:)
      real(wp), allocatable :: p(:,:),w(:,:),ds(:,:,:),dh(:,:,:),hh(:,:),g(:,:),oldg(:,:)
      real(wp), allocatable :: sp(:,:),sm(:,:),hp(:,:),hm(:,:),gp(:,:),gm(:,:),fd(:,:)
      real(wp), allocatable :: translation(:)
      real(wp) :: step,err_s,err_h,err_hh,err_g,previous_error
      integer :: n,nao,n3,i,j,a,atom,coord,level
      logical :: ok,failed,use_gradient_fd
      use_gradient_fd=.true.
      if(present(nodal)) use_gradient_fd=.not.nodal
      call init(env); call init(mol,at,xyz)
      call newXTBCalculator(env,mol,calc,method=1)
      call env%check(failed)
      call require(.not.failed,'GFN1 basis/parameter construction')
      n=mol%n; nao=calc%basis%nao; n3=3*n
      allocate(cn(n),dcn(3,n,n),hcn(3,n,3,n,n),se(5,n),dse(5,n))
      allocate(p(nao,nao),w(nao,nao),ds(nao,nao,n3),dh(nao,nao,n3),hh(n3,n3))
      allocate(g(3,n),oldg(3,n),sp(nao,nao),sm(nao,nao),hp(nao,nao),hm(nao,nao))
      allocate(gp(3,n),gm(3,n),fd(n3,n3),translation(n3))
      ! Independent symmetric test densities probe every AO matrix entry.
      ! They are held fixed, isolating geometric integrals from electronic response.
      do i=1,nao
         do j=1,nao
            p(j,i)=0.01_wp*cos(real(i+j,wp)); w(j,i)=0.003_wp*sin(real(i+j,wp))
         enddo
         p(i,i)=p(i,i)+0.3_wp
      enddo
      call getCoordinationNumberHessian(mol,cnType%exp,cn,dcn,hcn,ok)
      call require(ok,'GFN1 CN first/second derivatives')
      call getSelfEnergy(calc%xtbData%hamiltonian,calc%xtbData%nShell,at,cn=cn,selfEnergy=se,dSEdcn=dse)
      call buildGFN1IntegralResponse(mol,calc%basis,calc%xtbData%hamiltonian,calc%xtbData%nShell, &
         & 25.0_wp,se,dse,dcn,hcn,p,w,ds,dh,hh,ok)
      call require(ok,'GFN1 molecular integral response')
      do coord=1,n3
         a=mod(coord-1,3)+1; atom=(coord-1)/3+1
         g(a,atom)=sum(p*dh(:,:,coord))-sum(w*ds(:,:,coord))
      enddo
      call production(mol,calc,p,w,sp,hp,oldg)
      err_g=maxval(abs(g-oldg))
      print *, label,' NAO=',nao,' fixed-density Gradient error=',err_g
      ! The direct-prefactor kernel also retains derivatives at exact overlap
      ! nodes. Use an independent scalar-energy Hessian oracle for this case.
      call require(err_g<1.0e-12_wp,'new integral Gradient vs direct-prefactor Gradient, including nodes')
      call require(maxval(abs(hh-transpose(hh)))<1.0e-12_wp,'unsymmetrized geometric Hessian reciprocity')
      do a=1,3
         translation=0.0_wp; translation(a:n3:3)=1.0_wp
         call require(maxval(abs(matmul(hh,translation)))<1.0e-12_wp,'three translational null directions')
      enddo
      previous_error=huge(1.0_wp)
      do level=1,2
         step=4.0e-4_wp/real(2**(level-1),wp)
         err_s=0.0_wp; err_h=0.0_wp
         do coord=1,n3
            a=mod(coord-1,3)+1; atom=(coord-1)/3+1
            call moved%copy(mol); moved%xyz(a,atom)=mol%xyz(a,atom)+step
            call production(moved,calc,p,w,sp,hp,gp)
            moved%xyz(a,atom)=mol%xyz(a,atom)-step
            call production(moved,calc,p,w,sm,hm,gm)
            err_s=max(err_s,maxval(abs((sp-sm)/(2.0_wp*step)-ds(:,:,coord))))
            err_h=max(err_h,maxval(abs((hp-hm)/(2.0_wp*step)-dh(:,:,coord))))
            fd(:,coord)=reshape((gp-gm)/(2.0_wp*step),[n3])
         enddo
         if(use_gradient_fd) then
            err_hh=maxval(abs(fd-hh))
         else
            call energy_hessian(mol,calc,p,w,2.5_wp*step,fd)
            print *, 'independent fixed-density energy Hessian step=',2.5_wp*step
            err_hh=maxval(abs(fd-hh))
         endif
         print *, label,' step=',step,' dS=',err_s,' dH0=',err_h,' geometric Hessian=',err_hh
         call require(err_s<1.0e-7_wp.and.err_h<1.0e-7_wp,'all-coordinate first integral derivatives')
         call require(err_hh<1.0e-7_wp,'all-coordinate Hessian vs independent production finite difference')
         if(level==2) then
            call require(err_hh<0.35_wp*previous_error,'second-order finite-difference convergence')
            call require(err_hh<1.0e-8_wp,'smaller-step molecular geometric Hessian tolerance')
         endif
         previous_error=err_hh
      enddo
      call moved%copy(mol); moved%xyz(1,1)=ieee_value(0.0_wp,ieee_quiet_nan)
      call buildGFN1IntegralResponse(moved,calc%basis,calc%xtbData%hamiltonian,calc%xtbData%nShell, &
         & 25.0_wp,se,dse,dcn,hcn,p,w,ds,dh,hh,ok)
      call require(.not.ok,'nonfinite coordinates rejected without floating-point trap')
      p(1,2)=p(1,2)+0.1_wp
      call buildGFN1IntegralResponse(mol,calc%basis,calc%xtbData%hamiltonian,calc%xtbData%nShell, &
         & 25.0_wp,se,dse,dcn,hcn,p,w,ds,dh,hh,ok)
      call require(.not.ok,'nonsymmetric density rejected')
      p(1,2)=p(2,1)
      call buildGFN1IntegralResponse(mol,calc%basis,calc%xtbData%hamiltonian,calc%xtbData%nShell, &
         & 25.0_wp,se,dse,dcn,hcn,p,w,ds(:,:,1:n3-1),dh,hh,ok)
      call require(.not.ok,'inconsistent derivative shape rejected')
      call moved%copy(mol)
      moved%npbc=3
      call buildGFN1IntegralResponse(moved,calc%basis,calc%xtbData%hamiltonian,calc%xtbData%nShell, &
         & 25.0_wp,se,dse,dcn,hcn,p,w,ds,dh,hh,ok)
      call require(.not.ok,'periodic integral response defers to fallback')
      print *, 'PASS ',label,' molecular integral response'
   end subroutine

   subroutine energy_hessian(mol,calc,p,w,step,fd)
      type(TMolecule), intent(in) :: mol
      type(TxTBCalculator), intent(in) :: calc
      real(wp), intent(in) :: p(:,:),w(:,:),step
      real(wp), intent(out) :: fd(:,:)
      type(TMolecule) :: moved
      real(wp) :: energies(2,2),s(size(p,1),size(p,1)),h(size(p,1),size(p,1)),g(3,mol%n)
      integer :: ca,cb,ia,ib,aa,ab,signa,signb
      do ca=1,3*mol%n
         aa=mod(ca-1,3)+1; ia=(ca-1)/3+1
         do cb=1,ca
            ab=mod(cb-1,3)+1; ib=(cb-1)/3+1
            do signa=1,2
               do signb=1,2
                  call moved%copy(mol)
                  moved%xyz(aa,ia)=moved%xyz(aa,ia)+real(2*signa-3,wp)*step
                  moved%xyz(ab,ib)=moved%xyz(ab,ib)+real(2*signb-3,wp)*step
                  call production(moved,calc,p,w,s,h,g)
                  energies(signa,signb)=sum(p*h)-sum(w*s)
               enddo
            enddo
            fd(ca,cb)=(energies(2,2)-energies(2,1)-energies(1,2)+energies(1,1))/(4.0_wp*step**2)
            fd(cb,ca)=fd(ca,cb)
         enddo
      enddo
   end subroutine

   subroutine production(mol,calc,p,w,s,h,g)
      type(TMolecule), intent(in) :: mol
      type(TxTBCalculator), intent(in) :: calc
      real(wp), intent(in) :: p(:,:),w(:,:)
      real(wp), intent(out) :: s(:,:),h(:,:),g(:,:)
      real(wp) :: trans(3,1),cn(mol%n),dcn(3,mol%n,mol%n),strain(3,3,mol%n)
      real(wp) :: se(5,mol%n),dse(5,mol%n),ves(5,mol%n),dhdcn(mol%n),sigma(3,3)
      real(wp) :: packed(size(s,1)*(size(s,1)+1)/2)
      integer :: i,j,ij,nao
      nao=size(s,1); trans=0.0_wp
      call getCoordinationNumber(mol,trans,40.0_wp,cnType%exp,cn,dcn,strain)
      call getSelfEnergy(calc%xtbData%hamiltonian,calc%xtbData%nShell,mol%at,cn=cn,selfEnergy=se,dSEdcn=dse)
      call build_SH0_GFN1(calc%xtbData%nShell,calc%xtbData%hamiltonian,mol%n,mol%at, &
         & calc%basis%nbf,nao,mol%xyz,trans,se,25.0_wp,calc%basis%caoshell,calc%basis%saoshell, &
         & calc%basis%nprim,calc%basis%primcount,calc%basis%alp,calc%basis%cont,s,packed)
      ij=0
      do i=1,nao
         do j=1,i
            ij=ij+1; h(j,i)=packed(ij)*evtoau; h(i,j)=h(j,i)
         enddo
      enddo
      ves=0.0_wp; dhdcn=0.0_wp; g=0.0_wp; sigma=0.0_wp
      call build_dSH0_GFN1_noreset(calc%xtbData%nShell,calc%xtbData%hamiltonian,se,dse,25.0_wp, &
         & mol%n,nao,calc%basis%nbf,mol%at,mol%xyz,calc%basis%caoshell,calc%basis%saoshell, &
         & calc%basis%nprim,calc%basis%primcount,calc%basis%alp,calc%basis%cont, &
         & S=s,p=p,Pew=w,ves=ves,dhdcn=dhdcn,g=g,sigma=sigma,H0packed=packed,cartContract=.false., &
         & directPrefactor=.true.)
      g=g+reshape(matmul(reshape(dcn,[3*mol%n,mol%n]),dhdcn),shape(g))
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
