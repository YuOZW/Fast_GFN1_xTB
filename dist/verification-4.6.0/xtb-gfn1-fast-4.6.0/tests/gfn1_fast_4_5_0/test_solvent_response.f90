! Validate against the original TBorn functional and Gradient, not a
! finite difference of the newly assembled response itself.
program test_solvent_response
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : TSolvModel,init,newBornModel
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_cm5, only : calc_cm5
   use xtb_solv_response, only : buildGFN1SolventResponse
   implicit none
   type(TEnvironment) :: env
   type(TSolvModel) :: model
   type(TBorn) :: born
   integer :: atoms(3),modelIndex,kernelIndex,chargeIndex,level,x,a,iat,j
   real(wp) :: xyz(3,3),moved(3,3),q(3),qp(3),qm(3),energy,g(3,3),h(9,9),v(3),dv(3,9),amat(3,3)
   real(wp) :: referenceE,referenceG(3,3),referenceV(3),ep,em,gp(3,3),gm(3,3),vp(3),vm(3)
   real(wp) :: hfd(9,9),vfd(3,9),gfd(3,3),afd(3,3),step,he,ve,ge,previousH,previousV,translation(9)
   logical :: ok,failed
   call mctc_init('test_solvent_response',10,.true.);call init(env)
   atoms=[8,1,1];xyz(:,1)=0.0_wp
   xyz(:,2)=aatoau*[0.758602_wp,0.0_wp,0.504284_wp]
   xyz(:,3)=aatoau*[-0.3_wp,0.67_wp,0.62_wp]
   do modelIndex=1,2
      do kernelIndex=1,2
         call init(model,env,TSolvInput(solvent='water',alpb=modelIndex==2,kernel=kernelIndex),1)
         call env%check(failed);if(failed) error stop 'solvent parameters'
         call newBornModel(model,env,born,atoms)
         do chargeIndex=1,2
            q=[-0.5_wp,0.27_wp,0.23_wp]
            if(chargeIndex==2) q=[-0.5_wp,0.32_wp,0.33_wp]
            call parent(xyz,q,referenceE,referenceG,referenceV)
            call buildGFN1SolventResponse(born,atoms,xyz,q,energy,g,h,v,dv,amat,ok)
            if(.not.ok) error stop 'smooth solvent reference rejected'
            if(abs(energy-referenceE)>1.0e-12_wp) error stop 'original solvent energy parity'
            if(maxval(abs(g-referenceG))>1.0e-11_wp) error stop 'original solvent Gradient parity'
            if(maxval(abs(v-referenceV))>1.0e-12_wp) error stop 'original solvent potential parity'
            if(maxval(abs(amat-born%bornMat))>1.0e-12_wp) error stop 'original charge kernel parity'
            if(maxval(abs(h-transpose(h)))>1.0e-12_wp) error stop 'raw solvent reciprocity'
            do a=1,3
               translation=0.0_wp;translation(a:9:3)=1.0_wp
               if(maxval(abs(matmul(h,translation)))>1.0e-11_wp) error stop 'solvent translation'
               if(maxval(abs(matmul(dv,translation)))>1.0e-11_wp) error stop 'potential translation'
            enddo
            previousH=huge(1.0_wp);previousV=huge(1.0_wp)
            do level=1,3
               step=1.0e-4_wp/real(2**(level-1),wp)
               do x=1,9
                  a=mod(x-1,3)+1;iat=(x-1)/3+1
                  moved=xyz;moved(a,iat)=moved(a,iat)+step
                  call parent(moved,q,ep,gp,vp)
                  moved(a,iat)=xyz(a,iat)-step
                  call parent(moved,q,em,gm,vm)
                  hfd(:,x)=reshape(gp-gm,[9])/(2.0_wp*step)
                  vfd(:,x)=(vp-vm)/(2.0_wp*step)
                  gfd(a,iat)=(ep-em)/(2.0_wp*step)
               enddo
               he=maxval(abs(hfd-h));ve=maxval(abs(vfd-dv));ge=maxval(abs(gfd-g))
               print *, 'model/kernel/charge/step/H/V/G:',modelIndex,kernelIndex,chargeIndex,step,he,ve,ge
               if(he>2.0e-7_wp.or.ve>2.0e-7_wp.or.ge>1.0e-8_wp) error stop 'solvent geometry derivatives'
               if(level>1) then
                  if(he>0.4_wp*previousH+2.0e-10_wp) error stop 'H quadratic convergence'
                  if(ve>0.4_wp*previousV+2.0e-10_wp) error stop 'potential quadratic convergence'
               endif
               previousH=he;previousV=ve
            enddo
            step=1.0e-5_wp
            do j=1,3
               qp=q;qm=q;qp(j)=qp(j)+step;qm(j)=qm(j)-step
               call parent(xyz,qp,ep,gp,vp);call parent(xyz,qm,em,gm,vm)
               afd(:,j)=(vp-vm)/(2.0_wp*step)
               ! Maxwell identity links charge and coordinate derivatives.
               if(maxval(abs(reshape(gp-gm,[9])/(2.0_wp*step)-dv(j,:)))>1.0e-10_wp) &
                  & error stop 'charge/geometry mixed response'
               if(abs((ep-em)/(2.0_wp*step)-v(j))>1.0e-11_wp) error stop 'charge energy derivative'
            enddo
            if(maxval(abs(afd-amat))>1.0e-11_wp) error stop 'charge potential derivative'
         enddo
      enddo
   enddo
   ! Guarded discontinuity is an established feature of the original SASA.
   call init(model,env,TSolvInput(solvent='water',alpb=.false.,kernel=gbKernel%still),1)
   call newBornModel(model,env,born,atoms)
   moved=xyz;moved(2,3)=aatoau*0.65_wp
   call born%update(env,atoms,moved)
   call buildGFN1SolventResponse(born,atoms,moved,q,energy,g,h,v,dv,amat,ok,guardRadius=2.0e-4_wp)
   if(ok) error stop 'screened surface boundary accepted'
   call born%update(env,atoms,xyz);born%lsalt=.true.
   call buildGFN1SolventResponse(born,atoms,xyz,q,energy,g,h,v,dv,amat,ok)
   if(ok) error stop 'salt response accepted without derivative'
   born%lsalt=.false.
   call buildGFN1SolventResponse(born,atoms,xyz,q,energy,g,h(:8,:),v,dv,amat,ok)
   if(ok) error stop 'wrong Hessian shape accepted'
   print *, 'PASS fixed-bare-charge GBSA/ALPB Still/P16 response, CM5, mixed derivatives and guards'
contains
   subroutine parent(coordinates,bareCharges,e,gradient,potential)
      real(wp), intent(in) :: coordinates(3,3),bareCharges(3)
      real(wp), intent(out) :: e,gradient(3,3),potential(3)
      real(wp) :: cm(3),dc(3,3,3),charges(3)
      integer :: i
      call born%update(env,atoms,coordinates)
      call calc_cm5(3,atoms,coordinates,cm,dc);charges=bareCharges+cm
      call born%getEnergy(env,charges,bareCharges,e)
      gradient=0.0_wp
      call born%addGradient(env,atoms,coordinates,charges,bareCharges,gradient)
      potential=matmul(born%bornMat,charges)
      do i=1,3
         gradient=gradient+dc(:,:,i)*potential(i)
      enddo
      call env%check(failed);if(failed) error stop 'original solvent evaluation'
   end subroutine
end program
