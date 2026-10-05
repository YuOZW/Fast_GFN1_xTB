program test_surface_response
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment,init
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : TSolvModel,init,newBornModel
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_response, only : getSurfaceResponse
   implicit none
   type(TEnvironment) :: env
   type(TSolvModel) :: model
   type(TBorn) :: born
   real(wp) :: xyz(3,3),moved(3,3),surface(3),first(9,3),second(9,9,3),fp(9,3),fm(9,3),fd(9,9,3)
   real(wp) :: step,error,previous,translation(9)
   integer :: modelIndex,level,coordinate,a,atom,i
   logical :: ok,failed
   call mctc_init('test_surface_response',10,.true.);call init(env)
   xyz(:,1)=0.0_wp;xyz(:,2)=aatoau*[0.758602_wp,0.0_wp,0.504284_wp]
   xyz(:,3)=aatoau*[-0.3_wp,0.67_wp,0.62_wp]
   do modelIndex=1,2
      call init(model,env,TSolvInput(solvent='water',alpb=modelIndex==2,kernel=gbKernel%still),1)
      call env%check(failed);if(failed) error stop 'surface parameters'
      call newBornModel(model,env,born,[8,1,1]);call born%update(env,[8,1,1],xyz)
      call getSurfaceResponse(born,xyz,surface,first,second,ok)
      if(.not.ok) error stop 'smooth surface reference rejected'
      do i=1,3
         if(maxval(abs(second(:,:,i)-transpose(second(:,:,i))))>1.0e-12_wp) error stop 'surface reciprocity'
         do a=1,3
            translation=0.0_wp;translation(a:9:3)=1.0_wp
            if(maxval(abs(matmul(second(:,:,i),translation)))>1.0e-11_wp) error stop 'surface translation'
         enddo
      enddo
      previous=huge(1.0_wp)
      do level=1,2
         step=1.0e-4_wp/real(2**(level-1),wp)
         do coordinate=1,9
            a=mod(coordinate-1,3)+1;atom=(coordinate-1)/3+1
            moved=xyz;moved(a,atom)=moved(a,atom)+step
            call born%update(env,[8,1,1],moved);fp=reshape(born%dsdrt,[9,3])
            moved(a,atom)=xyz(a,atom)-step
            call born%update(env,[8,1,1],moved);fm=reshape(born%dsdrt,[9,3])
            fd(:,coordinate,:)=(fp-fm)/(2.0_wp*step)
         enddo
         error=maxval(abs(fd-second))
         print *, 'surface model/step/error:',modelIndex,step,error
         if(error>1.0e-7_wp*max(1.0_wp,maxval(abs(second)))) error stop 'all-coordinate surface Hessian'
         if(level==2.and.error>0.4_wp*previous+1.0e-8_wp) error stop 'surface quadratic convergence'
         previous=error
      enddo
      ! A fixed angular quadrature is not exactly rotation-invariant. Thus
      ! test translation and raw reciprocity, not a false continuum identity.
   enddo
   call init(model,env,TSolvInput(solvent='water',alpb=.false.,kernel=gbKernel%still),1)
   call newBornModel(model,env,born,[8,1,1])
   xyz(2,3)=aatoau*0.65_wp;call born%update(env,[8,1,1],xyz)
   call getSurfaceResponse(born,xyz,surface,first,second,ok,guardRadius=2.0e-4_wp)
   if(ok) error stop 'known surface screening displacement boundary accepted'
   print *, 'PASS angular surface Hessians and screened-branch guard'
end program
