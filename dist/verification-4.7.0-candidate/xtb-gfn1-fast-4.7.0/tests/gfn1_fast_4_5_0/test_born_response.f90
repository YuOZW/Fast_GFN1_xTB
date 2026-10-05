program test_born_response
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment,init
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : TSolvModel,init,newBornModel
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_response, only : getBornRadiusResponse
   use xtb_solv_cm5, only : calc_cm5
   implicit none
   type(TEnvironment) :: env
   type(TSolvModel) :: model
   type(TBorn) :: born
   integer :: atoms(2),cases(2,4),pair,modelIndex,geometry,coordinate,a,atom,level,i
   integer :: buried,partial,nonoverlap
   real(wp) :: xyz(3,2),moved(3,2),radii(2),first(6,2),second(6,6,2),fp(6,2),fm(6,2),fd(6,6,2)
   real(wp) :: cp(2),cm(2),c0(2),dc(3,2,2),dp(3,2,2),dm(3,2,2),weights(2),hc(6,6),fc(6,6)
   real(wp) :: distances(4),step,error,previous,unit(3),r
   logical :: ok,failed
   call mctc_init('test_born_response',10,.true.)
   call init(env)
   cases=reshape([8,1,6,6,17,7,53,8],[2,4])
   distances=[0.4_wp,2.3_wp,5.7_wp,10.0_wp];unit=[0.8_wp,0.48_wp,0.36_wp]
   buried=0;partial=0;nonoverlap=0
   do modelIndex=1,2
      call init(model,env,TSolvInput(solvent='water',alpb=modelIndex==2,kernel=gbKernel%still),1)
      call env%check(failed);if(failed) error stop 'solvent parameters'
      do pair=1,4
         atoms=cases(:,pair)
         call newBornModel(model,env,born,atoms)
         if(modelIndex==1.and.born%alpbet/=0.0_wp) error stop 'GBSA parameter fixture'
         if(modelIndex==2.and.born%alpbet<=0.0_wp) error stop 'ALPB parameter fixture'
         do geometry=1,4
            xyz(:,1)=0.0_wp;xyz(:,2)=distances(geometry)*unit
            call born%update(env,atoms,xyz)
            call getBornRadiusResponse(born,xyz,radii,first,second,ok)
            if(.not.ok) error stop 'Born response rejected ordinary reference'
            r=norm2(xyz(:,2))
            do i=1,2
               if(r>=born%vdwr(i)+born%rho(3-i)) then
                  nonoverlap=nonoverlap+1
               else if(r+born%rho(3-i)>born%vdwr(i)) then
                  partial=partial+1
               else
                  buried=buried+1
               endif
            enddo
            previous=huge(1.0_wp)
            do level=1,4
               step=4.0e-4_wp/real(2**(level-1),wp)
               do coordinate=1,6
                  a=mod(coordinate-1,3)+1;atom=(coordinate-1)/3+1
                  moved=xyz;moved(a,atom)=moved(a,atom)+step
                  call born%update(env,atoms,moved);fp=reshape(born%brdr,[6,2])
                  moved(a,atom)=xyz(a,atom)-step
                  call born%update(env,atoms,moved);fm=reshape(born%brdr,[6,2])
                  fd(:,coordinate,:)=(fp-fm)/(2.0_wp*step)
               enddo
               error=maxval(abs(fd-second))
               print *, 'Born model/pair/geometry/step/error:',modelIndex,pair,geometry,step,error
               if(level==4.and.error>2.0e-7_wp) error stop 'Born radius Hessian vs all Cartesian legacy derivative FD'
               if(level>1.and.error>0.4_wp*previous+2.0e-9_wp) error stop 'Born quadratic FD convergence'
               previous=error
            enddo
            call born%update(env,atoms,xyz)
            ! Same-type CM5 is identically zero; unequal elements test its
            ! actual correction with fixed atom weights, not solvent charges.
            weights=[0.37_wp,-0.29_wp]
            call calc_cm5(2,atoms,xyz,c0,dc,weights,hc,ok)
            if(.not.ok) error stop 'CM5 response rejected'
            previous=huge(1.0_wp)
            do level=1,4
            step=2.0e-4_wp/real(2**(level-1),wp)
            do coordinate=1,6
               a=mod(coordinate-1,3)+1;atom=(coordinate-1)/3+1
               moved=xyz;moved(a,atom)=moved(a,atom)+step
               call calc_cm5(2,atoms,moved,cp,dp)
               moved(a,atom)=xyz(a,atom)-step
               call calc_cm5(2,atoms,moved,cm,dm)
               fc(:,coordinate)=reshape((dp(:,:,1)*weights(1)+dp(:,:,2)*weights(2) &
                  & -dm(:,:,1)*weights(1)-dm(:,:,2)*weights(2))/(2.0_wp*step),[6])
            enddo
            error=maxval(abs(fc-hc))
            print *, 'CM5 pair/distance/step/error:',pair,distances(geometry),step,error
            if(level>1.and.error>0.4_wp*previous+2.0e-9_wp) error stop 'CM5 quadratic FD convergence'
            if(level==4.and.error>2.0e-7_wp) error stop 'CM5 weighted Hessian vs original derivative FD'
            previous=error
            enddo
            if(maxval(abs(hc-transpose(hc)))>1.0e-13_wp) error stop 'CM5 raw reciprocity'
            weights=0.37_wp
            call calc_cm5(2,atoms,xyz,c0,dc,weights,hc,ok)
            if(.not.ok.or.maxval(abs(hc))>1.0e-13_wp) error stop 'CM5 total-charge Hessian conservation'
         enddo
         xyz(:,1)=0.0_wp;xyz(:,2)=[born%lrcut,0.0_wp,0.0_wp]
         call born%update(env,atoms,xyz)
         call getBornRadiusResponse(born,xyz,radii,first,second,ok)
         if(ok) error stop 'Born cutoff boundary accepted'
         xyz(:,2)=[born%vdwr(1)+born%rho(2),0.0_wp,0.0_wp]
         call born%update(env,atoms,xyz)
         call getBornRadiusResponse(born,xyz,radii,first,second,ok)
         if(ok) error stop 'Born sphere-overlap boundary accepted'
      enddo
   enddo
   if(min(buried,partial,nonoverlap)<1) error stop 'descreening branch coverage incomplete'
   atoms=[6,6]
   call newBornModel(model,env,born,atoms)
   born%rho(2)=born%rho(1)+5.0e-9_wp
   xyz(:,1)=0.0_wp;xyz(:,2)=10.0_wp*unit
   call born%update(env,atoms,xyz)
   call getBornRadiusResponse(born,xyz,radii,first,second,ok)
   if(.not.ok) error stop 'near-equal reduced radii do not match original special branch'
   step=5.0e-5_wp
   do coordinate=1,6
      a=mod(coordinate-1,3)+1;atom=(coordinate-1)/3+1
      moved=xyz;moved(a,atom)=moved(a,atom)+step
      call born%update(env,atoms,moved);fp=reshape(born%brdr,[6,2])
      moved(a,atom)=xyz(a,atom)-step
      call born%update(env,atoms,moved);fm=reshape(born%brdr,[6,2])
      fd(:,coordinate,:)=(fp-fm)/(2.0_wp*step)
   enddo
   if(maxval(abs(fd-second))>2.0e-7_wp) error stop 'near-equal reduced-radius Hessian'
   atoms=[8,1]
   weights=[0.37_wp,-0.29_wp]
   call calc_cm5(2,atoms,xyz,c0,dc,hessianWeights=weights,ok=ok)
   if(ok) error stop 'CM5 incomplete optional response request accepted'
   call calc_cm5(2,atoms,xyz,c0,dc,weights,hc(:5,:5),ok)
   if(ok) error stop 'CM5 invalid Hessian dimensions accepted'
   print *, 'PASS Born and CM5 component response; buried/partial/nonoverlap=',buried,partial,nonoverlap
end program
