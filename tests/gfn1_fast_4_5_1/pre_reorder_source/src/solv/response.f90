! This file is part of xtb.
! SPDX-License-Identifier: LGPL-3.0-or-later
module xtb_solv_response
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_sasa, only : ah0,ah1,ah3,tolsesp
   use xtb_solv_cm5, only : calc_cm5
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_derivative, only : TDerivative,derivativeVariable,derivativeExp, &
      & operator(+),operator(-),operator(*),operator(/),operator(**)
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: getBornRadiusResponse,getSurfaceResponse
   public :: buildGFN1SolventResponse

   ! Scalar radial first and second derivatives. No displaced geometries.
   type :: TRadial
      real(wp) :: value=0.0_wp, first=0.0_wp, second=0.0_wp
   end type
   interface operator(+)
      module procedure add_rr,add_sr,add_rs
   end interface
   interface operator(-)
      module procedure sub_rr,sub_sr,sub_rs
   end interface
   interface operator(*)
      module procedure mul_rr,mul_sr,mul_rs
   end interface
   interface operator(/)
      module procedure div_rr,div_sr,div_rs
   end interface
contains

! Actual GFN1 solvent functional at fixed bare atomic charges q:
! Q=q+CM5(R), E=Q*A(R)*Q/2+SASA(R)+constant. Return the fixed-q
! geometry Hessian, charge potential V=A*Q, dV/dR, and charge kernel A.
! Streaming pair coefficients avoids an N*N*3N*3N Born-matrix tensor.
! Salt kernels are deliberately rejected pending their independent response.
subroutine buildGFN1SolventResponse(born,num,xyz,q,energy,gradient,hessian,potential, &
      & potentialDerivative,chargeKernel,ok,guardRadius)
   type(TBorn), intent(in) :: born
   integer, intent(in) :: num(:)
   real(wp), intent(in) :: xyz(:,:),q(:)
   real(wp), intent(out) :: energy,gradient(:,:),hessian(:,:),potential(:),potentialDerivative(:,:),chargeKernel(:,:)
   logical, intent(out) :: ok
   real(wp), intent(in), optional :: guardRadius
   real(wp), allocatable :: br(:),brG(:,:),brH(:,:,:),surface(:),surfaceG(:,:),surfaceH(:,:,:)
   real(wp), allocatable :: cm(:),dc(:,:,:),cmG(:,:),charges(:),db(:),dbJac(:,:),g(:),pairG(:),cross(:)
   real(wp), allocatable :: fixedPotentialDerivative(:,:),cmH(:,:),shapeG(:),shapeH(:,:)
   type(TDerivative) :: rj,bi,bj,ab,kernel
   real(wp) :: r,vec(3),unit(3),block(3,3),qq,factor,shapeFactor,totalCharge
   integer :: n,n3,k,i,j,ia,ja,a,b,allocationStatus
   logical :: good
   energy=0.0_wp;gradient=0.0_wp;hessian=0.0_wp;potential=0.0_wp
   potentialDerivative=0.0_wp;chargeKernel=0.0_wp;ok=.false.
   n=born%nat;n3=3*n
   if(n<1.or.any(shape(xyz)/=[3,n]).or.size(num)/=n.or.size(q)/=n) return
   if(any(shape(gradient)/=[3,n]).or.any(shape(hessian)/=[n3,n3]).or.size(potential)/=n) return
   if(any(shape(potentialDerivative)/=[n,n3]).or.any(shape(chargeKernel)/=[n,n])) return
   if(.not.born%geometryCacheValid.or.born%lsalt) return
   if(born%kernel/=gbKernel%still.and.born%kernel/=gbKernel%p16) return
   if(.not.allocated(born%at).or..not.allocated(born%bornMat).or..not.allocated(born%gamsasa)) return
   if(size(born%at)/=n.or.any(shape(born%bornMat)/=[n,n]).or.size(born%gamsasa)/=n) return
   if(any(num/=born%at).or.any(num<1).or.any(num>118)) return
   if(.not.all(ieee_is_finite(q)).or..not.all(ieee_is_finite(born%bornMat))) return
   if(.not.ieee_is_finite(born%keps).or..not.ieee_is_finite(born%alpbet)) return
   if(.not.ieee_is_finite(born%gshift).or..not.all(ieee_is_finite(born%gamsasa))) return
   if(born%alpbet<0.0_wp) return
   allocate(br(n),brG(n3,n),brH(n3,n3,n),surface(n),surfaceG(n3,n),surfaceH(n3,n3,n), &
      & cm(n),dc(3,n,n),cmG(n3,n),charges(n),db(n),dbJac(n,n),g(n3),pairG(n3),cross(n3), &
      & fixedPotentialDerivative(n,n3),cmH(n3,n3),shapeG(n3),shapeH(n3,n3),stat=allocationStatus)
   if(allocationStatus/=0) return
   call getBornRadiusResponse(born,xyz,br,brG,brH,good)
   if(.not.good) return
   call getSurfaceResponse(born,xyz,surface,surfaceG,surfaceH,good,guardRadius)
   if(.not.good) return
   call calc_cm5(n,num,xyz,cm,dc)
   if(.not.all(ieee_is_finite(cm)).or..not.all(ieee_is_finite(dc))) return
   cmG=reshape(dc,[n3,n]);charges=q+cm
   db=0.0_wp;dbJac=0.0_wp;g=0.0_wp;fixedPotentialDerivative=0.0_wp
   do k=1,n*(n-1)/2
      i=born%ppind(1,k);j=born%ppind(2,k);ia=3*i-2;ja=3*j-2
      vec=xyz(:,i)-xyz(:,j);r=norm2(vec)
      if(r<=1.0e-10_wp) return
      unit=vec/r;pairG=0.0_wp;pairG(ia:ia+2)=unit;pairG(ja:ja+2)=-unit
      rj=derivativeVariable(r,1);bi=derivativeVariable(br(i),2);bj=derivativeVariable(br(j),3)
      if(born%kernel==gbKernel%still) then
         ab=bi*bj
         kernel=born%keps*(rj*rj+ab*derivativeExp(-0.25_wp*rj*rj/ab))**(-0.5_wp)
      else
         ab=(bi*bj)**0.5_wp
         kernel=born%keps/(rj+ab*(ab/(ab+(1.028_wp/16.0_wp)*rj))**16)
      endif
      chargeKernel(i,j)=kernel%value;chargeKernel(j,i)=kernel%value
      cross=kernel%gradient(1)*pairG+kernel%gradient(2)*brG(:,i)+kernel%gradient(3)*brG(:,j)
      fixedPotentialDerivative(i,:)=fixedPotentialDerivative(i,:)+charges(j)*cross
      fixedPotentialDerivative(j,:)=fixedPotentialDerivative(j,:)+charges(i)*cross
      qq=charges(i)*charges(j)
      g=g+qq*kernel%gradient(1)*pairG
      db(i)=db(i)+qq*kernel%gradient(2);db(j)=db(j)+qq*kernel%gradient(3)
      dbJac(i,i)=dbJac(i,i)+qq*kernel%hessian(2,2)
      dbJac(j,j)=dbJac(j,j)+qq*kernel%hessian(3,3)
      dbJac(i,j)=dbJac(i,j)+qq*kernel%hessian(2,3)
      dbJac(j,i)=dbJac(j,i)+qq*kernel%hessian(3,2)
      do a=1,3
         do b=1,3
            block(a,b)=qq*(kernel%hessian(1,1)-kernel%gradient(1)/r)*unit(a)*unit(b)
            if(a==b) block(a,b)=block(a,b)+qq*kernel%gradient(1)/r
         enddo
      enddo
      hessian(ia:ia+2,ia:ia+2)=hessian(ia:ia+2,ia:ia+2)+block
      hessian(ja:ja+2,ja:ja+2)=hessian(ja:ja+2,ja:ja+2)+block
      hessian(ia:ia+2,ja:ja+2)=hessian(ia:ia+2,ja:ja+2)-block
      hessian(ja:ja+2,ia:ia+2)=hessian(ja:ja+2,ia:ia+2)-block
      cross=qq*(kernel%hessian(1,2)*brG(:,i)+kernel%hessian(1,3)*brG(:,j))
      do a=1,3
         hessian(ia+a-1,:)=hessian(ia+a-1,:)+unit(a)*cross
         hessian(:,ia+a-1)=hessian(:,ia+a-1)+unit(a)*cross
         hessian(ja+a-1,:)=hessian(ja+a-1,:)-unit(a)*cross
         hessian(:,ja+a-1)=hessian(:,ja+a-1)-unit(a)*cross
      enddo
   enddo
   do i=1,n
      kernel=born%keps/derivativeVariable(br(i),1)
      chargeKernel(i,i)=kernel%value
      fixedPotentialDerivative(i,:)=fixedPotentialDerivative(i,:)+charges(i)*kernel%gradient(1)*brG(:,i)
      db(i)=db(i)+0.5_wp*charges(i)**2*kernel%gradient(1)
      dbJac(i,i)=dbJac(i,i)+0.5_wp*charges(i)**2*kernel%hessian(1,1)
   enddo
   g=g+matmul(brG,db)
   hessian=hessian+matmul(brG,matmul(dbJac,transpose(brG)))
   do i=1,n
      hessian=hessian+db(i)*brH(:,:,i)
      factor=born%gamsasa(i)
      if(born%lhb) then
         if(.not.allocated(born%hbmag)) return
         if(size(born%hbmag)/=n) return
         if(.not.ieee_is_finite(born%hbmag(i))) return
         shapeFactor=2.0_wp*born%hbmag(i)/born%vdwsa(i)**2
         chargeKernel(i,i)=chargeKernel(i,i)+shapeFactor*surface(i)
         fixedPotentialDerivative(i,:)=fixedPotentialDerivative(i,:)+charges(i)*shapeFactor*surfaceG(:,i)
         factor=factor+0.5_wp*charges(i)**2*shapeFactor
      endif
      g=g+factor*surfaceG(:,i);hessian=hessian+factor*surfaceH(:,:,i)
   enddo
   if(born%alpbet>0.0_wp) then
      call shapeResponse(born,xyz,shapeFactor,shapeG,shapeH,good)
      if(.not.good) return
      totalCharge=sum(charges)
      chargeKernel=chargeKernel+shapeFactor
      g=g+0.5_wp*totalCharge**2*shapeG;hessian=hessian+0.5_wp*totalCharge**2*shapeH
      do i=1,n
         fixedPotentialDerivative(i,:)=fixedPotentialDerivative(i,:)+totalCharge*shapeG
      enddo
   endif
   if(maxval(abs(chargeKernel-born%bornMat))>2.0e-12_wp*max(1.0_wp,maxval(abs(born%bornMat)))) return
   potential=matmul(chargeKernel,charges)
   energy=0.5_wp*dot_product(charges,potential)+dot_product(surface,born%gamsasa)+born%gshift
   ! CM5 changes even while the bare q is fixed.
   g=g+matmul(cmG,potential)
   hessian=hessian+matmul(transpose(fixedPotentialDerivative),transpose(cmG)) &
      & +matmul(cmG,fixedPotentialDerivative)+matmul(cmG,matmul(chargeKernel,transpose(cmG)))
   call calc_cm5(n,num,xyz,cm,dc,potential,cmH,good)
   if(.not.good) return
   hessian=hessian+cmH
   potentialDerivative=fixedPotentialDerivative+matmul(chargeKernel,transpose(cmG))
   gradient=reshape(g,[3,n])
   ok=ieee_is_finite(energy).and.all(ieee_is_finite(gradient)).and.all(ieee_is_finite(hessian)) &
      & .and.all(ieee_is_finite(potential)).and.all(ieee_is_finite(potentialDerivative))
end subroutine

! ALPB coefficient keps*alpbet/aDet. Differentiate the six independent
! weighted inertia components, including the moving radius-weighted center.
subroutine shapeResponse(born,xyz,factor,first,second,ok)
   type(TBorn), intent(in) :: born
   real(wp), intent(in) :: xyz(:,:)
   real(wp), intent(out) :: factor,first(:),second(:,:)
   logical, intent(out) :: ok
   real(wp), allocatable :: weights(:),vectors(:,:),fieldG(:,:),fieldH(:,:,:)
   real(wp) :: total,center(3),inertia(6),base,term,identity
   type(TDerivative) :: t(6),det,f
   integer, parameter :: rows(6)=[1,2,3,1,1,2],cols(6)=[1,2,3,2,3,3]
   integer :: n,n3,i,j,c,d,e,a,b,x,y,allocationStatus
   factor=0.0_wp;first=0.0_wp;second=0.0_wp;ok=.false.;n=born%nat;n3=3*n
   if(.not.ieee_is_finite(born%aDet).or.born%aDet<=0.0_wp) return
   allocate(weights(n),vectors(3,n),fieldG(n3,6),fieldH(n3,n3,6),stat=allocationStatus)
   if(allocationStatus/=0) return
   ! Accumulate the moving center explicitly, as in getADet. Keeping the
   ! small moment construction separate from the derivative loops also
   ! avoids incorrect moments observed with the fused version under ifx /O3.
   center=0.0_wp;total=0.0_wp
   do i=1,n
      weights(i)=born%vdwr(i)**3;total=total+weights(i)
      center=center+xyz(:,i)*weights(i)
   enddo
   if(total<=0.0_wp.or..not.ieee_is_finite(total)) return
   center=center/total
   do i=1,n
      vectors(:,i)=xyz(:,i)-center
   enddo
   inertia=0.0_wp;fieldG=0.0_wp;fieldH=0.0_wp
   do e=1,6
      a=rows(e);b=cols(e);identity=0.0_wp
      if(a==b) identity=1.0_wp
      do i=1,n
         inertia(e)=inertia(e)+weights(i)*(identity*(sum(vectors(:,i)**2)+0.4_wp*born%vdwr(i)**2) &
            & -vectors(a,i)*vectors(b,i))
      enddo
      t(e)=derivativeVariable(inertia(e),e)
   enddo
   do e=1,6
      a=rows(e);b=cols(e);identity=0.0_wp
      if(a==b) identity=1.0_wp
      do i=1,n
         do c=1,3
            x=3*(i-1)+c;term=2.0_wp*vectors(c,i)*identity
            if(c==a) term=term-vectors(b,i)
            if(c==b) term=term-vectors(a,i)
            fieldG(x,e)=weights(i)*term
            do j=1,n
               base=-weights(i)*weights(j)/total
               if(i==j) base=base+weights(i)
               do d=1,3
                  y=3*(j-1)+d;term=0.0_wp
                  if(c==d) term=2.0_wp*identity
                  if(c==a.and.d==b) term=term-1.0_wp
                  if(c==b.and.d==a) term=term-1.0_wp
                  fieldH(x,y,e)=base*term
               enddo
            enddo
         enddo
      enddo
   enddo
   det=t(1)*t(2)*t(3)+2.0_wp*t(4)*t(5)*t(6)-t(1)*t(6)**2-t(2)*t(5)**2-t(3)*t(4)**2
   if(.not.ieee_is_finite(det%value).or.det%value<=0.0_wp) return
   f=(born%keps*born%alpbet*sqrt(0.4_wp*total))*det**(-1.0_wp/6.0_wp)
   factor=f%value;first=matmul(fieldG,f%gradient)
   second=matmul(fieldG,matmul(f%hessian,transpose(fieldG)))
   do e=1,6
      second=second+f%gradient(e)*fieldH(:,:,e)
   enddo
   if(abs(factor-born%keps*born%alpbet/born%aDet)>2.0e-12_wp*max(1.0_wp,abs(factor))) return
   ok=all(ieee_is_finite(first)).and.all(ieee_is_finite(second))
end subroutine

! Differentiate the actual screened angular quadrature, retaining its exact
! point weights and branch choices. guardRadius optionally rejects references
! whose screening/soft-sphere branches may change over a displacement window.
! first/second orders are coordinate, [coordinate,] surface atom.
subroutine getSurfaceResponse(born,xyz,surfaces,first,second,ok,guardRadius)
   type(TBorn), intent(in) :: born
   real(wp), intent(in) :: xyz(:,:)
   real(wp), intent(out) :: surfaces(:),first(:,:),second(:,:,:)
   logical, intent(out) :: ok
   real(wp), intent(in), optional :: guardRadius
   real(wp), allocatable :: pointG(:),pointH(:,:),factorG(:),factorH(:,:)
   real(wp) :: vec(3),unit(3),point(3),r,r2,u,f,f1,f2,weight,quadrature,block(3,3),radius,margin,inner
   integer :: n,n3,i,j,k,ip,a,b,ia,ja,allocationStatus
   logical :: buried
   surfaces=0.0_wp;first=0.0_wp;second=0.0_wp;ok=.false.
   n=born%nat;n3=3*n;radius=0.0_wp
   if(present(guardRadius)) radius=guardRadius
   if(.not.ieee_is_finite(radius).or.radius<0.0_wp) return
   if(n<1.or.any(shape(xyz)/=[3,n])) return
   if(size(surfaces)/=n.or.any(shape(first)/=[n3,n]).or.any(shape(second)/=[n3,n3,n])) return
   if(.not.born%geometryCacheValid) return
   if(.not.allocated(born%cachedXYZ).or..not.allocated(born%angGrid).or..not.allocated(born%angWeight)) return
   if(.not.allocated(born%nnsas).or..not.allocated(born%nnlists).or..not.allocated(born%vdwsa)) return
   if(.not.allocated(born%wrp).or..not.allocated(born%trj2).or..not.allocated(born%sasa)) return
   if(.not.allocated(born%dsdrt)) return
   if(any(shape(born%cachedXYZ)/=[3,n]).or.any(shape(born%angGrid)/=[3,size(born%angWeight)])) return
   if(size(born%nnsas)/=n.or.any(shape(born%nnlists)/=[n,n]).or.size(born%vdwsa)/=n) return
   if(size(born%wrp)/=n.or.any(shape(born%trj2)/=[2,n]).or.size(born%sasa)/=n) return
   if(any(shape(born%dsdrt)/=[3,n,n])) return
   if(.not.all(ieee_is_finite(xyz)).or..not.all(xyz==born%cachedXYZ)) return
   if(.not.all(ieee_is_finite(born%vdwsa)).or..not.all(ieee_is_finite(born%trj2))) return
   if(any(born%vdwsa<=0.0_wp).or.any(born%trj2<0.0_wp)) return
   if(any(born%nnsas<0).or.any(born%nnsas>n)) return
   if(.not.all(ieee_is_finite(born%wrp)).or..not.all(ieee_is_finite(born%angGrid))) return
   if(.not.all(ieee_is_finite(born%angWeight))) return
   if(.not.all(ieee_is_finite(born%sasa)).or..not.all(ieee_is_finite(born%dsdrt))) return
   if(.not.ieee_is_finite(born%srcut).or.born%srcut<=0.0_wp) return
   ! The pair neighbor list also has a branch independent of quadrature.
   do i=1,n
      do j=1,i-1
         r=norm2(xyz(:,i)-xyz(:,j))
         if(abs(r-born%srcut)<=radius+1.0e-10_wp*max(1.0_wp,born%srcut)) return
      enddo
   enddo
   allocate(pointG(n3),pointH(n3,n3),factorG(n3),factorH(n3,n3),stat=allocationStatus)
   if(allocationStatus/=0) return
   do i=1,n
      ia=3*i-2
      do ip=1,size(born%angWeight)
         point=xyz(:,i)+born%vdwsa(i)*born%angGrid(:,ip)
         ! A point deep inside any neighbour remains identically zero
         ! throughout the single-coordinate displacement window. Other
         ! neighbours' sphere/screening branches cannot affect it. Inspect
         ! burial first, independent of neighbour order, before rejecting
         ! irrelevant branch changes or building dense point Hessians.
         buried=.false.
         do k=1,born%nnsas(i)
            j=born%nnlists(k,i)
            if(j<1.or.j>n.or.j==i) return
            inner=sqrt(born%trj2(1,j))-radius-1.0e-10_wp*max(1.0_wp,born%vdwsa(j))
            if(inner<=0.0_wp) cycle
            r2=sum((point-xyz(:,j))**2)
            if(r2<inner*inner) then
               buried=.true.;exit
            endif
         enddo
         if(buried) cycle
         weight=1.0_wp;pointG=0.0_wp;pointH=0.0_wp;buried=.false.
         do k=1,born%nnsas(i)
            j=born%nnlists(k,i)
            if(j<1.or.j>n.or.j==i) return
            ja=3*j-2;vec=point-xyz(:,j);r2=sum(vec*vec)
            r=sqrt(r2)
            margin=radius+1.0e-10_wp*max(1.0_wp,r,born%vdwsa(j))
            if(abs(r-sqrt(born%trj2(1,j)))<=margin.or.abs(r-sqrt(born%trj2(2,j)))<=margin) return
            if(r2>=born%trj2(2,j)) cycle
            if(r2<=born%trj2(1,j)) then
               buried=.true.;exit
            endif
            if(r<=tiny(1.0_wp)) return
            u=r-born%vdwsa(j)
            f=ah0+(ah1+ah3*u*u)*u
            f1=ah1+3.0_wp*ah3*u*u;f2=6.0_wp*ah3*u
            unit=vec/r;factorG=0.0_wp;factorH=0.0_wp
            factorG(ia:ia+2)=f1*unit;factorG(ja:ja+2)=-f1*unit
            do a=1,3
               do b=1,3
                  block(a,b)=(f2-f1/r)*unit(a)*unit(b)
                  if(a==b) block(a,b)=block(a,b)+f1/r
               enddo
            enddo
            factorH(ia:ia+2,ia:ia+2)=block;factorH(ja:ja+2,ja:ja+2)=block
            factorH(ia:ia+2,ja:ja+2)=-block;factorH(ja:ja+2,ia:ia+2)=-block
            pointH=f*pointH+weight*factorH+spread(pointG,1,n3)*spread(factorG,2,n3) &
               & +spread(factorG,1,n3)*spread(pointG,2,n3)
            pointG=f*pointG+weight*factorG;weight=weight*f
         enddo
         if(buried) cycle
         margin=max(1.0e-12_wp,radius*sum(abs(pointG))+0.5_wp*radius*radius*sum(abs(pointH)))
         if(abs(weight-tolsesp)<=margin) return
         if(weight<=tolsesp) cycle
         quadrature=born%angWeight(ip)*born%wrp(i)
         surfaces(i)=surfaces(i)+quadrature*weight
         first(:,i)=first(:,i)+quadrature*pointG
         second(:,:,i)=second(:,:,i)+quadrature*pointH
      enddo
   enddo
   if(.not.all(ieee_is_finite(surfaces)).or..not.all(ieee_is_finite(first))) return
   if(.not.all(ieee_is_finite(second))) return
   if(maxval(abs(surfaces-born%sasa))>2.0e-12_wp*max(1.0_wp,maxval(abs(born%sasa)))) return
   if(maxval(abs(first-reshape(born%dsdrt,[n3,n])))>2.0e-11_wp*max(1.0_wp,maxval(abs(born%dsdrt)))) return
   ok=.true.
end subroutine

! The directed pair descreening integral matches compute_psi, including
! nonoverlapping, partially overlapping and fully buried sphere branches.
! The near-equality branch in compute_psi is handled by the pair caller.
pure function descreen(r,vdw,rho) result(p)
   type(TRadial), intent(in) :: r
   real(wp), intent(in) :: vdw,rho
   type(TRadial) :: p,ap,am
   ap=r+rho; am=r-rho
   if(r%value>=vdw+rho) then
      p=rho/(ap*am)+0.5_wp*radial_log(am/ap)/r
   else if(r%value+rho>vdw) then
      p=1.0_wp/vdw-1.0_wp/ap &
         & +(0.5_wp*am*(1.0_wp/ap-ap/(vdw*vdw))-radial_log(ap/vdw))/(2.0_wp*r)
   else
      p=TRadial()
   endif
end function

! Reconstruct the actual Born radii and their Cartesian derivatives from the
! cached TBorn geometry. Output Hessian order: coordinate, coordinate, atom.
! Reference and first-derivative parity with compute_bornr are mandatory.
! Boundary/invalid input failures return ok=false; callers must fall back.
subroutine getBornRadiusResponse(born,xyz,radii,first,second,ok)
   type(TBorn), intent(in) :: born
   real(wp), intent(in) :: xyz(:,:)
   real(wp), intent(out) :: radii(:),first(:,:),second(:,:,:)
   logical, intent(out) :: ok
   real(wp), allocatable :: psi(:),dpsi(:,:),hpsi(:,:,:),dr(:),hr(:,:)
   real(wp) :: vec(3),r,unit(3),block(3,3),scale,threshold
   type(TRadial) :: distance,p,t,arg,f,equalPair
   integer :: n,n3,i,j,k,a,b,aa,bb,direction,allocationStatus
   logical :: equalReduced
   radii=0.0_wp; first=0.0_wp; second=0.0_wp; ok=.false.
   n=born%nat; n3=3*n
   if(n<1.or.any(shape(xyz)/=[3,n])) return
   if(size(radii)/=n.or.any(shape(first)/=[n3,n]).or.any(shape(second)/=[n3,n3,n])) return
   if(.not.born%geometryCacheValid) return
   if(.not.allocated(born%cachedXYZ).or..not.allocated(born%vdwr).or..not.allocated(born%rho)) return
   if(.not.allocated(born%svdw).or..not.allocated(born%brad).or..not.allocated(born%brdr)) return
   if(.not.allocated(born%ppind)) return
   if(any(shape(born%cachedXYZ)/=[3,n]).or.size(born%vdwr)/=n.or.size(born%rho)/=n) return
   if(size(born%svdw)/=n.or.size(born%brad)/=n.or.any(shape(born%brdr)/=[3,n,n])) return
   if(any(shape(born%ppind)/=[2,n*(n-1)/2])) return
   if(.not.all(ieee_is_finite(xyz)).or..not.all(xyz==born%cachedXYZ)) return
   if(.not.all(ieee_is_finite(born%vdwr)).or..not.all(ieee_is_finite(born%rho))) return
   if(.not.all(ieee_is_finite(born%svdw)).or..not.ieee_is_finite(born%bornScale)) return
   if(.not.all(ieee_is_finite(born%brad)).or..not.all(ieee_is_finite(born%brdr))) return
   if(any(born%vdwr<=0.0_wp).or.any(born%rho<0.0_wp).or.any(born%svdw<=0.0_wp)) return
   if(born%bornScale<=0.0_wp.or..not.ieee_is_finite(born%lrcut).or.born%lrcut<=0.0_wp) return
   allocate(psi(n),dpsi(n3,n),hpsi(n3,n3,n),dr(n3),hr(n3,n3),stat=allocationStatus)
   if(allocationStatus/=0) return
   psi=0.0_wp; dpsi=0.0_wp; hpsi=0.0_wp
   do k=1,n*(n-1)/2
      i=born%ppind(1,k); j=born%ppind(2,k)
      if(i<1.or.i>n.or.j<1.or.j>n.or.i==j) return
      vec=xyz(:,i)-xyz(:,j); r=norm2(vec)
      threshold=1.0e-10_wp*max(1.0_wp,r,born%lrcut)
      if(r<=threshold.or.abs(r-born%lrcut)<=threshold) return
      if(r>=born%lrcut) cycle
      unit=vec/r; dr=0.0_wp; hr=0.0_wp
      dr(3*i-2:3*i)=unit; dr(3*j-2:3*j)=-unit
      do a=1,3
         do b=1,3
            block(a,b)=-unit(a)*unit(b)/r
            if(a==b) block(a,b)=block(a,b)+1.0_wp/r
         enddo
      enddo
      hr(3*i-2:3*i,3*i-2:3*i)=block
      hr(3*j-2:3*j,3*j-2:3*j)=block
      hr(3*i-2:3*i,3*j-2:3*j)=-block
      hr(3*j-2:3*j,3*i-2:3*i)=-block
      distance=TRadial(r,1.0_wp,0.0_wp)
      ! Preserve compute_psi's near-equal-rho branch byte-for-byte at the
      ! parameter level: its nonoverlapping pair uses rho(j) for both atoms.
      equalReduced=abs(born%rho(i)-born%rho(j))<1.0e-8_wp &
         & .and.r>=born%vdwr(i)+born%rho(j).and.r>=born%rho(i)+born%vdwr(j)
      if(equalReduced) equalPair=descreen(distance,born%vdwr(i),born%rho(j))
      do direction=1,2
         if(direction==1) then
            aa=i; bb=j
         else
            aa=j; bb=i
         endif
         if(abs(r-born%vdwr(aa)-born%rho(bb))<=threshold) return
         if(abs(r+born%rho(bb)-born%vdwr(aa))<=threshold) return
         if(equalReduced) then
            p=equalPair
         else
            p=descreen(distance,born%vdwr(aa),born%rho(bb))
         endif
         psi(aa)=psi(aa)+p%value
         dpsi(:,aa)=dpsi(:,aa)+p%first*dr
         hpsi(:,:,aa)=hpsi(:,:,aa)+p%first*hr &
            & +p%second*spread(dr,1,n3)*spread(dr,2,n3)
      enddo
   enddo
   do i=1,n
      t=0.5_wp*born%svdw(i)*TRadial(psi(i),1.0_wp,0.0_wp)
      arg=t*(1.0_wp+t*(4.85_wp*t-0.8_wp))
      f=1.0_wp/born%svdw(i)-radial_tanh(arg)/born%vdwr(i)
      if(f%value<=tiny(1.0_wp)) return
      f=born%bornScale/f
      radii(i)=f%value
      first(:,i)=f%first*dpsi(:,i)
      second(:,:,i)=f%first*hpsi(:,:,i) &
         & +f%second*spread(dpsi(:,i),1,n3)*spread(dpsi(:,i),2,n3)
   enddo
   if(.not.all(ieee_is_finite(radii)).or..not.all(ieee_is_finite(first))) return
   if(.not.all(ieee_is_finite(second))) return
   scale=max(1.0_wp,maxval(abs(born%brad)))
   if(maxval(abs(radii-born%brad))>2.0e-12_wp*scale) return
   scale=max(1.0_wp,maxval(abs(born%brdr)))
   if(maxval(abs(first-reshape(born%brdr,[n3,n])))>2.0e-11_wp*scale) return
   ok=.true.
end subroutine

pure function chain(a,v,d1,d2) result(b)
   type(TRadial), intent(in) :: a
   real(wp), intent(in) :: v,d1,d2
   type(TRadial) :: b
   b=TRadial(v,d1*a%first,d1*a%second+d2*a%first*a%first)
end function

pure function radial_log(a) result(b)
   type(TRadial), intent(in) :: a
   type(TRadial) :: b
   b=chain(a,log(a%value),1.0_wp/a%value,-1.0_wp/(a%value*a%value))
end function

pure function radial_tanh(a) result(b)
   type(TRadial), intent(in) :: a
   type(TRadial) :: b
   real(wp) :: v,d
   v=tanh(a%value); d=1.0_wp-v*v
   b=chain(a,v,d,-2.0_wp*v*d)
end function

pure function add_rr(a,b) result(c)
   type(TRadial), intent(in) :: a,b
   type(TRadial) :: c
   c=TRadial(a%value+b%value,a%first+b%first,a%second+b%second)
end function
pure function add_sr(a,b) result(c)
   real(wp), intent(in) :: a
   type(TRadial), intent(in) :: b
   type(TRadial) :: c
   c=b; c%value=a+b%value
end function
pure function add_rs(a,b) result(c)
   type(TRadial), intent(in) :: a
   real(wp), intent(in) :: b
   type(TRadial) :: c
   c=add_sr(b,a)
end function
pure function sub_rr(a,b) result(c)
   type(TRadial), intent(in) :: a,b
   type(TRadial) :: c
   c=TRadial(a%value-b%value,a%first-b%first,a%second-b%second)
end function
pure function sub_sr(a,b) result(c)
   real(wp), intent(in) :: a
   type(TRadial), intent(in) :: b
   type(TRadial) :: c
   c=TRadial(a-b%value,-b%first,-b%second)
end function
pure function sub_rs(a,b) result(c)
   type(TRadial), intent(in) :: a
   real(wp), intent(in) :: b
   type(TRadial) :: c
   c=a; c%value=a%value-b
end function
pure function mul_rr(a,b) result(c)
   type(TRadial), intent(in) :: a,b
   type(TRadial) :: c
   c=TRadial(a%value*b%value,a%first*b%value+a%value*b%first, &
      & a%second*b%value+2.0_wp*a%first*b%first+a%value*b%second)
end function
pure function mul_sr(a,b) result(c)
   real(wp), intent(in) :: a
   type(TRadial), intent(in) :: b
   type(TRadial) :: c
   c=TRadial(a*b%value,a*b%first,a*b%second)
end function
pure function mul_rs(a,b) result(c)
   type(TRadial), intent(in) :: a
   real(wp), intent(in) :: b
   type(TRadial) :: c
   c=mul_sr(b,a)
end function
pure function reciprocal(a) result(b)
   type(TRadial), intent(in) :: a
   type(TRadial) :: b
   real(wp) :: v
   v=1.0_wp/a%value; b=chain(a,v,-v*v,2.0_wp*v*v*v)
end function
pure function div_rr(a,b) result(c)
   type(TRadial), intent(in) :: a,b
   type(TRadial) :: c
   c=a*reciprocal(b)
end function
pure function div_sr(a,b) result(c)
   real(wp), intent(in) :: a
   type(TRadial), intent(in) :: b
   type(TRadial) :: c
   c=a*reciprocal(b)
end function
pure function div_rs(a,b) result(c)
   type(TRadial), intent(in) :: a
   real(wp), intent(in) :: b
   type(TRadial) :: c
   c=mul_sr(1.0_wp/b,a)
end function
end module
