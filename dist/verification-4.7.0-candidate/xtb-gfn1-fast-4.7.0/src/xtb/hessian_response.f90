! This file is part of xtb.
! SPDX-License-Identifier: LGPL-3.0-or-later
module xtb_xtb_hessian_response
   !$ use omp_lib, only : omp_get_num_threads
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : evtoau
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_basisset, only : TBasisset
   use xtb_type_wavefunction, only : TWavefunction
   use xtb_xtb_data, only : TxTBData
   use xtb_xtb_coulomb, only : TxTBCoulomb,init
   use xtb_coulomb_klopmanohno, only : TKlopmanOhno,init,gamAverage
   use xtb_xtb_hamiltonian, only : getSelfEnergy,build_SH0_GFN1,build_dSH0_GFN1_noreset
   use xtb_xtb_hessian_integrals, only : buildGFN1IntegralResponse
   use xtb_xtb_response, only : TGFN1ElectronicResponse,TGFN1CoupledResponse
   use xtb_disp_coordinationnumber, only : cnType,getCoordinationNumberHessian
   use xtb_xtb_repulsion, only : repulsionMolecularHessian
   use xtb_disp_dftd3, only : d3MolecularHessian
   use xtb_xtb_halogen, only : xbMolecularHessian
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_response, only : buildGFN1SolventResponse
   use xtb_mctc_blas_level3, only : mctc_gemm
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy,getGFN1FastPolicy,getGFN1HessianResources
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: buildGFN1ElectrostaticHessian
   public :: buildGFN1GasHessian
   public :: buildGFN1SolvatedHessian
contains

! Complete gas-phase GFN1 Hessian, including the halogen correction.
! Solvation, external potentials and constraints require additional responses.
! References for other models
! must be dispatched to the numerical fallback by the calculator.
subroutine buildGFN1GasHessian(env,mol,basis,data,wfn,temperature,accuracy,gradient,hessian,ok,halogenTerms,guardRadius)
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   type(TBasisset), intent(in) :: basis
   type(TxTBData), intent(in) :: data
   type(TWavefunction), intent(in) :: wfn
   real(wp), intent(in) :: temperature,accuracy
   real(wp), intent(out) :: gradient(:,:),hessian(:,:)
   logical, intent(out) :: ok
   integer, intent(out), optional :: halogenTerms
   real(wp), intent(in), optional :: guardRadius
   call buildGFN1MolecularHessian(env,mol,basis,data,wfn,temperature,accuracy,gradient,hessian,ok, &
      & halogenTerms,guardRadius=guardRadius)
end subroutine

! Complete Born-solvated GFN1 Hessian. The solvent must be updated at this
! geometry. CM5 and solvent-induced SCC response are part of the assembly.
subroutine buildGFN1SolvatedHessian(env,mol,basis,data,wfn,temperature,accuracy,born, &
      & gradient,hessian,ok,halogenTerms,guardRadius)
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   type(TBasisset), intent(in) :: basis
   type(TxTBData), intent(in) :: data
   type(TWavefunction), intent(in) :: wfn
   real(wp), intent(in) :: temperature,accuracy
   type(TBorn), intent(in) :: born
   real(wp), intent(out) :: gradient(:,:),hessian(:,:)
   logical, intent(out) :: ok
   integer, intent(out), optional :: halogenTerms
   real(wp), intent(in), optional :: guardRadius
   call buildGFN1MolecularHessian(env,mol,basis,data,wfn,temperature,accuracy,gradient,hessian,ok, &
      & halogenTerms,born,guardRadius)
end subroutine

subroutine buildGFN1MolecularHessian(env,mol,basis,data,wfn,temperature,accuracy,gradient,hessian,ok, &
      & halogenTerms,solvent,guardRadius)
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   type(TBasisset), intent(in) :: basis
   type(TxTBData), intent(in) :: data
   type(TWavefunction), intent(in) :: wfn
   real(wp), intent(in) :: temperature,accuracy
   real(wp), intent(out) :: gradient(:,:),hessian(:,:)
   logical, intent(out) :: ok
   integer, intent(out), optional :: halogenTerms
   type(TBorn), intent(in), optional :: solvent
   real(wp), intent(in), optional :: guardRadius
   real(wp), allocatable :: g(:,:),hh(:,:),cn(:),dcndr(:,:,:),hcn(:,:,:,:,:)
   real(wp) :: energy
   integer :: n,terms
   logical :: good
   gradient=0.0_wp; hessian=0.0_wp; ok=.false.; n=mol%n
   if(present(halogenTerms)) halogenTerms=0
   call buildGFN1ElectrostaticHessian(env,mol,basis,data,wfn,temperature,accuracy,gradient,hessian,good, &
      & solvent=solvent,guardRadius=guardRadius)
   if(.not.good) return
   allocate(g(3,n),hh(3*n,3*n),cn(n),dcndr(3,n,n),hcn(3,n,3,n,n))
   call repulsionMolecularHessian(mol,data%repulsion,40.0_wp,energy,g,hh,good,guardRadius)
   if(.not.good) return
   gradient=gradient+g; hessian=hessian+hh
   call getCoordinationNumberHessian(mol,cnType%exp,cn,dcndr,hcn,good,guardRadius=guardRadius)
   if(.not.good) return
   call d3MolecularHessian(mol,data%dispersion%dpar,4.0_wp,60.0_wp,cn,dcndr,hcn,energy,g,hh,good,guardRadius)
   if(.not.good) return
   gradient=gradient+g; hessian=hessian+hh
   if(allocated(data%halogen)) then
      ! Match SCF's exponent, active pair cutoff and nearest-neighbour choice.
      ! Nonsmooth neighbour/cutoff references reject the analytic assembly.
      call xbMolecularHessian(mol,data%halogen,12.0_wp,energy,g,hh,good,terms,guardRadius)
      if(.not.good) return
      gradient=gradient+g; hessian=hessian+hh
      if(present(halogenTerms)) halogenTerms=terms
   endif
   ok=all(ieee_is_finite(hessian)).and.all(ieee_is_finite(gradient))
end subroutine

! Molecular GFN1 electronic + isotropic Coulomb/third-order Hessian.
! Differentiate the analytic Gradient, including moving AO metric, canonical
! finite-temperature occupations and self-consistent shell charge response.
! This component excludes repulsion/dispersion/halogen/external terms. The
! optional Born model adds its full geometry and coupled charge response.
! A caller supplying a solvated or externally biased reference must instead
! use the complete response model or numerical fallback. The Hamiltonian
! residual below rejects references inconsistent with this gas-phase model.
! No symmetrization is applied: reciprocity is a validation of the assembly.
subroutine buildGFN1ElectrostaticHessian(env,mol,basis,data,wfn,temperature,accuracy, &
      & gradient,hessian,ok,dqsh,densityResponse,rcond,solvent,guardRadius)
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   type(TBasisset), intent(in) :: basis
   type(TxTBData), intent(in) :: data
   type(TWavefunction), intent(in) :: wfn
   real(wp), intent(in) :: temperature,accuracy
   real(wp), intent(out) :: gradient(:,:),hessian(:,:)
   logical, intent(out) :: ok
   real(wp), intent(out), optional :: dqsh(:,:),densityResponse(:,:,:),rcond
   type(TBorn), intent(in), optional :: solvent
   real(wp), intent(in), optional :: guardRadius
   type(TxTBCoulomb) :: ies
   type(TKlopmanOhno) :: coulomb
   type(TGFN1ElectronicResponse) :: electronic
   type(TGFN1CoupledResponse) :: coupled
   type(TGFN1FastPolicy) :: policy
   real(wp), allocatable :: cn(:),dcndr(:,:,:),hcn(:,:,:,:,:),se(:,:),dse(:,:),s(:,:),packed(:)
   real(wp), allocatable :: w(:,:),effectiveW(:,:),occupation(:,:),scaledC(:,:),ham(:,:),sc(:,:)
   real(wp), allocatable, target :: ds(:,:,:),dh0(:,:,:)
   real(wp), pointer :: flatS(:,:),flatH(:,:)
   real(wp), allocatable :: geometric(:,:),coulombHessian(:,:),coulombGradient(:,:)
   real(wp), allocatable :: v(:),atomicV(:),jacobian(:,:),dvGeometry(:,:)
   real(wp), allocatable :: solvG(:,:),solvH(:,:),solvV(:),solvDV(:,:),solvA(:,:)
   real(wp), allocatable :: rawPacked(:),parentG(:,:),parentCN(:),shellV(:,:)
   integer, allocatable :: idnum(:),sh2at(:)
   logical, allocatable :: columnsOK(:)
   real(wp) :: intcut,neglect,trans(3,1),condition,scale,solvEnergy,parentSigma(3,3)
   real(wp) :: stageWall(5),stageStart,stageEnd,stageCPU
   integer(int64) :: estimatedBytes
   integer :: n,nao,ns,n3,i,j,x,y,ish,allocationStatus,responseThreads,actualResponseThreads
   logical :: good,failed
   gradient=0.0_wp; hessian=0.0_wp; ok=.false.
   if(present(dqsh)) dqsh=0.0_wp
   if(present(densityResponse)) densityResponse=0.0_wp
   if(present(rcond)) rcond=0.0_wp
   n=mol%n; nao=basis%nao; ns=basis%nshell; n3=3*n
   call getGFN1FastPolicy(policy)
   call getGFN1HessianResources(n,nao,ns,present(solvent),responseThreads,estimatedBytes)
   call timing(stageCPU,stageStart)
   if(n<1.or.nao<1.or.ns<1.or.mol%npbc/=0.or.data%level/=1) return
   if(any(shape(gradient)/=[3,n]).or.any(shape(hessian)/=[n3,n3])) return
   if(present(dqsh)) then
      if(any(shape(dqsh)/=[ns,n3])) return
   endif
   if(present(densityResponse)) then
      if(any(shape(densityResponse)/=[nao,nao,n3])) return
   endif
   if(estimatedBytes>int(policy%analyticHessianMaxMiB,int64)*1048576_int64) return
   if(.not.ieee_is_finite(accuracy).or.accuracy<=0.0_wp) return
   if(.not.ieee_is_finite(temperature).or.temperature<0.0_wp) return
   if(.not.allocated(wfn%c).or..not.allocated(wfn%p).or..not.allocated(wfn%emo)) return
   if(.not.allocated(wfn%focc).or..not.allocated(wfn%focca).or..not.allocated(wfn%foccb)) return
   if(.not.allocated(wfn%q).or..not.allocated(wfn%qsh)) return
   if(any(shape(wfn%c)/=[nao,nao]).or.any(shape(wfn%p)/=[nao,nao])) return
   if(size(wfn%emo)/=nao.or.size(wfn%focc)/=nao.or.size(wfn%q)/=n.or.size(wfn%qsh)/=ns) return
   if(size(wfn%focca)/=nao.or.size(wfn%foccb)/=nao) return
   if(.not.all(ieee_is_finite(wfn%c)).or..not.all(ieee_is_finite(wfn%p))) return
   if(.not.all(ieee_is_finite(wfn%emo)).or..not.all(ieee_is_finite(wfn%focc))) return
   if(.not.all(ieee_is_finite(wfn%focca)).or..not.all(ieee_is_finite(wfn%foccb))) return
   if(.not.all(ieee_is_finite(wfn%q)).or..not.all(ieee_is_finite(wfn%qsh))) return
   if(maxval(abs(wfn%focc-wfn%focca-wfn%foccb))>1.0e-12_wp) return
   if(present(solvent)) then
      ! Complete the O(N^3) solvent geometry temporaries before allocating
      ! the two O(Nao^2*3N) AO derivative tensors. Only the compact completed
      ! solvent fields below survive into the electronic response stage.
      allocate(solvG(3,n),solvH(n3,n3),solvV(n),solvDV(n,n3),solvA(n,n),stat=allocationStatus)
      if(allocationStatus/=0) return
      call buildGFN1SolventResponse(solvent,mol%at,mol%xyz,wfn%q,solvEnergy,solvG,solvH, &
         & solvV,solvDV,solvA,good,guardRadius,threads=responseThreads)
      if(.not.good) return
   endif
   call timing(stageCPU,stageEnd);stageWall(1)=stageEnd-stageStart;stageStart=stageEnd
   intcut=max(20.0_wp,25.0_wp-10.0_wp*log10(accuracy)); neglect=1.0e-8_wp*accuracy
   allocate(cn(n),dcndr(3,n,n),hcn(3,n,3,n,n),se(5,n),dse(5,n),s(nao,nao),packed(nao*(nao+1)/2), &
      & w(nao,nao),effectiveW(nao,nao),occupation(nao,2),scaledC(nao,nao),ham(nao,nao),sc(nao,nao), &
      & ds(nao,nao,n3),dh0(nao,nao,n3),geometric(n3,n3),coulombHessian(n3,n3),coulombGradient(3,n), &
      & v(ns),atomicV(n),jacobian(ns,ns),dvGeometry(ns,n3), &
      & idnum(maxval(mol%id)),sh2at(ns),rawPacked(nao*(nao+1)/2),columnsOK(n3), &
      & parentG(3,n),parentCN(n),shellV(5,n),stat=allocationStatus)
   if(allocationStatus/=0) return
   call getCoordinationNumberHessian(mol,cnType%exp,cn,dcndr,hcn,good,guardRadius=guardRadius)
   if(.not.good) return
   call getSelfEnergy(data%hamiltonian,data%nShell,mol%at,cn=cn,selfEnergy=se,dSEdcn=dse)
   trans=0.0_wp
   call build_SH0_GFN1(data%nShell,data%hamiltonian,n,mol%at,basis%nbf,nao,mol%xyz,trans,se,intcut, &
      & basis%caoshell,basis%saoshell,basis%nprim,basis%primcount,basis%alp,basis%cont,s,packed)
   rawPacked=packed
   ! Match the actual SCC overlap sparsity and Hamiltonian at the reference.
   do i=1,nao
      do j=1,i-1
         if(abs(s(j,i))<neglect) then
            s(j,i)=0.0_wp; s(i,j)=0.0_wp
            packed(i*(i-1)/2+j)=0.0_wp
         endif
      enddo
   enddo
   do i=1,n; idnum(mol%id(i))=mol%at(i); enddo
   call init(ies,data%coulomb,data%nShell,mol%at)
   call init(coulomb,env,mol,gamAverage%harmonic,data%coulomb%shellHardness, &
      & data%coulomb%gExp,num=idnum,nshell=data%nShell)
   call env%check(failed)
   if(failed) return
   call coulomb%getCoulombMatrix(mol,ies%jmat)
   sh2at=basis%ash(1:ns)
   call ies%getResponseKernel(wfn%q,wfn%qsh,sh2at,jacobian,good)
   if(.not.good) return
   v=0.0_wp; atomicV=0.0_wp
   call ies%addShift(wfn%q,wfn%qsh,atomicV,v)
   v=v+atomicV(sh2at)
   call coulomb%getMolecularHessianResponse(mol,wfn%qsh,dvGeometry,coulombHessian,good,coulombGradient)
   if(.not.good) return
   if(present(solvent)) then
      coulombGradient=coulombGradient+solvG;coulombHessian=coulombHessian+solvH
      do ish=1,ns
         v(ish)=v(ish)+solvV(sh2at(ish))
         dvGeometry(ish,:)=dvGeometry(ish,:)+solvDV(sh2at(ish),:)
         do j=1,ns
            jacobian(ish,j)=jacobian(ish,j)+solvA(sh2at(ish),sh2at(j))
         enddo
      enddo
   endif
   do i=1,nao
      do j=1,i
         ham(j,i)=packed(i*(i-1)/2+j)*evtoau-0.5_wp*s(j,i)*(v(basis%ao2sh(i))+v(basis%ao2sh(j)))
         ham(i,j)=ham(j,i)
      enddo
   enddo
   sc=matmul(s,wfn%c); scaledC=sc
   do i=1,nao; scaledC(:,i)=scaledC(:,i)*wfn%emo(i)*evtoau; enddo
   scale=max(1.0_wp,maxval(abs(ham)))
   if(maxval(abs(matmul(ham,wfn%c)-scaledC))>5.0e-7_wp*scale) return
   occupation(:,1)=wfn%focca; occupation(:,2)=wfn%foccb
   call electronic%initialize(wfn%c,wfn%emo*evtoau,occupation,temperature,good)
   if(.not.good) return
   call coupled%initialize(electronic,basis%ao2sh,s,wfn%p,jacobian,good,condition,threads=responseThreads)
   if(.not.good) return
   if(present(rcond)) rcond=condition
   call timing(stageCPU,stageEnd);stageWall(2)=stageEnd-stageStart;stageStart=stageEnd
   scaledC=wfn%c
   do i=1,nao; scaledC(:,i)=scaledC(:,i)*wfn%focc(i)*wfn%emo(i)*evtoau; enddo
   w=matmul(scaledC,transpose(wfn%c)); effectiveW=w
   do j=1,nao
      do i=1,nao
         effectiveW(i,j)=effectiveW(i,j)+0.5_wp*wfn%p(i,j)*(v(basis%ao2sh(i))+v(basis%ao2sh(j)))
      enddo
   enddo
   call buildGFN1IntegralResponse(mol,basis,data%hamiltonian,data%nShell,intcut,se,dse,dcndr,hcn, &
      & wfn%p,effectiveW,ds,dh0,geometric,good,threads=responseThreads)
   if(.not.good) return
   call timing(stageCPU,stageEnd);stageWall(3)=stageEnd-stageStart;stageStart=stageEnd
   do x=1,n3
      gradient(mod(x-1,3)+1,(x-1)/3+1)=sum(wfn%p*dh0(:,:,x))-sum(effectiveW*ds(:,:,x))
      ! The fixed-charge overlap potential belongs to both the analytic
      ! Gradient and the external Hamiltonian's nuclear derivative.
      do j=1,nao
         do i=1,nao
            dh0(i,j,x)=dh0(i,j,x)-0.5_wp*ds(i,j,x)*(v(basis%ao2sh(i))+v(basis%ao2sh(j)))
         enddo
      enddo
   enddo
   ! Validate against the production derivative, including its treatment of
   ! small nonzero SCF-screened overlaps. Exact nodes retain the direct H0
   ! derivative; screened nonzero nodes can be nonvariational and must fall
   ! back rather than silently return the derivative of a different function.
   shellV=0.0_wp; ish=0
   do i=1,n
      do j=1,data%nShell(mol%at(i))
         ish=ish+1; shellV(j,i)=v(ish)
      enddo
   enddo
   parentG=0.0_wp; parentCN=0.0_wp; parentSigma=0.0_wp
   call build_dSH0_GFN1_noreset(data%nShell,data%hamiltonian,se,dse,intcut,n,nao,basis%nbf, &
      & mol%at,mol%xyz,basis%caoshell,basis%saoshell,basis%nprim,basis%primcount,basis%alp, &
      & basis%cont,S=s,p=wfn%p,Pew=w,ves=shellV,dhdcn=parentCN,g=parentG,sigma=parentSigma, &
      & H0packed=rawPacked,cartContract=.true.,directPrefactor=.true.)
   do i=1,n
      do j=1,n
         parentG(:,i)=parentG(:,i)+dcndr(:,i,j)*parentCN(j)
      enddo
   enddo
   if(.not.all(ieee_is_finite(parentG))) return
   if(maxval(abs(parentG-gradient))>1.0e-10_wp) return
   gradient=gradient+coulombGradient
   hessian=geometric+coulombHessian
   flatS(1:nao*nao,1:n3)=>ds
   flatH(1:nao*nao,1:n3)=>dh0
   call timing(stageCPU,stageEnd);stageWall(4)=stageEnd-stageStart;stageStart=stageEnd
   columnsOK=.false.
   actualResponseThreads=1
   !$omp parallel num_threads(responseThreads) default(none) &
   !$omp shared(coupled,basis,wfn,ds,dh0,s,dvGeometry,jacobian,hessian,dqsh,densityResponse, &
   !$omp& columnsOK,nao,ns,n3,flatS,flatH,actualResponseThreads)
   !$omp single
   !$ actualResponseThreads=omp_get_num_threads()
   !$omp end single
   call coordinateWorker()
   !$omp end parallel
   if(.not.all(columnsOK)) return
   call timing(stageCPU,stageEnd);stageWall(5)=stageEnd-stageStart
   ok=all(ieee_is_finite(hessian)).and.all(ieee_is_finite(gradient))
   if(policy%profile) then
      write(env%unit,'(3x,a,i0)') 'analytic response threads = ',actualResponseThreads
      write(env%unit,'(3x,a,i0)') 'analytic requested workers = ',responseThreads
      write(env%unit,'(3x,a,5(f12.6,1x))') 'analytic stage wall (solvent/setup/integrals/gradient/columns) = ',stageWall
   endif
contains
   recursive subroutine coordinateWorker()
      real(wp), allocatable :: localH(:,:),localP(:,:),localW(:,:),weightedW(:,:),localQ(:),localV(:)
      real(wp), allocatable :: blockP(:,:),blockW(:,:),blockQ(:,:),blockH(:,:)
      integer, parameter :: blockSize=4
      integer :: column,ia,ja,allocationStatus,firstColumn,lastColumn,blockIndex,slot,count
      logical :: valid
      allocate(localH(nao,nao),localP(nao,nao),localW(nao,nao),weightedW(nao,nao), &
         & localQ(ns),localV(ns),blockP(nao*nao,blockSize),blockW(nao*nao,blockSize), &
         & blockQ(ns,blockSize),blockH(n3,blockSize),stat=allocationStatus)
      !$omp do schedule(static)
      do blockIndex=1,(n3+blockSize-1)/blockSize
         if(allocationStatus/=0) cycle
         firstColumn=1+(blockIndex-1)*blockSize;lastColumn=min(n3,firstColumn+blockSize-1)
         count=lastColumn-firstColumn+1
         blockP=0.0_wp;blockW=0.0_wp;blockQ=0.0_wp
         do column=firstColumn,lastColumn
         slot=column-firstColumn+1
         localH=dh0(:,:,column)
         do ja=1,nao
            do ia=1,nao
               localH(ia,ja)=localH(ia,ja)-0.5_wp*s(ia,ja) &
                  & *(dvGeometry(basis%ao2sh(ia),column)+dvGeometry(basis%ao2sh(ja),column))
            enddo
         enddo
         call coupled%apply(localH,ds(:,:,column),localP,localW,localQ,valid)
         if(.not.valid) cycle
         if(present(dqsh)) dqsh(:,column)=localQ
         if(present(densityResponse)) densityResponse(:,:,column)=localP
         localV=dvGeometry(:,column)+matmul(jacobian,localQ)
         weightedW=localW
         do ja=1,nao
            do ia=1,nao
               weightedW(ia,ja)=weightedW(ia,ja)+0.5_wp*wfn%p(ia,ja) &
                  & *(localV(basis%ao2sh(ia))+localV(basis%ao2sh(ja)))
            enddo
         enddo
         blockP(:,slot)=reshape(localP,[nao*nao]);blockW(:,slot)=reshape(weightedW,[nao*nao])
         blockQ(:,slot)=localQ
         columnsOK(column)=.true.
         enddo
         ! Equivalent trace contractions, batched so BLAS reuses each AO
         ! derivative panel across several response columns. No NAO^2*3N
         ! response tensor or additional global density tensor is retained.
         call mctc_gemm(flatH,blockP(:,1:count),blockH(:,1:count),transa='t')
         call mctc_gemm(flatS,blockW(:,1:count),blockH(:,1:count),transa='t',alpha=-1.0_wp,beta=1.0_wp)
         call mctc_gemm(dvGeometry,blockQ(:,1:count),blockH(:,1:count),transa='t',beta=1.0_wp)
         hessian(:,firstColumn:lastColumn)=hessian(:,firstColumn:lastColumn)+blockH(:,1:count)
      enddo
      !$omp end do
   end subroutine coordinateWorker
end subroutine
end module
