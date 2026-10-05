program test_halogen_hessian
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment,init
   use xtb_type_molecule, only : TMolecule,init
   use xtb_xtb_calculator, only : TxTBCalculator,newXTBCalculator
   use xtb_xtb_data, only : THalogenData
   use xtb_xtb_halogen, only : xbpot,xbMolecularHessian
   implicit none
   type(TEnvironment) :: env
   type(TMolecule) :: water,mol,moved
   type(TxTBCalculator) :: calc
   type(THalogenData) :: data,custom,missing
   real(wp) :: xyz(3,3),g(3,3),oldg(3,3),gp(3,3),gm(3,3),h(9,9),fd(9,9),translation(9)
   real(wp) :: energy,oldenergy,step,error,previous,exponent,maximumG,maximumH,maximumE
   integer :: donors(4),acceptors(4),ix,ia,geometry,model,level,coordinate,a,atom,terms,count,list(3,1)
   logical :: ok,failed
   call mctc_init('test_halogen_hessian',10,.true.)
   call init(env)
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]; xyz(:,2)=[1.4_wp,0.0_wp,0.0_wp]; xyz(:,3)=[-0.3_wp,1.3_wp,0.0_wp]
   call init(water,[8,1,1],xyz)
   call newXTBCalculator(env,water,calc,method=1)
   call env%check(failed); call require(.not.failed,'GFN1 halogen parameters')
   data=calc%xtbData%halogen; donors=[17,35,53,85]; acceptors=[7,8,15,16]
   list(:,1)=[1,3,2]; maximumG=0.0_wp; maximumH=0.0_wp; maximumE=0.0_wp; count=0
   do model=1,2
      custom=data
      ! Test a nonzero custom Cl term and changed strengths for every donor.
      if(model==2) custom%bondStrength(donors)=0.02_wp
      exponent=12.0_wp
      if(model==2) exponent=10.5_wp
      do ix=1,4
         do ia=1,4
            do geometry=1,3
               xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]; xyz(:,2)=[-3.4_wp,0.0_wp,0.0_wp]
               select case(geometry)
               case(1); xyz(:,3)=[6.1_wp,0.0_wp,0.0_wp]
               case(2); xyz(:,3)=[4.8_wp,3.5_wp,0.7_wp]
               case(3); xyz(:,3)=[-5.9_wp,0.1_wp,0.3_wp]
               end select
               call init(mol,[donors(ix),6,acceptors(ia)],xyz)
               call xbMolecularHessian(mol,custom,exponent,energy,g,h,ok,terms)
               call require(ok,'analytic halogen triplet')
               oldg=0.0_wp; oldenergy=0.0_wp
               call xbpot(custom,3,mol%at,mol%xyz,list,1,exponent,oldenergy,oldg)
               maximumG=max(maximumG,maxval(abs(g-oldg))); maximumE=max(maximumE,abs(energy-oldenergy))
               call require(maxval(abs(g-oldg))<5.0e-13_wp,'Gradient vs original xbpot')
               call require(abs(energy-oldenergy)<5.0e-13_wp,'Energy vs original xbpot')
               call require(maxval(abs(h-transpose(h)))<5.0e-13_wp,'halogen Hessian reciprocity')
               do a=1,3
                  translation=0.0_wp; translation(a:9:3)=1.0_wp
                  call require(maxval(abs(matmul(h,translation)))<5.0e-13_wp,'halogen Hessian translation')
               enddo
               previous=huge(1.0_wp)
               do level=1,2
                  step=1.0e-3_wp/real(2**(level-1),wp)
                  do coordinate=1,9
                     a=mod(coordinate-1,3)+1; atom=(coordinate-1)/3+1
                     call moved%copy(mol); moved%xyz(a,atom)=mol%xyz(a,atom)+step
                     gp=0.0_wp; oldenergy=0.0_wp
                     call xbpot(custom,3,moved%at,moved%xyz,list,1,exponent,oldenergy,gp)
                     moved%xyz(a,atom)=mol%xyz(a,atom)-step
                     gm=0.0_wp; oldenergy=0.0_wp
                     call xbpot(custom,3,moved%at,moved%xyz,list,1,exponent,oldenergy,gm)
                     fd(:,coordinate)=reshape((gp-gm)/(2.0_wp*step),[9])
                  enddo
                  error=maxval(abs(h-fd)); maximumH=max(maximumH,error)
                  call require(error<5.0e-8_wp,'all-coordinate original Gradient finite differences')
                  if(level==2.and.previous>1.0e-12_wp) then
                     call require(error<0.35_wp*previous+1.0e-13_wp,'quadratic halogen Hessian FD convergence')
                  endif
                  previous=error
               enddo
               count=count+1
            enddo
         enddo
      enddo
   enddo
   print *, 'PASS halogen triplets=',count,' max Energy=',maximumE,' max Gradient=',maximumG,' max Hessian FD=',maximumH
   ! Explicitly exercise the production choice when the acceptor is nearest.
   xyz(:,1)=[0.0_wp,0.0_wp,0.0_wp]; xyz(:,2)=[-8.0_wp,0.0_wp,0.0_wp]; xyz(:,3)=[3.0_wp,0.2_wp,0.1_wp]
   call init(mol,[35,6,8],xyz)
   call xbMolecularHessian(mol,data,12.0_wp,energy,g,h,ok,terms)
   call require(ok.and.terms==0,'acceptor-neighbour identity')
   call require(abs(energy)+maxval(abs(g))+maxval(abs(h))==0.0_wp,'identity angle contribution is zero')
   ! Equal-distance competing neighbours and exact cutoff are nonsmooth.
   xyz(:,2)=[-3.0_wp,0.0_wp,0.0_wp]; xyz(:,3)=[3.0_wp,0.0_wp,0.0_wp]
   call init(mol,[35,6,8],xyz)
   call xbMolecularHessian(mol,data,12.0_wp,energy,g,h,ok)
   call require(.not.ok,'nearest-neighbour tie rejects analytic Hessian')
   xyz(:,3)=[20.0_wp,0.0_wp,0.0_wp]
   call init(mol,[35,6,8],xyz)
   call xbMolecularHessian(mol,data,12.0_wp,energy,g,h,ok)
   call require(.not.ok,'cutoff branch rejects analytic Hessian')
   mol%xyz(:,3)=mol%xyz(:,1)
   call xbMolecularHessian(mol,data,12.0_wp,energy,g,h,ok)
   call require(.not.ok,'coincident pair rejects analytic Hessian')
   call xbMolecularHessian(water,missing,12.0_wp,energy,g,h,ok)
   call require(.not.ok,'missing halogen data rejects')
   print *, 'PASS halogen neighbour identity and branch/invalid-input fallback guards'
contains
   subroutine require(condition,label)
      logical,intent(in) :: condition
      character(len=*),intent(in) :: label
      if(.not.condition) then
         print *, 'FAIL ',label
         error stop 1
      endif
   end subroutine
end program
