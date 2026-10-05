! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Extended tight binding calculator
module xtb_xtb_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_model, only : info,newSolvationModel,newBornModel
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_cpx, only: TCpcmx
   use xtb_type_basisset, only : TBasisset
   use xtb_type_calculator, only : TCalculator, numericalHessian => hessian
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_param, only : scc_parameter, TxTBParameter, chrg_parameter
   use xtb_type_pcem
   use xtb_type_solvation, only : TSolvation
   use xtb_type_restart, only : TRestart
   use xtb_xtb_data, only : TxTBData
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam
   use xtb_scf, only : scf
   use xtb_peeq, only : peeq
   use xtb_embedding, only : read_pcem
   use xtb_metadynamic
   use xtb_constrainpot
   use xtb_basis, only : newBasisset
   use xtb_mctc_systools, only : rdpath
   use xtb_readparam, only : readParam
   use xtb_paramset, only : use_parameterset
   use xtb_chargemodel, only : new_charge_model_2019
   use xtb_disp_ncoord, only : ncoord_erf
   use xtb_eeq, only : eeq_chrgeq
   use xtb_iniq, only : iniqcn
   use xtb_scc_core, only : iniqshell
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy,getGFN1FastPolicy,getGFN1HessianResources
   use xtb_xtb_hessian_response, only : buildGFN1GasHessian,buildGFN1SolvatedHessian
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none

   private

   public :: TxTBCalculator, newXTBCalculator, newWavefunction


   !> Calculator interface for xTB based methods
   type, extends(TCalculator) :: TxTBCalculator

      !> Tight binding basis set
      type(TBasisset), allocatable :: basis

      !> Parametrisation data base
      type(TxTBData), allocatable :: xtbData

      !> Electronic temperature
      real(wp) :: etemp

      !> Maximum number of cycles for SCC convergence
      integer :: maxiter

      !> External potential
      type(tb_pcem) :: pcem

   contains

      !> Perform xTB single point calculation
      procedure :: singlepoint

      !> Analytic molecular GFN1 response with numerical model fallback
      procedure :: hessian

      !> Write informative printout
      procedure :: writeInfo

   end type TxTBCalculator

   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'


contains


subroutine hessian(self,env,mol0,chk0,list,step,hess,dipgrad,polgrad)
   class(TxTBCalculator), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol0
   type(TRestart), intent(in) :: chk0
   integer, intent(in) :: list(:)
   real(wp), intent(in) :: step
   real(wp), intent(inout) :: hess(:,:),dipgrad(:,:)
   real(wp), intent(inout), optional :: polgrad(:,:)
   type(TGFN1FastPolicy) :: policy
   type(TBorn) :: born
   real(wp), allocatable :: candidate(:,:),gradient(:,:),translation(:)
   real(wp) :: t0,w0,t1,w1,scale
   integer(int64) :: n,nao,ns,estimatedBytes
   integer :: a,k,col,allocationStatus,halogenTerms,responseThreads
   logical :: used,ok,failed,solvated
   character(len=64) :: reason

   call getGFN1FastPolicy(policy)
   used=.false.; reason='disabled'; estimatedBytes=0_int64;solvated=.false.
   halogenTerms=0
   call timing(t0,w0)
   analytic: block
      if(.not.policy%analyticHessian) exit analytic
      reason='unsupported model'
      if(.not.allocated(self%xtbData).or..not.allocated(self%basis)) exit analytic
      reason='method or periodic model'
      if(self%xtbData%level/=1.or.mol0%npbc/=0) exit analytic
      reason='legacy Gradient derivative policy'
      if(.not.policy%gradientKernel.or..not.policy%gradientDirectH0) exit analytic
      reason='solvation model'
      if(allocated(self%solvation)) then
         associate(model=>self%solvation)
            if(model%cosmo.or.model%tmcosmo.or.allocated(model%cpxsolvent)) exit analytic
            if(model%ionStrength/=0.0_wp) exit analytic
            if(model%kernel/=gbKernel%still.and.model%kernel/=gbKernel%p16) exit analytic
            solvated=.true.
         end associate
      endif
      reason='point charges or electric field'
      if(self%pcem%n>0.or.any(set%efield/=0.0_wp)) exit analytic
      reason='polarizability derivatives'
      if(present(polgrad)) exit analytic
      reason='constraints'
      if(nconstr>0.or.potset%n>0.or.potset%pos%n>0.or.potset%dist%n>0 &
         & .or.potset%angle%n>0.or.potset%dihedral%n>0) exit analytic
      reason='fixed or frozen atoms'
      if(fixset%n>0.or.freezeset%n>0.or.shakeset%n>0) exit analytic
      reason='wall potential'
      if(number_walls>0) exit analytic
      reason='metadynamic bias'
      if(metaset%nstruc>0.or.rmsdset%nstruc>0) exit analytic
      n=int(mol0%n,int64); nao=int(self%basis%nao,int64); ns=int(self%basis%nshell,int64)
      reason='invalid dimensions'
      if(n<1.or.nao<1.or.ns<1.or.any(list<1).or.any(list>mol0%n)) exit analytic
      if(any(shape(hess)/=[3*mol0%n,3*mol0%n]).or.any(shape(dipgrad)/=[3,3*mol0%n])) exit analytic
      ! Conservative bound for first AO tensors, CN second derivatives and
      ! simultaneous response workspaces. Do not allocate beyond this budget.
      reason='memory budget'
      call getGFN1HessianResources(mol0%n,self%basis%nao,self%basis%nshell,solvated,responseThreads,estimatedBytes)
      if(estimatedBytes>int(policy%analyticHessianMaxMiB,int64)*1048576_int64) exit analytic
      allocate(candidate(3*n,3*n),gradient(3,n),translation(3*n),stat=allocationStatus)
      reason='allocation failed'
      if(allocationStatus/=0) exit analytic
      reason='reference or response rejected'
      if(.not.ieee_is_finite(step).or.step<=0.0_wp) exit analytic
      if(solvated) then
         reason='solvent response or branch rejected'
         if(.not.ieee_is_finite(step).or.step<=0.0_wp) exit analytic
         call newBornModel(self%solvation,env,born,mol0%at)
         call born%update(env,mol0%at,mol0%xyz)
         call env%check(failed)
         if(failed) exit analytic
         call buildGFN1SolvatedHessian(env,mol0,self%basis,self%xtbData,chk0%wfn,self%etemp, &
            & self%accuracy,born,gradient,candidate,ok,halogenTerms,guardRadius=step)
      else
         call buildGFN1GasHessian(env,mol0,self%basis,self%xtbData,chk0%wfn,self%etemp, &
            & self%accuracy,gradient,candidate,ok,halogenTerms,guardRadius=step)
      endif
      if(.not.ok) exit analytic
      if(.not.all(ieee_is_finite(candidate))) exit analytic
      scale=max(1.0_wp,maxval(abs(candidate)))
      reason='Hessian reciprocity rejected'
      if(maxval(abs(candidate-transpose(candidate)))>1.0e-7_wp*scale) exit analytic
      reason='Hessian translation rejected'
      do a=1,3
         translation=0.0_wp; translation(a:3*n:3)=1.0_wp
         if(maxval(abs(matmul(candidate,translation)))>1.0e-7_wp*scale) exit analytic
      enddo
      ! Match the base API's additive Hessian and assigned dipole columns.
      ! The energy/Gradient-only GFN1 SCF deliberately returns zero dipoles.
      do k=1,size(list)
         do a=1,3
            col=3*(list(k)-1)+a
            hess(:,col)=hess(:,col)+candidate(:,col)
            dipgrad(:,col)=0.0_wp
         enddo
      enddo
      used=.true.; reason='analytic gas GFN1'
      if(solvated) reason='analytic solvated GFN1'
   end block analytic
   call timing(t1,w1)
   if(policy%profile) then
      write(env%unit,'(/,1x,a,l1)') 'GFN1-fast 4.7.0 analytic Hessian: used=',used
      write(env%unit,'(3x,a,a)') 'dispatch reason = ',trim(reason)
      write(env%unit,'(3x,a,i0)') 'analytic halogen terms = ',halogenTerms
      write(env%unit,'(3x,a,i0)') 'estimated workspace bytes = ',estimatedBytes
      write(env%unit,'(3x,a,f12.6,a)') 'analytic attempt wall = ',w1-w0,' s'
   endif
   if(.not.used) call numericalHessian(self,env,mol0,chk0,list,step,hess,dipgrad,polgrad)
end subroutine hessian


subroutine newXTBCalculator(env, mol, calc, fname, method, accuracy)

   character(len=*), parameter :: source = 'xtb_calculator_newXTBCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   type(TxTBCalculator), intent(out) :: calc

   integer, intent(in), optional :: method

   real(wp), intent(in), optional :: accuracy

   character(len=*), intent(in), optional :: fname

   character(len=:), allocatable :: filename
   type(TxTBParameter) :: globpar
   integer :: ich
   logical :: exist, okbas
   logical :: exitRun

   if (present(fname)) then
      filename = fname
   else
      if (present(method)) then
         select case(method)
         case(0)
            call rdpath(env%xtbpath, 'param_gfn0-xtb.txt', filename, exist)
            if (.not.exist) filename = 'param_gfn0-xtb.txt'
         case(1)
            call rdpath(env%xtbpath, 'param_gfn1-xtb.txt', filename, exist)
            if (.not.exist) filename = 'param_gfn1-xtb.txt'
         case(2)
            call rdpath(env%xtbpath, 'param_gfn2-xtb.txt', filename, exist)
            if (.not.exist) filename = 'param_gfn2-xtb.txt'
         end select
      end if
   end if
   if (.not.allocated(filename)) then
      call env%error("No parameter file or parametrisation info provided", source)
      return
   end if

   if (present(accuracy)) then
      calc%accuracy = accuracy
   else
      calc%accuracy = 1.0_wp
   end if

   calc%etemp = set%etemp
   calc%maxiter = set%maxscciter

   !> Obtain the parameter file
   allocate(calc%xtbData)
   call open_file(ich, filename, 'r')
   exist = ich /= -1
   if (exist) then
      call readParam(env, ich, globpar, calc%xtbData, .true.)
      call close_file(ich)
   else ! no parameter file, check if we have one compiled into the code
      call use_parameterset(filename, globpar, calc%xtbData, exist)
      if (.not.exist) then
         call env%error('Parameter file '//filename//' not found!', source)
         return
      end if
   endif

   if (present(method)) then
      if (method /= calc%xtbData%level) then
         call env%error("Requested method does not match loaded method", source)
         return
      end if
   end if

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Could not load parameters", source)
      return
   end if

   !> set up the basis set for the tb-Hamiltonian
   allocate(calc%basis)
   call newBasisset(calc%xtbData, mol%n, mol%at, calc%basis, okbas)
   if (.not.okbas) then
      call env%error('basis set could not be setup completely', source)
      return
   end if

   !> check for external point charge field
   if (allocated(set%pcem_file)) then
      call open_file(ich, set%pcem_file, 'r')
      if (ich /= -1) then
         call read_pcem(ich, env, calc%pcem, calc%xtbData%coulomb)
         call close_file(ich)
      end if
   end if

end subroutine newXTBCalculator


subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)

   !> Source of the generated errors
   character(len=*), parameter :: source = 'xtb_calculator_singlepoint'

   !> Calculator instance
   class(TxTBCalculator), intent(inout) :: self

   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   !> Molecular structure data
   type(TMolecule), intent(inout) :: mol

   !> Wavefunction data
   type(TRestart), intent(inout) :: chk

   !> Print level for IO
   integer, intent(in) :: printlevel

   !> Restart from previous results
   logical, intent(in) :: restart

   !> Total energy
   real(wp), intent(out) :: energy

   !> Molecular gradient
   real(wp), intent(out) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(out) :: sigma(:, :)

   !> HOMO-LUMO gap
   real(wp), intent(out) :: hlgap

   !> Detailed results
   type(scc_results), intent(out) :: results

   class(TSolvation), allocatable :: solvation
   type(TBorn), allocatable :: gbsa
   integer :: i,ich
   integer :: mode_sp_run = 1
   real(wp) :: efix
   logical :: inmol
   logical, parameter :: ccm = .true.
   logical :: exitRun

   call mol%update

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   hlgap = 0.0_wp
   efix = 0.0_wp

   ! ------------------------------------------------------------------------
   !  actual calculation
   select case(self%xtbData%level)
   case(1, 2)
      if (allocated(self%solvation)) then
         call newSolvationModel(self%solvation, env, solvation, mol%at)
      end if
      call scf(env,mol,chk%wfn,self%basis,self%pcem,self%xtbData,solvation, &
         &   hlgap,self%etemp,self%maxiter,printlevel,restart,.true., &
         &   self%accuracy,energy,gradient,results)

   case(0)
      if (allocated(self%solvation)) then
         allocate(gbsa)
         call newBornModel(self%solvation, env, gbsa, mol%at)
      end if
      call peeq &
         & (env,mol,chk%wfn,self%basis,self%xtbData,gbsa,hlgap,self%etemp, &
         &  printlevel,.true.,ccm,self%accuracy,energy,gradient,sigma,results)

   end select

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Electronic structure method terminated", source)
      return
   end if

   ! ------------------------------------------------------------------------
   !  post processing of gradient and energy

   ! point charge embedding gradient file
   if (allocated(set%pcem_grad) .and. self%pcem%n > 0) then
      call open_file(ich,set%pcem_grad,'w')
      do i=1,self%pcem%n
         write(ich,'(3f12.8)')self%pcem%grd(1:3,i)
      enddo
      call close_file(ich)
   endif

   ! ------------------------------------------------------------------------
   !  various external potentials
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (rmsdset,mol%n,mol%at,mol%xyz,efix,gradient)

   ! ------------------------------------------------------------------------
   !  fixing of certain atoms
   energy = energy + efix
   results%e_total = energy
   results%gnorm = norm2(gradient)
   if (fixset%n.gt.0) then
      do i=1, fixset%n
         gradient(1:3,fixset%atoms(i))=0
      enddo
   endif

   ! save point charge gradients in results
   if (self%pcem%n > 0) then
      results%pcem = self%pcem
   endif

   if (printlevel.ge.2) then
      ! start with summary header
      if (.not.set%silent) then
         write(env%unit,'(9x,53(":"))')
         write(env%unit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,outfmt) "total energy      ", results%e_total,"Eh   "
      if (.not.set%silent.and.allocated(self%solvation)) then
         write(env%unit,outfmt) "total w/o Gsasa/hb", &
            &  results%e_total-results%g_sasa-results%g_hb-results%g_shift, "Eh   "
      endif
      write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
      write(env%unit,outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
      if (.not.set%silent) then
         if (set%verbose) then
            write(env%unit,'(9x,"::",49("."),"::")')
            write(env%unit,outfmt) "HOMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo),  "eV   "
            write(env%unit,outfmt) "LUMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo+1),"eV   "
         endif
         write(env%unit,'(9x,"::",49("."),"::")')
         if (self%xtbData%level.eq.2) call print_gfn2_results(env%unit,results,set%verbose,allocated(self%solvation))
         if (self%xtbData%level.eq.1) call print_gfn1_results(env%unit,results,set%verbose,allocated(self%solvation))
         if (self%xtbData%level.eq.0) call print_gfn0_results(env%unit,results,set%verbose,allocated(self%solvation))
         write(env%unit,outfmt) "add. restraining  ", efix,       "Eh   "
         write(env%unit,outfmt) "total charge      ", sum(chk%wfn%q), "e    "
         if (set%verbose) then
            write(env%unit,'(9x,"::",49("."),"::")')
            write(env%unit,outfmt) "atomisation energy", results%e_atom, "Eh   "
         endif
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(a)')
   endif

end subroutine singlepoint


subroutine print_gfn0_results(iunit,res,verbose,lsolv)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lsolv
   write(iunit,outfmt) "H0 energy         ", res%e_elec, "Eh   "
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
   write(iunit,outfmt) "electrostat energy", res%e_es,   "Eh   "
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   !write(iunit,outfmt) "   -> Gborn       ", res%g_born, "Eh   " ! not saved
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   !write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   " ! not saved
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   write(iunit,outfmt) "dispersion energy ", res%e_disp, "Eh   "
   write(iunit,outfmt) "short-range corr. ", res%e_xb,   "Eh   "
end subroutine print_gfn0_results

subroutine print_gfn1_results(iunit,res,verbose,lsolv)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lsolv
   write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
   write(iunit,outfmt) "-> electrostatic  ", res%e_es,   "Eh   "
   if (lsolv) then
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   write(iunit,outfmt) "   -> Gelec       ", res%g_born, "Eh   "
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   endif
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
   write(iunit,outfmt) "dispersion energy ", res%e_disp, "Eh   "
   write(iunit,outfmt) "halogen bond corr.", res%e_xb,   "Eh   "
end subroutine print_gfn1_results

subroutine print_gfn2_results(iunit,res,verbose,lsolv)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lsolv
   write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
   write(iunit,outfmt) "-> isotropic ES   ", res%e_es,   "Eh   "
   write(iunit,outfmt) "-> anisotropic ES ", res%e_aes,  "Eh   "
   write(iunit,outfmt) "-> anisotropic XC ", res%e_axc,  "Eh   "
   write(iunit,outfmt) "-> dispersion     ", res%e_disp, "Eh   "
   if (lsolv) then
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   write(iunit,outfmt) "   -> Gelec       ", res%g_born, "Eh   "
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   endif
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
end subroutine print_gfn2_results


subroutine writeInfo(self, unit, mol)

   !> Calculator instance
   class(TxTBCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   call self%xtbData%writeInfo(unit, mol%at)

   if (allocated(self%solvation)) then
      call info(self%solvation, unit)
   end if

end subroutine writeInfo

!---------------------------------------------
! Initialize new wavefunction
!---------------------------------------------
subroutine newWavefunction(env, mol, calc, chk)
   
   implicit none
   character(len=*), parameter :: source = 'xtb_calculator_newWavefunction'
      !! Name of error producer routine
   type(TEnvironment), intent(inout) :: env
      !! Calculation environment to handle I/O stream and error log
   type(TRestart), intent(inout) :: chk
      !! Restart data wrapper for wfn and nlist      
   type(TxTBCalculator), intent(in) :: calc
      !! Instance of xTB Calculator
   type(TMolecule), intent(in) :: mol
      !! Molecular structure data
   real(wp), allocatable :: cn(:)
      !! Coordination number
   type(chrg_parameter) :: chrgeq
      !! guess charges(gasteiger/goedecker/sad) 
   logical :: exitRun
      !! if it is recommended to terminate the run 

   associate(wfn => chk%wfn)
      allocate(cn(mol%n))
      call wfn%allocate(mol%n,calc%basis%nshell,calc%basis%nao)
      
      !> find partial charges
      if (mol%npbc > 0) then
         !! if periodic
         wfn%q = mol%chrg/real(mol%n,wp)
            !! evenly distribute charge with the equal partial charges 
      else
         if (set%guess_charges.eq.p_guess_gasteiger) then
            call iniqcn(mol%n,mol%at,mol%z,mol%xyz,nint(mol%chrg),1.0_wp, &
               & wfn%q,cn,calc%xtbData%level,.true.)
         else if (set%guess_charges.eq.p_guess_goedecker) then
            !! default
            
            call new_charge_model_2019(chrgeq,mol%n,mol%at)
               !! to get parametrized values for q (en,gam,kappa,alpha)

            call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
               !! to obtain CN
               !! (49) Extended Tight-Binding Quantum Chemistry Mehods 2020
            
            call eeq_chrgeq(mol,env,chrgeq,cn,wfn%q)
               !! to obtain partial charges q
               !! (47) Extended Tight-Binding Quantum Chemistry Mehods 2020

            call env%check(exitRun)
               !! to check status of environment
            if (exitRun) then
               call env%rescue("EEQ guess failed, falling back to SAD guess", source)
               wfn%q = mol%chrg/real(mol%n,wp)
            end if
         else
            wfn%q = mol%chrg/real(mol%n,wp)
         end if
      end if
      
      !> find shell charges
      call iniqshell(calc%xtbData,mol%n,mol%at,mol%z,calc%basis%nshell, &
         & wfn%q,wfn%qsh,calc%xtbData%level)
   
   end associate

end subroutine newWavefunction

end module xtb_xtb_calculator
