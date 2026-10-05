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

!> Implementation of an isotropic electrostatics container
module xtb_xtb_coulomb
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : blas_dot, blas_symv
   use xtb_xtb_data, only : TCoulombData
   use xtb_xtb_thirdorder, only : TThirdOrder, init
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   public :: TxTBCoulomb, init


   !> Isotropic electrostatics
   type :: TxTBCoulomb

      !> Third order electrostatics
      type(TThirdOrder) :: thirdOrder

      !> Coulomb matrix
      real(wp), allocatable :: jmat(:, :)

      !> Scratch memory for storing the shifts
      real(wp), allocatable :: shift(:)

   contains

      !> Add shifts from isotropic electrostatics
      procedure :: addShift

      !> Get energy from isotropic electrostatics
      procedure :: getEnergy

      !> Fixed-geometry derivative of the total shell potential, in Eh.
      procedure :: getResponseKernel

   end type TxTBCoulomb


   !> Initialize isotropic electrostatics
   interface init
      module procedure :: initCoulomb
   end interface init


contains

!> GFN1 Hessian response: d(shellShift)/d(q_shell), including atomic
!> shifts expanded to shells. The Coulomb matrix's lower triangle is the
!> authoritative data, exactly as in addShift. Geometry/solvent derivatives
!> are supplied separately by the future complete Hessian assembler.
pure subroutine getResponseKernel(self,qat,qsh,sh2at,kernel,ok)
   class(TxTBCoulomb), intent(in) :: self
   real(wp), intent(in) :: qat(:),qsh(:)
   integer, intent(in) :: sh2at(:)
   real(wp), intent(out) :: kernel(:,:)
   logical, intent(out) :: ok
   integer :: n,ns,i,j,atom
   kernel=0.0_wp; ok=.false.; n=size(qat); ns=size(qsh)
   if(n<1.or.ns<1.or.size(sh2at)/=ns) return
   if(size(kernel,1)/=ns.or.size(kernel,2)/=ns) return
   if(.not.allocated(self%jmat)) return
   if(size(self%jmat,1)/=ns.or.size(self%jmat,2)/=ns) return
   if(any(sh2at<1).or.any(sh2at>n)) return
   if(.not.all(ieee_is_finite(qat)).or..not.all(ieee_is_finite(qsh))) return
   if(allocated(self%thirdOrder%atomicGam)) then
      if(size(self%thirdOrder%atomicGam)/=n) return
   endif
   if(allocated(self%thirdOrder%shellGam)) then
      if(size(self%thirdOrder%shellGam)/=ns) return
   endif
   do j=1,ns
      do i=j,ns
         kernel(i,j)=self%jmat(i,j)
         if(allocated(self%thirdOrder%atomicGam)) then
            atom=sh2at(i)
            if(atom==sh2at(j)) kernel(i,j)=kernel(i,j)+2.0_wp*qat(atom)*self%thirdOrder%atomicGam(atom)
         endif
         kernel(j,i)=kernel(i,j)
      enddo
      if(allocated(self%thirdOrder%shellGam)) &
         & kernel(j,j)=kernel(j,j)+2.0_wp*qsh(j)*self%thirdOrder%shellGam(j)
   enddo
   ok=all(ieee_is_finite(kernel))
   if(.not.ok) kernel=0.0_wp
end subroutine


!> Initialize isotropic electrostatics from parametrisation data
subroutine initCoulomb(self, input, nshell, num)

   !> Instance of the isotropic electrostatics
   type(TxTBCoulomb), intent(out) :: self

   !> Parametrisation data for coulombic interactions
   type(TCoulombData), intent(in) :: input

   !> Number of shells for each species
   integer, intent(in) :: nshell(:)

   !> Atomic numbers of each element
   integer, intent(in) :: num(:)

   integer :: nsh

   call init(self%thirdOrder, input, nshell, num)

   nsh = sum(nshell(num))
   allocate(self%jmat(nsh, nsh))
   allocate(self%shift(nsh))

end subroutine initCoulomb


!> Add shifts from isotropic electrostatics
subroutine addShift(self, qat, qsh, atomicShift, shellShift)

   !> Instance of the isotropic electrostatics
   class(TxTBCoulomb), intent(inout) :: self

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Atomic potential shift
   real(wp), intent(inout) :: atomicShift(:)

   !> Shell-resolved potential shift
   real(wp), intent(inout) :: shellShift(:)

   integer :: nsh

   nsh = size(shellShift)

   call self%thirdOrder%addShift(qat, qsh, atomicShift, shellShift)

   call blas_symv('l', nsh, 1.0_wp, self%jmat, nsh, qsh, 1, 1.0_wp, shellShift, 1)

end subroutine addShift


!> Get energy from isotropic electrostatics
pure subroutine getEnergy(self, qat, qsh, energy)

   !> Instance of the isotropic electrostatics
   class(TxTBCoulomb), intent(inout) :: self

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Third order contribution to the energy
   real(wp), intent(out) :: energy

   integer :: nsh
   real(wp) :: eThird

   nsh = size(qsh)

   call self%thirdOrder%getEnergy(qat, qsh, eThird)

   call blas_symv('l', nsh, 1.0_wp, self%jmat, nsh, qsh, 1, 0.0_wp, self%shift, 1)
   energy = 0.5_wp * blas_dot(nsh, self%shift, 1, qsh, 1) + eThird

end subroutine getEnergy


end module xtb_xtb_coulomb
