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

!> Implementation of the repulsion energy used in the xTB Hamiltonian
module xtb_xtb_repulsion
   use xtb_mctc_accuracy, only : wp
   use xtb_type_identitymap, only : TIdentityMap
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_neighbourlist, only : TNeighbourlist
   use xtb_xtb_data
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   public :: repulsionEnGrad
   public :: repulsionMolecularHessian


   interface repulsionEnGrad
      module procedure :: repulsionEnGrad_latp
      module procedure :: repulsionEnGrad_neighs
   end interface repulsionEnGrad


contains

! Molecular second derivative of the same radial repulsion used by
! repulsionEnGrad_latp. Outputs are fresh contributions, in atomic units.
pure subroutine repulsionMolecularHessian(mol,repData,cutoff,energy,gradient,hessian,ok,guardRadius)
   type(TMolecule), intent(in) :: mol
   type(TRepulsionData), intent(in) :: repData
   real(wp), intent(in) :: cutoff
   real(wp), intent(out) :: energy,gradient(:,:),hessian(:,:)
   logical, intent(out) :: ok
   real(wp), intent(in), optional :: guardRadius
   real(wp) :: vec(3),r2,r,power,alpha,zeff,kexp,value,amplitude,radial,anisotropic,hh(3,3),dg(3),window
   integer :: n,iat,jat,iz,jz,a,b,ci,cj,di,dj
   energy=0.0_wp; gradient=0.0_wp; hessian=0.0_wp; ok=.false.; n=mol%n
   if(n<1.or.mol%npbc/=0) return
   if(any(shape(gradient)/=[3,n]).or.any(shape(hessian)/=[3*n,3*n])) return
   if(.not.ieee_is_finite(cutoff).or.cutoff<=0.0_wp) return
   window=0.0_wp
   if(present(guardRadius)) window=guardRadius
   if(.not.ieee_is_finite(window).or.window<0.0_wp) return
   if(.not.allocated(repData%alpha).or..not.allocated(repData%zeff)) return
   if(any(mol%at<1).or.any(mol%at>min(size(repData%alpha),size(repData%zeff)))) return
   if(.not.all(ieee_is_finite(mol%xyz))) return
   if(.not.all(ieee_is_finite(repData%alpha(mol%at))).or.any(repData%alpha(mol%at)<=0.0_wp)) return
   if(.not.all(ieee_is_finite(repData%zeff(mol%at)))) return
   if(.not.ieee_is_finite(repData%kExp).or..not.ieee_is_finite(repData%kExpLight) &
      & .or..not.ieee_is_finite(repData%rExp)) return
   do iat=1,n
      iz=mol%at(iat)
      do jat=1,iat-1
         jz=mol%at(jat); vec=mol%xyz(:,iat)-mol%xyz(:,jat); r2=sum(vec**2)
         if(r2<1.0e-8_wp) return
         if(abs(sqrt(r2)-cutoff)<=window+1.0e-10_wp*max(1.0_wp,cutoff)) return
         if(r2>cutoff**2) cycle
         r=sqrt(r2); alpha=sqrt(repData%alpha(iz)*repData%alpha(jz)); zeff=repData%zeff(iz)*repData%zeff(jz)
         kexp=repData%kExp
         if(iz<=2.and.jz<=2) kexp=repData%kExpLight
         power=r**kexp; value=zeff*exp(-alpha*power)/r**repData%rExp
         amplitude=alpha*power*kexp+repData%rExp
         radial=-amplitude*value/r2
         anisotropic=(amplitude**2+2.0_wp*amplitude-alpha*power*kexp**2)*value/r2**2
         dg=radial*vec; hh=anisotropic*spread(vec,2,3)*spread(vec,1,3)
         do a=1,3; hh(a,a)=hh(a,a)+radial; enddo
         energy=energy+value; gradient(:,iat)=gradient(:,iat)+dg; gradient(:,jat)=gradient(:,jat)-dg
         do a=1,3
            ci=3*(iat-1)+a; cj=3*(jat-1)+a
            do b=1,3
               di=3*(iat-1)+b; dj=3*(jat-1)+b
               hessian(ci,di)=hessian(ci,di)+hh(a,b); hessian(cj,dj)=hessian(cj,dj)+hh(a,b)
               hessian(ci,dj)=hessian(ci,dj)-hh(a,b); hessian(cj,di)=hessian(cj,di)-hh(a,b)
            enddo
         enddo
      enddo
   enddo
   ok=ieee_is_finite(energy).and.all(ieee_is_finite(gradient)).and.all(ieee_is_finite(hessian))
end subroutine repulsionMolecularHessian


!> Lattice point based implementation of the repulsion energy
subroutine repulsionEnGrad_latp(mol, repData, trans, cutoff, energy, gradient, &
      & sigma)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Repulsion parametrisation
   type(TRepulsionData), intent(in) :: repData

   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(inout) :: sigma(:, :)

   !> Repulsion energy
   real(wp), intent(inout) :: energy

   integer  :: iat, jat, iZp, jZp, itr, k, l
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, zeff, kExp, cutoff2
   real(wp) :: r1, r2, rij(3), dS(3, 3), dG(3), dE
   real(wp), allocatable :: energies(:)

   cutoff2 = cutoff**2

   allocate(energies(mol%n))

   !$acc enter data create(energies, rij, dS, dG) copyin(gradient, sigma, mol, mol%at, &
   !$acc& mol%xyz, repData, repData%alpha, repData%zeff, trans)

   !$acc kernels default(present)
   energies(:) = 0.0_wp
   !$acc end kernels

#ifdef XTB_GPU
   !$acc parallel default(present)
   !$acc loop gang collapse(2) private(iZp, jZp, alpha, zeff, kExp)
#else
   !$omp parallel do default(none) reduction(+:energies, gradient, sigma) &
   !$omp shared(repData, mol, cutoff2, trans) &
   !$omp private(iat, jat, itr, iZp, jZp, r2, rij, r1, alpha, zeff, kexp, &
   !$omp& t16, t26, t27, dE, dG, dS, k, l)
#endif
   do iAt = 1, mol%n
      do jAt = 1, mol%n
         if (jAt > iAt) cycle
         iZp = mol%at(iAt)
         jZp = mol%at(jAt)
         alpha = sqrt(repData%alpha(iZp)*repData%alpha(jZp))
         zeff = repData%zeff(iZp)*repData%zeff(jZp)
         if (iZp > 2 .or. jZp > 2) then
            kExp = repData%kExp
         else
            kExp = repData%kExpLight
         end if
         !$acc loop vector private(rij, r2, r1, t16, t26, t27, dE, dG, dS)
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iAt) - mol%xyz(:, jAt) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-8_wp) cycle
            r1 = sqrt(r2)

            t16 = r1**kExp
            t26 = exp(-alpha*t16)
            t27 = r1**repData%rExp
            dE = zeff * t26/t27
            dG = -(alpha*t16*kExp + repData%rExp) * dE * rij/r2
            dS = spread(dG, 1, 3) * spread(rij, 2, 3)
            !$acc atomic
            energies(iAt) = energies(iAt) + 0.5_wp * dE
            if (iAt /= jAt) then
               !$acc atomic
               energies(jAt) = energies(jAt) + 0.5_wp * dE
               !$acc loop seq
               do k = 1, 3
                  !$acc atomic
                  gradient(k, iAt) = gradient(k, iAt) + dG(k)
                  !$acc atomic
                  gradient(k, jAt) = gradient(k, jAt) - dG(k)
               end do
               !$acc loop seq collapse(2)
               do k = 1, 3
                  do l = 1, 3
                     !$acc atomic
                     sigma(l, k) = sigma(l, k) + dS(l, k)
                  end do
               end do
            else
               !$acc loop seq collapse(2)
               do k = 1, 3
                  do l = 1, 3
                     !$acc atomic
                     sigma(l, k) = sigma(l, k) + 0.5_wp * dS(l, k)
                  end do
               end do
            endif
         enddo
      enddo
   enddo
#ifdef XTB_GPU
   !$acc end parallel

   !$acc exit data copyout(energies, gradient, sigma) delete(mol, mol%at, &
   !$acc& mol%xyz, repData, repData%alpha, repData%zeff, rij, dG, dS, trans)
#endif

   energy = energy + sum(energies)

end subroutine repulsionEnGrad_latp


!> Lattice point based implementation of the repulsion energy
subroutine repulsionEnGrad_neighs(mol, repData, neighs, neighList, energy, &
      & gradient, sigma)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Repulsion parametrisation
   type(TRepulsionData), intent(in) :: repData

   !> Number of neighbours for each atom
   integer, intent(in) :: neighs(:)

   !> Neighbourlist
   class(TNeighbourList), intent(in) :: neighList

   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(inout) :: sigma(:, :)

   !> Repulsion energy
   real(wp), intent(inout) :: energy

   integer  :: iat, jat, iZp, jZp, ij, img
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, zeff, kExp, cutoff2
   real(wp) :: r1, r2, rij(3), dS(3, 3), dG(3), dE
   real(wp), allocatable :: energies(:)

   allocate(energies(len(mol)))
   energies(:) = 0.0_wp

   !$omp parallel do default(none) reduction(+:energies, gradient, sigma) &
   !$omp shared(repData, mol, neighs, neighList) &
   !$omp private(iat, jat, ij, img, iZp, jZp, r2, rij, r1, alpha, zeff, kexp, &
   !$omp& t16, t26, t27, dE, dG, dS)
   do iAt = 1, len(mol)
      iZp = mol%at(iAt)
      do ij = 1, neighs(iAt)
         img = neighList%ineigh(ij, iAt)
         jAt = neighList%image(img)
         jZp = mol%at(jAt)
         rij = neighList%coords(:, iAt) - neighList%coords(:, img)
         r2 = neighList%dist2(ij, iAt)
         r1 = sqrt(r2)

         alpha = sqrt(repData%alpha(iZp)*repData%alpha(jZp))
         zeff = repData%zeff(iZp)*repData%zeff(jZp)
         if (iZp > 2 .or. jZp > 2) then
            kExp = repData%kExp
         else
            kExp = repData%kExpLight
         end if
         t16 = r1**kExp
         t26 = exp(-alpha*t16)
         t27 = r1**repData%rExp
         dE = zeff * t26/t27
         dG = -(alpha*t16*kExp + repData%rExp) * dE * rij/r2
         dS = spread(dG, 1, 3) * spread(rij, 2, 3)
         energies(iAt) = energies(iAt) + 0.5_wp * dE
         sigma = sigma + 0.5_wp * dS
         if (iAt /= jAt) then
            energies(jAt) = energies(jAt) + 0.5_wp * dE
            sigma = sigma + 0.5_wp * dS
            gradient(:, iAt) = gradient(:, iAt) + dG
            gradient(:, jAt) = gradient(:, jAt) - dG
         endif
      enddo
   enddo
   !$omp end parallel do

   energy = energy + sum(energies)

end subroutine repulsionEnGrad_neighs


end module xtb_xtb_repulsion
