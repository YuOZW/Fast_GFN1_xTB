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

!> Implementation of the xTB core Hamiltonian
module xtb_xtb_hamiltonian
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_convert, only : evtoau
   use xtb_xtb_data, only : THamiltonianData
   use xtb_intgrad
   use xtb_lin
   use xtb_scc_core, only : shellPoly, h0scal
   use xtb_grad_core, only : dshellPoly
   implicit none
   private

   public :: getSelfEnergy, build_SH0_GFN1, build_SDQH0, build_dSDQH0, build_dSDQH0_noreset
   public :: build_dSH0_GFN1_noreset
   public :: count_dpint, count_qpint


   interface getSelfEnergy
      module procedure :: getSelfEnergyFlat
      module procedure :: getSelfEnergy2D
   end interface getSelfEnergy


   integer,private, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer,private, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

contains


subroutine getSelfEnergyFlat(hData, nShell, at, cn, qat, selfEnergy, dSEdcn, dSEdq)
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: at(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfEnergy(:)
   real(wp), intent(out), optional :: dSEdcn(:)
   real(wp), intent(out), optional :: dSEdq(:)

   integer :: ind, iAt, iZp, iSh

   selfEnergy(:) = 0.0_wp
   if (present(dSEdcn)) dSEdcn(:) = 0.0_wp
   if (present(dSEdq)) dSEdq(:) = 0.0_wp
   ind = 0
   do iAt = 1, size(cn)
      iZp = at(iAt)
      do iSh = 1, nShell(iZp)
         selfEnergy(ind+iSh) = hData%selfEnergy(iSh, iZp)
      end do
      ind = ind + nShell(iZp)
   end do
   if (present(dSEdcn) .and. present(cn)) then
      ind = 0
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(ind+iSh) = selfEnergy(ind+iSh) &
               & - hData%kCN(iSh, iZp) * cn(iAt)
            dSEdcn(ind+iSh) = -hData%kCN(iSh, iZp)
         end do
         ind = ind + nShell(iZp)
      end do
   end if
   if (present(dSEdq) .and. present(qat)) then
      ind = 0
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(ind+iSh) = selfEnergy(ind+iSh) &
               & - hData%kQShell(iSh,iZp)*qat(iAt) - hData%kQAtom(iZp)*qat(iAt)**2
            dSEdq(ind+iSh) = -hData%kQShell(iSh,iZp) - hData%kQAtom(iZp)*2*qat(iAt)
         end do
         ind = ind + nShell(iZp)
      end do
   end if

end subroutine getSelfEnergyFlat


subroutine getSelfEnergy2D(hData, nShell, at, cn, qat, selfEnergy, dSEdcn, dSEdq)
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: at(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfEnergy(:, :)
   real(wp), intent(out), optional :: dSEdcn(:, :)
   real(wp), intent(out), optional :: dSEdq(:, :)

   integer :: iAt, iZp, iSh

   selfEnergy(:, :) = 0.0_wp
   if (present(dSEdcn)) dSEdcn(:, :) = 0.0_wp
   if (present(dSEdq)) dSEdq(:, :) = 0.0_wp
   do iAt = 1, size(cn)
      iZp = at(iAt)
      do iSh = 1, nShell(iZp)
         selfEnergy(iSh, iAt) = hData%selfEnergy(iSh, iZp)
      end do
   end do
   if (present(dSEdcn) .and. present(cn)) then
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(iSh, iAt) = selfEnergy(iSh, iAt) &
               & - hData%kCN(iSh, iZp) * cn(iAt)
            dSEdcn(iSh, iAt) = -hData%kCN(iSh, iZp)
         end do
      end do
   end if
   if (present(dSEdq) .and. present(qat)) then
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(iSh, iAt) = selfEnergy(iSh, iAt) &
               & - hData%kQShell(iSh,iZp)*qat(iAt) - hData%kQAtom(iZp)*qat(iAt)**2
            dSEdq(iSh, iAt) = -hData%kQShell(iSh,iZp) &
               & - hData%kQAtom(iZp)*2*qat(iAt)
         end do
      end do
   end if

end subroutine getSelfEnergy2D


!> GFN1-only overlap/core-Hamiltonian builder.
!>
!> GFN1 has no anisotropic multipole SCC term and Energy/Gradient/Hessian do
!> not require dipole/quadrupole moment integrals.  Avoid computing and storing
!> those nine NAO x NAO tensors and evaluate only the overlap recurrence needed
!> by S and H0.  The algebra, screening threshold, shell polynomial and H0
!> accumulation order are the same as build_SDQH0.
subroutine build_SH0_GFN1(nShell, hData, nat, at, nbf, nao, xyz, trans, selfEnergy, &
      & intcut, caoshell, saoshell, nprim, primcount, alp, cont, sint, H0)
   implicit none
   integer, intent(in) :: nShell(:), nat, nao, nbf
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: at(nat)
   real(wp), intent(in) :: xyz(3,nat), trans(:, :), selfEnergy(:, :)
   real(wp), intent(in) :: intcut
   integer, intent(in) :: caoshell(:,:), saoshell(:,:), nprim(:), primcount(:)
   real(wp), intent(in) :: alp(:), cont(:)
   real(wp), intent(out) :: sint(nao,nao), H0(:)

   integer :: iat,jat,izp,jzp,ish,jsh,ishtyp,jshtyp,icao,jcao,naoi,naoj
   integer :: ii,jj,iao,jao,ij,il,jl,itr
   real(wp) :: ra(3),rb(3),ss(6,6)
   real(wp) :: zi,zj,zetaij,km,hii,hjj,hav,shpoly
   real(wp), parameter :: origin(3) = 0.0_wp

   H0(:) = 0.0_wp
   sint(:,:) = 0.0_wp

   !$omp parallel do if(nao >= 64) default(none) schedule(dynamic) &
   !$omp shared(nat,xyz,at,nShell,hData,selfEnergy,caoshell,saoshell,nprim, &
   !$omp& primcount,alp,cont,intcut,trans,sint,H0) &
   !$omp private(iat,jat,izp,jzp,ra,rb,ish,jsh,ishtyp,jshtyp,icao,jcao, &
   !$omp& naoi,naoj,ii,jj,iao,jao,ij,il,jl,itr,ss,zi,zj,zetaij, &
   !$omp& km,hii,hjj,hav,shpoly)
   do iat = 1, nat
      ra = xyz(:,iat)
      izp = at(iat)
      do jat = 1, iat-1
         jzp = at(jat)
         do ish = 1, nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            do jsh = 1, nShell(jzp)
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               il = ishtyp+1
               jl = jshtyp+1
               hii = selfEnergy(ish,iat)
               hjj = selfEnergy(jsh,jat)
               zi = hData%slaterExponent(ish,izp)
               zj = hData%slaterExponent(jsh,jzp)
               zetaij = (2*sqrt(zi*zj)/(zi+zj))**hData%wExp
               call h0scal(hData,il,jl,izp,jzp,hData%valenceShell(ish,izp).ne.0, &
                  & hData%valenceShell(jsh,jzp).ne.0,km)
               hav = 0.5_wp*km*(hii+hjj)*zetaij
               do itr = 1, size(trans,dim=2)
                  rb = xyz(:,jat) + trans(:,itr)
                  shpoly = shellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp), &
                     & hData%atomicRad(izp),hData%atomicRad(jzp),ra,rb)
                  ss = 0.0_wp
                  call get_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,rb, &
                     & origin,intcut,nprim,primcount,alp,cont,ss)
                  call dtrf2(ss,ishtyp,jshtyp)
                  do ii = 1,llao2(ishtyp)
                     iao = ii+saoshell(ish,iat)
                     do jj = 1,llao2(jshtyp)
                        jao = jj+saoshell(jsh,jat)
                        ij = lin(iao,jao)
                        H0(ij) = H0(ij) + hav*shpoly*ss(jj,ii)
                        sint(jao,iao) = sint(jao,iao) + ss(jj,ii)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

   !$omp parallel do if(nao >= 64) default(none) shared(nao,sint) private(iao,jao)
   do iao = 1, nao
      do jao = 1, iao-1
         sint(iao,jao) = sint(jao,iao)
      end do
   end do

   !$omp parallel do if(nao >= 64) default(none) shared(nat,nao,at,nShell,hData,saoshell,selfEnergy,sint,H0) &
   !$omp private(iat,izp,ish,ishtyp,iao,ii,ij)
   do iat = 1, nat
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1,llao2(ishtyp)
            ii = iao+saoshell(ish,iat)
            ij = ii*(1+ii)/2
            sint(ii,ii) = 1.0_wp + sint(ii,ii)
            H0(ij) = H0(ij) + selfEnergy(ish,iat)
         end do
      end do
   end do
end subroutine build_SH0_GFN1


!> Computes the dipole and quadrupole integrals and performs screening to
!  determine, which contribute to potential
subroutine build_SDQH0(nShell, hData, nat, at, nbf, nao, xyz, trans, selfEnergy, &
      & intcut, caoshell, saoshell, nprim, primcount, alp, cont, &
      & sint, dpint, qpint, H0)
   implicit none
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   !> # of atoms
   integer, intent(in)  :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)  :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(nat)
   !> Cartesian coordinates
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)    :: trans(:, :)
   real(wp), intent(in) :: selfEnergy(:, :)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)  :: intcut
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)  :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   !> Overlap integral matrix
   real(wp),intent(out) :: sint(nao,nao)
   !> Dipole integral matrix
   real(wp),intent(out) :: dpint(3,nao,nao)
   !> Quadrupole integral matrix
   real(wp),intent(out) :: qpint(6,nao,nao)
   !> Core Hamiltonian
   real(wp),intent(out) :: H0(:)


   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rab2,ab,est

   real(wp)  ra(3),rb(3),f1,f2,point(3)
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,jshmax
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer :: il, jl, itr
   real(wp) :: zi, zj, zetaij, km, hii, hjj, hav, shpoly
   integer itt(0:3)
   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)


   ! integrals
   H0(:) = 0.0_wp
   sint = 0.0_wp
   dpint = 0.0_wp
   qpint = 0.0_wp
   ! --- Aufpunkt for moment operator
   point = 0.0_wp

   !$omp parallel do default(none) schedule(dynamic) &
   !$omp shared(nat, xyz, at, nShell, hData, selfEnergy, caoshell, saoshell, &
   !$omp& nprim, primcount, alp, cont, intcut, trans, point) &
   !$omp private (iat,jat,izp,ci,ra,rb,saw, &
   !$omp& rab2,jzp,ish,ishtyp,icao,naoi,iptyp, &
   !$omp& jsh,jshmax,jshtyp,jcao,naoj,jptyp,ss,dd,qq,shpoly, &
   !$omp& est,alpi,alpj,ab,iprim,jprim,ip,jp,il,jl,hii,hjj,km,zi,zj,zetaij,hav, &
   !$omp& mli,mlj,tmp,tmp1,tmp2,iao,jao,ii,jj,k,ij,itr) &
   !$omp shared(sint,dpint,qpint,H0)
   do iat = 1, nat
      ra(1:3) = xyz(1:3,iat)
      izp = at(iat)
      do jat = 1, iat-1
         jzp = at(jat)
         do ish = 1, nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            do jsh = 1, nShell(jzp)
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, izp)
               zj = hData%slaterExponent(jsh, jzp)
               zetaij = (2 * sqrt(zi*zj)/(zi+zj))**hData%wExp
               call h0scal(hData,il,jl,izp,jzp,hData%valenceShell(ish, izp).ne.0, &
                  & hData%valenceShell(jsh, jzp).ne.0,km)

               hav = 0.5_wp * km * (hii + hjj) * zetaij

               do itr = 1, size(trans, dim=2)
                  rb(1:3) = xyz(1:3,jat) + trans(:, itr)
                  rab2 = sum( (rb-ra)**2 )

                  ! distance dependent polynomial
                  shpoly=shellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp),&
                     &             hData%atomicRad(izp),hData%atomicRad(jzp),ra,rb)

                  ss = 0.0_wp
                  dd = 0.0_wp
                  qq = 0.0_wp
                  call get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,rb,point, &
                     &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
                  !transform from CAO to SAO
                  call dtrf2(ss,ishtyp,jshtyp)
                  do k = 1,3
                     tmp(1:6,1:6) = dd(k,1:6,1:6)
                     call dtrf2(tmp,ishtyp,jshtyp)
                     dd(k,1:6,1:6) = tmp(1:6,1:6)
                  enddo
                  do k = 1,6
                     tmp(1:6,1:6) = qq(k,1:6,1:6)
                     call dtrf2(tmp,ishtyp,jshtyp)
                     qq(k,1:6,1:6) = tmp(1:6,1:6)
                  enddo
                  do ii = 1,llao2(ishtyp)
                     iao = ii+saoshell(ish,iat)
                     do jj = 1,llao2(jshtyp)
                        jao = jj+saoshell(jsh,jat)
                        ij = lin(iao, jao)
                        H0(ij) = H0(ij) + hav * shpoly * ss(jj, ii)
                        !sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                        sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                        !dpint(:, iao, jao) = dpint(:, iao, jao) + dd(:, jj, ii)
                        dpint(:, jao, iao) = dpint(:, jao, iao) + dd(:, jj, ii)
                        !qpint(:, iao, jao) = qpint(:, iao, jao) + qq(:, jj, ii)
                        qpint(:, jao, iao) = qpint(:, jao, iao) + qq(:, jj, ii)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp parallel do default(none) shared(nao, sint, dpint, qpint) private(iao, jao)
   do iao = 1, nao
      do jao = 1, iao - 1
         sint(iao, jao) = sint(jao, iao)
         dpint(:, iao, jao) = dpint(:, jao, iao)
         qpint(:, iao, jao) = qpint(:, jao, iao)
      end do
   end do

   ! diagonal elements
   !$omp parallel do default(none) schedule(dynamic) &
   !$omp shared(H0, sint, dpint, qpint) &
   !$omp shared(nat, xyz, at, nShell, hData, saoshell, selfEnergy, caoshell, &
   !$omp& point, intcut, nprim, primcount, alp, cont) &
   !$omp private(iat, ra, izp, ish, ishtyp, iao, i, ii, icao, naoi, iptyp, &
   !$omp& jsh, jshtyp, jcao, ss, dd, qq, k, tmp, jao, jj, naoj, jptyp)
   do iat = 1, nat
      ra = xyz(:, iat)
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)
            ii = i*(1+i)/2
            sint(i,i) = 1.0_wp + sint(i,i)
            H0(ii) = H0(ii) + selfEnergy(ish, iat)
         end do

         icao = caoshell(ish,iat)
         naoi = llao(ishtyp)
         iptyp = itt(ishtyp)
         do jsh = 1, ish
            jshtyp = hData%angShell(jsh,izp)
            jcao = caoshell(jsh,iat)
            naoj = llao(jshtyp)
            jptyp = itt(jshtyp)
            ss = 0.0_wp
            dd = 0.0_wp
            qq = 0.0_wp
            call get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,ra,point, &
               &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
            !transform from CAO to SAO
            !call dtrf2(ss,ishtyp,jshtyp)
            do k = 1,3
               tmp(1:6, 1:6) = dd(k,1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               dd(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do k = 1,6
               tmp(1:6, 1:6) = qq(k, 1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               qq(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do ii = 1, llao2(ishtyp)
               iao = ii + saoshell(ish,iat)
               do jj = 1, llao2(jshtyp)
                  jao = jj + saoshell(jsh,iat)
                  if (jao > iao .and. ish ==  jsh) cycle
                  dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dd(1:3, jj, ii)
                  if (iao /= jao) then
                     dpint(1:3, jao, iao) = dpint(1:3, jao, iao) + dd(1:3, jj, ii)
                  end if
                  qpint(1:6, iao, jao) = qpint(1:6, iao, jao) + qq(1:6, jj, ii)
                  if (jao /= iao) then
                     qpint(1:6, jao, iao) = qpint(1:6, jao, iao) + qq(1:6, jj, ii)
                  end if
               end do
            end do
         end do
      end do
   end do

end subroutine build_SDQH0


!> Computes the gradient of the dipole/qpole integral contribution
subroutine build_dSDQH0(nShell, hData, selfEnergy, dSEdcn, intcut, nat, nao, nbf, &
      & at, xyz, trans, caoshell, saoshell, nprim, primcount, alp, cont, &
      & p, Pew, ves, vs, vd, vq, dhdcn, g, sigma)
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   real(wp), intent(in) :: selfEnergy(:, :)
   real(wp), intent(in) :: dSEdcn(:, :)
   !> # of atoms
   integer, intent(in)    :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)    :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)    :: nbf
   !> Atomic numbers of atoms
   integer, intent(in)    :: at(nat)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)    :: intcut
   real(wp),intent(in)    :: ves(:, :)
   real(wp),intent(in)    :: vs(nat)
   real(wp),intent(in)    :: vd(3,nat)
   real(wp),intent(in)    :: vq(6,nat)
   !> Cartesian coordinates
   real(wp),intent(in)    :: xyz(3,nat)
   real(wp),intent(in)    :: trans(:, :)
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)    :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)    :: saoshell(:,:)
   integer, intent(in)    :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO, dimension: nbf
   integer, intent(in)    :: primcount(:)
   real(wp),intent(in)    :: alp(:)
   real(wp),intent(in)    :: cont(:)
   !> Density matrix
   real(wp),intent(in) :: p(:, :)
   !> Energy weighted density matrix
   real(wp),intent(in) :: Pew(:, :)
   real(wp),intent(inout) :: g(:, :)
   real(wp),intent(inout) :: sigma(:, :)
   real(wp),intent(inout) :: dhdcn(:)

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cj,alpi,rij2,ab,est
   real(wp) f1,f2,point(3),tmp(6,6),rij(3),ri(3),rj(3)
   real(wp) stmp,ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,jshmax
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,ixyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3,10),dum(10),sdq(10,6,6),sdqg(3,19,6,6)
   integer :: il, jl, itr
   real(wp) :: zi, zj, zetaij, km, hii, hjj, hav, shpoly, dshpoly(3)
   real(wp) :: Pij, Hij, HPij, g_xyz(3)
   real(wp), parameter :: rthr = 1600.0_wp

   thr2 = intcut
   point = 0.0_wp
   ! call timing(t1,t3)
   !$omp parallel do default(none) schedule(dynamic) &
   !$omp shared(nat, at, xyz, trans, nShell, hData, selfEnergy, dSEdcn, P, Pew, &
   !$omp& ves, vs, vd, vq, intcut, nprim, primcount, caoshell, saoshell, alp, cont) &
   !$omp private(iat,jat,ixyz,izp,ci,rij2,jzp,ish,ishtyp, &
   !$omp& icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp& sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,ri,rj,rij,km,shpoly,dshpoly, &
   !$omp& mli,mlj,dum,dumdum,tmp,stmp,dtmp,qtmp,il,jl,zi,zj,zetaij,hii,hjj,hav, &
   !$omp& iao,jao,ii,jj,k,pij,hij,hpij,g_xyz,itr) &
   !$omp reduction(+:g,sigma,dhdcn)
   do iat = 1,nat
      ri = xyz(:,iat)
      izp = at(iat)
      do jat = 1,iat-1
         !           if (jat.eq.iat) cycle
         jzp = at(jat)

         do ish = 1,nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(jzp)
            !              if(iat.eq.jat) jshmax = ish
            do jsh = 1,jshmax ! jshells
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, izp)
               zj = hData%slaterExponent(jsh, jzp)
               zetaij = (2 * sqrt(zi*zj)/(zi+zj))**hData%wExp
               call h0scal(hData,il,jl,izp,jzp,hData%valenceShell(ish, izp).ne.0, &
                  & hData%valenceShell(jsh, jzp).ne.0,km)

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj) * zetaij * evtoau

               do itr = 1, size(trans, dim=2)
                  rj = xyz(:,jat) + trans(:, itr)
                  rij = ri - rj
                  rij2 =  sum( rij**2 )

                  if (rij2 > rthr) cycle

                  ! distance dependent polynomial
                  call dshellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp),&
                     & hData%atomicRad(izp),hData%atomicRad(jzp),rij2,ri,rj,&
                     & shpoly,dshpoly)

                  sdqg = 0;sdq = 0
                  call get_grad_multiint(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
                     &                   intcut,nprim,primcount,alp,cont,sdq,sdqg)
                  tmp(1:6,1:6) = sdq(1,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  sdq(1,1:6,1:6) = tmp(1:6,1:6)
                  do k = 1,19 ! 1 S, 2-4 D, 5-10 Q, 11-13 D, 14-19 Q
                     do ixyz = 1,3
                        ! transform from CAO to SAO
                        !                        call dtrf2(sdqg(ixyz,k,1:6,1:6),ishtyp,jshtyp)
                        tmp(1:6,1:6) = sdqg(ixyz,k,1:6,1:6)
                        call dtrf2(tmp,ishtyp,jshtyp)
                        sdqg(ixyz,k,1:6,1:6) = tmp(1:6,1:6)
                     enddo
                  enddo
                  g_xyz(:) = 0.0_wp
                  do ii = 1,llao2(ishtyp)
                     iao = ii+saoshell(ish,iat)
                     do jj = 1,llao2(jshtyp)
                        jao = jj+saoshell(jsh,jat)
                        Pij = p(jao,iao)

                        ! Hamiltonian element without overlap
                        Hij  = hav * shpoly
                        HPij = Hij * Pij

                        g_xyz(:) = g_xyz + 2*HPij*sdq(1,jj,ii)*dshpoly/shpoly

                        do ixyz = 1,3
                           stmp = sdqg(ixyz,1,jj,ii)*(2*HPij - 2*Pew(jao, iao) &
                              & -Pij*(ves(ish,iat)+ves(jsh,jat)) &
                              & +Pij*(vs(iat)+vs(jat)))
                           dtmp = Pij*sum(sdqg(ixyz,11:13,jj,ii)*vd(1:3,iat) &
                              & +sdqg(ixyz, 2:4, jj,ii)*vd(1:3,jat) )
                           qtmp = Pij*sum( sdqg(ixyz,14:19,jj,ii)*vq(1:6,iat) &
                              & +sdqg(ixyz, 5:10,jj,ii)*vq(1:6,jat) )
                           g_xyz(ixyz) = g_xyz(ixyz)+stmp+dtmp+qtmp

                        enddo ! ixyz

                        ! Hamiltonian without Hav
                        HPij = km * zetaij * shpoly * Pij * sdq(1,jj,ii) * evtoau
                        ! save dE/dCN for CNi
                        dhdcn(iat) = dhdcn(iat) + HPij*dSEdcn(ish, iat)
                        ! save dE/dCN for CNj
                        dhdcn(jat) = dhdcn(jat) + HPij*dSEdcn(jsh, jat)
                     enddo
                  enddo
                  g(:,iat) = g(:,iat)+g_xyz
                  g(:,jat) = g(:,jat)-g_xyz
                  sigma(:, :) = sigma + spread(g_xyz, 1, 3) * spread(rij, 2, 3)
               enddo ! lattice translations
            enddo ! jsh : loop over shells on jat
         enddo  ! ish : loop over shells on iat
      enddo ! jat
   enddo  ! iat

   ! diagonal contributions
   !$omp parallel do default(none) schedule(dynamic) reduction(+:dhdcn) &
   !$omp shared(nat, at, nshell, hData, saoshell, P, dSEdcn) &
   !$omp private(iat, izp, ish, ishtyp, iao, i, Pij)
   do iat = 1, nat
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)

            Pij = P(i,i)
            ! save dE/dCN for CNi
            dhdcn(iat) = dhdcn(iat) + Pij*dSEdcn(ish, iat)*evtoau
         end do
      end do
   end do

end subroutine build_dSDQH0


!> GFN1-fast 2.2.0: overlap-only electronic gradient kernel for native GFN1.
!!
!! GFN1 has no anisotropic dipole/quadrupole SCC term.  The generic xTB
!! electronic gradient path therefore spends most of its time evaluating 19
!! derivative channels although only the overlap derivative contributes.  This
!! routine preserves the algebra of build_dSDQH0_noreset for the GFN1 case
!! (vs=vd=vq=0) while evaluating only dS/dR. The optional packed H0 path
!! reconstructs exactly the screened off-atom H0/S ratio on demand. The
!! Cartesian contraction transforms weights once and contracts primitives
!! directly; both improvements retain independent runtime fallbacks. Optional
!! directPrefactor forms the shell H0 coefficient from parameters and removes
!! divisions by overlap, shell polynomial and average self energy.
subroutine build_dSH0_GFN1_noreset(nShell, hData, selfEnergy, dSEdcn, intcut, &
      & nat, nao, nbf, at, xyz, caoshell, saoshell, nprim, primcount, &
      & alp, cont, H0, S, p, Pew, ves, dhdcn, g, sigma, H0packed, cartContract,directPrefactor)
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   real(wp), intent(in) :: selfEnergy(:, :)
   real(wp), intent(in) :: dSEdcn(:, :)
   integer, intent(in) :: nat, nao, nbf
   integer, intent(in) :: at(nat)
   real(wp), intent(in) :: intcut
   real(wp), intent(in) :: ves(:, :)
   real(wp), intent(in) :: xyz(3,nat)
   integer, intent(in) :: caoshell(:,:), saoshell(:,:), nprim(:), primcount(:)
   real(wp), intent(in) :: alp(:), cont(:)
   real(wp), intent(in), optional :: H0(:, :), H0packed(:)
   logical, intent(in), optional :: cartContract,directPrefactor
   real(wp), intent(in) :: S(:, :), p(:, :), Pew(:, :)
   real(wp), intent(inout) :: g(:, :), sigma(:, :), dhdcn(:)

   integer :: iat,jat,izp,jzp,ish,jsh,ishtyp,jshtyp,icao,jcao,naoi,naoj
   integer :: ii,jj,iao,jao,ixyz,il,jl,i
   real(wp) :: ri(3),rj(3),rij(3),rij2,point(3)
   real(wp) :: shpoly,dshpoly(3),hii,hjj,hav
   real(wp) :: sraw(6,6),sgraw(3,6,6),tmp(6,6)
   real(wp) :: Pij,Hij,HPij,g_xyz(3),dCN,weights(6,6),overlapGrad(3)
   real(wp) :: zi,zj,zeta,km,shellScale,directH0,polyForce(3)
   logical :: contractionActive,directActive
   real(wp), parameter :: rthr = 1600.0_wp

   ! nao/nbf are part of the stable interface and document the matrix sizes.
   ! They are intentionally not otherwise needed by the shell-pair kernel.
   point = 0.0_wp
   contractionActive=.false.
   if(present(cartContract)) contractionActive=cartContract
   directActive=.false.
   if(present(directPrefactor)) directActive=directPrefactor

   !$omp parallel do default(none) schedule(dynamic) &
   !$omp shared(nat,at,xyz,nShell,hData,selfEnergy,dSEdcn,P,Pew,H0,H0packed,S,ves,contractionActive,directActive, &
   !$omp& intcut,nprim,primcount,caoshell,saoshell,alp,cont,point) &
   !$omp private(iat,jat,izp,jzp,ri,rj,rij,rij2,ish,jsh,ishtyp,jshtyp, &
   !$omp& icao,jcao,naoi,naoj,il,jl,hii,hjj,hav,shpoly,dshpoly,sraw,sgraw, &
   !$omp& tmp,ixyz,ii,jj,iao,jao,Pij,Hij,HPij,g_xyz,dCN,weights,overlapGrad,zi,zj,zeta,km,shellScale,directH0,polyForce) &
   !$omp reduction(+:g,sigma,dhdcn)
   do iat = 1, nat
      ri = xyz(:,iat)
      izp = at(iat)
      do jat = 1, iat-1
         jzp = at(jat)
         rj = xyz(:,jat)
         rij = ri-rj
         rij2 = sum(rij**2)
         if (rij2 > rthr) cycle

         do ish = 1, nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            do jsh = 1, nShell(jzp)
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               il = ishtyp+1
               jl = jshtyp+1
               hii = selfEnergy(ish,iat)
               hjj = selfEnergy(jsh,jat)
               hav = 0.5_wp*(hii+hjj)

               call dshellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp), &
                  & hData%atomicRad(izp),hData%atomicRad(jzp),rij2,ri,rj, &
                  & shpoly,dshpoly)
               if(directActive) then
                  zi=hData%slaterExponent(ish,izp); zj=hData%slaterExponent(jsh,jzp)
                  zeta=(2.0_wp*sqrt(zi*zj)/(zi+zj))**hData%wExp
                  call h0scal(hData,il,jl,izp,jzp,hData%valenceShell(ish,izp)/=0, &
                     & hData%valenceShell(jsh,jzp)/=0,km)
                  shellScale=km*zeta*evtoau
                  directH0=hav*shellScale
               endif

               sraw = 0.0_wp
               sgraw = 0.0_wp
               if(.not.contractionActive) then
                  call get_grad_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
                     & point,intcut,nprim,primcount,alp,cont,sraw,sgraw)
                  do ixyz = 1, 3
                     tmp = sgraw(ixyz,:,:)
                     call dtrf2(tmp,ishtyp,jshtyp)
                     sgraw(ixyz,:,:) = tmp
                  end do
               endif

               g_xyz = 0.0_wp
               dCN = 0.0_wp
               weights=0.0_wp
               do ii = 1, llao2(ishtyp)
                  iao = ii+saoshell(ish,iat)
                  do jj = 1, llao2(jshtyp)
                     jao = jj+saoshell(jsh,jat)
                     Pij = p(jao,iao)
                     if(directActive) then
                        Hij=directH0*shpoly
                        if(present(H0packed)) then
                           ! SCF clips small nonzero overlaps outside matlist2.
                           ! Their nonzero raw H0 distinguishes screening from
                           ! a true overlap node, where the derivative survives.
                           if(S(jao,iao)==0.0_wp) then
                              if(H0packed(lin(iao,jao))/=0.0_wp) Hij=0.0_wp
                           endif
                        endif
                     else if(present(H0packed)) then
                        Hij=0.0_wp
                        ! S was set to exact zero outside matlist2. Keep
                        ! the derivative/Pew terms there, but never divide
                        ! by a screened overlap or use its packed H0 value.
                        if(S(jao,iao)/=0.0_wp) Hij=H0packed(lin(iao,jao))*evtoau/S(jao,iao)
                     else
                        Hij = H0(jao,iao)
                     endif
                     HPij = Hij*Pij
                     if(directActive) then
                        polyForce=2.0_wp*Pij*directH0*S(jao,iao)*dshpoly
                        dCN=dCN+Pij*shellScale*shpoly*S(jao,iao)
                     else
                        polyForce=2.0_wp*HPij*S(jao,iao)*dshpoly/shpoly
                        dCN=dCN+HPij/hav*S(jao,iao)
                     endif

                     if(contractionActive) then
                        weights(jj,ii)=2.0_wp*HPij - 2.0_wp*Pew(jao,iao) &
                           & - Pij*(ves(ish,iat)+ves(jsh,jat))
                        g_xyz=g_xyz+polyForce
                     else
                        g_xyz = g_xyz + polyForce &
                           & + sgraw(:,jj,ii)*(2.0_wp*HPij - 2.0_wp*Pew(jao,iao) &
                           & - Pij*(ves(ish,iat)+ves(jsh,jat)))
                     endif

                  end do
               end do

               if(contractionActive) then
                  call dtrf2_weights(weights,ishtyp,jshtyp)
                  call get_grad_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
                     & point,intcut,nprim,primcount,alp,cont,sraw,sgraw,weights,overlapGrad)
                  g_xyz=g_xyz+overlapGrad
               endif

               dhdcn(iat) = dhdcn(iat) + dCN*dSEdcn(ish,iat)
               dhdcn(jat) = dhdcn(jat) + dCN*dSEdcn(jsh,jat)
               g(:,iat) = g(:,iat)+g_xyz
               g(:,jat) = g(:,jat)-g_xyz
               sigma(:,:) = sigma + spread(g_xyz,1,3)*spread(rij,2,3)
            end do
         end do
      end do
   end do

   !$omp parallel do default(none) schedule(dynamic) reduction(+:dhdcn) &
   !$omp shared(nat,at,nShell,hData,saoshell,P,dSEdcn) &
   !$omp private(iat,izp,ish,ishtyp,iao,i,Pij)
   do iat = 1, nat
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)
            Pij = P(i,i)
            dhdcn(iat) = dhdcn(iat) + Pij*dSEdcn(ish,iat)*evtoau
         end do
      end do
   end do
end subroutine build_dSH0_GFN1_noreset


!> Computes the gradient of the dipole/qpole integral contribution
subroutine build_dSDQH0_noreset(nShell, hData, selfEnergy, dSEdcn, intcut, &
      & nat, nao, nbf, at, xyz, caoshell, saoshell, nprim, primcount, &
      & alp, cont, H0, S, p, Pew, ves, vs, vd, vq, dhdcn, g, sigma)
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   real(wp), intent(in) :: selfEnergy(:, :)
   real(wp), intent(in) :: dSEdcn(:, :)
   !> # of atoms
   integer, intent(in)    :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)    :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)    :: nbf
   !> Atomic numbers of atoms
   integer, intent(in)    :: at(nat)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)    :: intcut
   real(wp),intent(in)    :: ves(:, :)
   real(wp),intent(in)    :: vs(nat)
   real(wp),intent(in)    :: vd(3,nat)
   real(wp),intent(in)    :: vq(6,nat)
   !> Cartesian coordinates
   real(wp),intent(in)    :: xyz(3,nat)
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)    :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)    :: saoshell(:,:)
   integer, intent(in)    :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO, dimension: nbf
   integer, intent(in)    :: primcount(:)
   real(wp),intent(in)    :: alp(:)
   real(wp),intent(in)    :: cont(:)
   !> Molecular Hamiltonian
   real(wp),intent(in) :: H0(:, :)
   !> Molecular Overlap
   real(wp),intent(in) :: S(:, :)
   !> Density matrix
   real(wp),intent(in) :: p(:, :)
   !> Energy weighted density matrix
   real(wp),intent(in) :: Pew(:, :)
   real(wp),intent(inout) :: g(:, :)
   real(wp),intent(inout) :: sigma(:, :)
   real(wp),intent(inout) :: dhdcn(:)

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cj,alpi,rij2,ab,est
   real(wp) f1,f2,point(3),tmp(6,6),rij(3),ri(3),rj(3)
   real(wp) ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,jshmax
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,ixyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3,10),dum(10),sdq(10,6,6),sdqg(3,19,6,6)
   integer :: il, jl, itr
   real(wp) :: zi, zj, zetaij, km, hii, hjj, hav, shpoly, dshpoly(3), dCN
   real(wp) :: Pij, Hij, HPij, g_xyz(3)
   real(wp), parameter :: rthr = 1600.0_wp

   thr2 = intcut
   point = 0.0_wp
   ! call timing(t1,t3)
   !$omp parallel do default(none) schedule(dynamic) &
   !$omp shared(nat, at, xyz, nShell, hData, selfEnergy, dSEdcn, P, Pew, &
   !$omp& H0, S, ves, vs, vd, vq, intcut, nprim, primcount, caoshell, saoshell, &
   !$omp& alp, cont) &
   !$omp private(iat,jat,ixyz,izp,ci,rij2,jzp,ish,ishtyp,ij, &
   !$omp& icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp,dCN, &
   !$omp& sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,ri,rj,rij,km,shpoly,dshpoly, &
   !$omp& mli,mlj,dum,dumdum,tmp,dtmp,qtmp,il,jl,zi,zj,zetaij,hii,hjj,hav, &
   !$omp& iao,jao,ii,jj,k,pij,hij,hpij,g_xyz,itr) &
   !$omp reduction(+:g,sigma,dhdcn)
   do iat = 1,nat
      ri = xyz(:,iat)
      izp = at(iat)
      do jat = 1,iat-1
         !           if (jat.eq.iat) cycle
         jzp = at(jat)

         rj = xyz(:,jat)
         rij = ri - rj
         rij2 =  sum( rij**2 )

         if (rij2 > rthr) cycle
         do ish = 1,nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(jzp)
            !              if(iat.eq.jat) jshmax = ish
            do jsh = 1,jshmax ! jshells
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! distance dependent polynomial
               call dshellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp),&
                  & hData%atomicRad(izp),hData%atomicRad(jzp),rij2,ri,rj,&
                  & shpoly,dshpoly)

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * (hii + hjj)

               sdqg = 0;sdq = 0
               call get_grad_multiint(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
                  &                   intcut,nprim,primcount,alp,cont,sdq,sdqg)
               do k = 1,19 ! 1 S, 2-4 D, 5-10 Q, 11-13 D, 14-19 Q
                  do ixyz = 1,3
                     ! transform from CAO to SAO
                     !                        call dtrf2(sdqg(ixyz,k,1:6,1:6),ishtyp,jshtyp)
                     tmp(1:6,1:6) = sdqg(ixyz,k,1:6,1:6)
                     call dtrf2(tmp,ishtyp,jshtyp)
                     sdqg(ixyz,k,1:6,1:6) = tmp(1:6,1:6)
                  enddo
               enddo
               g_xyz(:) = 0.0_wp
               dCN = 0.0_wp
               do ii = 1,llao2(ishtyp)
                  iao = ii+saoshell(ish,iat)
                  do jj = 1,llao2(jshtyp)
                     jao = jj+saoshell(jsh,jat)
                     Pij = p(jao,iao)

                     ! Hamiltonian element without overlap
                     Hij  = H0(jao, iao)
                     HPij = Hij * Pij

                     g_xyz(:) = g_xyz + 2*HPij*S(jao,iao)*dshpoly/shpoly &
                        & + sdqg(:,1,jj,ii)*(2*HPij - 2*Pew(jao, iao) &
                        & - Pij*(ves(ish,iat)+ves(jsh,jat)) &
                        & + Pij*(vs(iat)+vs(jat)))

                     do ixyz = 1,3
                        dtmp = Pij*sum(sdqg(ixyz,11:13,jj,ii)*vd(1:3,iat) &
                           & +sdqg(ixyz, 2:4, jj,ii)*vd(1:3,jat) )
                        qtmp = Pij*sum( sdqg(ixyz,14:19,jj,ii)*vq(1:6,iat) &
                           & +sdqg(ixyz, 5:10,jj,ii)*vq(1:6,jat) )
                        g_xyz(ixyz) = g_xyz(ixyz)+dtmp+qtmp

                     enddo ! ixyz

                     ! Hamiltonian without Hav
                     dCN = dCN + HPij / hav * S(jao, iao)
                  enddo
               enddo
               ! save dE/dCN for CNi
               dhdcn(iat) = dhdcn(iat) + dCN*dSEdcn(ish, iat)
               ! save dE/dCN for CNj
               dhdcn(jat) = dhdcn(jat) + dCN*dSEdcn(jsh, jat)
               g(:,iat) = g(:,iat)+g_xyz
               g(:,jat) = g(:,jat)-g_xyz
               sigma(:, :) = sigma + spread(g_xyz, 1, 3) * spread(rij, 2, 3)
            enddo ! jsh : loop over shells on jat
         enddo  ! ish : loop over shells on iat
      enddo ! jat
   enddo  ! iat

   ! diagonal contributions
   !$omp parallel do default(none) schedule(dynamic) reduction(+:dhdcn) &
   !$omp shared(nat, at, nshell, hData, saoshell, P, dSEdcn) &
   !$omp private(iat, izp, ish, ishtyp, iao, i, Pij)
   do iat = 1, nat
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)

            Pij = P(i,i)
            ! save dE/dCN for CNi
            dhdcn(iat) = dhdcn(iat) + Pij*dSEdcn(ish, iat)*evtoau
         end do
      end do
   end do

end subroutine build_dSDQH0_noreset


!> Count number of significant dipole integrals
subroutine count_dpint(ndp, dpint, thr)

   !> Number of significant dipole integrals
   integer, intent(out) :: ndp

   !> Dipole integrals
   real(wp), intent(in) :: dpint(:, :, :)

   !> Neglect threshold for dipole integrals
   real(wp), intent(in) :: thr

   integer :: i, j
   real(wp) :: tmp1, thr2

   ndp = 0
   thr2 = (thr*1.0e-2_wp)-thr*1.0e-12_wp

   do i = 1, size(dpint, dim=3)
      do j = 1, i
         tmp1 = sum(dpint(1:3, j, i)**2)
         if (tmp1 > thr2) ndp = ndp + 1
      enddo
   enddo

end subroutine count_dpint


!> Count number of significant quadrupole integrals
subroutine count_qpint(nqp, qpint, thr)

   !> Number of significant quadrupole integrals
   integer, intent(out) :: nqp

   !> Quadrupole integrals
   real(wp), intent(in) :: qpint(:, :, :)

   !> Neglect threshold for quadrupole integrals
   real(wp), intent(in) :: thr

   integer :: i, j
   real(wp) :: tmp2, thr2

   nqp = 0
   thr2 = (thr*1.0e-2_wp)-thr*1.0e-12_wp

   do i = 1, size(qpint, dim=3)
      do j = 1, i
         tmp2 = sum(qpint(1:3, j, i)**2) + 2*sum(qpint(4:6, j, i)**2)
         if (tmp2 > thr2) nqp = nqp + 1
      enddo
   enddo

end subroutine count_qpint


end module xtb_xtb_hamiltonian
