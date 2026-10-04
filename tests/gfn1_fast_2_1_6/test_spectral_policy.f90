program test_gfn1_fast_216_spectral_policy
   use iso_fortran_env, only : real64
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy, gfn1SpectralPartialUnsafe
   implicit none
   type(TGFN1FastPolicy) :: policy
   real(real64), parameter :: zero = 0.0_real64

   ! Default zero-temperature gap gate: unsafe strictly below 0.10 eV.
   call assert_true(gfn1SpectralPartialUnsafe(policy,0.05_real64,zero,zero), 1)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.10_real64,zero,zero), 2)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.20_real64,zero,zero), 3)

   ! At 300 K the 8 kBT criterion (~0.2068 eV) is stricter than 0.10 eV.
   call assert_true(gfn1SpectralPartialUnsafe(policy,0.15_real64,300.0_real64,zero), 4)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.25_real64,300.0_real64,zero), 5)

   ! Fractional occupations disable partial solvers only at finite T.
   call assert_true(gfn1SpectralPartialUnsafe(policy,0.25_real64,300.0_real64,1.0e-5_real64), 6)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.25_real64,300.0_real64,1.0e-6_real64), 7)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.25_real64,zero,1.0e-2_real64), 8)

   ! Policy overrides are honored by the shared decision function.
   policy%partialMinGapEV = 0.30_real64
   call assert_true(gfn1SpectralPartialUnsafe(policy,0.25_real64,zero,zero), 9)

   print *, 'PASS spectral policy 2.1.6'
contains
   subroutine assert_true(value, code)
      logical, intent(in) :: value
      integer, intent(in) :: code
      if (.not.value) then
         print *, "FAILED assertion", code
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_false(value, code)
      logical, intent(in) :: value
      integer, intent(in) :: code
      if (value) then
         print *, "FAILED assertion", code
         error stop 1
      end if
   end subroutine assert_false
end program test_gfn1_fast_216_spectral_policy
