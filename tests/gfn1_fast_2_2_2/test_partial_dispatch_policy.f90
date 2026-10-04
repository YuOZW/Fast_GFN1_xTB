program test_gfn1_fast_222_partial_dispatch_policy
   use iso_fortran_env, only : real64
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy, gfn1SpectralPartialUnsafe, &
      & gfn1PartialSpectrumReady, gfn1PartialSpectrumReadyThermal, gfn1ThermalRootGuard
   implicit none
   type(TGFN1FastPolicy) :: policy
   real(real64), parameter :: zero = 0.0_real64

   ! Legacy zero-temperature certification remains unchanged.
   call assert_false(gfn1PartialSpectrumReady(policy,0,.false.,.false.), 1)
   call assert_true (gfn1PartialSpectrumReady(policy,1,.false.,.false.), 2)

   ! 2.2.2: finite-T partial solving waits for a more stable full-spectrum history.
   call assert_false(gfn1PartialSpectrumReadyThermal(policy,1,300.0_real64,.false.,.false.), 3)
   call assert_false(gfn1PartialSpectrumReadyThermal(policy,3,300.0_real64,.false.,.false.), 4)
   call assert_true (gfn1PartialSpectrumReadyThermal(policy,4,300.0_real64,.false.,.false.), 5)
   call assert_false(gfn1PartialSpectrumReadyThermal(policy,4,300.0_real64,.true., .false.), 6)
   call assert_false(gfn1PartialSpectrumReadyThermal(policy,4,300.0_real64,.false.,.true. ), 7)
   call assert_true (gfn1PartialSpectrumReadyThermal(policy,1,0.0_real64,.false.,.false.), 8)

   ! Thermal partial path no longer rejects harmless fractional occupation.
   policy = TGFN1FastPolicy()
   call assert_true (gfn1SpectralPartialUnsafe(policy,0.05_real64,zero,zero), 9)
   call assert_true (gfn1SpectralPartialUnsafe(policy,0.15_real64,300.0_real64,zero), 10)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.25_real64,300.0_real64,zero), 11)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.25_real64,300.0_real64,1.0e-3_real64), 12)

   ! Disabling thermal-aware partial restores the conservative 2.1.5 rule.
   policy%thermalPartial = .false.
   call assert_true(gfn1SpectralPartialUnsafe(policy,0.25_real64,300.0_real64,1.0e-3_real64), 13)

   policy = TGFN1FastPolicy()
   call assert_equal(gfn1ThermalRootGuard(policy,300.0_real64),8,14)
   call assert_equal(gfn1ThermalRootGuard(policy,0.0_real64),0,15)
   call assert_equal(gfn1ThermalRootGuard(policy,300.0_real64,16),16,16)

   print *, 'PASS thermal partial dispatcher policy 2.2.2'
contains
   subroutine assert_true(value, code)
      logical, intent(in) :: value
      integer, intent(in) :: code
      if (.not.value) then
         print *, 'FAILED assertion', code
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_false(value, code)
      logical, intent(in) :: value
      integer, intent(in) :: code
      if (value) then
         print *, 'FAILED assertion', code
         error stop 1
      end if
   end subroutine assert_false

   subroutine assert_equal(value, expected, code)
      integer, intent(in) :: value, expected, code
      if (value /= expected) then
         print *, 'FAILED assertion', code, value, expected
         error stop 1
      end if
   end subroutine assert_equal
end program test_gfn1_fast_222_partial_dispatch_policy
