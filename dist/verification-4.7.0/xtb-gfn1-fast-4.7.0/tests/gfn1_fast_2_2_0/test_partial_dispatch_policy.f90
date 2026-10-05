program test_gfn1_fast_220_partial_dispatch_policy
   use iso_fortran_env, only : real64
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy, gfn1SpectralPartialUnsafe, &
      & gfn1PartialSpectrumReady
   implicit none
   type(TGFN1FastPolicy) :: policy
   real(real64), parameter :: zero = 0.0_real64

   ! 2.2.0: no speculative partial solve before one exact full spectrum.
   call assert_false(gfn1PartialSpectrumReady(policy,0,.false.,.false.), 1)
   call assert_true (gfn1PartialSpectrumReady(policy,1,.false.,.false.), 2)
   call assert_false(gfn1PartialSpectrumReady(policy,1,.true., .false.), 3)
   call assert_false(gfn1PartialSpectrumReady(policy,1,.false.,.true. ), 4)

   policy%partialWarmupFullSolves = 2
   call assert_false(gfn1PartialSpectrumReady(policy,1,.false.,.false.), 5)
   call assert_true (gfn1PartialSpectrumReady(policy,2,.false.,.false.), 6)

   ! Retain the 2.1.5/2.1.6 spectral safety semantics.
   policy = TGFN1FastPolicy()
   call assert_true (gfn1SpectralPartialUnsafe(policy,0.05_real64,zero,zero), 7)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.10_real64,zero,zero), 8)
   call assert_true (gfn1SpectralPartialUnsafe(policy,0.15_real64,300.0_real64,zero), 9)
   call assert_false(gfn1SpectralPartialUnsafe(policy,0.25_real64,300.0_real64,zero), 10)
   call assert_true (gfn1SpectralPartialUnsafe(policy,0.25_real64,300.0_real64,1.0e-5_real64), 11)

   print *, 'PASS partial dispatcher policy 2.2.0'
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
end program test_gfn1_fast_220_partial_dispatch_policy
