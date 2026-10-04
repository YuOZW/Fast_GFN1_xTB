! This file is part of xtb.
!
! GFN1-fast adaptive execution policy, introduced in 2.0.2.
module xtb_gfn1_fast_policy
   use iso_fortran_env, only : int64, real64
   implicit none
   private

   integer, parameter, public :: gfn1_full_eig_backend_auto = 0
   integer, parameter, public :: gfn1_full_eig_backend_syevd = 1
   integer, parameter, public :: gfn1_full_eig_backend_syevr = 2
   integer, parameter, public :: gfn1_fast_default_full_eig_autotune_min_nao = 768
   integer, parameter, public :: gfn1_fast_default_compact_density_min_nao = 128
   integer, parameter, public :: gfn1_fast_default_partial_min_nao = 128
   integer, parameter, public :: gfn1_fast_default_partial_warmup_full_solves = 1
   integer, parameter, public :: gfn1_fast_default_partial_probe_min_nao = 768
   integer, parameter, public :: gfn1_fast_default_thermal_warmup_full_solves = 4
   integer, parameter, public :: gfn1_fast_default_thermal_guard_roots = 8
   integer, parameter, public :: gfn1_fast_default_davidson_min_nao = 1024
   integer, parameter, public :: gfn1_fast_default_davidson_max_root_percent = 45
   integer, parameter, public :: gfn1_fast_default_sparse_max_pair_percent = 10
   integer, parameter, public :: gfn1_fast_default_hess_static_max_nao = 128
   integer, parameter, public :: gfn1_fast_default_h1_omp_min_nmat = 262144
   integer, parameter, public :: gfn1_fast_default_solvent_cache_min_nao = 128
   integer, parameter, public :: gfn1_fast_default_hess_parallel_min_nao = 24
   integer, parameter, public :: gfn1_fast_default_hess_min_jobs_per_thread = 2
   real(real64), parameter, public :: gfn1_fast_default_partial_min_gap_ev = 0.10_real64
   real(real64), parameter, public :: gfn1_fast_default_partial_min_gap_kbt = 8.0_real64
   real(real64), parameter, public :: gfn1_fast_default_partial_max_fractionality = 1.0e-6_real64
   real(real64), parameter, public :: gfn1_fast_default_thermal_tail_kbt = 50.0_real64
   real(real64), parameter, public :: gfn1_fast_default_partial_autotune_max_ratio = 0.95_real64
   real(real64), parameter, public :: gfn1_fast_default_full_eig_autotune_max_ratio = 0.98_real64

   type, public :: TGFN1FastPolicy
      integer :: fullEigAutotuneMinNao = gfn1_fast_default_full_eig_autotune_min_nao
      integer :: compactDensityMinNao = gfn1_fast_default_compact_density_min_nao
      integer :: partialMinNao = gfn1_fast_default_partial_min_nao
      integer :: partialWarmupFullSolves = gfn1_fast_default_partial_warmup_full_solves
      integer :: partialProbeMinNao = gfn1_fast_default_partial_probe_min_nao
      integer :: thermalWarmupFullSolves = gfn1_fast_default_thermal_warmup_full_solves
      integer :: thermalGuardRoots = gfn1_fast_default_thermal_guard_roots
      integer :: davidsonMinNao = gfn1_fast_default_davidson_min_nao
      integer :: davidsonMaxRootPercent = gfn1_fast_default_davidson_max_root_percent
      integer :: sparseMaxPairPercent = gfn1_fast_default_sparse_max_pair_percent
      integer :: hessianStaticMaxNao = gfn1_fast_default_hess_static_max_nao
      integer :: h1OMPMinNmat = gfn1_fast_default_h1_omp_min_nmat
      integer :: solventCacheMinNao = gfn1_fast_default_solvent_cache_min_nao
      integer :: hessianParallelMinNao = gfn1_fast_default_hess_parallel_min_nao
      integer :: hessianMinJobsPerThread = gfn1_fast_default_hess_min_jobs_per_thread
      real(real64) :: partialMinGapEV = gfn1_fast_default_partial_min_gap_ev
      real(real64) :: partialMinGapKBT = gfn1_fast_default_partial_min_gap_kbt
      real(real64) :: partialMaxFractionality = gfn1_fast_default_partial_max_fractionality
      real(real64) :: thermalTailKBT = gfn1_fast_default_thermal_tail_kbt
      real(real64) :: partialAutotuneMaxRatio = gfn1_fast_default_partial_autotune_max_ratio
      real(real64) :: fullEigAutotuneMaxRatio = gfn1_fast_default_full_eig_autotune_max_ratio
      integer :: fullEigBackend = gfn1_full_eig_backend_syevd
      logical :: profile = .false.
      logical :: davidsonAudit = .false.
      ! GFN1-fast 2.2.0: use the GFN1-specific overlap-only electronic
      ! gradient kernel by default. Set XTB_GFN1_FAST_DISABLE_GRADIENT_KERNEL
      ! to a true value to recover the generic multipole-derivative path.
      logical :: gradientKernel = .true.
      ! GFN1-fast 2.2.2: finite-temperature partial diagonalization is allowed
      ! only when the computed root window extends beyond the exact fermismear
      ! zero-occupation cutoff. This replaces the overly conservative 2.1.5
      ! rule that disabled partial solves for any fractional occupation.
      logical :: thermalPartial = .true.
      ! GFN1-fast 2.2.4: retain indexed partial diagonalization only when a
      ! measured accepted partial path is materially faster than the most
      ! recent full factorized solve at the same matrix size.
      logical :: partialAutotune = .true.
      ! GFN1-fast 2.2.6: choose the fastest exact full-spectrum standard
      ! eigensolver backend after the cached-overlap SYGST reduction. AUTO
      ! audits SYEVD versus full-spectrum SYEVR once on the same matrix.
      logical :: fullEigAutotune = .true.
      ! GFN1-fast 2.2.6: on intermediate native-GFN1 SCC steps,
      ! build only the upper triangle of P with a weighted SYRK and fuse the
      ! Mulliken/H0 contractions.  The final verification step uses the legacy
      ! full-density path, preserving the converged Energy/Gradient reference.
      logical :: compactDensitySCC = .true.
   end type TGFN1FastPolicy

   public :: getGFN1FastPolicy
   public :: useGFN1PartialEigensolver
   public :: gfn1PartialSpectrumReady
   public :: gfn1PartialSpectrumReadyThermal
   public :: gfn1ThermalRootGuard
   public :: gfn1PartialAutotuneKeep
   public :: useGFN1Davidson
   public :: useGFN1SparseDavidson
   public :: useGFN1StaticHessianSchedule
   public :: useGFN1ParallelHessian
   public :: gfn1SpectralPartialUnsafe

   type(TGFN1FastPolicy), save :: cachedPolicy
   logical, save :: policyInitialized = .false.

contains

subroutine getGFN1FastPolicy(policy)
   type(TGFN1FastPolicy), intent(out) :: policy

   ! Read environment overrides once per process.  In Hessian calculations
   ! this avoids repeating getenv/string parsing for every displaced SCC.
   if (.not.policyInitialized) then
      !$omp critical(xtb_gfn1_fast_policy_init)
      if (.not.policyInitialized) then
         cachedPolicy = TGFN1FastPolicy()
         call readPositiveIntegerEnv('XTB_GFN1_FAST_FULL_EIG_AUTOTUNE_MIN_NAO', cachedPolicy%fullEigAutotuneMinNao)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_COMPACT_DENSITY_MIN_NAO', cachedPolicy%compactDensityMinNao)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_PARTIAL_MIN_NAO', cachedPolicy%partialMinNao)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_PARTIAL_WARMUP_FULL_SOLVES', cachedPolicy%partialWarmupFullSolves)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_PARTIAL_PROBE_MIN_NAO', cachedPolicy%partialProbeMinNao)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_THERMAL_WARMUP_FULL_SOLVES', cachedPolicy%thermalWarmupFullSolves)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_THERMAL_GUARD_ROOTS', cachedPolicy%thermalGuardRoots)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_DAVIDSON_MIN_NAO', cachedPolicy%davidsonMinNao)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_DAVIDSON_MAX_ROOT_PERCENT', cachedPolicy%davidsonMaxRootPercent)
         cachedPolicy%davidsonMaxRootPercent = min(100,max(1,cachedPolicy%davidsonMaxRootPercent))
         call readPositiveIntegerEnv('XTB_GFN1_FAST_SPARSE_MAX_PAIR_PERCENT', cachedPolicy%sparseMaxPairPercent)
         cachedPolicy%sparseMaxPairPercent = min(100,max(1,cachedPolicy%sparseMaxPairPercent))
         call readPositiveIntegerEnv('XTB_GFN1_FAST_HESS_STATIC_MAX_NAO', cachedPolicy%hessianStaticMaxNao)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_H1_OMP_MIN_NMAT', cachedPolicy%h1OMPMinNmat)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_SOLVENT_MIN_NAO', cachedPolicy%solventCacheMinNao)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_HESS_PARALLEL_MIN_NAO', cachedPolicy%hessianParallelMinNao)
         call readPositiveIntegerEnv('XTB_GFN1_FAST_HESS_MIN_JOBS_PER_THREAD', cachedPolicy%hessianMinJobsPerThread)
         call readPositiveRealEnv('XTB_GFN1_FAST_PARTIAL_MIN_GAP_EV', cachedPolicy%partialMinGapEV)
         call readPositiveRealEnv('XTB_GFN1_FAST_PARTIAL_MIN_GAP_KBT', cachedPolicy%partialMinGapKBT)
         call readPositiveRealEnv('XTB_GFN1_FAST_PARTIAL_MAX_FRACTIONALITY', cachedPolicy%partialMaxFractionality)
         call readPositiveRealEnv('XTB_GFN1_FAST_THERMAL_TAIL_KBT', cachedPolicy%thermalTailKBT)
         call readPositiveRealEnv('XTB_GFN1_FAST_PARTIAL_AUTOTUNE_MAX_RATIO', cachedPolicy%partialAutotuneMaxRatio)
         call readPositiveRealEnv('XTB_GFN1_FAST_FULL_EIG_AUTOTUNE_MAX_RATIO', cachedPolicy%fullEigAutotuneMaxRatio)
         cachedPolicy%partialAutotuneMaxRatio = min(2.0_real64,max(0.10_real64,cachedPolicy%partialAutotuneMaxRatio))
         cachedPolicy%fullEigAutotuneMaxRatio = min(2.0_real64,max(0.10_real64,cachedPolicy%fullEigAutotuneMaxRatio))
         call readFullEigBackendEnv('XTB_GFN1_FAST_FULL_EIG_BACKEND', cachedPolicy%fullEigBackend)
         ! The original fermismear routine sets occupations to exact zero only
         ! at 50 kBT. Never permit an override below that reference cutoff.
         cachedPolicy%thermalTailKBT = max(gfn1_fast_default_thermal_tail_kbt, cachedPolicy%thermalTailKBT)
         cachedPolicy%profile = readLogicalEnv('XTB_GFN1_FAST_PROFILE')
         cachedPolicy%davidsonAudit = readLogicalEnv('XTB_GFN1_FAST_DAVIDSON_AUDIT')
         cachedPolicy%gradientKernel = .not.readLogicalEnv('XTB_GFN1_FAST_DISABLE_GRADIENT_KERNEL')
         cachedPolicy%thermalPartial = .not.readLogicalEnv('XTB_GFN1_FAST_DISABLE_THERMAL_PARTIAL')
         cachedPolicy%partialAutotune = .not.readLogicalEnv('XTB_GFN1_FAST_DISABLE_PARTIAL_AUTOTUNE')
         cachedPolicy%fullEigAutotune = .not.readLogicalEnv('XTB_GFN1_FAST_DISABLE_FULL_EIG_AUTOTUNE')
         cachedPolicy%compactDensitySCC = .not.readLogicalEnv('XTB_GFN1_FAST_DISABLE_COMPACT_DENSITY')
         policyInitialized = .true.
      end if
      !$omp end critical(xtb_gfn1_fast_policy_init)
   end if
   policy = cachedPolicy
end subroutine getGFN1FastPolicy

pure logical function useGFN1PartialEigensolver(policy, nao, nroots) result(usePartial)
   type(TGFN1FastPolicy), intent(in) :: policy
   integer, intent(in) :: nao, nroots

   ! The partial eigensolver has fixed setup/validation costs.  Below the
   ! threshold the factorized full eigensolver is faster for the small
   ! matrices seen in the 2.0.0 benchmark.  For larger matrices retain the
   ! established 0.7-0.9 guard that rejects partial solves spanning more than
   ! 75 percent of the spectrum.
   usePartial = nao >= policy%partialMinNao .and. nroots > 0 .and. &
      & nroots < nao .and. 4*nroots <= 3*nao
   ! GFN1-fast 2.2.4: on the 594-AO benchmark an accepted SYEVR partial
   ! probe cost about twice the cached-factor full solve. Avoid paying this
   ! known fixed-cost probe in the medium-size regime. Larger matrices still
   ! use the 2.2.3 runtime autotuner, so the crossover remains hardware-aware.
   if (policy%partialAutotune .and. nao < policy%partialProbeMinNao) usePartial = .false.
end function useGFN1PartialEigensolver


pure logical function gfn1PartialSpectrumReady(policy, certifiedFullSolves, &
      & spectralDisabled, pathDisabled) result(ready)
   type(TGFN1FastPolicy), intent(in) :: policy
   integer, intent(in) :: certifiedFullSolves
   logical, intent(in) :: spectralDisabled, pathDisabled

   ! GFN1-fast 2.1.7: partial diagonalization is speculative until at least
   ! one exact full solve has certified the current SCC spectrum.  This avoids
   ! paying partial validation followed by an immediate exact fallback during
   ! the unstable first iterations.  A non-spectral partial failure also
   ! disables the path for the rest of this SCC.
   ready = certifiedFullSolves >= policy%partialWarmupFullSolves .and. &
      & .not.spectralDisabled .and. .not.pathDisabled
end function gfn1PartialSpectrumReady



pure logical function gfn1PartialSpectrumReadyThermal(policy, certifiedFullSolves, &
      & electronicTemperature, spectralDisabled, pathDisabled) result(ready)
   type(TGFN1FastPolicy), intent(in) :: policy
   integer, intent(in) :: certifiedFullSolves
   real(real64), intent(in) :: electronicTemperature
   logical, intent(in) :: spectralDisabled, pathDisabled
   integer :: requiredFullSolves

   requiredFullSolves = policy%partialWarmupFullSolves
   if (policy%thermalPartial .and. electronicTemperature > 0.1_real64) then
      requiredFullSolves = max(requiredFullSolves, policy%thermalWarmupFullSolves)
   end if
   ready = certifiedFullSolves >= requiredFullSolves .and. &
      & .not.spectralDisabled .and. .not.pathDisabled
end function gfn1PartialSpectrumReadyThermal

pure integer function gfn1ThermalRootGuard(policy, electronicTemperature, currentGuard) result(guard)
   type(TGFN1FastPolicy), intent(in) :: policy
   real(real64), intent(in) :: electronicTemperature
   integer, intent(in), optional :: currentGuard

   guard = 0
   if (policy%thermalPartial .and. electronicTemperature > 0.1_real64) then
      guard = policy%thermalGuardRoots
      if (present(currentGuard)) guard = max(guard,currentGuard)
   end if
end function gfn1ThermalRootGuard


pure logical function gfn1PartialAutotuneKeep(policy, partialTime, fullTime) result(keepPartial)
   type(TGFN1FastPolicy), intent(in) :: policy
   real(real64), intent(in) :: partialTime, fullTime

   ! If no trustworthy full timing exists yet, do not disable a numerically
   ! valid partial solve. Otherwise require a configurable speed margin.
   keepPartial = .true.
   if (.not.policy%partialAutotune) return
   if (partialTime <= 0.0_real64 .or. fullTime <= 0.0_real64) return
   keepPartial = partialTime <= policy%partialAutotuneMaxRatio*fullTime
end function gfn1PartialAutotuneKeep

pure logical function useGFN1Davidson(policy, nao, nroots) result(useDavidson)
   type(TGFN1FastPolicy), intent(in) :: policy
   integer, intent(in) :: nao, nroots

   ! Block-Davidson is beneficial only when the requested invariant subspace
   ! is a minority of the full AO spectrum.  A conservative default keeps the
   ! mature factorized LAPACK path for small/medium matrices and for systems
   ! whose occupied space is too large for an iterative lowest-roots method.
   useDavidson = nao >= policy%davidsonMinNao .and. nroots > 0 .and. &
      & nroots < nao .and. 100*nroots <= policy%davidsonMaxRootPercent*nao
end function useGFN1Davidson

pure logical function useGFN1SparseDavidson(policy, nao, nmat) result(useSparse)
   type(TGFN1FastPolicy), intent(in) :: policy
   integer, intent(in) :: nao, nmat
   integer(int64) :: lhs, rhs

   ! The pair-list kernel wins only while the screened upper triangle is
   ! sufficiently sparse.  Above this threshold optimized dense BLAS is
   ! faster, so skip Davidson and retain the factorized LAPACK path.
   if (nao <= 0 .or. nmat <= 0) then
      useSparse = .false.
      return
   end if
   lhs = 200_int64*int(nmat,int64)
   rhs = int(policy%sparseMaxPairPercent,int64)*int(nao,int64)*int(nao+1,int64)
   useSparse = lhs <= rhs
end function useGFN1SparseDavidson


pure logical function gfn1SpectralPartialUnsafe(policy, gapEV, electronicTemperature, &
      & maxFractionality) result(unsafe)
   type(TGFN1FastPolicy), intent(in) :: policy
   real(real64), intent(in) :: gapEV, electronicTemperature, maxFractionality
   real(real64) :: gapThreshold
   ! Keep the 2.1.5 dispatcher semantics in one testable location.  The
   ! constants reproduce xtb_mctc_constants:kB * xtb_mctc_convert:autoev.
   real(real64), parameter :: kb_au = 3.166808578545117e-6_real64
   real(real64), parameter :: hartree_to_ev = 27.211386245988_real64

   gapThreshold = policy%partialMinGapEV
   if (electronicTemperature > 0.1_real64) gapThreshold = max(gapThreshold, &
      & policy%partialMinGapKBT*kb_au*hartree_to_ev*electronicTemperature)

   unsafe = gapEV < gapThreshold
   ! 2.2.2: fractional occupation by itself is not unsafe.  A partial solve is
   ! accepted only after the explicit 50-kBT Fermi-tail coverage gate in SCC.
   ! Retain the 2.1.5 behavior when thermal-aware partial solving is disabled.
   if (electronicTemperature > 0.1_real64 .and. .not.policy%thermalPartial) unsafe = unsafe .or. &
      & maxFractionality > policy%partialMaxFractionality
end function gfn1SpectralPartialUnsafe

pure logical function useGFN1StaticHessianSchedule(policy, nao) result(useStatic)
   type(TGFN1FastPolicy), intent(in) :: policy
   integer, intent(in) :: nao

   ! Small systems have short, nearly uniform SCC displacement jobs; dynamic
   ! scheduling overhead is then measurable. Larger systems retain dynamic
   ! scheduling because SCC iteration counts can vary substantially.
   useStatic = nao > 0 .and. nao <= policy%hessianStaticMaxNao
end function useGFN1StaticHessianSchedule


pure logical function useGFN1ParallelHessian(policy, nao, njobs, nthreads) result(useParallel)
   type(TGFN1FastPolicy), intent(in) :: policy
   integer, intent(in) :: nao, njobs, nthreads

   ! Creating an OpenMP team is a measurable fixed cost for tiny Hessians.
   ! Keep water-sized and similarly tiny systems serial, while allowing the
   ! outer displacement parallelism as soon as there is enough work to give
   ! every worker at least a few independent Cartesian jobs.
   useParallel = nthreads > 1 .and. nao >= policy%hessianParallelMinNao .and. &
      & njobs >= policy%hessianMinJobsPerThread*nthreads
end function useGFN1ParallelHessian

subroutine readPositiveIntegerEnv(name, value)
   character(len=*), intent(in) :: name
   integer, intent(inout) :: value
   character(len=64) :: buffer
   integer :: stat, parsed, ios

   buffer = ''
   call get_environment_variable(name, buffer, status=stat)
   if (stat /= 0 .or. len_trim(buffer) == 0) return
   read(buffer, *, iostat=ios) parsed
   if (ios == 0 .and. parsed > 0) value = parsed
end subroutine readPositiveIntegerEnv

subroutine readPositiveRealEnv(name, value)
   character(len=*), intent(in) :: name
   real(real64), intent(inout) :: value
   character(len=64) :: buffer
   integer :: stat, ios
   real(real64) :: parsed

   buffer = ''
   call get_environment_variable(name, buffer, status=stat)
   if (stat /= 0 .or. len_trim(buffer) == 0) return
   read(buffer, *, iostat=ios) parsed
   if (ios == 0 .and. parsed >= 0.0_real64) value = parsed
end subroutine readPositiveRealEnv

subroutine readFullEigBackendEnv(name, value)
   character(len=*), intent(in) :: name
   integer, intent(inout) :: value
   character(len=64) :: buffer
   integer :: stat

   buffer = ''
   call get_environment_variable(name, buffer, status=stat)
   if (stat /= 0 .or. len_trim(buffer) == 0) return
   select case (lowerCase(trim(adjustl(buffer))))
   case ('auto')
      value = gfn1_full_eig_backend_auto
   case ('syevd','dsyevd')
      value = gfn1_full_eig_backend_syevd
   case ('syevr','dsyevr')
      value = gfn1_full_eig_backend_syevr
   case default
      value = gfn1_full_eig_backend_auto
   end select
end subroutine readFullEigBackendEnv

logical function readLogicalEnv(name) result(enabled)
   character(len=*), intent(in) :: name
   character(len=64) :: buffer
   integer :: stat

   enabled = .false.
   buffer = ''
   call get_environment_variable(name, buffer, status=stat)
   if (stat /= 0 .or. len_trim(buffer) == 0) return
   select case (lowerCase(trim(adjustl(buffer))))
   case ('1','true','yes','on','enable','enabled')
      enabled = .true.
   case default
      enabled = .false.
   end select
end function readLogicalEnv

pure function lowerCase(input) result(output)
   character(len=*), intent(in) :: input
   character(len=len(input)) :: output
   integer :: i, c

   output = input
   do i = 1, len(input)
      c = iachar(input(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) output(i:i) = achar(c + 32)
   end do
end function lowerCase

end module xtb_gfn1_fast_policy
