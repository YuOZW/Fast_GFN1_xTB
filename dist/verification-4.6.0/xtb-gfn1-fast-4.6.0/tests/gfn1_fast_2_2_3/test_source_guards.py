from pathlib import Path

root = Path(__file__).resolve().parents[2]
scc = (root / "src/scc_core.f90").read_text()
policy = (root / "src/gfn1_fast_policy.f90").read_text()
main = (root / "src/prog/main.F90").read_text()
scf = (root / "src/scf_module.F90").read_text()
ham = (root / "src/xtb/hamiltonian.f90").read_text()

# Keep the 2.1.6 noreset gradient correctness fix.
start = ham.index("subroutine build_dSDQH0_noreset")
end = ham.index("end subroutine build_dSDQH0_noreset", start)
noreset = ham[start:end]
assert "stmp" not in noreset
assert "g_xyz(ixyz) = g_xyz(ixyz)+dtmp+qtmp" in noreset

# Strict-Debug findings retained in the 2.2.3 source baseline.
assert "er = 0.0_wp" in main
assert "exist = .false." in main
assert "sigma = 0.0_wp" in scf

# 2.1.7 certified-spectrum dispatcher retained in 2.2.3.
assert "gfn1PartialSpectrumReady" in policy
assert "partialWarmupFullSolves" in policy
assert "XTB_GFN1_FAST_PARTIAL_WARMUP_FULL_SOLVES" in policy
assert "partialSpectrumReady" in scc
assert "certifiedFullSolves" in scc
assert "partialPathDisabled = .true." in scc
assert "profPartialDeferred" in scc
assert "profFallbackSolver" in scc
assert "profFallbackSpectral" in scc
assert "profFallbackOccupation" in scc
assert "profFallbackDegeneracy" in scc
assert "profFallbackFermi" in scc
assert "GFN1-fast 2.2.3 SCC profile" in scc
assert "GFN1-fast 2.2.3 gradient profile" in scf

# Both Davidson and exact indexed partial paths must require certification.
assert scc.count("partialSpectrumReady .and.") >= 2

print("PASS source guards 2.2.3")

# 2.2.3 GFN1 overlap-only electronic gradient kernel.
assert "build_dSH0_GFN1_noreset" in ham
start = ham.index("subroutine build_dSH0_GFN1_noreset")
end = ham.index("end subroutine build_dSH0_GFN1_noreset", start)
gfn1_grad = ham[start:end]
assert "get_grad_overlap" in gfn1_grad
assert "get_grad_multiint" not in gfn1_grad
assert "vd(" not in gfn1_grad and "vq(" not in gfn1_grad
assert "2.0_wp*HPij - 2.0_wp*Pew" in gfn1_grad
assert "XTB_GFN1_FAST_DISABLE_GRADIENT_KERNEL" in policy
assert "fastPolicy%gradientKernel" in scf
assert "GFN1 overlap kernel" in scf
assert "GFN1-fast 2.2.3 gradient profile" in scf
# Generic fallback remains available.
assert "call build_dSDQH0_noreset" in scf

print("PASS GFN1 overlap-only gradient guards 2.2.3")


# 2.2.3 thermal-aware partial eigensolver.
assert "thermalPartial" in policy
assert "thermalWarmupFullSolves" in policy
assert "thermalGuardRoots" in policy
assert "thermalTailKBT" in policy
assert "XTB_GFN1_FAST_DISABLE_THERMAL_PARTIAL" in policy
assert "XTB_GFN1_FAST_THERMAL_WARMUP_FULL_SOLVES" in policy
assert "XTB_GFN1_FAST_THERMAL_GUARD_ROOTS" in policy
assert "gfn1PartialSpectrumReadyThermal" in policy
assert "gfn1ThermalRootGuard" in policy
assert ".not.fastPolicy%thermalPartial" in scc
assert "thermalGuardRoots = gfn1ThermalRootGuard" in scc
assert "fastPolicy%thermalTailKBT" in scc
assert "thermalFermiFallbacks" in scc
assert "thermal retries=" in scc
assert "GFN1-fast 2.2.3 SCC profile" in scc
assert "GFN1-fast 2.2.3 gradient profile" in scf
print("PASS thermal-aware partial guards 2.2.3")


# 2.2.3 runtime autotuning and true SYEVR workspace reuse.
eigensolve = (root / "src/mctc/lapack/eigensolve.F90").read_text()
assert "partialAutotune" in policy
assert "XTB_GFN1_FAST_DISABLE_PARTIAL_AUTOTUNE" in policy
assert "XTB_GFN1_FAST_PARTIAL_AUTOTUNE_MAX_RATIO" in policy
assert "gfn1PartialAutotuneKeep" in policy
assert "partial autotune=" in scc
assert "partial path    :" in scc
assert "full solve core :" in scc
assert "lapack_syevr('v', 'i'" in eigensolve
assert "spartial_work" in eigensolve and "dpartial_work" in eigensolve
assert "partial_iwork" in eigensolve and "partial_isuppz" in eigensolve
# The exact indexed path must no longer allocate/deallocate temporary SYEVX work arrays per SCC step.
sp = eigensolve[eigensolve.index("subroutine mctc_ssyevr_factorized"):eigensolve.index("end subroutine mctc_ssyevr_factorized")]
dp = eigensolve[eigensolve.index("subroutine mctc_dsyevr_factorized"):eigensolve.index("end subroutine mctc_dsyevr_factorized")]
assert "lapack_syevx" not in sp and "lapack_syevx" not in dp
assert "allocate(work" not in sp and "allocate(work" not in dp
print("PASS partial autotune/SYEVR guards 2.2.3")
