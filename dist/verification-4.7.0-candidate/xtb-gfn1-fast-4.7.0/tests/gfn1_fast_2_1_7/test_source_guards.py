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

# Strict-Debug findings promoted into the 2.1.7 source baseline.
assert "er = 0.0_wp" in main
assert "exist = .false." in main
assert "sigma = 0.0_wp" in scf

# 2.1.7 certified-spectrum dispatcher and one-failure auto-disable.
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
assert "GFN1-fast 2.1.7 SCC profile" in scc
assert "GFN1-fast 2.1.7 gradient profile" in scf

# Both Davidson and exact indexed partial paths must require certification.
assert scc.count("partialSpectrumReady .and.") >= 2

print("PASS source guards 2.1.7")
