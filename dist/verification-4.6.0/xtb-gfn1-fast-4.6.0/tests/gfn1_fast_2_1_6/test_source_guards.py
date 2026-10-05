from pathlib import Path

root = Path(__file__).resolve().parents[2]
ham = (root / "src/xtb/hamiltonian.f90").read_text()
scf = (root / "src/scf_module.F90").read_text()
policy = (root / "src/gfn1_fast_policy.f90").read_text()

start = ham.index("subroutine build_dSDQH0_noreset")
end = ham.index("end subroutine build_dSDQH0_noreset", start)
noreset = ham[start:end]
assert "stmp" not in noreset, "2.1.6 noreset gradient path must not reference uninitialized stmp"
assert "g_xyz(ixyz) = g_xyz(ixyz)+dtmp+qtmp" in noreset

# Ensure the periodic/general routine still retains its explicitly assigned overlap term.
general = ham[ham.index("subroutine build_dSDQH0("):start]
assert "stmp =" in general and "+stmp+dtmp+qtmp" in general

# Profiling labels are part of the 2.1.6 regression surface.
assert "GFN1-fast 2.1.6 SCC profile" in (root / "src/scc_core.f90").read_text()
assert "GFN1-fast 2.1.6 gradient profile" in scf
assert "gfn1SpectralPartialUnsafe" in policy

# Reproducible full-Debug build guards added after re-evaluating the 6.7.1 build.
cmake_top = (root / "CMakeLists.txt").read_text()
cmake_opts = (root / "cmake/CMakeLists.txt").read_text()
debug_script = (root / "tests/gfn1_fast_2_1_6/run_debug_build.sh").read_text()
assert 'option(WITH_TESTS "Build xtb test suite" TRUE)' in cmake_opts
assert 'if(WITH_TESTS AND NOT TARGET "test-drive::test-drive")' in cmake_top
assert 'if(WITH_TESTS)\n  add_subdirectory("test")' in cmake_top
assert '-DWITH_TBLITE=OFF' in debug_script
assert '-DWITH_CPCMX=OFF' in debug_script
assert '-DWITH_TESTS=OFF' in debug_script
assert '-DWITH_JSON=OFF' in debug_script
assert 'MCTCLIB_FIND_METHOD=cmake;pkgconf;subproject' in debug_script
assert 'finit-real=snan' in debug_script and 'finit-integer=-999999' in debug_script
prep_script = (root / 'tests/gfn1_fast_2_1_6/prepare_debug_deps.sh').read_text()
assert '77f65c6f2cf6330d05d0757ca173da097096780e' in prep_script

print("PASS source guards 2.1.6")
