# GFN1-fast 4.7.0: analytic Hessian OpenMP

The 4.7.0 implementation parallelizes susceptibility
columns, molecular coordinate responses, AO first/second derivative assembly
and solvent SASA second derivatives. Immutable reference caches are shared;
scratch is created inside reentrant workers and results use disjoint columns
or per-worker geometric matrices. Numerical fallback and all 4.6.0 guards
remain in place. Release/Debug correctness checks have passed; the frozen
performance comparison and source distribution are recorded in RESULTS.md.

Coordinate trace contractions use four-column BLAS batches rather than
forming a full AO density-response tensor. Born-radius pair Hessians update
only the affected six coordinates; nonlinear radius derivatives retain the
complete dense second derivative. SASA screening factors scatter their six
coordinate derivatives without dense zero Hessians or outer-product
temporaries. Coupled electronic applications reuse C^T dS C in each worker;
the zero overlap perturbation in susceptibility assembly skips this transform.

The central native/CLI resource policy accounts for per-worker AO and
geometric scratch. It reduces workers to respect the unchanged 1 GiB budget,
uses the serial path for small systems or nested OpenMP, and preserves the
numerical fallback if even the serial workspace exceeds the cap. Current
defaults are 64 AO minimum and eight workers maximum. Counts also respect the
OpenMP runtime, number of shells/coordinates and the workspace budget. The
estimated cap covers analytic scratch, not all external process allocations;
Windows process memory is measured separately in the benchmark records.

Environment controls are `XTB_GFN1_FAST_DISABLE_HESSIAN_OPENMP=1`,
`XTB_GFN1_FAST_HESSIAN_OMP_MAX_THREADS` and
`XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO`. These do not disable existing SCF/Gradient
OpenMP. MKL remains sequential to avoid nested BLAS teams.

`run_windows_parallel.ps1` exercises actual 2/4-worker teams, compares full
CLI matrices and the explicit serial switch, and covers gas, zero/high
temperature and open-shell cases in Release and Debug. It also checks eight
workers, an explicit two-worker maximum, a two-worker runtime limit and four
numerical fallbacks. These total 19 comparisons per build. The native solvent
runner accepts `-Threads 4` or `-Threads 8`; forcing minimum AO to one exercises
teams in all 38 small molecular Gradient-FD cases. Water has fewer than eight
shells, so its selected team is limited accordingly. The gas response test
also verifies optional density and shell-charge response tensors against
nonlinear SCC finite differences.

`run_windows_resources.ps1` checks concurrent first policy calls, the worker
cap, reduction under a 900 MiB budget, nested teams and the small-system
crossover. Its checked standalone program is compiled both with and without
OpenMP. `run_windows_calculator.ps1 -Threads 4` checks a real nested analytic
call and partial additive Hessians, retaining polarizability fallback.

`run_windows_probe.ps1` runs single profiled taxol stage/resource diagnostics.
It checks actual runtime team size and complete 339-by-339 matrix parity at
1/4/8 threads. These single-run timings are not repeated benchmarks. The large
native test independently differentiates converged nonlinear SCC gradients
for all 339 raw columns in GBSA and ALPB. At eight requested threads it passes
with maximum errors 9.42e-9 and 1.60e-8 Eh/bohr^2, respectively, at 1e-5 bohr.

`run_windows_benchmark.ps1` uses the frozen 4.6.0 and 4.7.0 executables, one
warm-up and five randomized measured repeats for gas/GBSA/ALPB. Every complete
matrix, original Energy and Gradient norm is compared; profiled team/resource
audits are separate. It records both serial algorithm improvement and parallel
scaling, without inferring a multithreaded stock-xTB comparison. Earlier 4.6.0
snapshots and ZIP remain unchanged.

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildType Release -BuildName build-gfn1-fast-current-windows-ifx
.\tests\gfn1_fast_4_7_0\run_windows_parallel.ps1
.\tests\gfn1_fast_4_7_0\run_windows_resources.ps1
$env:XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO = '1'
.\tests\gfn1_fast_4_5_0\run_windows_molecular.ps1 -Threads 4
.\tests\gfn1_fast_4_1_0\run_windows_molecular.ps1 -Threads 8
```

Use `-BuildType Debug` on the build and native solvent runners, and
`-LibraryBuildType Debug` on the gas response runner, for checked libraries.
For production use the default threshold and set the runtime explicitly:

```powershell
Remove-Item Env:XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO -ErrorAction SilentlyContinue
$env:OMP_NUM_THREADS = '8'
$env:MKL_NUM_THREADS = '1'
.\build-gfn1-fast-current-windows-ifx-release\fast-gfn1-xtb.exe input.xyz --gfn 1 --hess --norestart
```
