# GFN1-fast 2.4.1 / 4.4.0 Gradient comparison

Same Windows native ifx 2025.2 / oneMKL sequential environment, one thread,
profiling off, one warmup and seven randomized interleaved process timings
per mode. All timed Energy/Gradient pairs pass the established stock parity
checks. A separate profile run verifies use of the direct-prefactor path.
Raw trials and executable hashes: `GRADIENT_BENCHMARK_20261005.json`.

| System | Stock / s | Previous fast Gradient / s | Corrected Gradient / s | Stock speedup |
|---|---:|---:|---:|---:|
| Disilane | 0.036165 | 0.034965 | 0.035439 | 1.02× |
| Taxol gas, 350 AO | 0.329856 | 0.222324 | 0.221678 | 1.49× |
| Taxol ALPB water | 0.329134 | 0.221912 | 0.224156 | 1.47× |

Compared with the previous fast policy, the corrected path is 0.29% faster
for taxol gas, 1.01% slower for taxol ALPB and 1.36% slower for the small
disilane case. These small differences establish no additional speedup.
The prefactor change is adopted for derivative correctness at exact nodes;
the complete fast implementation retains its measured stock speedup.
This comparison makes no claim about the prior WSL run or missing 594-AO input.

## Large taxol full Hessian: completed

The complete 339×339 Hessian of taxol (113 atoms, 350 AO) uses the preserved
4.4 executable, one warmup and five randomized profile-off stock/analytic
trials, plus a separate profile audit. Every matrix comparison passes with
maximum difference `2.991e-7 Eh/bohr²`. The audit reports actual analytic
dispatch; the executable hash remained unchanged during all measurements.

| Mode | Median process wall / s | Five timed samples / s |
|---|---:|---|
| Stock numerical | 250.238602 | 249.0892, 250.2386, 250.6061, 253.3931, 248.7179 |
| Analytic 4.4 | 30.070810 | 30.0708, 31.1309, 29.7632, 30.8039, 29.8760 |

The measured speedup over stock is **8.32×** for this one-thread gas-phase
case. Hessian SCC accuracy is `1e-7`, displacement step `5e-4 bohr`.
Energy and Gradient norm differences in every timed pair satisfy `1e-10`.
No componentwise stock Gradient equality is claimed for untested nodal
geometries. Full trials, matrix hashes and policy/profile evidence are in
`LARGE_BENCHMARK_20261005.json`; the completed log is also preserved here.

During the separate profile audit, a live process sample observed private
memory 841,154,560 bytes and peak working set 831,291,392 bytes. This is a
single process observation, not a complete peak-allocation measurement.
The calculator's conservative preallocation estimate was 872,726,600 bytes.

The stock numerical Hessian uses the documented one-thread Windows ifx
compatibility repair: only its Hessian OpenMP directives are suppressed;
the original numerical differentiation formula and calculations remain.
No multithreaded stock crossover or solvent speedup follows from this result.
Keep the implementation available for validated gas-phase cases; broader
screening consistency and performance-based default dispatch still require
their own evidence.
