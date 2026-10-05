"""Run native GFN1 Energy/Gradient and compare production/fallback paths."""

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
NUMBER = r"[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?"


def run_case(exe, output, name, xyz, options=(), overrides=None):
    work = output / name
    work.mkdir(parents=True, exist_ok=True)
    (work / "input.xyz").write_text(xyz, encoding="ascii")
    env = {k: v for k, v in os.environ.items() if not k.startswith("XTB_GFN1_FAST_")}
    # Exercise the normal Windows HOME-missing startup, including Debug's
    # allocatable checks, without changing the parent process environment.
    env.pop("HOME", None)
    env.pop("XTBHOME", None)
    # Python normalizes Windows environment keys. Launchers containing both
    # PATH and Path can lose setvars' PATH during that normalization.
    runtime_dirs = [str(Path(env[key]) / "bin") for key in ("CMPLR_ROOT", "MKLROOT") if key in env]
    env["PATH"] = os.pathsep.join([*runtime_dirs, env.get("PATH", "")])
    env.update({key: "1" for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                                   "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS")})
    env.update(OMP_STACKSIZE="64M", KMP_STACKSIZE="64M")
    env.update(XTBPATH=str(ROOT), XTB_GFN1_FAST_PROFILE="1")
    env.update(overrides or {})
    command = [str(exe), "input.xyz", "--gfn", "1", "--grad", "--norestart", *options]
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=work, env=env, text=True,
                               encoding="utf-8", errors="replace", capture_output=True, timeout=180)
    elapsed = time.perf_counter() - started
    log = completed.stdout + "\n" + completed.stderr
    (work / "run.log").write_text(log, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"{name}: exit {completed.returncode}; see {work / 'run.log'}")
    if "convergence criteria satisfied" not in log:
        raise AssertionError(f"{name}: SCC convergence not reported")
    energy = float(re.search(r"TOTAL ENERGY\s+(" + NUMBER + ")", log).group(1))
    norm = float(re.search(r"GRADIENT NORM\s+(" + NUMBER + ")", log).group(1))
    gradient_lines = (work / "gradient").read_text().splitlines()
    # Turbomole gradient file: header, N coordinate records, N gradient records.
    nat = int(xyz.splitlines()[0])
    numeric_lines = [line.split() for line in gradient_lines
                     if re.fullmatch(r"\s*" + NUMBER + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s*", line)]
    if len(numeric_lines) != nat:
        raise AssertionError(f"{name}: expected {nat} gradient records, got {len(numeric_lines)}")
    gradient = [float(value.replace("D", "E").replace("d", "e"))
                for line in numeric_lines for value in line]
    if not all(math.isfinite(value) for value in [energy, norm, *gradient]):
        raise AssertionError(f"{name}: nonfinite Energy/Gradient")
    if abs(math.sqrt(sum(value * value for value in gradient)) - norm) > 2e-10:
        raise AssertionError(f"{name}: gradient file/norm disagreement")
    # Profile-off runs still report the normal basis summary.
    nao_match = re.search(r"NAO=(\d+)", log) or re.search(r"#\s+atomic orbitals\s+(\d+)", log)
    nao = int(nao_match.group(1))
    compact = re.search(r"compact density SCC=([TF]), min NAO=\d+, iterations=(\d+)", log)
    result = {"name": name, "energy_Eh": energy, "gradient_norm_Eh_bohr": norm,
              "nao": nao, "compact_iterations": int(compact.group(2)) if compact else None,
              "wall_seconds": elapsed,
              "iterations": int(re.search(r"satisfied after (\d+) iterations", log).group(1)),
              "log": str(work / "run.log")}
    print(f"PASS {name}: NAO={nao}, E={energy:.12f}, |G|={norm:.12f}", flush=True)
    return result, gradient, log


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--small-only", action="store_true", help="Debug startup/ALPB checks on water only")
    args = parser.parse_args()
    exe, output = args.exe.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    water = "3\nwater\nO 0 0 0\nH 0.758602 0 0.504284\nH -0.758602 0 0.504284\n"
    taxol = (ROOT / "assets/inputs/xyz/taxol.xyz").read_text()
    cases = [
        ("water_gas", water, (), {}),
        ("water_alpb", water, ("--alpb", "water"), {}),
        ("taxol_fast", taxol, (), {}),
        ("taxol_legacy_density", taxol, (), {"XTB_GFN1_FAST_DISABLE_COMPACT_DENSITY": "1"}),
        ("taxol_generic_gradient", taxol, (), {"XTB_GFN1_FAST_DISABLE_GRADIENT_KERNEL": "1"}),
        ("taxol_alpb_fast", taxol, ("--alpb", "water"), {}),
        ("taxol_alpb_legacy_density", taxol, ("--alpb", "water"),
         {"XTB_GFN1_FAST_DISABLE_COMPACT_DENSITY": "1"}),
    ]
    if args.small_only:
        cases = cases[:2]
    results = {}
    for name, xyz, options, overrides in cases:
        results[name] = run_case(exe, output, name, xyz, options, overrides)
    if args.small_only:
        summary = {"executable": str(exe), "cases": [r[0] for r in results.values()],
                   "comparisons": [], "threads": 1, "scope": "water startup and ALPB"}
        (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return
    for name in ("taxol_fast", "taxol_alpb_fast"):
        if results[name][0]["compact_iterations"] < 1:
            raise AssertionError(f"{name}: compact density path not exercised")
    for name in ("taxol_legacy_density", "taxol_alpb_legacy_density"):
        if results[name][0]["compact_iterations"] != 0:
            raise AssertionError(f"{name}: legacy density override not effective")
    if "used=F" not in results["taxol_generic_gradient"][2]:
        raise AssertionError("generic gradient fallback not exercised")
    comparisons = []
    for fast, reference in (("taxol_fast", "taxol_legacy_density"),
                            ("taxol_fast", "taxol_generic_gradient"),
                            ("taxol_alpb_fast", "taxol_alpb_legacy_density")):
        delta_e = abs(results[fast][0]["energy_Eh"] - results[reference][0]["energy_Eh"])
        delta_g = max(abs(a - b) for a, b in zip(results[fast][1], results[reference][1]))
        if delta_e > 1e-11 or delta_g > 1e-10:
            raise AssertionError(f"{fast}/{reference}: dE={delta_e}, max dG={delta_g}")
        comparisons.append({"fast": fast, "reference": reference,
                            "abs_delta_energy_Eh": delta_e, "max_abs_delta_gradient_Eh_bohr": delta_g})
        print(f"PASS parity {fast}/{reference}: dE={delta_e:.3e}, max dG={delta_g:.3e}", flush=True)
    summary = {"executable": str(exe), "cases": [r[0] for r in results.values()],
               "comparisons": comparisons, "threads": 1}
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
