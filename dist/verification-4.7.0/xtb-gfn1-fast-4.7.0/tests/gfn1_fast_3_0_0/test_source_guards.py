"""Architectural guards complement numerical tests; these do not prove parity."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[2]
operator = (root/'src/gfn1_fermi_operator.f90').read_text().lower()
scc = (root/'src/scc_core.f90').read_text()
policy = (root/'src/gfn1_fast_policy.f90').read_text()

# A density-only implementation must not call or import an eigensolver.
assert not re.search(r'\buse\s+\w*(?:eigensolve|stdeigval|geneigval)', operator)
calls = set(re.findall(r'\bcall\s+(lapack_\w+)', operator))
assert calls <= {'lapack_potrf', 'lapack_sygst', 'lapack_potrs', 'lapack_sytrf'}
assert re.search(r'logical\s*::\s*fermiOperatorSCC\s*=\s*\.false\.', policy)
for setting in ('ENABLE_FERMI_OPERATOR', 'DISABLE_FERMI_OPERATOR',
                'DISABLE_FERMI_OPERATOR_AUTOTUNE', 'FERMI_OPERATOR_MAX_STEPS',
                'FERMI_OPERATOR_AUDIT_PERIOD'):
    assert 'XTB_GFN1_FAST_' + setting in policy
assert re.search(r'if\(foeEnabled\.and\.\.not\.foeDisabled\.and\.\.not\.lastdiag', scc)
assert 'call auditFermiDensity(' in scc
assert 'call validateFermiDensity(' in scc
assert 'audit electron count error: FOE=' in scc
print('PASS diagonalization-free architecture, opt-in policy and final full-solve guards 3.0.0')
