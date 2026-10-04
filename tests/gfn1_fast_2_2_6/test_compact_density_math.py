import math
import random

random.seed(226)
nao = 12
nocc = 5
nshell = 4
ao2sh = [i % nshell for i in range(nao)]
focc = [2.0, 1.999, 1.7, 0.25, 0.01]

# Symmetric overlap and arbitrary MO coefficient matrix.  The algebraic
# identity P = (C sqrt(f))(C sqrt(f))^T does not require orthonormal C.
S = [[0.0]*nao for _ in range(nao)]
for i in range(nao):
    for j in range(i+1):
        v = (1.0 if i == j else 0.03*random.uniform(-1,1))
        S[i][j] = S[j][i] = v
C = [[random.uniform(-0.5,0.5) for _ in range(nocc)] for _ in range(nao)]

P_legacy = [[0.0]*nao for _ in range(nao)]
for i in range(nao):
    for j in range(nao):
        P_legacy[i][j] = sum(C[i][m]*focc[m]*C[j][m] for m in range(nocc))
B = [[C[i][m]*math.sqrt(max(0.0,focc[m])) for m in range(nocc)] for i in range(nao)]
P_compact = [[0.0]*nao for _ in range(nao)]
for i in range(nao):
    for j in range(i,nao):
        P_compact[i][j] = sum(B[i][m]*B[j][m] for m in range(nocc))

max_p = max(abs(P_legacy[i][j]-P_compact[i][j]) for i in range(nao) for j in range(i,nao))
assert max_p < 2e-15, max_p

# Check the fused upper-triangle Mulliken and H0 contractions against the
# original two independent traversals.
H0 = []
for i in range(nao):
    for j in range(i+1):
        H0.append(random.uniform(-3.0, 1.0))

q_ref = [0.0]*nshell
for i in range(nao):
    ii = ao2sh[i]
    for j in range(i):
        jj = ao2sh[j]
        ps = P_legacy[j][i]*S[j][i]
        q_ref[ii] += ps
        q_ref[jj] += ps
    q_ref[ii] += P_legacy[i][i]*S[i][i]

k = 0
h_ref = 0.0
for i in range(nao):
    for j in range(i):
        h_ref += P_legacy[j][i]*H0[k]
        k += 1
    h_ref += 0.5*P_legacy[i][i]*H0[k]
    k += 1

q_fused = [0.0]*nshell
k = 0
h_fused = 0.0
for i in range(nao):
    ii = ao2sh[i]
    for j in range(i):
        jj = ao2sh[j]
        ps = P_compact[j][i]*S[j][i]
        q_fused[ii] += ps
        q_fused[jj] += ps
        h_fused += P_compact[j][i]*H0[k]
        k += 1
    ps = P_compact[i][i]*S[i][i]
    q_fused[ii] += ps
    h_fused += 0.5*P_compact[i][i]*H0[k]
    k += 1

assert max(abs(a-b) for a,b in zip(q_ref,q_fused)) < 5e-15
assert abs(h_ref-h_fused) < 5e-14
print('PASS compact weighted-SYRK and fused contraction math 2.2.6')
