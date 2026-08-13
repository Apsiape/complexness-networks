import sys
from cases import *

# =====================================================================
# KILLER A : two fully-covered, PAIR-ALIGNED codewords whose UNMEASURED
#            supports INTERSECT in one holder.
#   holders 0,1,2 kept ; 3,4,5,6 measured ; aux pairs (3,4),(5,6)
#   c1 = {0,1,3,4,5,6}  wt6, target {0,1}, D1 = all 4 measured (2 aux pairs)
#   c2 = {1,2,3,4}      wt4, D2 = {3,4} = ONE whole aux pair, R2 = {1,2}
#   R1 = {0,1} and R2 = {1,2} share holder 1.
# =====================================================================
SA = Setup(H=7, auxpairs=[(3, 4), (5, 6)],
           gens=[[0, 1, 3, 4, 5, 6], [1, 2, 3, 4]], qvals=[0, 0])
vA, _ = run_case("KILLER A : intersecting unmeasured supports (11 rebits)",
                 SA, dense=True, dense_limit=2)

x = dict(zip(SA.U, (1, 1, 1, 1)))
pB, outB = engineB(SA, x)
nrm = normalise(outB, len(SA.keep))
print("\n  output state on kept sites %s, branch x=(+,+,+,+):" % SA.keep)
for S_, v in sorted(nrm.items(), key=lambda kv: (len(kv[0]), sorted(kv[0]))):
    print("     J_%-22s coeff = %+.6f  (x 2^%d = %+.0f)"
          % (sorted(S_), v, len(SA.keep), v * 2 ** len(SA.keep)))

# rank / local-equivalence obstruction on the kept RESOURCE rebits {0,1,2}
auxs = set(SA.aux.values())
res3 = sym_ptrace(nrm, sorted(auxs))
res3 = sym_scale(res3, 1.0 / sym_trace(res3, 3))
M = sym_to_dense(res3, 3, sites_order=[0, 1, 2])
ev = np.linalg.eigvalsh(M)
print("\n  reduced state on kept resource rebits {0,1,2}:")
print("     eigenvalues =", np.round(ev, 8), " rank =", int(np.sum(ev > 1e-9)))
tgt = np.kron(0.25 * (np.eye(4) - jstring_dense([0, 1], 2)), np.eye(2) / 2)
print("     rank of  Omega_{0,1} (x) I/2  =", int(np.linalg.matrix_rank(tgt)))
print("     -> product form claim requires rank 4; actual rank %d  => IMPOSSIBLE by ANY"
      % int(np.sum(ev > 1e-9)))
print("        local (or even global) unitary/orthogonal correction.")
m01 = sym_ptrace(nrm, sorted(auxs) + [2]); m01 = sym_scale(m01, 1.0 / sym_trace(m01, 2))
print("     BUT target-pair MARGINAL on {0,1} =", {tuple(sorted(k)): round(v, 6) for k, v in m01.items()},
      " -> = Omega up to the prescribed flip:", abs(abs(m01.get(frozenset([0, 1]), 0)) - 0.25) < 1e-12)

# =====================================================================
# KILLER B : a fully-covered codeword whose measured support is NOT a
#            union of aux pairs (one holder from each of two aux pairs)
#   c2 = {1,2,3,5}
# =====================================================================
SB = Setup(H=7, auxpairs=[(3, 4), (5, 6)],
           gens=[[0, 1, 3, 4, 5, 6], [1, 2, 3, 5]], qvals=[0, 0])
run_case("KILLER B : measured support NOT aux-pair-aligned (11 rebits)",
         SB, dense=True, dense_limit=2)
pB, outB = engineB(SB, x)
nrm = normalise(outB, len(SB.keep))
print("\n  output state, branch x=(+,+,+,+)   (aux sites are %s):" % sorted(SB.aux.values()))
for S_, v in sorted(nrm.items(), key=lambda kv: (len(kv[0]), sorted(kv[0]))):
    tag = "  <-- RESOURCE-AUX CROSS TERM" if (set(S_) & set(SB.aux.values())
                                              and set(S_) - set(SB.aux.values())) else ""
    print("     J_%-24s coeff = %+.6f%s" % (sorted(S_), v, tag))
# are the auxiliaries returned in kind?
for t, (u, v_) in enumerate(SB.auxpairs):
    a, b = SB.aux[u], SB.aux[v_]
    others = [s for s in SB.keep if s not in (a, b)]
    red = sym_ptrace(nrm, others); red = sym_scale(red, 1.0 / sym_trace(red, 2))
    isOm = abs(red.get(frozenset([a, b]), 0.0) + 0.25) < 1e-12
    print("     aux pair %s marginal == Omega ? %s   (coeff of K = %+.4f)"
          % ((a, b), isOm, red.get(frozenset([a, b]), 0.0) * 4))
