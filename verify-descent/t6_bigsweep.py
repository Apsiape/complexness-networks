import sys, itertools, random
sys.path.insert(0, r"C:\Users\PC\AppData\Local\Temp\claude\c--Infanox-finite-contact\5e738f3c-f46b-4ebc-9d44-8215213f7a82\scratchpad\fc_ref")
from cases import *
from t4_general import predict_general, structure_flags, omega_product_after_flips

random.seed(7)
H = 8                      # holders 0,1 target ; 2,3,4,5 measured ; 6,7 spare
auxpairs = [(2, 3), (4, 5)]
c1 = [0, 1, 2, 3, 4, 5]
n_setup = 0; mism = 0; iff_viol = 0; marg_bad = 0; psd_bad = 0; prod_yes = 0; nbr = 0
seen = set()
trials = 0
while n_setup < 120 and trials < 4000:
    trials += 1
    k = random.choice([2, 3])
    gens = [c1]
    for _ in range(k - 1):
        while True:
            w = [h for h in range(H) if random.random() < 0.45]
            if len(w) % 2 == 0 and w: break
        gens.append(sorted(w))
    key = tuple(map(tuple, gens))
    if key in seen: continue
    seen.add(key)
    q = [random.randint(0, 1) for _ in gens]
    try:
        S = Setup(H=H, auxpairs=auxpairs, gens=gens, qvals=q)
    except AssertionError:
        continue
    n_setup += 1
    H2, H3, why, Rs = structure_flags(S, {u: 1 for u in S.U})
    allprod = True
    for bvals in itertools.product([1, -1], repeat=4):
        xx = dict(zip(S.U, bvals))
        pB, outB = engineB(S, xx)
        pg = predict_general(S, xx)
        if sym_add(outB, sym_scale(pg, -1.0)): mism += 1
        if abs(pB) < 1e-13: continue
        nbr += 1
        nrm = normalise(outB, len(S.keep))
        if omega_product_after_flips(S, outB) is not True: allprod = False
        M = sym_to_dense(nrm, len(S.keep), sites_order=S.keep)
        if np.linalg.eigvalsh(M).min() < -1e-10: psd_bad += 1
        others = [s for s in S.keep if s not in (0, 1)]
        m = sym_ptrace(nrm, others); m = sym_scale(m, 1.0 / sym_trace(m, 2))
        if abs(abs(m.get(frozenset([0, 1]), 0.0) * 4) - 1.0) > 1e-10: marg_bad += 1
    if allprod: prod_yes += 1
    if (H2 and H3) != allprod: iff_viol += 1

print("RANDOM SWEEP, H=8 holders (+4 aux = 12 rebits), codes of dim 2-3, random signs")
print("  setups                                  : %d  (branches: %d)" % (n_setup, nbr))
print("  general closed-form mismatches          : %d   <-- must be 0" % mism)
print("  violations of  (H2 and H3) <=> product  : %d   <-- must be 0" % iff_viol)
print("  setups achieving the Omega-product      : %d / %d" % (prod_yes, n_setup))
print("  positivity failures                     : %d" % psd_bad)
print("  target-pair marginal not sharp +-Omega  : %d   <-- must be 0" % marg_bad)
