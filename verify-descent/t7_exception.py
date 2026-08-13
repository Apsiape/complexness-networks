import sys, itertools, random
from cases import *
from t4_general import predict_general, structure_flags, omega_product_after_flips

random.seed(7)
H = 8; auxpairs = [(2, 3), (4, 5)]; c1 = [0, 1, 2, 3, 4, 5]
seen = set(); n_setup = 0; trials = 0; found = []
while n_setup < 120 and trials < 4000:
    trials += 1
    k = random.choice([2, 3]); gens = [c1]
    for _ in range(k - 1):
        while True:
            w = [h for h in range(H) if random.random() < 0.45]
            if len(w) % 2 == 0 and w: break
        gens.append(sorted(w))
    key = tuple(map(tuple, gens))
    if key in seen: continue
    seen.add(key)
    q = [random.randint(0, 1) for _ in gens]
    try: S = Setup(H=H, auxpairs=auxpairs, gens=gens, qvals=q)
    except AssertionError: continue
    n_setup += 1
    H2, H3, why, Rs = structure_flags(S, {u: 1 for u in S.U})
    allprod = True; zero = 0
    for bvals in itertools.product([1, -1], repeat=4):
        pB, outB = engineB(S, dict(zip(S.U, bvals)))
        if abs(pB) < 1e-13: zero += 1; continue
        if omega_product_after_flips(S, outB) is not True: allprod = False
    if (H2 and H3) != allprod:
        found.append((gens, q, H2, H3, why, allprod, zero, Rs))

for gens, q, H2, H3, why, allprod, zero, Rs in found:
    print("EXCEPTION")
    print("  gens =", gens, " q =", q)
    S = Setup(H=H, auxpairs=auxpairs, gens=gens, qvals=q)
    print("  words:", [(sorted(w), int(s)) for w, b, s in S.words])
    print("  R-image (all codewords) =", sorted([sorted(r) for r in Rs], key=lambda z: (len(z), z)))
    print("  H2(aligned)=%s  H3(disjoint-pairs)=%s (%s)   product-on-all-branches=%s  zerobranches=%d"
          % (H2, H3, why, allprod, zero))
    # which words actually contribute on a surviving branch?
    for bvals in itertools.product([1, -1], repeat=4):
        pB, outB = engineB(S, dict(zip(S.U, bvals)))
        if abs(pB) < 1e-13: continue
        nrm = normalise(outB, len(S.keep))
        auxs = set(S.aux.values())
        Vpart = sorted([sorted(k - auxs) for k in support_group(nrm)], key=lambda z: (len(z), z))
        print("  branch x=%s p=%.4f  V-strings actually present: %s"
              % (bvals, pB, sorted(set(map(tuple, Vpart)))))
        break
    print("  -> the FULL R-image is not pair-generated, but the codewords whose R would")
    print("     break it are killed on every surviving branch (zero-probability sector),")
    print("     so the state that actually occurs IS a pair product.")
