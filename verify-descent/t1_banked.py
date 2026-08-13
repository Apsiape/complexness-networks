import sys
sys.path.insert(0, r"C:\Users\PC\AppData\Local\Temp\claude\c--Infanox-finite-contact\5e738f3c-f46b-4ebc-9d44-8215213f7a82\scratchpad\fc_ref")
from cases import *

# ---- banked 1: |c|=4, one auxiliary pair.  target {0,1}, measured {2,3}
S1 = Setup(H=4, auxpairs=[(2, 3)], gens=[[0, 1, 2, 3]], qvals=[0])
run_case("BANKED-1  |c|=4 descent (r should be 1, flip at h=-1)", S1)

# ---- banked 2: |c|=6, two auxiliary pairs. target {0,1}, measured {2,3,4,5}
S2 = Setup(H=6, auxpairs=[(2, 3), (4, 5)], gens=[[0, 1, 2, 3, 4, 5]], qvals=[0])
run_case("BANKED-2  |c|=6 descent (r should be 2, flip at h=+1)", S2)

# report the flip rule explicitly for these two
for nm, S in (("|c|=4", S1), ("|c|=6", S2)):
    print("\nflip table for %s :" % nm)
    for bvals in itertools.product([1, -1], repeat=len(S.U)):
        x = dict(zip(S.U, bvals))
        pred, diag, al = predict(S, x)
        d = [d for d in diag if len(d['word'])][0]
        print("   x=%s  h=%+d  r=%d  coeff(J_target)=%+d  -> flip needed: %s"
              % (bvals, d['h'], d['r'], d['coeff'], d['coeff'] > 0))
