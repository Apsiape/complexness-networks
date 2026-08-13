# Descent-accounting verification (Corollary 7.3, v1.0.2)

Two independent engines for the declared-auxiliary descent of
frame-code states: `lib.py` (dense real-matrix) and `proto.py`
(symbolic J-string coefficient algebra). `t1_banked.py` reproduces
the |c|=4, |c|=6, and overlapping-code instances exhaustively;
`t3_killers.py` exhibits the two counterexamples that force the
v1.0.2 scoping (intersecting output pairs: rank-obstructed
correlated output; non-aligned covered word: auxiliaries entangled
with the resource); `t6_bigsweep.py` runs the alignment/
disjointness biconditional over 124 exhaustive + 150 randomized
setups and the unconditional target-pair check over all branches.
`independent_killer_a_check.py` is a third, separately written
reproduction of the rank-obstructed counterexample. NumPy only.
