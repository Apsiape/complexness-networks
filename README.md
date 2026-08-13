# The Resource Theory of Complexness References in Real Quantum Networks

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21907062.svg)](https://doi.org/10.5281/zenodo.21907062)

Real quantum theory with independent sources cannot reproduce certain
complex network correlations (Renou et al., Nature 600, 625 (2021));
the known remedy is a shared two-rebit state Ω = ¼(I − J⊗J), the
"reference frame of complexness" of Weilenmann, Gisin, and Sekatski
(arXiv:2502.20102). This paper develops the resource theory in which
Ω is the unit.

**Main results**

- **Frame locking**: two exact edge self-tests at one source do not fix
  a single frame; for commuting induced complex structures the missing
  datum is one binary observable (the relative frame parity), of whose
  +1 sector Ω is the normalized projector.
- **A graded layer**: the partial references η_t = ¼(I − tJ⊗J) are
  useful above an explicit threshold yet admit no exact finite-copy
  conversion to Ω under real separable operations with charge-free
  ancillas, while distilling Ω under rLOCC at rate at least one.
- **A classification**: every resource on m separated holders carries a
  transpose-charge code H(σ) ≤ F₂^m that such operations can never
  enlarge; frame-code states realize every even-weight code, with a
  maximal reference one-shot obtainable between holders a,b iff
  e_a + e_b is a codeword, and rate zero under real separable operations with charge-free ancillas otherwise.
- **Visibility and activation**: which charges can influence network
  statistics obeys a parity law on the source–party incidence graph;
  for frame-code resources, used charges are reducible to the weight-2
  unit by partition coarsening or by activation consuming declared Ω
  auxiliaries — the frame-sector transcription of the unlocking and
  activation of multipartite bound entanglement.
- **Bound complexness** is partition-relative and exists in the
  activatable form (τ₄ = (1/16)(I + J⊗⁴) with an explicit four-source
  network); all statements are relative to a declared
  separable-source-independence axiom.

## Contents

- `paper/complexness-networks.tex` (and compiled PDF) — the paper.
- `lean/` — a Lean 4 / Mathlib project machine-verifying the algebraic
  cores: transpose-parity preservation under real conjugation (the
  single-map germ of the span-monotonicity theorem), the determinant
  lemma AᵀJA = (det A)·J behind the one-copy gap, the absence of real
  product vectors in supp Ω behind the finite-copy no-go, and (from a
  companion manuscript in preparation) the radical theorem, the
  twelve-ray frame identities, and the 4×4 Pfaffian parity.

## Verifying the Lean artifact

With Lean 4.30.0 and Lake installed (e.g. via elan):

```
cd lean
lake exe cache get
lake build
lake env lean axiom_check.lean
```

The final command prints the axiom dependencies of the seven public
theorems; each depends only on Lean's standard three
(`propext`, `Classical.choice`, `Quot.sound`). No `sorry`, no custom
axioms. The combinatorial network arguments of the paper's later
sections are **not** formalized; the artifact's scope is exactly the
list above.

## Citation

See `CITATION.cff`. Archived at Zenodo: doi:10.5281/zenodo.21907062 (concept, resolves to latest); this version (v1.0.3): doi:10.5281/zenodo.21911918.

## License

Code (`lean/`) under the MIT License; paper text under CC BY 4.0.

## Changelog

- **v1.0.3** (2026-08-12): repairs the v1.0.2 verification artifact
  and sharpens the Corollary 7.3 footnote. The shipped sweep had a
  machine-local import path (unrunnable from a clean clone) and
  reported one biconditional violation without diagnosis; the
  repaired suite (all scripts self-contained, exit codes armed)
  identifies it as a sign-killed exception class: the formal
  alignment/disjointness conditions can fail while the sign
  character kills precisely the offending codewords on every
  surviving branch, so the realized state is the product. The
  footnote now states: sufficiency unconditional; necessity at the
  level of the surviving (positive-probability) sector, where the
  equivalence holds with zero violations (120 randomized setups,
  1,784 surviving branches; diagnostic in t7_exception.py, missing
  module t4_general.py now shipped). Target-pair exactness
  unaffected.

- **v1.0.2** (2026-08-12): scopes Corollary 7.3's general clause
  and strengthens its aligned case. Continued adversarial
  verification found explicit counterexamples to the
  arbitrary-frame-code decoupled-product claim (one rank-obstructed
  against any correction; one entangling the auxiliaries with the
  resource); the corrected statement gives the exact
  alignment/disjointness characterization (necessary and
  sufficient), the closed correction rule flip iff (-1)^r h s = +1,
  bonus references from overlapping covered codewords, and the
  unconditional guarantee that the target-pair reference is exact
  on every branch. Verification engines ship in `verify-descent/`.
  Also records wgs2025's journal publication (PRL 135, 180201).
  No other statement is affected.

- **v1.0.1** (2026-08-12): corrects the accessibility corollary (Corollary 7.3 in the compiled numbering). Version 1.0.0
  stated the declared-auxiliary descent accounting as "unit-neutral
  at best"; the paper's own Theorem 6.3 protocol is nondemolition on
  the auxiliary pair's frame character, so the descent of a used
  weight-2k codeword consumes k-1 auxiliaries and yields k pairwise
  references, each auxiliary returned in kind (net +1, catalytic).
  Found in internal adversarial review; verified by explicit
  computation over all measurement branches at |c| = 4 and |c| = 6.
  No other statement is affected.
- **v1.0.0** (2026-08-12): initial release.
