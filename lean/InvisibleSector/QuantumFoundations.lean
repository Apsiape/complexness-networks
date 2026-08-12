/-
QuantumFoundations.lean — machine-checked results accompanying
the companion manuscript "Dynamics as the invisible sector"
(in preparation).

Three targets:
1. Radical theorem: a real matrix X has tr(Xρ) = 0 for every real
   PSD ρ iff X is skew-symmetric (generalized to any finite index).
2. The twelve rays (e_i ± e_j)/√2 form a tight frame on ℝ⁴
   (Σ qqᵀ = 3I) with the reconstruction identity
   Σ (qᵀρq) qqᵀ = ρ + (tr ρ / 2) I for symmetric ρ.
3. Pfaffian parity: pf₄(−J) = pf₄(J) for every 4×4 matrix.
-/
import Mathlib

open Matrix Finset

namespace QuantumFoundations

variable {n : Type*} [Fintype n] [DecidableEq n]

/-! ## Helpers -/

/-- Standard basis vector as an explicit if-function. -/
def e (i : n) : n → ℝ := fun t => if t = i then 1 else 0

lemma trace_mul_outer (X : Matrix n n ℝ) (v : n → ℝ) :
    (X * vecMulVec v v).trace = v ⬝ᵥ X.mulVec v := by
  simp only [Matrix.trace, Matrix.diag, Matrix.mul_apply,
    Matrix.vecMulVec_apply, dotProduct, Matrix.mulVec, Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ =>
    Finset.sum_congr rfl fun j _ => by ring

lemma posSemidef_outer (v : n → ℝ) : (vecMulVec v v).PosSemidef := by
  simpa using Matrix.posSemidef_vecMulVec_self_star v

lemma dot_e_mulVec (X : Matrix n n ℝ) (i j : n) :
    e i ⬝ᵥ X.mulVec (e j) = X i j := by
  simp [e, dotProduct, Matrix.mulVec, mul_ite, ite_mul,
    Finset.sum_ite_eq', Finset.mul_sum]

/-- Quadratic form of `e_i + c e_j`, computed symbolically — this is
what keeps the twelve-ray computations small. -/
lemma quad_pair (X : Matrix n n ℝ) (i j : n) (c : ℝ) :
    (e i + c • e j) ⬝ᵥ X.mulVec (e i + c • e j) =
      X i i + c * (X i j + X j i) + c * c * X j j := by
  simp only [Matrix.mulVec_add, Matrix.mulVec_smul, dotProduct_add,
    add_dotProduct, smul_dotProduct, dotProduct_smul, dot_e_mulVec,
    smul_eq_mul]
  ring

/-! ## Theorem 1: the radical of single-time receivers -/

/-- A real matrix annihilates every real density (PSD suffices; the
trace-one normalization is irrelevant to a linear condition) iff it
is skew-symmetric. The invisible sector of single-time real
statistics is exactly the skew part. -/
theorem annihilates_real_densities_iff_skew (X : Matrix n n ℝ) :
    (∀ ρ : Matrix n n ℝ, ρ.PosSemidef → (X * ρ).trace = 0) ↔
      Xᵀ = -X := by
  constructor
  · intro h
    have hq : ∀ v : n → ℝ, v ⬝ᵥ X.mulVec v = 0 := fun v => by
      have := h (vecMulVec v v) (posSemidef_outer v)
      rwa [trace_mul_outer] at this
    have hdiag : ∀ i, X i i = 0 := fun i => by
      have := hq (e i)
      rwa [dot_e_mulVec] at this
    have hoff : ∀ i j, X i i + 1 * (X i j + X j i) + 1 * 1 * X j j = 0 := by
      intro i j
      have h0 := hq (e i + (1 : ℝ) • e j)
      rwa [quad_pair] at h0
    ext i j
    have h1 := hoff i j
    have h2 := hdiag i
    have h3 := hdiag j
    simp only [Matrix.transpose_apply, Matrix.neg_apply]
    linarith
  · intro hX ρ hρ
    have hsym : ρᵀ = ρ := by
      ext i j
      simpa using congrFun (congrFun hρ.1 i) j
    have h1 : (X * ρ).trace = -(X * ρ).trace := by
      conv_lhs => rw [← Matrix.trace_transpose, Matrix.transpose_mul,
        hX, hsym]
      rw [Matrix.mul_neg, Matrix.trace_neg, Matrix.trace_mul_comm]
    linarith

/-- Specialization to `Fin n`. -/
theorem radical_fin (m : ℕ) (X : Matrix (Fin m) (Fin m) ℝ) :
    (∀ ρ : Matrix (Fin m) (Fin m) ℝ, ρ.PosSemidef →
      (X * ρ).trace = 0) ↔ Xᵀ = -X :=
  annihilates_real_densities_iff_skew X

/-! ## Theorem 2: the twelve-ray tight frame on ℝ⁴ -/

abbrev M4 := Matrix (Fin 4) (Fin 4) ℝ

/-- The six unordered pairs of distinct indices in `Fin 4`. -/
def pairIdx : Fin 6 → Fin 4 × Fin 4
  | 0 => (0, 1)
  | 1 => (0, 2)
  | 2 => (0, 3)
  | 3 => (1, 2)
  | 4 => (1, 3)
  | 5 => (2, 3)

def sgn : Fin 2 → ℝ
  | 0 => 1
  | 1 => -1

/-- Unnormalized ray `e_i ± e_j`. -/
def uray (k : Fin 6) (s : Fin 2) : Fin 4 → ℝ :=
  e (pairIdx k).1 + sgn s • e (pairIdx k).2

lemma uray_apply (k : Fin 6) (s : Fin 2) (t : Fin 4) :
    uray k s t = (if t = (pairIdx k).1 then (1 : ℝ) else 0) +
      sgn s * (if t = (pairIdx k).2 then (1 : ℝ) else 0) := by
  simp [uray, e]

/-- Normalized ray `(e_i ± e_j)/√2`. -/
noncomputable def ray12 (k : Fin 6) (s : Fin 2) : Fin 4 → ℝ :=
  (Real.sqrt 2)⁻¹ • uray k s

lemma vecMulVec_smul_smul (c : ℝ) (v : n → ℝ) :
    vecMulVec (c • v) (c • v) = (c * c) • vecMulVec v v := by
  ext i j
  simp only [Matrix.vecMulVec_apply, Pi.smul_apply, smul_eq_mul,
    Matrix.smul_apply]
  ring

lemma dot_smul_quad (c : ℝ) (u : n → ℝ) (M : Matrix n n ℝ) :
    (c • u) ⬝ᵥ M.mulVec (c • u) = (c * c) * (u ⬝ᵥ M.mulVec u) := by
  simp only [dotProduct, Matrix.mulVec, Pi.smul_apply, smul_eq_mul,
    Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ =>
    Finset.sum_congr rfl fun j _ => by ring

lemma inv_sqrt_two_mul_self :
    (Real.sqrt 2)⁻¹ * (Real.sqrt 2)⁻¹ = (2 : ℝ)⁻¹ := by
  rw [← mul_inv, Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 2)]

set_option maxHeartbeats 1600000
set_option maxRecDepth 8192

/-- `fin_cases` produces `⟨k, h⟩`-form indices; these `rfl` lemmas
normalize them to numeral form so `Fin.reduceEq` can fire and atoms
unify. -/
lemma fin4_mk0 (h : 0 < 4) : (⟨0, h⟩ : Fin 4) = 0 := rfl
lemma fin4_mk1 (h : 1 < 4) : (⟨1, h⟩ : Fin 4) = 1 := rfl
lemma fin4_mk2 (h : 2 < 4) : (⟨2, h⟩ : Fin 4) = 2 := rfl
lemma fin4_mk3 (h : 3 < 4) : (⟨3, h⟩ : Fin 4) = 3 := rfl

/-- Unnormalized tightness: `Σ u uᵀ = 6 I`. -/
theorem uray_tight :
    (∑ k : Fin 6, ∑ s : Fin 2, vecMulVec (uray k s) (uray k s)) =
      (6 : ℝ) • (1 : M4) := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    norm_num [fin4_mk0, fin4_mk1, fin4_mk2, fin4_mk3, Fin.ext_iff,
      Fin.val_ofNat, Matrix.sum_apply, Fin.sum_univ_six,
      Fin.sum_univ_two, Matrix.add_apply, Matrix.vecMulVec_apply,
      uray_apply, pairIdx, sgn, Matrix.smul_apply, Matrix.one_apply,
      smul_eq_mul]

/-- Tight frame: `Σ q qᵀ = 3 I` for the twelve normalized rays. -/
theorem rays12_tight :
    (∑ k : Fin 6, ∑ s : Fin 2,
      vecMulVec (ray12 k s) (ray12 k s)) = (3 : ℝ) • (1 : M4) := by
  have hv : ∀ (k : Fin 6) (s : Fin 2),
      vecMulVec (ray12 k s) (ray12 k s) =
        ((2 : ℝ)⁻¹) • vecMulVec (uray k s) (uray k s) := by
    intro k s
    rw [ray12, vecMulVec_smul_smul, inv_sqrt_two_mul_self]
  simp only [hv, ← Finset.smul_sum]
  rw [uray_tight, smul_smul]
  norm_num

/-- Unnormalized reconstruction: for symmetric ρ,
`Σ (uᵀρu) u uᵀ = 4 ρ + 2 (tr ρ) I`. -/
theorem uray_reconstruction (ρ : M4) (hρ : ρ.IsSymm) :
    (∑ k : Fin 6, ∑ s : Fin 2,
      ((uray k s) ⬝ᵥ ρ.mulVec (uray k s)) •
        vecMulVec (uray k s) (uray k s)) =
      (4 : ℝ) • ρ + (2 * ρ.trace) • (1 : M4) := by
  have hsym : ∀ i j, ρ j i = ρ i j := fun i j => by
    simpa using congrFun (congrFun hρ i) j
  have h10 : ρ 1 0 = ρ 0 1 := hsym 0 1
  have h20 : ρ 2 0 = ρ 0 2 := hsym 0 2
  have h30 : ρ 3 0 = ρ 0 3 := hsym 0 3
  have h21 : ρ 2 1 = ρ 1 2 := hsym 1 2
  have h31 : ρ 3 1 = ρ 1 3 := hsym 1 3
  have h32 : ρ 3 2 = ρ 2 3 := hsym 2 3
  have hq : ∀ (k : Fin 6) (s : Fin 2),
      (uray k s) ⬝ᵥ ρ.mulVec (uray k s) =
        ρ (pairIdx k).1 (pairIdx k).1 +
          sgn s * (ρ (pairIdx k).1 (pairIdx k).2 +
            ρ (pairIdx k).2 (pairIdx k).1) +
          sgn s * sgn s * ρ (pairIdx k).2 (pairIdx k).2 :=
    fun k s => quad_pair ρ (pairIdx k).1 (pairIdx k).2 (sgn s)
  ext a b
  fin_cases a <;> fin_cases b <;>
    (norm_num [fin4_mk0, fin4_mk1, fin4_mk2, fin4_mk3, Fin.ext_iff,
      Fin.val_ofNat, Matrix.sum_apply, Fin.sum_univ_six,
      Fin.sum_univ_two, Fin.sum_univ_four, Matrix.add_apply,
      Matrix.smul_apply, smul_eq_mul, hq, Matrix.vecMulVec_apply,
      uray_apply, pairIdx, sgn, Matrix.one_apply, Matrix.trace,
      Matrix.diag];
     try simp only [h10, h20, h30, h21, h31, h32];
     try ring)

/-- Reconstruction for the normalized rays: for symmetric ρ,
`Σ (qᵀρq) q qᵀ = ρ + (tr ρ / 2) I`. -/
theorem rays12_reconstruction (ρ : M4) (hρ : ρ.IsSymm) :
    (∑ k : Fin 6, ∑ s : Fin 2,
      ((ray12 k s) ⬝ᵥ ρ.mulVec (ray12 k s)) •
        vecMulVec (ray12 k s) (ray12 k s)) =
      ρ + (ρ.trace / 2) • (1 : M4) := by
  have hterm : ∀ (k : Fin 6) (s : Fin 2),
      ((ray12 k s) ⬝ᵥ ρ.mulVec (ray12 k s)) •
          vecMulVec (ray12 k s) (ray12 k s) =
        ((4 : ℝ)⁻¹) • (((uray k s) ⬝ᵥ ρ.mulVec (uray k s)) •
          vecMulVec (uray k s) (uray k s)) := by
    intro k s
    rw [ray12, dot_smul_quad, inv_sqrt_two_mul_self,
      vecMulVec_smul_smul, inv_sqrt_two_mul_self, smul_smul,
      smul_smul]
    congr 1
    ring
  simp only [hterm, ← Finset.smul_sum]
  rw [uray_reconstruction ρ hρ, smul_add, smul_smul, smul_smul]
  rw [show ((4 : ℝ)⁻¹ * 4) = 1 by norm_num, one_smul,
    show ((4 : ℝ)⁻¹ * (2 * ρ.trace)) = ρ.trace / 2 by ring]

/-- Density-matrix specialization: trace one makes the correction
term exactly `I/2`. -/
theorem rays12_reconstruction_trace_one (ρ : M4) (hρ : ρ.IsSymm)
    (htr : ρ.trace = 1) :
    (∑ k : Fin 6, ∑ s : Fin 2,
      ((ray12 k s) ⬝ᵥ ρ.mulVec (ray12 k s)) •
        vecMulVec (ray12 k s) (ray12 k s)) =
      ρ + ((2 : ℝ)⁻¹) • (1 : M4) := by
  rw [rays12_reconstruction ρ hρ, htr]
  norm_num

/-! ## Theorem 3: Pfaffian parity in dimension 4 -/

/-- The 4×4 Pfaffian formula (evaluated on the upper triangle). -/
def pf4 (J : M4) : ℝ :=
  J 0 1 * J 2 3 - J 0 2 * J 1 3 + J 0 3 * J 1 2

/-- Parity: `pf₄(−J) = pf₄(J)` for EVERY 4×4 matrix — the parity is
a property of the Pfaffian polynomial (each monomial is quadratic in
the entries); antisymmetry is not needed. -/
theorem pf4_neg (J : M4) : pf4 (-J) = pf4 J := by
  simp only [pf4, Matrix.neg_apply]
  ring

/-- The skew-specialized statement, an immediate corollary. -/
theorem pf4_neg_skew (J : M4) (_ : Jᵀ = -J) : pf4 (-J) = pf4 J :=
  pf4_neg J

end QuantumFoundations
