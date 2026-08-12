/-
ComplexnessNetworks.lean — machine-checked cores for the paper
"The Resource Theory of Complexness References in Real Quantum
Networks."

Three targets (the algebraic hearts of Sections 4-6):
1. Transpose-parity preservation: real conjugation maps preserve
   the symmetric/skew parity of an operator — the single-map germ
   of the charge-span monotonicity theorem.
2. The determinant lemma: Aᵀ J A = (det A) • J for real 2×2 A —
   the engine of the one-copy approximation gap.
3. No real product vector lies in span{Φ⁻, Ψ⁺} = supp Ω — the
   engine of the no-exact-finite-copy-conversion theorem.
-/
import Mathlib

open Matrix

namespace ComplexnessNetworks

/-! ## 1. Transpose-parity preservation -/

/-- A real conjugation map `X ↦ K X Kᵀ` preserves transpose parity:
if `Xᵀ = ε • X` then the image satisfies the same relation with the
same ε. Applied to each Kraus term of a real CP map and each local
tensor factor, this is the mechanism of charge-sector preservation. -/
theorem conj_preserves_parity {n m : Type*} [Fintype n] [Fintype m]
    (K : Matrix m n ℝ) (X : Matrix n n ℝ) (ε : ℝ)
    (h : Xᵀ = ε • X) :
    (K * X * Kᵀ)ᵀ = ε • (K * X * Kᵀ) := by
  rw [Matrix.transpose_mul, Matrix.transpose_mul,
    Matrix.transpose_transpose, ← Matrix.mul_assoc, h,
    Matrix.mul_smul, Matrix.smul_mul]

/-- Symmetric case (charge 0). -/
theorem conj_preserves_symm {n m : Type*} [Fintype n] [Fintype m]
    (K : Matrix m n ℝ) (X : Matrix n n ℝ) (h : Xᵀ = X) :
    (K * X * Kᵀ)ᵀ = K * X * Kᵀ := by
  have := conj_preserves_parity K X 1 (by simpa using h)
  simpa using this

/-- Skew case (charge 1). -/
theorem conj_preserves_skew {n m : Type*} [Fintype n] [Fintype m]
    (K : Matrix m n ℝ) (X : Matrix n n ℝ) (h : Xᵀ = -X) :
    (K * X * Kᵀ)ᵀ = -(K * X * Kᵀ) := by
  have := conj_preserves_parity K X (-1)
    (by simpa [neg_smul, one_smul] using h)
  simpa [neg_smul, one_smul] using this

/-! ## 2. The determinant lemma -/

abbrev M2 := Matrix (Fin 2) (Fin 2) ℝ

/-- The standard complex structure on ℝ². -/
def J2 : M2 := !![0, -1; 1, 0]

/-- `Aᵀ J A = (det A) • J` for every real 2×2 matrix. Consequence
used in the paper: local real Kraus factors scale the `J ⊗ J`
moment by `det A · det B`, bounding any separable branch's overlap
with the reference and yielding the one-copy gap `(1−t)/2`. -/
theorem transpose_J_conj (A : M2) :
    Aᵀ * J2 * A = A.det • J2 := by
  have hdet : A.det = A 0 0 * A 1 1 - A 0 1 * A 1 0 :=
    Matrix.det_fin_two A
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [J2, Matrix.mul_apply, Matrix.transpose_apply,
      Matrix.smul_apply, Fin.sum_univ_two, hdet] <;>
    ring

/-! ## 3. No real product vector in the support of Ω -/

/-- Unnormalized Φ⁻ component function on index pairs. -/
def phiM : Fin 2 → Fin 2 → ℝ := fun i j =>
  if i = 0 ∧ j = 0 then 1 else if i = 1 ∧ j = 1 then -1 else 0

/-- Unnormalized Ψ⁺ component function on index pairs. -/
def psiP : Fin 2 → Fin 2 → ℝ := fun i j =>
  if i = 0 ∧ j = 1 then 1 else if i = 1 ∧ j = 0 then 1 else 0

/-- The support of Ω = ¼(I − J⊗J) is span{Φ⁻, Ψ⁺}, and it contains
no nonzero real product vector: if `v ⊗ w = a·Φ⁻ + b·Ψ⁺` entrywise
then `a = b = 0`. This is the engine of the theorem that no finite
number of copies of a partial reference converts exactly to Ω under
stochastic real separable operations. -/
theorem no_real_product_vector (v w : Fin 2 → ℝ) (a b : ℝ)
    (h : ∀ i j, v i * w j = a * phiM i j + b * psiP i j) :
    a = 0 ∧ b = 0 := by
  have h00 : v 0 * w 0 = a := by simpa [phiM, psiP] using h 0 0
  have h11 : v 1 * w 1 = -a := by simpa [phiM, psiP] using h 1 1
  have h01 : v 0 * w 1 = b := by simpa [phiM, psiP] using h 0 1
  have h10 : v 1 * w 0 = b := by simpa [phiM, psiP] using h 1 0
  have key : (v 0 * w 0) * (v 1 * w 1) = (v 0 * w 1) * (v 1 * w 0) := by
    ring
  rw [h00, h11, h01, h10] at key
  constructor
  · nlinarith [sq_nonneg a, sq_nonneg b]
  · nlinarith [sq_nonneg a, sq_nonneg b]

/-- Packaging: the determinant of the matrix `[[a, b], [b, -a]]`
associated to `a·Φ⁻ + b·Ψ⁺` is `-(a² + b²)`, so it vanishes only at
`a = b = 0` — the coordinate-free reason for the theorem above. -/
theorem assoc_det (a b : ℝ) :
    (!![a, b; b, -a] : M2).det = -(a^2 + b^2) := by
  simp [Matrix.det_fin_two]
  ring

end ComplexnessNetworks
