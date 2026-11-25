/-
  prize_work/flagship-erdos-86-f7/lean/Stars.lean

  **Star capacity, and what it gives without any assumption:**

      f(4) <= 24,  f(5) <= 56,  f(6) <= 133,  f(7) <= 306,  f(8) <= 691,

  where f(n) is the largest number of edges of a C4-free subgraph of the cube Q_n.

  A "full three-star" is a vertex v together with three edges of H at v.  Its three directions
  span one 3-dimensional face.  Two different vertices of the same parity in one 3-face are two
  flips apart, so if both carried a full star in that face, the four edges between them would
  be a 4-cycle.  Hence a C4-free H has at most one star per 3-face and parity side:

      sum over the vertices v of one parity side of C(deg v, 3)  <=  C(n,3) * 2^(n-3).

  With C(d,3) >= a d - c (a tangent line) and "sum of degrees = 2 e(H)" this bounds e(H).

  What is and is not claimed.  The star-capacity inequality is Lemma 2.2 of Zhou, Ma, Yang, Ren,
  Yang, Li, "The Four-Cycle Extremal Number of the Seven-Cube", arXiv:2611.27284 (20 Nov 2026),
  which proves f(7) = 304; counting of this kind is classical.  No new theorem about f(7) is
  claimed here.  The file is a machine-checked proof of the inequality for every n and of the
  five numerical bounds above.  f(4) = 24 and f(5) = 56 are the known exact values, so the bound
  is sharp there; the true values f(6) = 132 and f(7) = 304 are smaller than what this file
  proves.  For n = 8 the lab found no published upper bound to compare with (OEIS A245762 lists
  only a(8) >= 680; averaging from f(7) = 304 gives 694), so "f(8) <= 691" is to be read as
  "elementary, machine-checked, novelty not established" until a referee has checked the
  literature (Brass-Harborth-Nienborg 1995 was not available to the lab).

  The file stands alone (the lab's Lean gate checks single files); section 1 is copied word for
  word from Cube.lean (lemma F1).

  MAIN RESULTS
    * `root_unique`    -- two same-parity vertices with a full star in the same 3-face are equal.
    * `roots_card`     -- for a 3-set T of directions and a parity p: 8 * #roots <= 2^n.
    * `star_side`      -- 8 * (sum over one parity side of C(deg v, 3)) <= C(n,3) * 2^n.
    * `star_capacity`  -- both sides together.
    * `edge_bound`     -- the general edge bound from a tangent line a d <= C(d,3) + c.
    * `bound_4_24`, `bound_5_56`, `bound_6_133`, `bound_7_306`, `bound_8_691`.

  Written and checked by the Vector session of 2026-11-25
  (cd prize_work/lean ; python check_claim.py ../flagship-erdos-86-f7/lean/Stars.lean).
-/

import Mathlib

open Finset
open scoped symmDiff

namespace Stars

/-! ## 1. The cube, word for word as in Cube.lean -/

/-- Vertex set of the `n`-cube `Q_n`: bit vectors of length `n`. -/
abbrev V (n : ℕ) := Fin n → Bool

/-- The coordinate set on which a vertex is `true`. -/
def supp {n : ℕ} (x : V n) : Finset (Fin n) := univ.filter fun i => x i = true

theorem mem_supp {n : ℕ} (x : V n) (i : Fin n) : i ∈ supp x ↔ x i = true := by
  simp [supp]

theorem supp_injective {n : ℕ} : Function.Injective (supp : V n → Finset (Fin n)) := by
  intro x y h
  funext i
  have hxy : x i = true ↔ y i = true := by rw [← mem_supp x i, h, mem_supp y i]
  cases hx : x i <;> cases hy : y i <;> simp_all

/-! ## 2. Flipping one coordinate -/

/-- Flip one coordinate of a vertex. -/
def flip {n : ℕ} (x : V n) (i : Fin n) : V n := fun k => if k = i then !x k else x k

theorem flip_self {n : ℕ} (x : V n) (i : Fin n) : flip x i i = !x i := by
  simp [flip]

theorem flip_of_ne {n : ℕ} {x : V n} {i j : Fin n} (h : j ≠ i) : flip x i j = x j := by
  simp [flip, h]

theorem flip_flip_self {n : ℕ} (x : V n) (i : Fin n) : flip (flip x i) i = x := by
  funext k
  by_cases hk : k = i
  · subst hk; simp [flip]
  · rw [flip_of_ne hk, flip_of_ne hk]

theorem flip_comm {n : ℕ} (x : V n) (i j : Fin n) : flip (flip x i) j = flip (flip x j) i := by
  by_cases h : i = j
  · subst h; rfl
  · funext k
    by_cases hk : k = i
    · subst hk
      rw [flip_self, flip_of_ne h]
      simp [flip, h]
    · by_cases hj : k = j
      · subst hj
        rw [flip_of_ne (Ne.symm h), flip_self]
        simp [flip, hk]
      · rw [flip_of_ne hj, flip_of_ne hk, flip_of_ne hk, flip_of_ne hj]

/-- Flipping one coordinate is injective in the coordinate: two one-step neighbours of `x`
coincide only if the flipped coordinate does. -/
theorem flip_inj {n : ℕ} {x : V n} {i j : Fin n} (h : flip x i = flip x j) : i = j := by
  by_contra hne
  have hi := congr_fun h i
  simp only [flip, if_neg hne] at hi
  cases hxi : x i <;> simp_all

theorem supp_flip {n : ℕ} (x : V n) (i : Fin n) : supp (flip x i) = supp x ∆ {i} := by
  ext k
  rw [mem_supp, mem_symmDiff, mem_supp, mem_singleton]
  by_cases hk : k = i
  · subst hk; simp [flip_self]
  · rw [flip_of_ne hk]; simp [hk]

theorem supp_flip_flip {n : ℕ} (x : V n) (i j : Fin n) :
    supp (flip (flip x i) j) = supp x ∆ ({i} ∆ {j}) := by
  rw [supp_flip, supp_flip, symmDiff_assoc]

/-- A flipped vertex is a different vertex. -/
theorem flip_ne {n : ℕ} (x : V n) (i : Fin n) : flip x i ≠ x := by
  intro h
  have := congr_fun h i
  simp only [flip, if_pos] at this
  cases hx : x i <;> simp_all

/-! ## 3. Hamming distance and adjacency -/

/-- Hamming distance: the number of coordinates where two vertices differ. -/
def hamDist {n : ℕ} (x y : V n) : ℕ := (univ.filter fun i => x i ≠ y i).card

theorem hamDist_eq_card_symmDiff {n : ℕ} (x y : V n) :
    hamDist x y = (supp x ∆ supp y).card := by
  unfold hamDist
  congr 1
  ext i
  rw [mem_filter, mem_symmDiff, mem_supp, mem_supp]
  cases x i <;> cases y i <;> simp

/-- Adjacency in the cube: Hamming distance one. -/
def Adj {n : ℕ} (x y : V n) : Prop := hamDist x y = 1

theorem adj_irrefl {n : ℕ} (x : V n) : ¬ Adj x x := by
  unfold Adj hamDist
  simp

theorem hamDist_comm {n : ℕ} (x y : V n) : hamDist x y = hamDist y x := by
  rw [hamDist_eq_card_symmDiff, hamDist_eq_card_symmDiff]
  exact congrArg Finset.card (symmDiff_comm _ _)

theorem adj_comm {n : ℕ} {x y : V n} (h : Adj x y) : Adj y x := by
  unfold Adj at *
  rw [hamDist_comm]
  exact h

theorem hamDist_self {n : ℕ} (x : V n) : hamDist x x = 0 := by
  unfold hamDist; simp

theorem hamDist_flip {n : ℕ} (x : V n) (i : Fin n) : hamDist x (flip x i) = 1 := by
  rw [hamDist_eq_card_symmDiff, supp_flip, symmDiff_symmDiff_cancel_left, card_singleton]

theorem adj_flip {n : ℕ} (x : V n) (i : Fin n) : Adj x (flip x i) := hamDist_flip x i

theorem adj_flip' {n : ℕ} (x : V n) (i : Fin n) : Adj (flip x i) x := adj_comm (adj_flip x i)

/-- Adjacent vertices are different. -/
theorem adj_ne {n : ℕ} {x y : V n} (h : Adj x y) : x ≠ y := by
  rintro rfl; exact adj_irrefl x h

/-- **Adjacency is exactly a single coordinate flip.** -/
theorem adj_iff_flip {n : ℕ} {x y : V n} : Adj x y ↔ ∃ i, y = flip x i := by
  constructor
  · intro h
    have hc : (supp x ∆ supp y).card = 1 := by
      rw [← hamDist_eq_card_symmDiff]; exact h
    obtain ⟨i, hi⟩ := card_eq_one.mp hc
    refine ⟨i, supp_injective ?_⟩
    show supp y = supp (flip x i)
    rw [supp_flip]
    have hcon := congrArg (fun t => supp x ∆ t) hi
    rw [symmDiff_symmDiff_cancel_left] at hcon
    exact hcon
  · rintro ⟨i, rfl⟩
    exact adj_flip x i

/-- Symmetric difference of two distinct singletons is the unordered pair. -/
theorem symmDiff_singletons {n : ℕ} {i j : Fin n} (h : i ≠ j) :
    ({i} ∆ {j} : Finset (Fin n)) = {i, j} := by
  ext k
  rw [mem_symmDiff, mem_singleton, mem_insert, mem_singleton]
  constructor
  · rintro (⟨rfl, _⟩ | ⟨rfl, _⟩) <;> simp
  · rintro (rfl | rfl)
    · exact Or.inl ⟨rfl, h⟩
    · exact Or.inr ⟨rfl, Ne.symm h⟩

/-- Two flips of one vertex are at Hamming distance two. -/
theorem hamDist_two_flips {n : ℕ} (x : V n) {i j : Fin n} (hij : i ≠ j) :
    hamDist x (flip (flip x i) j) = 2 := by
  have hset : supp x ∆ supp (flip (flip x i) j) = ({i} ∆ {j} : Finset (Fin n)) := by
    rw [supp_flip_flip]
    ext k
    simp only [mem_symmDiff]
    tauto
  rw [hamDist_eq_card_symmDiff, hset, symmDiff_singletons hij, card_pair hij]

/-- The vertex two flips away differs from the original one. -/
theorem two_flips_ne {n : ℕ} (x : V n) {i j : Fin n} (hij : i ≠ j) :
    x ≠ flip (flip x i) j := by
  intro h
  have h1 : hamDist x (flip (flip x i) j) = 2 := hamDist_two_flips x hij
  rw [← h, hamDist_self] at h1
  exact absurd h1 (by decide)

/-- Two different single flips of one vertex differ. -/
theorem flip_ne_flip {n : ℕ} (x : V n) {i j : Fin n} (hij : i ≠ j) :
    flip x i ≠ flip x j := by
  intro h
  exact hij (flip_inj h)

/-- `e` is an edge of `Q_n`: an unordered pair of adjacent vertices. -/
def IsEdge {n : ℕ} (e : Finset (V n)) : Prop := ∃ x y : V n, Adj x y ∧ e = {x, y}

/-- The edge set `H` contains a 4-cycle. -/
def HasC4 {n : ℕ} (H : Finset (Finset (V n))) : Prop :=
  ∃ a b c d : V n, a ≠ c ∧ b ≠ d ∧
    ({a, b} : Finset (V n)) ∈ H ∧ ({b, c} : Finset (V n)) ∈ H ∧
    ({c, d} : Finset (V n)) ∈ H ∧ ({d, a} : Finset (V n)) ∈ H

/-- "f(n) ≤ M": every C4-free set of edges of `Q_n` has at most `M` edges. -/
def Bound (n M : ℕ) : Prop :=
  ∀ H : Finset (Finset (V n)), (∀ e ∈ H, IsEdge e) → ¬ HasC4 H → H.card ≤ M

/-! ## 2. Directions at a vertex, parity, flipping a set of coordinates -/

/-- The directions of the edges of `H` at the vertex `v`; its size is the degree of `v`. -/
def nbr {n : ℕ} (H : Finset (Finset (V n))) (v : V n) : Finset (Fin n) :=
  univ.filter fun i => ({v, flip v i} : Finset (V n)) ∈ H

theorem mem_nbr {n : ℕ} {H : Finset (Finset (V n))} {v : V n} {i : Fin n} :
    i ∈ nbr H v ↔ ({v, flip v i} : Finset (V n)) ∈ H := by
  simp [nbr]

/-- Flip every coordinate of `x` that lies in `D`. -/
def flipSet {n : ℕ} (x : V n) (D : Finset (Fin n)) : V n := fun k => if k ∈ D then !x k else x k

/-- Parity of the Hamming weight, as an element of `ZMod 2`. -/
def par {n : ℕ} (x : V n) : ZMod 2 := ∑ k, if x k = true then 1 else 0

/-- Flipping the coordinates in `D` changes the parity by `|D|`. -/
theorem par_flipSet {n : ℕ} (x : V n) (D : Finset (Fin n)) :
    par (flipSet x D) = par x + (D.card : ZMod 2) := by
  unfold par flipSet
  have h : ∀ k ∈ (univ : Finset (Fin n)),
      (if (if k ∈ D then !x k else x k) = true then (1 : ZMod 2) else 0)
        = (if x k = true then 1 else 0) + (if k ∈ D then 1 else 0) := by
    intro k _
    by_cases hk : k ∈ D
    · simp only [hk, ↓reduceIte]
      cases x k <;> decide
    · simp only [hk, ↓reduceIte]
      cases x k <;> decide
  rw [Finset.sum_congr rfl h, Finset.sum_add_distrib]
  congr 1
  simp

theorem flipSet_empty {n : ℕ} (x : V n) : flipSet x ∅ = x := by
  funext k
  simp [flipSet]

theorem flipSet_pair {n : ℕ} (x : V n) {a b : Fin n} (hab : a ≠ b) :
    flipSet x {a, b} = flip (flip x a) b := by
  funext k
  unfold flipSet flip
  by_cases hka : k = a
  · subst hka
    simp [hab]
  · by_cases hkb : k = b
    · subst hkb
      simp [hka]
    · simp [hka, hkb]

/-! ## 3. One star per 3-face and parity side -/

/-- Two vertices two flips apart cannot both carry edges of `H` in those two directions. -/
theorem c4_of_two {n : ℕ} {H : Finset (Finset (V n))} {r : V n} {a b : Fin n} (hab : a ≠ b)
    (ha : a ∈ nbr H r) (hb : b ∈ nbr H r)
    (ha' : a ∈ nbr H (flip (flip r a) b)) (hb' : b ∈ nbr H (flip (flip r a) b)) : HasC4 H := by
  rw [mem_nbr] at ha hb ha' hb'
  rw [flip_flip_self] at hb'
  have h3 : flip (flip (flip r a) b) a = flip r b := by
    rw [flip_comm (flip r a) b a, flip_flip_self]
  rw [h3] at ha'
  refine ⟨r, flip r a, flip (flip r a) b, flip r b, two_flips_ne r hab, flip_ne_flip r hab,
    ha, ?_, ha', ?_⟩
  · rw [pair_comm]
    exact hb'
  · rw [pair_comm]
    exact hb

/-- **One root per 3-face and side.** In a C4-free `H`, two vertices of the same parity that
agree outside a 3-set `T` of directions and both have edges of `H` in all directions of `T`
are equal. -/
theorem root_unique {n : ℕ} {H : Finset (Finset (V n))} (hC : ¬ HasC4 H) {T : Finset (Fin n)}
    (hT : T.card = 3) {r r' : V n} (hr : T ⊆ nbr H r) (hr' : T ⊆ nbr H r')
    (hpar : par r = par r') (hoff : ∀ k, k ∉ T → r k = r' k) : r = r' := by
  set D := T.filter (fun m => r m ≠ r' m) with hD
  have hDT : D ⊆ T := filter_subset _ _
  have hr'eq : r' = flipSet r D := by
    funext m
    unfold flipSet
    by_cases hm : m ∈ D
    · have hne : r m ≠ r' m := (mem_filter.mp hm).2
      simp only [hm, ↓reduceIte]
      revert hne
      cases r m <;> cases r' m <;> decide
    · simp only [hm, ↓reduceIte]
      by_cases hmT : m ∈ T
      · by_contra hne
        exact hm (mem_filter.mpr ⟨hmT, fun h => hne h.symm⟩)
      · exact (hoff m hmT).symm
  have hcard : D.card ≤ 3 := hT ▸ card_le_card hDT
  have hparD : (D.card : ZMod 2) = 0 := by
    have h := par_flipSet r D
    rw [← hr'eq, ← hpar] at h
    linear_combination -h
  have h4 : D.card = 0 ∨ D.card = 1 ∨ D.card = 2 ∨ D.card = 3 := by omega
  rcases h4 with hc | hc | hc | hc
  · rw [hr'eq, card_eq_zero.mp hc, flipSet_empty]
  · rw [hc] at hparD
    exfalso
    revert hparD
    decide
  · obtain ⟨a, b, hab, hDab⟩ := card_eq_two.mp hc
    exfalso
    apply hC
    have ha : a ∈ T := hDT (by rw [hDab]; simp)
    have hb : b ∈ T := hDT (by rw [hDab]; simp)
    have hr'2 : r' = flip (flip r a) b := by rw [hr'eq, hDab, flipSet_pair r hab]
    exact c4_of_two hab (hr ha) (hr hb) (hr'2 ▸ hr' ha) (hr'2 ▸ hr' hb)
  · rw [hc] at hparD
    exfalso
    revert hparD
    decide

/-- The vertices of parity `p` that carry edges of `H` in all directions of `T`. -/
def roots {n : ℕ} (H : Finset (Finset (V n))) (p : ZMod 2) (T : Finset (Fin n)) :
    Finset (V n) :=
  univ.filter fun v => par v = p ∧ T ⊆ nbr H v

/-- For a 3-set `T`: the 3-faces of direction `T` are disjoint 8-sets of vertices and each
holds at most one root, so `8 * #roots ≤ 2^n`. -/
theorem roots_card {n : ℕ} {H : Finset (Finset (V n))} (hC : ¬ HasC4 H) (p : ZMod 2)
    (T : Finset (Fin n)) (hT : T.card = 3) : (roots H p T).card * 8 ≤ 2 ^ n := by
  have hinj : Set.InjOn (fun q : V n × Finset (Fin n) => flipSet q.1 q.2)
      (↑((roots H p T) ×ˢ T.powerset) : Set (V n × Finset (Fin n))) := by
    rintro ⟨r, S⟩ h1 ⟨r', S'⟩ h2 heq
    simp only [coe_product, Set.mem_prod, mem_coe, roots, mem_filter, mem_univ, true_and,
      mem_powerset] at h1 h2
    obtain ⟨⟨hp, hr⟩, hS⟩ := h1
    obtain ⟨⟨hp', hr'⟩, hS'⟩ := h2
    have heq' : ∀ k, flipSet r S k = flipSet r' S' k := fun k => congr_fun heq k
    have hrr : r = r' := by
      refine root_unique hC hT hr hr' (hp.trans hp'.symm) ?_
      intro k hk
      have h := heq' k
      have hkS : k ∉ S := fun h' => hk (hS h')
      have hkS' : k ∉ S' := fun h' => hk (hS' h')
      simpa [flipSet, hkS, hkS'] using h
    subst hrr
    have hSS : S = S' := by
      ext k
      have h := heq' k
      unfold flipSet at h
      by_cases h1 : k ∈ S <;> by_cases h2 : k ∈ S'
      · exact ⟨fun _ => h2, fun _ => h1⟩
      · simp only [h1, h2, ↓reduceIte] at h
        exfalso
        revert h
        cases r k <;> decide
      · simp only [h1, h2, ↓reduceIte] at h
        exfalso
        revert h
        cases r k <;> decide
      · exact ⟨fun h' => absurd h' h1, fun h' => absurd h' h2⟩
    rw [hSS]
  have h := Finset.card_le_card_of_injOn (t := (univ : Finset (V n)))
    (fun q : V n × Finset (Fin n) => flipSet q.1 q.2) (by intro q _; simp) hinj
  rw [card_product, card_powerset, hT, card_univ] at h
  have hV : Fintype.card (V n) = 2 ^ n := by simp
  rw [hV] at h
  exact h

/-- **Star capacity on one parity side.** -/
theorem star_side {n : ℕ} {H : Finset (Finset (V n))} (hC : ¬ HasC4 H) (p : ZMod 2) :
    (∑ v ∈ univ.filter (fun v : V n => par v = p), (nbr H v).card.choose 3) * 8
      ≤ n.choose 3 * 2 ^ n := by
  have h1 : ∀ v : V n, (nbr H v).card.choose 3
      = ((powersetCard 3 (univ : Finset (Fin n))).filter (fun T => T ⊆ nbr H v)).card := by
    intro v
    rw [← card_powersetCard 3 (nbr H v)]
    congr 1
    ext T
    simp only [mem_powersetCard, mem_filter, subset_univ, true_and]
    exact and_comm
  have h2 : ∑ v ∈ univ.filter (fun v : V n => par v = p), (nbr H v).card.choose 3
      = ∑ T ∈ powersetCard 3 (univ : Finset (Fin n)), (roots H p T).card := by
    simp_rw [h1, card_filter]
    rw [sum_comm]
    refine sum_congr rfl fun T _ => ?_
    rw [← card_filter, filter_filter]
    rfl
  rw [h2, sum_mul]
  calc ∑ T ∈ powersetCard 3 (univ : Finset (Fin n)), (roots H p T).card * 8
      ≤ ∑ _T ∈ powersetCard 3 (univ : Finset (Fin n)), 2 ^ n :=
        sum_le_sum fun T hT => roots_card hC p T (mem_powersetCard.mp hT).2
    _ = n.choose 3 * 2 ^ n := by
        rw [sum_const, card_powersetCard, card_univ, Fintype.card_fin, smul_eq_mul]

/-- **Star capacity**, both parity sides together:
`8 * Σ_v C(deg v, 3) ≤ 2 * C(n,3) * 2^n`, i.e. `Σ_v C(deg v, 3) ≤ 2 * C(n,3) * 2^(n-3)`. -/
theorem star_capacity {n : ℕ} {H : Finset (Finset (V n))} (hC : ¬ HasC4 H) :
    (∑ v : V n, (nbr H v).card.choose 3) * 8 ≤ 2 * (n.choose 3 * 2 ^ n) := by
  rw [← Finset.sum_fiberwise (univ : Finset (V n)) par (fun v => (nbr H v).card.choose 3),
    sum_mul]
  calc ∑ p : ZMod 2, (∑ v ∈ univ.filter (fun v : V n => par v = p), (nbr H v).card.choose 3) * 8
      ≤ ∑ _p : ZMod 2, n.choose 3 * 2 ^ n := sum_le_sum fun p _ => star_side hC p
    _ = 2 * (n.choose 3 * 2 ^ n) := by
        rw [sum_const, card_univ, ZMod.card, smul_eq_mul]

/-! ## 4. Degrees add up to twice the number of edges -/

/-- The edges of `H` through `v` are given by the directions in `nbr H v`. -/
theorem deg_le {n : ℕ} {H : Finset (Finset (V n))} (hH : ∀ e ∈ H, IsEdge e) (v : V n) :
    (H.filter fun e => v ∈ e).card ≤ (nbr H v).card := by
  refine le_trans (card_le_card ?_)
    (card_image_le (s := nbr H v) (f := fun i => ({v, flip v i} : Finset (V n))))
  intro e he
  obtain ⟨heH, hve⟩ := mem_filter.mp he
  obtain ⟨x, y, hxy, rfl⟩ := hH e heH
  obtain ⟨i, rfl⟩ := adj_iff_flip.mp hxy
  rw [mem_image]
  simp only [mem_insert, mem_singleton] at hve
  rcases hve with hv | hv
  · subst hv
    exact ⟨i, mem_nbr.mpr heH, rfl⟩
  · subst hv
    refine ⟨i, ?_, ?_⟩
    · rw [mem_nbr, flip_flip_self, pair_comm]
      exact heH
    · rw [flip_flip_self, pair_comm]

/-- Every edge has two ends. -/
theorem sum_deg {n : ℕ} {H : Finset (Finset (V n))} (hH : ∀ e ∈ H, IsEdge e) :
    ∑ v : V n, (H.filter fun e => v ∈ e).card = 2 * H.card := by
  simp_rw [card_filter]
  rw [sum_comm]
  have h : ∀ e ∈ H, ∑ v : V n, (if v ∈ e then 1 else 0) = 2 := by
    intro e he
    obtain ⟨x, y, hxy, rfl⟩ := hH e he
    rw [← card_filter, filter_mem_eq_inter, univ_inter, card_pair (adj_ne hxy)]
  rw [sum_congr rfl h, sum_const, smul_eq_mul, mul_comm]

/-- `2 e(H) ≤ Σ_v deg v`. -/
theorem two_mul_card_le {n : ℕ} {H : Finset (Finset (V n))} (hH : ∀ e ∈ H, IsEdge e) :
    2 * H.card ≤ ∑ v : V n, (nbr H v).card := by
  rw [← sum_deg hH]
  exact sum_le_sum fun v _ => deg_le hH v

/-! ## 5. The edge bound -/

/-- **Edge bound from star capacity.** If `a d ≤ C(d,3) + c` for all degrees `d ≤ n`, then a
C4-free `H` in `Q_n` has `16 a e(H) ≤ 2 C(n,3) 2^n + 8 c 2^n`. -/
theorem edge_bound {n : ℕ} (a c : ℕ) (hconv : ∀ d ≤ n, a * d ≤ d.choose 3 + c)
    {H : Finset (Finset (V n))} (hH : ∀ e ∈ H, IsEdge e) (hC : ¬ HasC4 H) :
    16 * (a * H.card) ≤ 2 * (n.choose 3 * 2 ^ n) + 8 * (c * 2 ^ n) := by
  have hdeg := two_mul_card_le hH
  have hsum : ∑ v : V n, a * (nbr H v).card ≤ ∑ v : V n, ((nbr H v).card.choose 3 + c) :=
    sum_le_sum fun v _ => hconv _ (by simpa using card_le_univ (nbr H v))
  rw [← mul_sum, sum_add_distrib, sum_const, card_univ, smul_eq_mul] at hsum
  have hV : Fintype.card (V n) = 2 ^ n := by simp
  rw [hV] at hsum
  have hcap := star_capacity hC
  have h1 : a * (2 * H.card) ≤ a * ∑ v : V n, (nbr H v).card := Nat.mul_le_mul_left a hdeg
  nlinarith [h1, hsum, hcap]

/-- **f(8) ≤ 691**, with no assumption. -/
theorem bound_8_691 : Bound 8 691 := by
  intro H hH hC
  have h := edge_bound (n := 8) 10 40 (by intro d hd; interval_cases d <;> decide) hH hC
  have h56 : Nat.choose 8 3 = 56 := by decide
  rw [h56] at h
  omega

/-- **f(7) ≤ 306**, with no assumption (the true value is 304, arXiv:2611.27284). -/
theorem bound_7_306 : Bound 7 306 := by
  intro H hH hC
  have h := edge_bound (n := 7) 6 20 (by intro d hd; interval_cases d <;> decide) hH hC
  have h35 : Nat.choose 7 3 = 35 := by decide
  rw [h35] at h
  omega

/-- **f(6) ≤ 133**, with no assumption (the true value is 132, Harborth-Nienborg 1994). -/
theorem bound_6_133 : Bound 6 133 := by
  intro H hH hC
  have h := edge_bound (n := 6) 6 20 (by intro d hd; interval_cases d <;> decide) hH hC
  have h20 : Nat.choose 6 3 = 20 := by decide
  rw [h20] at h
  omega

/-- **f(5) ≤ 56**, with no assumption; 56 is the true value. -/
theorem bound_5_56 : Bound 5 56 := by
  intro H hH hC
  have h := edge_bound (n := 5) 3 8 (by intro d hd; interval_cases d <;> decide) hH hC
  have h10 : Nat.choose 5 3 = 10 := by decide
  rw [h10] at h
  omega

/-- **f(4) ≤ 24**, with no assumption; 24 is the true value. -/
theorem bound_4_24 : Bound 4 24 := by
  intro H hH hC
  have h := edge_bound (n := 4) 1 2 (by intro d hd; interval_cases d <;> decide) hH hC
  have h4 : Nat.choose 4 3 = 4 := by decide
  rw [h4] at h
  omega

end Stars
