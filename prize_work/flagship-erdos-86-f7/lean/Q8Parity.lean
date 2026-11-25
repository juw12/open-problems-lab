/-
  prize_work/flagship-erdos-86-f7/lean/Q8Parity.lean

  **f(8) <= 690**: a C4-free subgraph of the 8-cube Q_8 has at most 690 edges
  (`Q8Parity.bound_8_690`).  Flagship lemma F5; the hand proof is ../f8_upper_bound.md.

  Sections 1-4 are copied word for word from Stars.lean (lemma F4: the cube, stars, one root per
  3-face and parity side, degree sums); the file stands alone because the lab's Lean gate checks
  single files.  Sections 5-8 are new.

  THE ARGUMENT.  Let H be C4-free in Q_8 with at least 691 edges.
    1. One parity side P (128 vertices) has degree sum >= 691.  Star capacity (`star_side`) says
       sum over P of C(deg v, 3) <= 1792.  Since C(d,3) >= 10 d - 40 with equality only for
       d = 5, 6 and a gap of at least 4 otherwise, every vertex of P has degree 5 or 6 and the
       sum is exactly 1790.
    2. For a 3-set T of directions, `cov H p T` is the set of vertices lying in a 3-face of
       direction T that has a root (a vertex of P with edges of H in all directions of T); it has
       8 * #roots elements.  `W` is the set of pairs (T, x) with x NOT covered: it has
       56 * 256 - 8 * 1790 = 16 elements (two root-less 3-faces of 8 vertices each).
    3. Facet parity.  Fix a facet {x : x c = b}.  The pairs of W with c not in T and x in the
       facet number 128 * (number of T) - 8 * (number of roots inside the facet, over all T that
       avoid c).  That number of roots is the sum over the P-vertices v of the facet of
       C(|nbr v minus c|, 3), and |nbr v minus c| is 4, 5 or 6, so every term is even
       (4, 10, 20).  Hence 16 divides the number of such pairs (`dvd_card_Wcb`).
    4. Take (T1, x1) in W.  For c outside T1 the facet {x : x c = x1 c} holds at least one pair,
       so by step 3 it holds all 16.  So every pair (T, x) of W has T = T1 and x = x1 outside
       T1: at most 2^3 = 8 pairs.  But W has 16.  Contradiction (`no_W`).

  What is and is not claimed.  This is a machine-checked proof that f(8) <= 690.  The best
  bound the lab found in print is 694 (averaging from f(7) = 304, arXiv:2611.27284), and the
  best lower bound is 682 (Minamoto, arXiv:2603.29127).  Whether 690 is new has NOT been
  established; that is ledger lemma F6.

  Written and checked by the Vector session of 2026-11-25
  (cd prize_work/lean ; python check_claim.py ../flagship-erdos-86-f7/lean/Q8Parity.lean).
-/

import Mathlib

open Finset
open scoped symmDiff

namespace Q8Parity

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

/-! ## 5. Sizes of a parity side and of a facet of `Q_8` -/

theorem flipSet_singleton {n : ℕ} (x : V n) (i : Fin n) : flipSet x {i} = flip x i := by
  funext k
  simp [flipSet, flip]

/-- Flipping one coordinate changes the parity. -/
theorem par_flip {n : ℕ} (x : V n) (i : Fin n) : par (flip x i) = par x + 1 := by
  rw [← flipSet_singleton, par_flipSet]
  simp

theorem flip_injective {n : ℕ} (i : Fin n) {v w : V n} (h : flip v i = flip w i) : v = w := by
  calc v = flip (flip v i) i := (flip_flip_self v i).symm
    _ = flip (flip w i) i := by rw [h]
    _ = w := flip_flip_self w i

theorem card_V8 : Fintype.card (V 8) = 256 := by
  simp

/-- Each parity side of `Q_8` has 128 vertices. -/
theorem card_side (p : ZMod 2) : (univ.filter fun v : V 8 => par v = p).card = 128 := by
  have hbij : ∀ q : ZMod 2, (univ.filter fun v : V 8 => par v = q).card
      ≤ (univ.filter fun v : V 8 => par v = q + 1).card := by
    intro q
    refine card_le_card_of_injOn (fun v => flip v 0) ?_ ?_
    · intro v hv
      have hv' : par v = q := by simpa using hv
      have h : par (flip v 0) = q + 1 := by rw [par_flip, hv']
      simpa using h
    · intro v _ w _ h
      exact flip_injective 0 h
  have h0 := hbij 0
  have h1 := hbij 1
  have h11 : (1 : ZMod 2) + 1 = 0 := by decide
  rw [zero_add] at h0
  rw [h11] at h1
  have hsum := card_filter_add_card_filter_not (s := (univ : Finset (V 8)))
    (fun v => par v = 0)
  have hneg : (univ.filter fun v : V 8 => ¬ par v = 0) = univ.filter fun v : V 8 => par v = 1 := by
    apply filter_congr
    intro v _
    have h : ∀ z : ZMod 2, ¬ z = 0 ↔ z = 1 := by decide
    exact h (par v)
  rw [hneg, card_univ, card_V8] at hsum
  have hp : p = 0 ∨ p = 1 := by
    have h : ∀ z : ZMod 2, z = 0 ∨ z = 1 := by decide
    exact h p
  rcases hp with rfl | rfl <;> omega

/-- Each facet `{x : x c = b}` of `Q_8` has 128 vertices. -/
theorem card_facet (c : Fin 8) (b : Bool) : (univ.filter fun x : V 8 => x c = b).card = 128 := by
  have hbij : ∀ a : Bool, (univ.filter fun x : V 8 => x c = a).card
      ≤ (univ.filter fun x : V 8 => x c = !a).card := by
    intro a
    refine card_le_card_of_injOn (fun v => flip v c) ?_ ?_
    · intro v hv
      have hv' : v c = a := by simpa using hv
      have h : flip v c c = !a := by rw [flip_self, hv']
      simpa using h
    · intro v _ w _ h
      exact flip_injective c h
  have ht := hbij true
  have hf := hbij false
  simp only [Bool.not_true, Bool.not_false] at ht hf
  have hsum := card_filter_add_card_filter_not (s := (univ : Finset (V 8)))
    (fun x => x c = true)
  have hneg : (univ.filter fun x : V 8 => ¬ x c = true) = univ.filter fun x : V 8 => x c = false := by
    apply filter_congr
    intro x _
    simp
  rw [hneg, card_univ, card_V8] at hsum
  cases b <;> omega

/-! ## 6. Covered vertices: the 3-faces that have a root -/

/-- The map (root, set of flips inside `T`) -> vertex is injective: one root per 3-face. -/
theorem flipSet_injOn {n : ℕ} {H : Finset (Finset (V n))} (hC : ¬ HasC4 H) (p : ZMod 2)
    {T : Finset (Fin n)} (hT : T.card = 3) :
    Set.InjOn (fun q : V n × Finset (Fin n) => flipSet q.1 q.2)
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

/-- The vertices that lie in a 3-face of direction `T` with a root of parity `p`. -/
def cov {n : ℕ} (H : Finset (Finset (V n))) (p : ZMod 2) (T : Finset (Fin n)) : Finset (V n) :=
  ((roots H p T) ×ˢ T.powerset).image fun q => flipSet q.1 q.2

/-- The roots of direction `T` inside the facet `{x : x c = b}`. -/
def rootsIn {n : ℕ} (H : Finset (Finset (V n))) (p : ZMod 2) (T : Finset (Fin n)) (c : Fin n)
    (b : Bool) : Finset (V n) :=
  (roots H p T).filter fun r => r c = b

theorem card_cov {n : ℕ} {H : Finset (Finset (V n))} (hC : ¬ HasC4 H) (p : ZMod 2)
    {T : Finset (Fin n)} (hT : T.card = 3) : (cov H p T).card = (roots H p T).card * 8 := by
  calc (cov H p T).card = ((roots H p T) ×ˢ T.powerset).card :=
        card_image_of_injOn (flipSet_injOn hC p hT)
    _ = (roots H p T).card * 2 ^ 3 := by rw [card_product, card_powerset, hT]
    _ = (roots H p T).card * 8 := by norm_num

/-- Inside a facet whose coordinate `c` is not a direction of `T`: 8 covered vertices per root
of the facet. -/
theorem card_cov_facet {n : ℕ} {H : Finset (Finset (V n))} (hC : ¬ HasC4 H) (p : ZMod 2)
    {T : Finset (Fin n)} (hT : T.card = 3) {c : Fin n} (hc : c ∉ T) (b : Bool) :
    ((cov H p T).filter fun x => x c = b).card = (rootsIn H p T c b).card * 8 := by
  have hset : (cov H p T).filter (fun x => x c = b)
      = ((rootsIn H p T c b) ×ˢ T.powerset).image fun q => flipSet q.1 q.2 := by
    ext x
    constructor
    · intro hx
      obtain ⟨hx1, hxc⟩ := mem_filter.mp hx
      have hx1' : x ∈ ((roots H p T) ×ˢ T.powerset).image (fun q => flipSet q.1 q.2) := hx1
      obtain ⟨⟨r, S⟩, hq, rfl⟩ := mem_image.mp hx1'
      obtain ⟨hr, hS⟩ := mem_product.mp hq
      have hcS : c ∉ S := fun h => hc (mem_powerset.mp hS h)
      have hrc : r c = b := by simpa [flipSet, hcS] using hxc
      exact mem_image.mpr ⟨(r, S), mem_product.mpr ⟨mem_filter.mpr ⟨hr, hrc⟩, hS⟩, rfl⟩
    · intro hx
      obtain ⟨⟨r, S⟩, hq, rfl⟩ := mem_image.mp hx
      obtain ⟨hr, hS⟩ := mem_product.mp hq
      obtain ⟨hr1, hrc⟩ := mem_filter.mp hr
      have hcS : c ∉ S := fun h => hc (mem_powerset.mp hS h)
      have hmem : flipSet r S ∈ cov H p T :=
        mem_image.mpr ⟨(r, S), mem_product.mpr ⟨hr1, hS⟩, rfl⟩
      refine mem_filter.mpr ⟨hmem, ?_⟩
      simpa [flipSet, hcS] using hrc
  have hsub : (↑((rootsIn H p T c b) ×ˢ T.powerset) : Set (V n × Finset (Fin n)))
      ⊆ ↑((roots H p T) ×ˢ T.powerset) := by
    intro q hq
    have hq' := mem_product.mp (mem_coe.mp hq)
    exact mem_coe.mpr (mem_product.mpr ⟨(mem_filter.mp hq'.1).1, hq'.2⟩)
  calc ((cov H p T).filter fun x => x c = b).card
      = (((rootsIn H p T c b) ×ˢ T.powerset).image fun q => flipSet q.1 q.2).card := by rw [hset]
    _ = ((rootsIn H p T c b) ×ˢ T.powerset).card :=
        card_image_of_injOn ((flipSet_injOn hC p hT).mono hsub)
    _ = (rootsIn H p T c b).card * 2 ^ 3 := by rw [card_product, card_powerset, hT]
    _ = (rootsIn H p T c b).card * 8 := by norm_num

/-! ## 7. Double counting stars and roots -/

/-- Stars on one side, counted by root and by direction. -/
theorem sum_choose_eq_sum_roots {n : ℕ} (H : Finset (Finset (V n))) (p : ZMod 2) :
    ∑ v ∈ univ.filter (fun v : V n => par v = p), (nbr H v).card.choose 3
      = ∑ T ∈ powersetCard 3 (univ : Finset (Fin n)), (roots H p T).card := by
  have h1 : ∀ v : V n, (nbr H v).card.choose 3
      = ((powersetCard 3 (univ : Finset (Fin n))).filter (fun T => T ⊆ nbr H v)).card := by
    intro v
    rw [← card_powersetCard 3 (nbr H v)]
    congr 1
    ext T
    simp only [mem_powersetCard, mem_filter, subset_univ, true_and]
    exact and_comm
  simp_rw [h1, card_filter]
  rw [sum_comm]
  refine sum_congr rfl fun T _ => ?_
  rw [← card_filter, filter_filter]
  rfl

/-- The same inside a facet `{x : x c = b}`, for the directions that avoid `c`. -/
theorem sum_roots_facet {n : ℕ} (H : Finset (Finset (V n))) (p : ZMod 2) (c : Fin n) (b : Bool) :
    ∑ v ∈ univ.filter (fun v : V n => par v = p ∧ v c = b), ((nbr H v).erase c).card.choose 3
      = ∑ T ∈ (powersetCard 3 (univ : Finset (Fin n))).filter (fun T => c ∉ T),
          (rootsIn H p T c b).card := by
  have h1 : ∀ v : V n, ((nbr H v).erase c).card.choose 3
      = (((powersetCard 3 (univ : Finset (Fin n))).filter (fun T => c ∉ T)).filter
          (fun T => T ⊆ nbr H v)).card := by
    intro v
    rw [← card_powersetCard 3 ((nbr H v).erase c)]
    congr 1
    ext T
    simp only [mem_powersetCard, mem_filter, subset_univ, true_and, subset_erase]
    tauto
  simp_rw [h1, card_filter]
  rw [sum_comm]
  refine sum_congr rfl fun T _ => ?_
  rw [← card_filter]
  congr 1
  ext v
  simp only [rootsIn, roots, mem_filter, mem_univ, true_and]
  tauto

/-- Counting a filtered product row by row. -/
theorem card_filter_product {α β : Type*} (A : Finset α) (B : Finset β) (P : α → β → Prop)
    [∀ a b, Decidable (P a b)] :
    ((A ×ˢ B).filter fun q => P q.1 q.2).card = ∑ a ∈ A, (B.filter (P a)).card := by
  rw [card_filter, sum_product]
  refine sum_congr rfl fun a _ => ?_
  rw [card_filter]

/-! ## 8. The pairs (direction, vertex) that are not covered -/

/-- The pairs `(T, x)`: `T` a 3-set of directions, `x` a vertex in a root-less 3-face of
direction `T`. -/
def W (H : Finset (Finset (V 8))) (p : ZMod 2) : Finset (Finset (Fin 8) × V 8) :=
  ((powersetCard 3 (univ : Finset (Fin 8))) ×ˢ (univ : Finset (V 8))).filter
    fun q => q.2 ∉ cov H p q.1

/-- The pairs of `W` whose direction avoids `c` and whose vertex lies in the facet `x c = b`. -/
def Wcb (H : Finset (Finset (V 8))) (p : ZMod 2) (c : Fin 8) (b : Bool) :
    Finset (Finset (Fin 8) × V 8) :=
  (W H p).filter fun q => c ∉ q.1 ∧ q.2 c = b

theorem mem_W {H : Finset (Finset (V 8))} {p : ZMod 2} {T : Finset (Fin 8)} {x : V 8} :
    (T, x) ∈ W H p ↔ T.card = 3 ∧ x ∉ cov H p T := by
  simp only [W, mem_filter, mem_product, mem_powersetCard, subset_univ, true_and, mem_univ,
    and_true]

theorem mem_Wcb {H : Finset (Finset (V 8))} {p : ZMod 2} {c : Fin 8} {b : Bool}
    {T : Finset (Fin 8)} {x : V 8} :
    (T, x) ∈ Wcb H p c b ↔ (T, x) ∈ W H p ∧ c ∉ T ∧ x c = b := by
  simp only [Wcb, mem_filter]

theorem Wcb_subset (H : Finset (Finset (V 8))) (p : ZMod 2) (c : Fin 8) (b : Bool) :
    Wcb H p c b ⊆ W H p := by
  unfold Wcb
  exact filter_subset _ _

theorem card_uncovered {H : Finset (Finset (V 8))} (hC : ¬ HasC4 H) (p : ZMod 2)
    {T : Finset (Fin 8)} (hT : T.card = 3) :
    (univ.filter fun x : V 8 => x ∉ cov H p T).card + (roots H p T).card * 8 = 256 := by
  have h := card_filter_add_card_filter_not (s := (univ : Finset (V 8)))
    (fun x => x ∈ cov H p T)
  rw [filter_mem_eq_inter, univ_inter, card_cov hC p hT, card_univ, card_V8] at h
  omega

/-- If the stars on side `p` number 1790, exactly 16 pairs are not covered. -/
theorem card_W {H : Finset (Finset (V 8))} (hC : ¬ HasC4 H) (p : ZMod 2)
    (hR : ∑ T ∈ powersetCard 3 (univ : Finset (Fin 8)), (roots H p T).card = 1790) :
    (W H p).card = 16 := by
  have h0 : (W H p).card = ∑ T ∈ powersetCard 3 (univ : Finset (Fin 8)),
      (univ.filter fun x : V 8 => x ∉ cov H p T).card :=
    card_filter_product _ _ (fun T x => x ∉ cov H p T)
  have h : ∑ T ∈ powersetCard 3 (univ : Finset (Fin 8)),
      ((univ.filter fun x : V 8 => x ∉ cov H p T).card + (roots H p T).card * 8)
      = ∑ _T ∈ powersetCard 3 (univ : Finset (Fin 8)), 256 :=
    sum_congr rfl fun T hT => card_uncovered hC p (mem_powersetCard.mp hT).2
  rw [sum_add_distrib, ← sum_mul, hR, sum_const, card_powersetCard, card_univ, Fintype.card_fin,
    smul_eq_mul] at h
  have h56 : Nat.choose 8 3 = 56 := by decide
  rw [h56] at h
  omega

/-- **Facet parity.** If all degrees on side `p` are 5 or 6, the number of uncovered pairs in
any facet is a multiple of 16. -/
theorem dvd_card_Wcb {H : Finset (Finset (V 8))} (hC : ¬ HasC4 H) (p : ZMod 2)
    (h56 : ∀ v : V 8, par v = p → (nbr H v).card = 5 ∨ (nbr H v).card = 6) (c : Fin 8)
    (b : Bool) : 16 ∣ (Wcb H p c b).card := by
  have hset : Wcb H p c b
      = (((powersetCard 3 (univ : Finset (Fin 8))).filter fun T => c ∉ T)
          ×ˢ (univ.filter fun x : V 8 => x c = b)).filter (fun q => q.2 ∉ cov H p q.1) := by
    ext ⟨T, x⟩
    simp only [Wcb, W, mem_filter, mem_product, mem_univ, and_true, true_and]
    tauto
  have h0 : (Wcb H p c b).card
      = ∑ T ∈ (powersetCard 3 (univ : Finset (Fin 8))).filter (fun T => c ∉ T),
          ((univ.filter fun x : V 8 => x c = b).filter fun x => x ∉ cov H p T).card := by
    rw [hset]
    exact card_filter_product _ _ (fun T x => x ∉ cov H p T)
  have hT : ∀ T ∈ (powersetCard 3 (univ : Finset (Fin 8))).filter (fun T => c ∉ T),
      ((univ.filter fun x : V 8 => x c = b).filter fun x => x ∉ cov H p T).card
        + (rootsIn H p T c b).card * 8 = 128 := by
    intro T hT
    obtain ⟨hT3, hcT⟩ := mem_filter.mp hT
    have hT3' : T.card = 3 := (mem_powersetCard.mp hT3).2
    have h := card_filter_add_card_filter_not
      (s := (univ.filter fun x : V 8 => x c = b)) (fun x => x ∈ cov H p T)
    rw [card_facet c b] at h
    have h2 : ((univ.filter fun x : V 8 => x c = b).filter fun x => x ∈ cov H p T).card
        = (rootsIn H p T c b).card * 8 := by
      have hs : (univ.filter fun x : V 8 => x c = b).filter (fun x => x ∈ cov H p T)
          = (cov H p T).filter fun x => x c = b := by
        ext x
        simp only [mem_filter, mem_univ, true_and]
        exact and_comm
      rw [hs]
      exact card_cov_facet hC p hT3' hcT b
    omega
  have hsum : ∑ T ∈ (powersetCard 3 (univ : Finset (Fin 8))).filter (fun T => c ∉ T),
      (((univ.filter fun x : V 8 => x c = b).filter fun x => x ∉ cov H p T).card
        + (rootsIn H p T c b).card * 8)
      = ∑ _T ∈ (powersetCard 3 (univ : Finset (Fin 8))).filter (fun T => c ∉ T), 128 :=
    sum_congr rfl hT
  rw [sum_add_distrib, ← sum_mul, sum_const, smul_eq_mul, ← h0] at hsum
  have heven : 2 ∣ ∑ T ∈ (powersetCard 3 (univ : Finset (Fin 8))).filter (fun T => c ∉ T),
      (rootsIn H p T c b).card := by
    rw [← sum_roots_facet H p c b]
    refine dvd_sum fun v hv => ?_
    have hv' := (mem_filter.mp hv).2
    have hd := h56 v hv'.1
    have he : ((nbr H v).erase c).card = (nbr H v).card
        ∨ ((nbr H v).erase c).card + 1 = (nbr H v).card := by
      by_cases hcv : c ∈ nbr H v
      · right
        rw [card_erase_of_mem hcv]
        have hpos := card_pos.mpr ⟨c, hcv⟩
        omega
      · left
        rw [erase_eq_self.mpr hcv]
    have h456 : ((nbr H v).erase c).card = 4 ∨ ((nbr H v).erase c).card = 5
        ∨ ((nbr H v).erase c).card = 6 := by omega
    rcases h456 with h | h | h <;> rw [h] <;> decide
  obtain ⟨k, hk⟩ := heven
  omega

-- `W` is a closed term over `Q_8`; keep the elaborator from unfolding it into 14336 pairs.
attribute [irreducible] W Wcb

/-- 16 uncovered pairs cannot sit in every facet in multiples of 16. -/
theorem no_W {H : Finset (Finset (V 8))} (p : ZMod 2) (hW : (W H p).card = 16)
    (hdvd : ∀ c b, 16 ∣ (Wcb H p c b).card) : False := by
  have hne : (W H p).Nonempty := card_pos.mp (by omega)
  obtain ⟨⟨T1, x1⟩, hq1⟩ := hne
  have hT1 : T1.card = 3 := (mem_W.mp hq1).1
  have hall : ∀ c, c ∉ T1 → ∀ T x, (T, x) ∈ W H p → c ∉ T ∧ x c = x1 c := by
    intro c hc
    have hsub : Wcb H p c (x1 c) ⊆ W H p := Wcb_subset H p c (x1 c)
    have hmem : (T1, x1) ∈ Wcb H p c (x1 c) := mem_Wcb.mpr ⟨hq1, hc, rfl⟩
    have hpos : 0 < (Wcb H p c (x1 c)).card := card_pos.mpr ⟨_, hmem⟩
    have hle := card_le_card hsub
    have h16 : (Wcb H p c (x1 c)).card = 16 := by
      obtain ⟨k, hk⟩ := hdvd c (x1 c)
      omega
    have heq : Wcb H p c (x1 c) = W H p := eq_of_subset_of_card_le hsub (by omega)
    intro T x hq
    rw [← heq] at hq
    exact (mem_Wcb.mp hq).2
  have hTq : ∀ T x, (T, x) ∈ W H p → T = T1 := by
    intro T x hq
    have hq3 : T.card = 3 := (mem_W.mp hq).1
    apply eq_of_subset_of_card_le
    · intro k hk
      by_contra hkT
      exact (hall k hkT T x hq).1 hk
    · omega
  have hinj : Set.InjOn (fun q : Finset (Fin 8) × V 8 => T1.filter fun k => q.2 k = true)
      (↑(W H p) : Set (Finset (Fin 8) × V 8)) := by
    rintro ⟨T, x⟩ hq ⟨T', x'⟩ hq' h
    have hqW : (T, x) ∈ W H p := hq
    have hqW' : (T', x') ∈ W H p := hq'
    have h' : T1.filter (fun k => x k = true) = T1.filter (fun k => x' k = true) := h
    have e1 : T = T' := (hTq T x hqW).trans (hTq T' x' hqW').symm
    have e2 : x = x' := by
      funext k
      by_cases hk : k ∈ T1
      · have hmem : k ∈ T1.filter (fun k => x k = true)
            ↔ k ∈ T1.filter (fun k => x' k = true) := by rw [h']
        have hiff : x k = true ↔ x' k = true := by simpa [hk] using hmem
        exact Bool.eq_iff_iff.mpr hiff
      · rw [(hall k hk T x hqW).2, (hall k hk T' x' hqW').2]
    rw [e1, e2]
  have hcard := card_le_card_of_injOn (t := T1.powerset)
    (fun q : Finset (Fin 8) × V 8 => T1.filter fun k => q.2 k = true)
    (fun q _ => mem_powerset.mpr (filter_subset _ _)) hinj
  rw [card_powerset, hT1, hW] at hcard
  norm_num at hcard

/-! ## 9. The bound -/

/-- **f(8) ≤ 690**: every C4-free set of edges of `Q_8` has at most 690 edges. -/
theorem bound_8_690 : Bound 8 690 := by
  intro H hH hC
  by_contra hcon
  have hcard : 691 ≤ H.card := by omega
  -- one parity side has degree sum at least 691
  have hdeg := two_mul_card_le hH
  rw [← Finset.sum_fiberwise (univ : Finset (V 8)) par (fun v => (nbr H v).card)] at hdeg
  have huniv : (univ : Finset (ZMod 2)) = {0, 1} := by decide
  rw [huniv, sum_pair (by decide)] at hdeg
  have hp : ∃ p : ZMod 2,
      691 ≤ ∑ v ∈ univ.filter (fun v : V 8 => par v = p), (nbr H v).card := by
    by_contra hno
    have h0 : ¬ 691 ≤ ∑ v ∈ univ.filter (fun v : V 8 => par v = 0), (nbr H v).card :=
      fun h => hno ⟨0, h⟩
    have h1 : ¬ 691 ≤ ∑ v ∈ univ.filter (fun v : V 8 => par v = 1), (nbr H v).card :=
      fun h => hno ⟨1, h⟩
    omega
  obtain ⟨p, hp⟩ := hp
  have hPcard := card_side p
  have hside := star_side hC p
  have h56c : Nat.choose 8 3 = 56 := by decide
  rw [h56c] at hside
  have hd8 : ∀ v : V 8, (nbr H v).card ≤ 8 := fun v => by
    simpa using card_le_univ (nbr H v)
  have hw : ∀ d : ℕ, d ≤ 8 → 10 * d ≤ d.choose 3 + 40 := by
    intro d hd
    interval_cases d <;> decide
  have hw4 : ∀ d : ℕ, d ≤ 8 → d = 5 ∨ d = 6 ∨ 10 * d + 4 ≤ d.choose 3 + 40 := by
    intro d hd
    interval_cases d <;> decide
  have e3 : ∑ v ∈ univ.filter (fun v : V 8 => par v = p), 10 * (nbr H v).card
      = 10 * ∑ v ∈ univ.filter (fun v : V 8 => par v = p), (nbr H v).card := by
    rw [mul_sum]
  have e4 : ∑ v ∈ univ.filter (fun v : V 8 => par v = p), ((nbr H v).card.choose 3 + 40)
      = ∑ v ∈ univ.filter (fun v : V 8 => par v = p), (nbr H v).card.choose 3 + 128 * 40 := by
    rw [sum_add_distrib, sum_const, hPcard, smul_eq_mul]
  -- every vertex of the side has degree 5 or 6
  have h56deg : ∀ v ∈ univ.filter (fun v : V 8 => par v = p),
      (nbr H v).card = 5 ∨ (nbr H v).card = 6 := by
    intro v0 hv0
    by_contra hbad
    have h4 : 10 * (nbr H v0).card + 4 ≤ (nbr H v0).card.choose 3 + 40 := by
      rcases hw4 _ (hd8 v0) with h | h | h
      · exact absurd (Or.inl h) hbad
      · exact absurd (Or.inr h) hbad
      · exact h
    have hrest : ∑ v ∈ (univ.filter (fun v : V 8 => par v = p)).erase v0, 10 * (nbr H v).card
        ≤ ∑ v ∈ (univ.filter (fun v : V 8 => par v = p)).erase v0,
            ((nbr H v).card.choose 3 + 40) :=
      sum_le_sum fun v _ => hw _ (hd8 v)
    have e1 := add_sum_erase (univ.filter (fun v : V 8 => par v = p))
      (fun v => 10 * (nbr H v).card) hv0
    have e2 := add_sum_erase (univ.filter (fun v : V 8 => par v = p))
      (fun v => (nbr H v).card.choose 3 + 40) hv0
    beta_reduce at e1 e2
    omega
  have hCeq : ∀ v ∈ univ.filter (fun v : V 8 => par v = p),
      (nbr H v).card.choose 3 + 40 = 10 * (nbr H v).card := by
    intro v hv
    rcases h56deg v hv with h | h <;> rw [h] <;> decide
  have hS : ∑ v ∈ univ.filter (fun v : V 8 => par v = p), (nbr H v).card.choose 3 = 1790 := by
    have e : ∑ v ∈ univ.filter (fun v : V 8 => par v = p), ((nbr H v).card.choose 3 + 40)
        = ∑ v ∈ univ.filter (fun v : V 8 => par v = p), 10 * (nbr H v).card :=
      sum_congr rfl hCeq
    omega
  have hR : ∑ T ∈ powersetCard 3 (univ : Finset (Fin 8)), (roots H p T).card = 1790 := by
    rw [← sum_choose_eq_sum_roots H p]
    exact hS
  have h56all : ∀ v : V 8, par v = p → (nbr H v).card = 5 ∨ (nbr H v).card = 6 :=
    fun v hv => h56deg v (mem_filter.mpr ⟨mem_univ v, hv⟩)
  exact no_W p (card_W hC p hR) (fun c b => dvd_card_Wcb hC p h56all c b)

end Q8Parity
