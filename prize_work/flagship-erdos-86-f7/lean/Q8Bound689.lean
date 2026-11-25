import Mathlib

open Finset
open scoped symmDiff

namespace Q8Bound689

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

/-! ## 9. The general facet-parity rule, and the e = 690 slack classification (F9) -/

/-- `C(d,3)` is odd exactly at `d = 3` and `d = 7` for `d ≤ 8`: the parity of a star count
sees only those two local degrees.  This is why the F5 rule ("even in every facet", valid
when all local degrees are 4, 5 or 6) has a general form. -/
theorem choose3_mod_two (d : ℕ) (hd : d ≤ 8) :
    ((d.choose 3 : ℕ) : ZMod 2) = (if d ∈ ({3, 7} : Finset ℕ) then 1 else 0) := by
  interval_cases d <;> decide

/-- Casting a sum of naturals into `ZMod 2`, term by term. -/
theorem nat_sum_cast {α : Type*} (s : Finset α) (f : α → ℕ) :
    ((∑ x ∈ s, f x : ℕ) : ZMod 2) = ∑ x ∈ s, ((f x : ℕ) : ZMod 2) := by
  push_cast
  rfl

/-- The local degree in a facet is at most `n`. -/
theorem erase_le {n : ℕ} (H : Finset (Finset (V n))) (v : V n) (c : Fin n) :
    ((nbr H v).erase c).card ≤ n :=
  calc ((nbr H v).erase c).card ≤ (nbr H v).card := card_le_card (erase_subset _ _)
    _ ≤ Fintype.card (Fin n) := card_le_univ _
    _ = n := by simp

/-- **The general facet-parity rule (F9, the rule the job names).**  Let `H` be C4-free in
`Q_8`, `p` a parity side and `K = {x : x c = b}` a facet.  The 3-faces of `K` that carry a
`p`-root are counted by `sum over v in K of C(l_v, 3)`, where `l_v` is the number of directions
present at `v` inside `K`; modulo 2 that count is the number of `p`-vertices of `K` whose local
degree is 3 or 7.  Hence (a facet of `Q_8` has exactly 560 three-faces)

    #(root-less 3-faces of K)  ==  #{v in K : l_v = 3 or 7}   (mod 2).

`dvd_card_Wcb` is the special case in which every local degree is 4, 5 or 6, so that the
right-hand side vanishes. -/
theorem parity_facet {H : Finset (Finset (V 8))} (_hC : ¬ HasC4 H) (p : ZMod 2) (c : Fin 8)
    (b : Bool) :
    ((∑ T ∈ (powersetCard 3 (univ : Finset (Fin 8))).filter (fun T => c ∉ T),
        (rootsIn H p T c b).card : ℕ) : ZMod 2)
      = ((univ.filter fun v : V 8 => par v = p ∧ v c = b
          ∧ ((nbr H v).erase c).card ∈ ({3, 7} : Finset ℕ)).card : ZMod 2) := by
  have h := sum_roots_facet (n := 8) H p c b
  have hle : ∀ v ∈ (univ.filter fun v : V 8 => par v = p ∧ v c = b),
      ((nbr H v).erase c).card ≤ 8 := fun v _ => by
    simpa using erase_le H v c
  have hset : (univ.filter fun v : V 8 => par v = p ∧ v c = b
        ∧ ((nbr H v).erase c).card ∈ ({3, 7} : Finset ℕ))
      = (univ.filter fun v : V 8 => par v = p ∧ v c = b).filter
        (fun v => ((nbr H v).erase c).card ∈ ({3, 7} : Finset ℕ)) := by
    ext v
    simp only [mem_filter, mem_univ, true_and]
    tauto
  calc ((∑ T ∈ (powersetCard 3 (univ : Finset (Fin 8))).filter (fun T => c ∉ T),
        (rootsIn H p T c b).card : ℕ) : ZMod 2)
      = ((∑ v ∈ (univ.filter fun v : V 8 => par v = p ∧ v c = b),
          ((nbr H v).erase c).card.choose 3 : ℕ) : ZMod 2) := by rw [h]
    _ = ∑ v ∈ (univ.filter fun v : V 8 => par v = p ∧ v c = b),
          (((nbr H v).erase c).card.choose 3 : ZMod 2) := nat_sum_cast _ _
    _ = ∑ v ∈ (univ.filter fun v : V 8 => par v = p ∧ v c = b),
          (if ((nbr H v).erase c).card ∈ ({3, 7} : Finset ℕ) then 1 else 0) :=
        sum_congr rfl fun v hv => choose3_mod_two _ (hle v hv)
    _ = ((univ.filter fun v : V 8 => par v = p ∧ v c = b
          ∧ ((nbr H v).erase c).card ∈ ({3, 7} : Finset ℕ)).card : ZMod 2) := by
        rw [hset, ← sum_boole]

/-- **At 691 the general rule is F5's rule.**  If every vertex of the side `p` has degree 5 or
6, every local degree in a facet is 4, 5 or 6, so no `p`-vertex of a facet has local degree 3
or 7 and the star count of every facet is even. -/
theorem special_empty_of_deg {H : Finset (Finset (V 8))} (p : ZMod 2)
    (h56 : ∀ v : V 8, par v = p → (nbr H v).card = 5 ∨ (nbr H v).card = 6) (c : Fin 8)
    (b : Bool) :
    (univ.filter fun v : V 8 => par v = p ∧ v c = b
        ∧ ((nbr H v).erase c).card ∈ ({3, 7} : Finset ℕ)).card = 0 := by
  rw [card_eq_zero, filter_eq_empty_iff]
  intro v _
  rintro ⟨hp, _, hsp⟩
  have hd := h56 v hp
  have hloc : ((nbr H v).erase c).card = (nbr H v).card
      ∨ ((nbr H v).erase c).card + 1 = (nbr H v).card := by
    by_cases hcv : c ∈ nbr H v
    · right
      rw [card_erase_of_mem hcv]
      have : 0 < (nbr H v).card := card_pos.mpr ⟨c, hcv⟩
      omega
    · left
      rw [erase_eq_self.mpr hcv]
  rcases hd with h | h <;> rw [h] at hloc <;> simp only [mem_insert, mem_singleton] at hsp <;>
    omega

/-! ## 10. The e = 690 slack classification -/

/-- Pointwise, `C(d,3) + 40 - 10 d ≥ 0` for `d ≤ 8`; equality holds at `d = 5` and `d = 6`. -/
theorem w_nonneg (d : ℕ) (hd : d ≤ 8) : 0 ≤ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 := by
  interval_cases d <;> decide

/-- Pointwise, `C(d,3) + 40 - 10 d ≥ 4` away from `d = 5, 6` (`d ≤ 8`). -/
theorem w_ge_four (d : ℕ) (hd : d ≤ 8) (h5 : d ≠ 5) (h6 : d ≠ 6) :
    4 ≤ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 := by
  interval_cases d <;> simp_all [Nat.choose]

/-- The exact values of the slack weight, for the case list. -/
theorem w_values : ∀ d ∈ ({0, 1, 2, 3, 4, 5, 6, 7, 8} : Finset ℕ),
    ((d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 = 40 ∨ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 = 30
      ∨ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 = 20
      ∨ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 = 11
      ∨ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 = 4
      ∨ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 = 0
      ∨ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 = 5
      ∨ (d.choose 3 : ℤ) - 10 * (d : ℤ) + 40 = 16) := by
  intro d hd
  simp only [mem_insert, mem_singleton] at hd
  rcases hd with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl <;> decide

/-- **The slack arithmetic at `e = 690`.**  On one parity side of `Q_8` (128 vertices) with
degree sum 690 and star sum at most 1792, the total slack `sum_v (C(d_v,3) - 10 d_v + 40)` is at
most 12, and at most three vertices of the side have degree other than 5 or 6.  Since a
non-zero slack is at least 4, this is exactly the F9 case list: no special vertex and 12
root-less faces, or one, two or three of degree 4, or one or two of degree 7, or one of each,
or one of degree 3. -/
theorem slack_690 {s : Finset ℕ} (d : ℕ → ℕ) (hcard : s.card = 128)
    (hsum : ∑ v ∈ s, (d v : ℤ) = 690) (hch : ∑ v ∈ s, ((d v).choose 3 : ℤ) ≤ 1792)
    (hd : ∀ v ∈ s, d v ≤ 8) :
    (s.filter fun v => d v ≠ 5 ∧ d v ≠ 6).card ≤ 3 ∧
      ∑ v ∈ s, ((d v).choose 3 - 10 * (d v : ℤ) + 40) ≤ 12 := by
  have h1 : ∑ v ∈ s, (10 * (d v : ℤ)) = 10 * ∑ v ∈ s, (d v : ℤ) := by rw [Finset.mul_sum]
  have h2 : ∑ _v ∈ s, (40 : ℤ) = 40 * (s.card : ℤ) := by
    rw [Finset.sum_const, nsmul_eq_mul]; ring
  have hsplit : ∑ v ∈ s, ((d v).choose 3 - 10 * (d v : ℤ) + 40)
      = (∑ v ∈ s, ((d v).choose 3 : ℤ)) - 10 * (∑ v ∈ s, (d v : ℤ)) + 40 * (s.card : ℤ) := by
    rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, h1, h2]
  constructor
  · set E := s.filter fun v => d v ≠ 5 ∧ d v ≠ 6 with hE
    have hsub : E ⊆ s := filter_subset _ _
    have h4 : 4 * (E.card : ℤ) ≤ ∑ v ∈ E, ((d v).choose 3 - 10 * (d v : ℤ) + 40) := by
      have hpt : ∀ v ∈ E, (4 : ℤ) ≤ (d v).choose 3 - 10 * (d v : ℤ) + 40 := by
        intro v hv
        obtain ⟨hv5, hv6⟩ := (mem_filter.mp hv).2
        exact w_ge_four (d v) (hd v (hsub hv)) hv5 hv6
      calc 4 * (E.card : ℤ) = ∑ _v ∈ E, (4 : ℤ) := by
            rw [Finset.sum_const, nsmul_eq_mul]; ring
        _ ≤ _ := sum_le_sum hpt
    have hnonneg : ∀ v ∈ s, (0 : ℤ) ≤ (d v).choose 3 - 10 * (d v : ℤ) + 40 :=
      fun v hv => w_nonneg (d v) (hd v hv)
    have hle : ∑ v ∈ E, ((d v).choose 3 - 10 * (d v : ℤ) + 40)
        ≤ ∑ v ∈ s, ((d v).choose 3 - 10 * (d v : ℤ) + 40) :=
      sum_le_sum_of_subset_of_nonneg hsub (fun v hv _ => hnonneg v hv)
    have h12 : ∑ v ∈ s, ((d v).choose 3 - 10 * (d v : ℤ) + 40) ≤ 12 := by
      rw [hsplit, hsum, hcard]; push_cast; linarith
    have hstep : (4 : ℤ) * (E.card : ℤ) ≤ 12 := le_trans h4 (le_trans hle h12)
    have h3 : (E.card : ℤ) ≤ 3 := by linarith
    omega
  · rw [hsplit, hsum, hcard]
    push_cast
    linarith

/-! ## 11. F9: `Bound 8 689` from the sub-5-cube inequality of Thomason-Wagner

Thomason-Wagner (Discrete Math. 309 (2009), Table 1 on p. 1732) state that the star vectors
`s_r = ∑_v C(d_v,r)` of the C4-free subgraphs of `Q_5` all satisfy

    12 * s1 - 5 * s2 + 11 * s4 - 20 * s5 ≤ 800                                (*)

(their table's 33 extreme points; equality at the 56-edge point `(112,144,80,16,0)`).
Since `s_r = ∑_v C(d_v,r)`, (*) is the single-vertex inequality `∑_v g(d_v) ≤ 800` with
`g d = 12 d - 5 C(d,2) + 11 C(d,4) - 20 C(d,5)`.

`Q_8` has `C(8,5) * 2^3 = 448` sub-5-cubes, and the 5-face of direction set `J` through `v`
contains `C(8-r,5-r)` of the `r`-subsets of the directions at `v`; so summing (*) over all
448 faces counts `h(d) = ∑_J g(|S_d ∩ J|) = 420 d - 100 C(d,2) + 44 C(d,4) - 20 C(d,5)`
at a vertex of degree `d`.  `h(d) - 260 d = 0,160,220,180,84,0,0,140,440 ≥ 0` for
`d = 0..8`, so `520 e ≤ 800 * 448 = 358400` and `e ≤ 689`.

Nothing below assumes the table is right: `bound_8_689_of_q5` takes (*) — face by face —
as its hypothesis, and proves `Bound 8 689` from it.  Section 12 then replaces the two
hypotheses by one statement about `Q_5` alone (`Q5Ineq`), so that

    theorem bound_8_689_of_Q5Ineq (h : Q5Ineq) : Bound 8 689

is all that is left.  `Q5Ineq` itself is *not* proved in Lean; it is checked by exact solver
runs (`../sources/q5_runs/`, a solver verdict) and stated in the Thomason-Wagner table. -/

/-- The 5-face of `Q_8` with direction set `J` through `x`: the vertices that agree with `x`
on every coordinate outside `J`. -/
def face5 (J : Finset (Fin 8)) (x : V 8) : Finset (V 8) :=
  univ.filter fun v => ∀ c, c ∉ J → v c = x c

/-- The degree of `v` inside its 5-face of direction set `J`: the present directions of `H`
at `v` that lie in `J`. -/
def dj (H : Finset (Finset (V 8))) (J : Finset (Fin 8)) (v : V 8) : ℕ := (nbr H v ∩ J).card

/-- The base point of the 5-face of `v`: `v` with the coordinates in `J` set to `false`. -/
def base (J : Finset (Fin 8)) (v : V 8) : V 8 := fun c => if c ∈ J then false else v c

/-- A face through a base point is a fibre of `base`. -/
theorem face5_eq (J : Finset (Fin 8)) {x : V 8} (hx : ∀ c ∈ supp x, c ∉ J) :
    face5 J x = univ.filter fun v => base J v = x := by
  ext v
  simp only [face5, mem_filter, mem_univ, true_and]
  constructor
  · intro h
    funext c
    by_cases hc : c ∈ J
    · have hxc : x c = false := by
        cases hb : x c
        · rfl
        · exact absurd hc (hx c ((mem_supp x c).mpr hb))
      show (if c ∈ J then false else v c) = x c
      rw [if_pos hc, hxc]
    · show (if c ∈ J then false else v c) = x c
      rw [if_neg hc]
      exact h c hc
  · intro h c hc
    have hcv : base J v c = x c := congr_fun h c
    have hb : base J v c = v c := by
      show (if c ∈ J then false else v c) = v c
      rw [if_neg hc]
    rw [← hb]
    exact hcv

/-- There are at most `2^3 = 8` base points for a direction set of size 5. -/
theorem card_X_le (J : Finset (Fin 8)) (hJ : J.card = 5) :
    (univ.filter fun x : V 8 => ∀ c ∈ supp x, c ∉ J).card ≤ 8 := by
  have hsub : (univ.filter fun x : V 8 => ∀ c ∈ supp x, c ∉ J).image supp ⊆ (Jᶜ).powerset := by
    intro S hS
    obtain ⟨x, hx, rfl⟩ := mem_image.mp hS
    rw [mem_powerset]
    intro c hc
    rw [mem_compl]
    exact (mem_filter.mp hx).2 c hc
  have h := card_le_card hsub
  rw [card_image_of_injective _ supp_injective, card_powerset, card_compl, Fintype.card_fin,
    hJ] at h
  exact le_trans h (by norm_num)

/-- Every vertex lies in exactly one 5-face of a given direction set: summing over the base
points and then over each face is summing over all vertices. -/
theorem sum_over_faces (J : Finset (Fin 8)) (F : V 8 → ℤ) :
    (∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J), ∑ v ∈ face5 J x, F v)
      = ∑ v : V 8, F v := by
  have hmaps : ∀ v ∈ (univ : Finset (V 8)),
      base J v ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J) := by
    intro v _
    rw [mem_filter]
    refine ⟨mem_univ _, ?_⟩
    intro c hc hcJ
    have h1 : base J v c = true := (mem_supp _ c).mp hc
    have h2 : base J v c = false := by
      show (if c ∈ J then false else v c) = false
      rw [if_pos hcJ]
    rw [h2] at h1
    exact Bool.noConfusion h1
  calc (∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J), ∑ v ∈ face5 J x, F v)
      = ∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J),
          ∑ v ∈ (univ : Finset (V 8)).filter (fun v => base J v = x), F v :=
        sum_congr rfl fun x hx => by rw [face5_eq J (mem_filter.mp hx).2]
    _ = ∑ v : V 8, F v := Finset.sum_fiberwise_of_maps_to hmaps F

/-- The Thomason-Wagner weight `g(d) = 12 d - 5 C(d,2) + 11 C(d,4) - 20 C(d,5)`; for
`d = 0..8` it takes the values `0, 12, 19, 21, 29, 45, 42, -56, -394`.  A C4-free subgraph of
`Q_5` with degree sequence `d` satisfies `12 s1 - 5 s2 + 11 s4 - 20 s5 = ∑_v g(d_v)`. -/
def g (d : ℕ) : ℤ :=
  12 * (d : ℤ) - 5 * ((d.choose 2 : ℕ) : ℤ) + 11 * ((d.choose 4 : ℕ) : ℤ) -
    20 * ((d.choose 5 : ℕ) : ℤ)

/-- **F9, from the sub-5-cube inequality.**  If every C4-free subgraph lying in a 5-face
`face5 J x` of `Q_8` (`|J| = 5`) satisfies Thomason-Wagner's `12 s1 - 5 s2 + 11 s4 - 20 s5
≤ 800` — written here in the equivalent single-vertex form `∑_v g(d_v) ≤ 800`, `d_v` the
local degree `dj H J v` — then every C4-free subgraph of `Q_8` has at most 689 edges. -/
theorem bound_8_689_of_q5
    (h5 : ∀ (H : Finset (Finset (V 8))) (J : Finset (Fin 8)) (x : V 8), J.card = 5 →
      (∀ c ∈ supp x, c ∉ J) → ¬ HasC4 H →
      (∑ v ∈ face5 J x, g (dj H J v)) ≤ 800)
    (hvert : ∀ S : Finset (Fin 8),
      (260 : ℤ) * ((S.card : ℕ) : ℤ)
        ≤ ∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)), g ((S ∩ J).card)) :
    Bound 8 689 := by
  intro H hH hC
  -- summing the hypothesis over the 448 faces
  have hfaces : (∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)),
      ∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J),
        ∑ v ∈ face5 J x, g (dj H J v)) ≤ 800 * 448 := by
    have h1 : ∀ J ∈ powersetCard 5 (univ : Finset (Fin 8)),
        (∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J),
          ∑ v ∈ face5 J x, g (dj H J v)) ≤ 800 * 8 := by
      intro J hJ
      have hc : J.card = 5 := (mem_powersetCard.mp hJ).2
      calc (∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J),
            ∑ v ∈ face5 J x, g (dj H J v))
          ≤ ∑ _x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J), (800 : ℤ) :=
            sum_le_sum fun x hx => h5 H J x hc (mem_filter.mp hx).2 hC
        _ = 800 * (((univ.filter fun x : V 8 => ∀ c ∈ supp x, c ∉ J).card : ℕ) : ℤ) := by
            rw [sum_const, nsmul_eq_mul]
            ring
        _ ≤ 800 * 8 := by
            have hcard := card_X_le J hc
            have : (((univ.filter fun x : V 8 => ∀ c ∈ supp x, c ∉ J).card : ℕ) : ℤ)
                ≤ ((8 : ℕ) : ℤ) := by exact_mod_cast hcard
            nlinarith
    calc (∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)),
          ∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J),
            ∑ v ∈ face5 J x, g (dj H J v))
        ≤ ∑ _J ∈ powersetCard 5 (univ : Finset (Fin 8)), (800 * 8 : ℤ) := sum_le_sum h1
      _ = 800 * 448 := by
          rw [sum_const, nsmul_eq_mul, card_powersetCard, card_univ, Fintype.card_fin]
          rw [show Nat.choose 8 5 = 56 from by decide]
          norm_num
  -- the same sum, reindexed: every vertex once per direction set
  have hEq : (∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)),
      ∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J),
        ∑ v ∈ face5 J x, g (dj H J v))
      = ∑ v : V 8, ∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)), g (dj H J v) := by
    calc (∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)),
          ∑ x ∈ univ.filter (fun x : V 8 => ∀ c ∈ supp x, c ∉ J),
            ∑ v ∈ face5 J x, g (dj H J v))
        = ∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)), ∑ v : V 8, g (dj H J v) :=
          sum_congr rfl fun J _ => sum_over_faces J (fun v => g (dj H J v))
      _ = ∑ v : V 8, ∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)), g (dj H J v) :=
          Finset.sum_comm
  -- and bounded below, vertex by vertex
  have hlow : (260 : ℤ) * (((∑ v : V 8, (nbr H v).card) : ℕ) : ℤ)
      ≤ ∑ v : V 8, ∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)), g (dj H J v) := by
    push_cast
    rw [Finset.mul_sum]
    exact sum_le_sum fun v _ => hvert (nbr H v)
  have hB : (260 : ℤ) * (((∑ v : V 8, (nbr H v).card) : ℕ) : ℤ) ≤ 800 * 448 :=
    le_trans hlow (by rw [← hEq]; exact hfaces)
  have hdeg := two_mul_card_le hH
  have hcast : (((2 * H.card : ℕ)) : ℤ) ≤ (((∑ v : V 8, (nbr H v).card) : ℕ) : ℤ) := by
    exact_mod_cast hdeg
  have hkey : (520 : ℤ) * ((H.card : ℕ) : ℤ) ≤ 800 * 448 := by
    have h1 : (260 : ℤ) * (((2 * H.card : ℕ)) : ℤ)
        ≤ (260 : ℤ) * (((∑ v : V 8, (nbr H v).card) : ℕ) : ℤ) := by
      have hnn : (0 : ℤ) ≤ 260 := by norm_num
      exact mul_le_mul_of_nonneg_left hcast hnn
    have h2 : (260 : ℤ) * (((2 * H.card : ℕ)) : ℤ) = 520 * ((H.card : ℕ) : ℤ) := by
      push_cast
      ring
    linarith
  omega

/-! ## 12. The two hypotheses, replaced by one statement about `Q_5`

`hvert` is a finite fact about the 256 subsets of the 8 directions; the kernel checks it.
`h5` follows from the inequality (*) for `Q_5` itself, because a 5-face of `Q_8` is a copy of
`Q_5`: `emb ι x` puts `Q_5` onto the face through `x`, `pull ι x H` is the part of `H` in the
face, and it has no 4-cycle if `H` has none. -/

/-- `260 d ≤ h(d)` at every vertex, written for the set `S` of its present directions. -/
theorem vertex_ineq : ∀ S : Finset (Fin 8),
    (260 : ℤ) * ((S.card : ℕ) : ℤ)
      ≤ ∑ J ∈ powersetCard 5 (univ : Finset (Fin 8)), g ((S ∩ J).card) := by
  decide +kernel

/-- (*) for `Q_5`: every C4-free edge set `K` of `Q_5` has `∑_u g(deg u) ≤ 800`. -/
def Q5Ineq : Prop :=
  ∀ K : Finset (Finset (V 5)), ¬ HasC4 K → (∑ u : V 5, g ((nbr K u).card)) ≤ 800

/-- Put the vertex `u` of `Q_5` into the 5-face of `Q_8` through `x`: coordinate `i` of `u`
goes to coordinate `ι i`, every other coordinate is the one of `x`. -/
def emb (ι : Fin 5 → Fin 8) (x : V 8) (u : V 5) : V 8 :=
  fun c => x c || decide (∃ i, ι i = c ∧ u i = true)

/-- A base point is `false` on the directions of its face. -/
theorem base_false {ι : Fin 5 → Fin 8} {J : Finset (Fin 8)} (hJ : ∀ c, c ∈ J ↔ ∃ i, ι i = c)
    {x : V 8} (hx : ∀ c ∈ supp x, c ∉ J) (i : Fin 5) : x (ι i) = false := by
  cases hb : x (ι i)
  · rfl
  · exact absurd ((hJ (ι i)).mpr ⟨i, rfl⟩) (hx (ι i) ((mem_supp x (ι i)).mpr hb))

theorem emb_on {ι : Fin 5 → Fin 8} (hι : Function.Injective ι) {x : V 8}
    (hx0 : ∀ i, x (ι i) = false) (u : V 5) (i : Fin 5) : emb ι x u (ι i) = u i := by
  show (x (ι i) || decide (∃ j, ι j = ι i ∧ u j = true)) = u i
  rw [hx0 i, Bool.false_or]
  cases hu : u i
  · rw [decide_eq_false_iff_not]
    rintro ⟨j, hj, huj⟩
    rw [hι hj, hu] at huj
    exact Bool.noConfusion huj
  · rw [decide_eq_true_iff]
    exact ⟨i, rfl, hu⟩

theorem emb_off {ι : Fin 5 → Fin 8} (x : V 8) (u : V 5) {c : Fin 8} (hc : ∀ i, ι i ≠ c) :
    emb ι x u c = x c := by
  show (x c || decide (∃ i, ι i = c ∧ u i = true)) = x c
  have h : decide (∃ i, ι i = c ∧ u i = true) = false := by
    rw [decide_eq_false_iff_not]
    rintro ⟨i, hi, _⟩
    exact hc i hi
  rw [h, Bool.or_false]

theorem emb_injective {ι : Fin 5 → Fin 8} (hι : Function.Injective ι) {x : V 8}
    (hx0 : ∀ i, x (ι i) = false) : Function.Injective (emb ι x) := by
  intro u u' h
  funext i
  have hi : emb ι x u (ι i) = emb ι x u' (ι i) := congr_fun h (ι i)
  rw [emb_on hι hx0 u i, emb_on hι hx0 u' i] at hi
  exact hi

/-- The embedding carries the edge of direction `i` to the edge of direction `ι i`. -/
theorem emb_flip {ι : Fin 5 → Fin 8} (hι : Function.Injective ι) {x : V 8}
    (hx0 : ∀ i, x (ι i) = false) (u : V 5) (i : Fin 5) :
    emb ι x (flip u i) = flip (emb ι x u) (ι i) := by
  funext c
  by_cases hc : ∃ j, ι j = c
  · obtain ⟨j, rfl⟩ := hc
    rw [emb_on hι hx0 (flip u i) j]
    by_cases hji : j = i
    · subst hji
      rw [flip_self, flip_self, emb_on hι hx0 u j]
    · rw [flip_of_ne hji, flip_of_ne (fun h => hji (hι h)), emb_on hι hx0 u j]
  · have hc' : ∀ j, ι j ≠ c := fun j hj => hc ⟨j, hj⟩
    rw [emb_off x (flip u i) hc', flip_of_ne (fun h => hc ⟨i, h.symm⟩), emb_off x u hc']

/-- The face through a base point is the image of `Q_5`. -/
theorem face5_eq_image {ι : Fin 5 → Fin 8} (hι : Function.Injective ι) {J : Finset (Fin 8)}
    (hJ : ∀ c, c ∈ J ↔ ∃ i, ι i = c) {x : V 8} (hx : ∀ c ∈ supp x, c ∉ J) :
    face5 J x = univ.image (emb ι x) := by
  have hx0 : ∀ i, x (ι i) = false := base_false hJ hx
  ext v
  rw [mem_image]
  simp only [face5, mem_filter, mem_univ, true_and]
  constructor
  · intro hv
    refine ⟨fun i => v (ι i), ?_⟩
    funext c
    by_cases hc : ∃ j, ι j = c
    · obtain ⟨j, rfl⟩ := hc
      exact emb_on hι hx0 (fun i => v (ι i)) j
    · have hc' : ∀ j, ι j ≠ c := fun j hj => hc ⟨j, hj⟩
      rw [emb_off x (fun i => v (ι i)) hc']
      exact (hv c (fun hcJ => hc ((hJ c).mp hcJ))).symm
  · rintro ⟨u, rfl⟩ c hc
    exact emb_off x u (fun j hj => hc ((hJ c).mpr ⟨j, hj⟩))

theorem sum_face5 {ι : Fin 5 → Fin 8} (hι : Function.Injective ι) {J : Finset (Fin 8)}
    (hJ : ∀ c, c ∈ J ↔ ∃ i, ι i = c) {x : V 8} (hx : ∀ c ∈ supp x, c ∉ J) (F : V 8 → ℤ) :
    ∑ v ∈ face5 J x, F v = ∑ u : V 5, F (emb ι x u) := by
  rw [face5_eq_image hι hJ hx, sum_image]
  intro a _ b _ h
  exact emb_injective hι (base_false hJ hx) h

/-- The part of `H` inside the face through `x`, read as an edge set of `Q_5`. -/
def pull (ι : Fin 5 → Fin 8) (x : V 8) (H : Finset (Finset (V 8))) : Finset (Finset (V 5)) :=
  univ.filter fun e => e.image (emb ι x) ∈ H

theorem pair_mem_pull {ι : Fin 5 → Fin 8} {x : V 8} {H : Finset (Finset (V 8))} (a b : V 5) :
    ({a, b} : Finset (V 5)) ∈ pull ι x H ↔ ({emb ι x a, emb ι x b} : Finset (V 8)) ∈ H := by
  simp only [pull, mem_filter, mem_univ, true_and, image_insert, image_singleton]

/-- A 4-cycle in the face is a 4-cycle of `H`. -/
theorem pull_noC4 {ι : Fin 5 → Fin 8} (hι : Function.Injective ι) {x : V 8}
    (hx0 : ∀ i, x (ι i) = false) {H : Finset (Finset (V 8))} (hC : ¬ HasC4 H) :
    ¬ HasC4 (pull ι x H) := by
  rintro ⟨a, b, c, d, hac, hbd, h1, h2, h3, h4⟩
  have hinj := emb_injective hι hx0
  exact hC ⟨emb ι x a, emb ι x b, emb ι x c, emb ι x d,
    fun h => hac (hinj h), fun h => hbd (hinj h),
    (pair_mem_pull a b).mp h1, (pair_mem_pull b c).mp h2,
    (pair_mem_pull c d).mp h3, (pair_mem_pull d a).mp h4⟩

/-- The degree in `Q_5` is the degree inside the face. -/
theorem card_nbr_pull {ι : Fin 5 → Fin 8} (hι : Function.Injective ι) {J : Finset (Fin 8)}
    (hJ : ∀ c, c ∈ J ↔ ∃ i, ι i = c) {x : V 8} (hx0 : ∀ i, x (ι i) = false)
    (H : Finset (Finset (V 8))) (u : V 5) :
    (nbr (pull ι x H) u).card = dj H J (emb ι x u) := by
  have hmap : (nbr (pull ι x H) u).image ι = nbr H (emb ι x u) ∩ J := by
    ext c
    rw [mem_image, mem_inter]
    constructor
    · rintro ⟨i, hi, rfl⟩
      rw [mem_nbr, pair_mem_pull, emb_flip hι hx0 u i] at hi
      exact ⟨mem_nbr.mpr hi, (hJ (ι i)).mpr ⟨i, rfl⟩⟩
    · rintro ⟨hc, hcJ⟩
      obtain ⟨i, rfl⟩ := (hJ c).mp hcJ
      refine ⟨i, ?_, rfl⟩
      rw [mem_nbr, pair_mem_pull, emb_flip hι hx0 u i]
      exact mem_nbr.mp hc
  show (nbr (pull ι x H) u).card = (nbr H (emb ι x u) ∩ J).card
  rw [← hmap, card_image_of_injective _ hι]

/-- (*) on one face, from (*) for `Q_5`. -/
theorem face_ineq_of_Q5Ineq (h : Q5Ineq) {ι : Fin 5 → Fin 8} (hι : Function.Injective ι)
    {J : Finset (Fin 8)} (hJ : ∀ c, c ∈ J ↔ ∃ i, ι i = c) (H : Finset (Finset (V 8)))
    {x : V 8} (hx : ∀ c ∈ supp x, c ∉ J) (hC : ¬ HasC4 H) :
    (∑ v ∈ face5 J x, g (dj H J v)) ≤ 800 := by
  have hx0 : ∀ i, x (ι i) = false := base_false hJ hx
  calc (∑ v ∈ face5 J x, g (dj H J v))
      = ∑ u : V 5, g (dj H J (emb ι x u)) := sum_face5 hι hJ hx (fun v => g (dj H J v))
    _ = ∑ u : V 5, g ((nbr (pull ι x H) u).card) :=
        sum_congr rfl fun u _ => by rw [card_nbr_pull hι hJ hx0 H u]
    _ ≤ 800 := h (pull ι x H) (pull_noC4 hι hx0 hC)

/-- A set of 5 directions is the range of an injective map from `Fin 5`. -/
theorem directions_of_card {J : Finset (Fin 8)} (hJ5 : J.card = 5) (c : Fin 8) :
    c ∈ J ↔ ∃ i, J.orderEmbOfFin hJ5 i = c := by
  constructor
  · intro hc
    have hmem : c ∈ Set.range (J.orderEmbOfFin hJ5) := by
      rw [Finset.range_orderEmbOfFin]
      exact Finset.mem_coe.mpr hc
    exact hmem
  · rintro ⟨i, rfl⟩
    exact Finset.orderEmbOfFin_mem J hJ5 i

/-- **`f(8) ≤ 689`, given the inequality (*) for `Q_5` and nothing else.** -/
theorem bound_8_689_of_Q5Ineq (h : Q5Ineq) : Bound 8 689 :=
  bound_8_689_of_q5
    (fun H J _x hJ5 hx hC =>
      face_ineq_of_Q5Ineq h (J.orderEmbOfFin hJ5).injective (directions_of_card hJ5) H hx hC)
    vertex_ineq

end Q8Bound689
