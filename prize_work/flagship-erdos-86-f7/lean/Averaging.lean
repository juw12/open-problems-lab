/-
  prize_work/flagship-erdos-86-f7/lean/Averaging.lean

  Lemma **F2** of the Erdős #86 flagship ledger (prize_work/flagship-erdos-86-f7/ledger.json):

      Each edge of Q_n lies in exactly n-1 of the 2n facets (sub-cubes Q_{n-1}), hence for
      every C4-free G, (n-1) e(G) <= 2n f(n-1).  Corollary with hypothesis f(6) = 132:
      f(7) <= 308.

  No proof of Erdős #86 is claimed here and no value of f(7) is claimed: the last theorem says
  only "IF every C4-free subgraph of Q_6 has at most 132 edges THEN every C4-free subgraph of
  Q_7 has at most 308 edges".  The hypothesis f(6) = 132 is not proved in this file.

  The file stands alone (the lab's Lean gate checks single files).  `V`, `hamDist`, `Adj`,
  `IsEdge` and `HasC4` are word for word the definitions of Cube.lean (lemma F1), where
  `HasC4` is shown to be the same as "contains a full 2-face".

  MAIN RESULTS
    * `card_facets_of_edge` -- an edge of Q_n lies in exactly n - 1 facets; there are 2n facets
      (`card_facets`).  A facet is a pair (i, b): the vertices with coordinate i equal to b.
    * `sum_card_facets` -- double count: the facets' edge counts add up to (n - 1) e(H).
    * `averaging` -- if every facet holds at most M edges of H then (n - 1) e(H) <= 2 n M.
    * `facet_card_le` -- a facet of Q_{n+1} is a copy of Q_n: for C4-free H, the edges of H in
      one facet form a C4-free edge set of Q_n of the same size, so there are at most f(n).
    * `averaging_c4free` -- f(n) <= M gives n e(H) <= 2 (n + 1) M for C4-free H in Q_{n+1}.
    * `bound_7_of_bound_6` -- f(6) <= 132 gives f(7) <= 308.

  Written and checked by the Vector session of 2026-11-25
  (cd prize_work/lean ; python check_claim.py ../flagship-erdos-86-f7/lean/Averaging.lean).
-/

import Mathlib

open Finset

namespace Averaging

/-! ## 1. The cube, as in Cube.lean -/

/-- Vertex set of the `n`-cube `Q_n`: bit vectors of length `n`. -/
abbrev V (n : ℕ) := Fin n → Bool

/-- Hamming distance: the number of coordinates where two vertices differ. -/
def hamDist {n : ℕ} (x y : V n) : ℕ := (univ.filter fun i => x i ≠ y i).card

/-- Adjacency in the cube: Hamming distance one. -/
def Adj {n : ℕ} (x y : V n) : Prop := hamDist x y = 1

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

/-! ## 2. Facets and the double count -/

/-- The edge `e` lies in the facet "coordinate `i` equals `b`". -/
def InFacet {n : ℕ} (i : Fin n) (b : Bool) (e : Finset (V n)) : Prop := ∀ x ∈ e, x i = b

instance {n : ℕ} (i : Fin n) (b : Bool) : DecidablePred (InFacet (n := n) i b) := by
  intro e
  unfold InFacet
  infer_instance

/-- `Q_n` has `2n` facets. -/
theorem card_facets (n : ℕ) : Fintype.card (Fin n × Bool) = 2 * n := by
  rw [Fintype.card_prod, Fintype.card_fin, Fintype.card_bool, mul_comm]

/-- **An edge of `Q_n` lies in exactly `n - 1` facets**: its two ends differ in one coordinate
`k`; the facets through the edge are "coordinate `i` = the common value", one for each `i ≠ k`. -/
theorem card_facets_of_edge {n : ℕ} {e : Finset (V n)} (he : IsEdge e) :
    (univ.filter fun p : Fin n × Bool => InFacet p.1 p.2 e).card = n - 1 := by
  obtain ⟨x, y, hxy, rfl⟩ := he
  unfold Adj hamDist at hxy
  obtain ⟨k, hk⟩ := card_eq_one.mp hxy
  have hdiff : ∀ i, x i ≠ y i ↔ i = k := by
    intro i
    have h := Finset.ext_iff.mp hk i
    simpa using h
  have hset : (univ.filter fun p : Fin n × Bool => InFacet p.1 p.2 ({x, y} : Finset (V n)))
      = (univ.erase k).image (fun i => (i, x i)) := by
    ext ⟨i, b⟩
    simp only [mem_filter, mem_univ, true_and, and_true, InFacet, mem_insert, mem_singleton,
      forall_eq_or_imp, forall_eq, mem_image, mem_erase, Prod.mk.injEq]
    constructor
    · rintro ⟨hx, hy⟩
      refine ⟨i, ?_, rfl, hx⟩
      intro hik
      exact (hdiff i).mpr hik (hx.trans hy.symm)
    · rintro ⟨j, hj, rfl, rfl⟩
      refine ⟨rfl, ?_⟩
      by_contra hne
      exact hj ((hdiff j).mp (fun h => hne h.symm))
  rw [hset, card_image_of_injective _ (fun i j h => (Prod.mk.inj h).1),
    card_erase_of_mem (mem_univ k), card_univ, Fintype.card_fin]

/-- **Double count.** Adding the number of edges of `H` in each of the `2n` facets counts every
edge `n - 1` times. -/
theorem sum_card_facets {n : ℕ} (H : Finset (Finset (V n))) (hH : ∀ e ∈ H, IsEdge e) :
    ∑ p : Fin n × Bool, (H.filter (InFacet p.1 p.2)).card = (n - 1) * H.card := by
  have h1 : ∀ p : Fin n × Bool, (H.filter (InFacet p.1 p.2)).card
      = ∑ e ∈ H, if InFacet p.1 p.2 e then 1 else 0 := fun p => card_filter _ _
  simp only [h1]
  rw [sum_comm]
  have h2 : ∀ e ∈ H, (∑ p : Fin n × Bool, if InFacet p.1 p.2 e then 1 else 0) = n - 1 := by
    intro e he
    rw [← card_filter]
    exact card_facets_of_edge (hH e he)
  rw [sum_congr rfl h2, sum_const, smul_eq_mul, mul_comm]

/-- **Facet averaging.** If every facet holds at most `M` edges of `H`, then
`(n - 1) e(H) ≤ 2 n M`. -/
theorem averaging {n : ℕ} (H : Finset (Finset (V n))) (hH : ∀ e ∈ H, IsEdge e) (M : ℕ)
    (hM : ∀ (i : Fin n) (b : Bool), (H.filter (InFacet i b)).card ≤ M) :
    (n - 1) * H.card ≤ 2 * n * M := by
  rw [← sum_card_facets H hH]
  calc ∑ p : Fin n × Bool, (H.filter (InFacet p.1 p.2)).card
      ≤ ∑ _p : Fin n × Bool, M := sum_le_sum (fun p _ => hM p.1 p.2)
    _ = 2 * n * M := by
      rw [sum_const, card_univ, card_facets, smul_eq_mul]

/-! ## 3. A facet of `Q_{n+1}` is a copy of `Q_n` -/

/-- Forget coordinate `i`. -/
def drop {n : ℕ} (i : Fin (n + 1)) (x : V (n + 1)) : V n := fun j => x (i.succAbove j)

/-- Put the value `b` in as coordinate `i`. -/
def ins {n : ℕ} (i : Fin (n + 1)) (b : Bool) (x : V n) : V (n + 1) :=
  Fin.insertNth (α := fun _ => Bool) i b x

theorem ins_same {n : ℕ} (i : Fin (n + 1)) (b : Bool) (x : V n) : ins i b x i = b := by
  simp [ins]

theorem ins_succAbove {n : ℕ} (i : Fin (n + 1)) (b : Bool) (x : V n) (j : Fin n) :
    ins i b x (i.succAbove j) = x j := by
  simp [ins]

theorem drop_ins {n : ℕ} (i : Fin (n + 1)) (b : Bool) (x : V n) : drop i (ins i b x) = x := by
  funext j
  simp [drop, ins]

theorem ins_drop {n : ℕ} {i : Fin (n + 1)} {b : Bool} {x : V (n + 1)} (hx : x i = b) :
    ins i b (drop i x) = x := by
  funext j
  by_cases hj : j = i
  · rw [hj, ins_same, hx]
  · obtain ⟨j', rfl⟩ := Fin.exists_succAbove_eq hj
    rw [ins_succAbove]
    rfl

theorem ins_injective {n : ℕ} (i : Fin (n + 1)) (b : Bool) : Function.Injective (ins (n := n) i b) := by
  intro x y h
  rw [← drop_ins i b x, h, drop_ins]

/-- Putting a coordinate in does not change the Hamming distance. -/
theorem hamDist_ins {n : ℕ} (i : Fin (n + 1)) (b : Bool) (x y : V n) :
    hamDist (ins i b x) (ins i b y) = hamDist x y := by
  unfold hamDist
  have hset : (univ.filter fun j => ins i b x j ≠ ins i b y j)
      = (univ.filter fun j => x j ≠ y j).image i.succAbove := by
    ext j
    simp only [mem_filter, mem_univ, true_and, mem_image]
    by_cases hj : j = i
    · rw [hj, ins_same, ins_same]
      constructor
      · intro h
        exact absurd rfl h
      · rintro ⟨j', _, hj'⟩
        exact absurd hj' (Fin.succAbove_ne i j')
    · obtain ⟨j', rfl⟩ := Fin.exists_succAbove_eq hj
      rw [ins_succAbove, ins_succAbove]
      constructor
      · intro h
        exact ⟨j', h, rfl⟩
      · rintro ⟨j'', h, hj''⟩
        rw [Fin.succAbove_right_injective hj''] at h
        exact h
  rw [hset, card_image_of_injective _ Fin.succAbove_right_injective]

/-- An edge that lies in the facet `(i, b)` is given back by forgetting and re-inserting `i`. -/
theorem image_ins_drop {n : ℕ} {i : Fin (n + 1)} {b : Bool} {e : Finset (V (n + 1))}
    (he : InFacet i b e) : (e.image (drop i)).image (ins i b) = e := by
  ext x
  simp only [mem_image, exists_exists_and_eq_and]
  constructor
  · rintro ⟨y, hy, rfl⟩
    rw [ins_drop (he y hy)]
    exact hy
  · intro hx
    exact ⟨x, hx, ins_drop (he x hx)⟩

/-- The edges of `H` inside the facet `(i, b)`, read as edges of `Q_n`. -/
def pull {n : ℕ} (i : Fin (n + 1)) (b : Bool) (H : Finset (Finset (V (n + 1)))) :
    Finset (Finset (V n)) :=
  (H.filter (InFacet i b)).image (fun e => e.image (drop i))

/-- Reading a facet as `Q_n` loses no edge. -/
theorem card_pull {n : ℕ} (i : Fin (n + 1)) (b : Bool) (H : Finset (Finset (V (n + 1)))) :
    (pull i b H).card = (H.filter (InFacet i b)).card := by
  unfold pull
  apply card_image_of_injOn
  intro e1 he1 e2 he2 h
  have h1 : InFacet i b e1 := (mem_filter.mp (mem_coe.mp he1)).2
  have h2 : InFacet i b e2 := (mem_filter.mp (mem_coe.mp he2)).2
  have h' : e1.image (drop i) = e2.image (drop i) := h
  rw [← image_ins_drop h1, h', image_ins_drop h2]

/-- A pair in the pulled-back set comes from the pair of the re-inserted vertices in `H`. -/
theorem pair_mem_of_mem_pull {n : ℕ} {i : Fin (n + 1)} {b : Bool}
    {H : Finset (Finset (V (n + 1)))} {u v : V n}
    (h : ({u, v} : Finset (V n)) ∈ pull i b H) :
    ({ins i b u, ins i b v} : Finset (V (n + 1))) ∈ H := by
  unfold pull at h
  obtain ⟨e, he, heq⟩ := mem_image.mp h
  rw [mem_filter] at he
  have hform : e = {ins i b u, ins i b v} := by
    rw [← image_ins_drop he.2, heq, image_insert, image_singleton]
  rw [← hform]
  exact he.1

/-- The pulled-back set consists of edges of `Q_n`. -/
theorem pull_isEdge {n : ℕ} (i : Fin (n + 1)) (b : Bool) (H : Finset (Finset (V (n + 1))))
    (hH : ∀ e ∈ H, IsEdge e) : ∀ e ∈ pull i b H, IsEdge e := by
  intro e' he'
  unfold pull at he'
  obtain ⟨e, he, rfl⟩ := mem_image.mp he'
  rw [mem_filter] at he
  obtain ⟨x, y, hxy, rfl⟩ := hH e he.1
  have hx : x i = b := he.2 x (by simp)
  have hy : y i = b := he.2 y (by simp)
  refine ⟨drop i x, drop i y, ?_, ?_⟩
  · unfold Adj
    rw [← hamDist_ins i b, ins_drop hx, ins_drop hy]
    exact hxy
  · rw [image_insert, image_singleton]

/-- A 4-cycle in the pulled-back set is a 4-cycle in `H`. -/
theorem hasC4_of_pull {n : ℕ} {i : Fin (n + 1)} {b : Bool} {H : Finset (Finset (V (n + 1)))}
    (h : HasC4 (pull i b H)) : HasC4 H := by
  obtain ⟨p, q, r, s, hpr, hqs, h1, h2, h3, h4⟩ := h
  exact ⟨ins i b p, ins i b q, ins i b r, ins i b s,
    fun e => hpr (ins_injective i b e), fun e => hqs (ins_injective i b e),
    pair_mem_of_mem_pull h1, pair_mem_of_mem_pull h2, pair_mem_of_mem_pull h3,
    pair_mem_of_mem_pull h4⟩

/-- **A facet holds at most `f(n)` edges of a C4-free subgraph of `Q_{n+1}`.** -/
theorem facet_card_le {n M : ℕ} (hB : Bound n M) (H : Finset (Finset (V (n + 1))))
    (hH : ∀ e ∈ H, IsEdge e) (hC : ¬ HasC4 H) (i : Fin (n + 1)) (b : Bool) :
    (H.filter (InFacet i b)).card ≤ M := by
  rw [← card_pull i b H]
  exact hB (pull i b H) (pull_isEdge i b H hH) (fun h => hC (hasC4_of_pull h))

/-! ## 4. F2 -/

/-- **F2.** If every C4-free subgraph of `Q_n` has at most `M` edges, then every C4-free
subgraph `H` of `Q_{n+1}` has `n · e(H) ≤ 2 (n + 1) M`. -/
theorem averaging_c4free {n M : ℕ} (hB : Bound n M) (H : Finset (Finset (V (n + 1))))
    (hH : ∀ e ∈ H, IsEdge e) (hC : ¬ HasC4 H) : n * H.card ≤ 2 * (n + 1) * M := by
  have h := averaging H hH M (facet_card_le hB H hH hC)
  rwa [Nat.add_sub_cancel] at h

/-- **F2, the corollary.** IF every C4-free subgraph of `Q_6` has at most 132 edges, THEN every
C4-free subgraph of `Q_7` has at most 308 edges. -/
theorem bound_7_of_bound_6 (h : Bound 6 132) : Bound 7 308 := by
  intro H hH hC
  have h6 := averaging_c4free h H hH hC
  omega

end Averaging
