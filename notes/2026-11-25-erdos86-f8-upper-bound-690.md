---
Target: the lab's public notes repo; later, if no earlier source turns up, a comment on erdosproblems.com problem 86
Date: 2026-11-25
Status: proof checked by Lean; no referee and no human reader yet; the bound 690 is NOT the best that follows
  from published work (see the second update)
Updated: 2026-11-25, same day: Chung (1992) was read. Step 1 and the bound 691 are credited to her Lemmas 4-5;
  the bound 690 is not in her paper. See "What is not known".
Updated: 2026-11-25, afternoon: Thomason and Wagner (2009) was read. A table printed in their paper gives
  f(8) <= 689 by averaging over sub-5-cubes; the number is not stated there. So 690 is weaker than a consequence
  of data published in 2009. See "Thomason-Wagner's table gives 689". The proof below is unchanged.
Updated: 2026-11-25, evening: the one inequality that 689 needs from that table was checked directly, without
  the table, by an exact solver (CP-SAT: no C4-free subgraph of Q5 above 800). A solver verdict, not a proof
  that Lean or a person has checked. Same section.
Updated: 2026-11-25, night: the step from that one inequality to f(8) <= 689 is now checked by Lean
  (`Q8Bound689.bound_8_689_of_Q5Ineq`). So 689 rests on one statement about the 5-cube, and that statement
  rests on solver runs, not on Lean. A second kind of run (SAT solver, each "no" with a proof file that a
  separate small checker verifies) was started. Same section.
Flagship: F2, F4, F5
Scope: one finite statement about the 8-cube. Claimed: f(8) <= 690, with a Lean 4 proof and no computer search.
  Not claimed: that 690 is new or the best known bound (it is not, see the second update), that it is the true
  value, or anything about the asymptotic question of problem 86.
AI-assisted: yes (statement, hand proof and Lean proof were written by an AI agent)
---

> **Status (published 2026-11-25).** This note comes from an AI-agent research lab run by
> Juwen Lie. The proof was written by an AI agent and checked by the Lean 4 kernel with
> Mathlib (standard axioms only); the Lean files are in this repository under `prize_work/`.
> The kernel checks the proof, not whether the formal statement says what this note says:
> please read the definitions. **No human mathematician has reviewed it**. The bound is
> not the best that follows from published work: a table of Thomason and Wagner (2009) gives
> 689 (section "Thomason-Wagner's table gives 689"; Lean checks the step from one inequality of
> that table to 689, not the inequality). It has not been posted to the problem's own site.
> "I" below means the lab. Nothing here settles the open problem itself. Corrections and
> references are welcome: please open an issue.

# f(8) <= 690: a C4-free subgraph of the 8-cube has at most 690 edges (Lean-checked)

f(n) is the largest number of edges in a subgraph of the n-cube Q_n with no 4-cycle. Erdos problem 86 asks
how f(n) grows; the small values are f(4) = 24, f(5) = 56, f(6) = 132 and, since 20 Nov 2026, f(7) = 304
(Zhou, Ma, Yang, Ren, Yang, Li, arXiv:2611.27284). For n = 8 the lab knows of:

| Bound | Source |
|---|---|
| f(8) >= 682 | odd-square set from the fully frustrated ground state of Marinari, Parisi and Ritort; edge list and proof of the odd-square optimum in Minamoto, arXiv:2603.29127 v5 |
| f(8) <= 694 | averaging over the 16 facets from f(7) = 304: 7 e <= 16 * 304 |
| f(8) <= 691 | star capacity (Chung 1992, Lemmas 4-5; per parity class it is Lemma 2.2 of arXiv:2611.27284; step 1 below) with C(d,3) >= 10 d - 40; `Stars.bound_8_691` |
| **f(8) <= 690** | **this note; `Q8Parity.bound_8_690`**; no computer search, checked by Lean |
| f(8) <= 689 | follows from one inequality (*) about C4-free subgraphs of Q5, read off Table 1 of Thomason-Wagner (2009), by averaging over the 448 sub-5-cubes; the number is not printed there. The step from (*) to 689 is checked by Lean (`Q8Bound689.bound_8_689_of_Q5Ineq`). (*) itself is not: it rests on their table (a computer enumeration, not re-computed here) or on the lab's exact-solver runs; section "Thomason-Wagner's table gives 689" |

So 682 <= f(8) <= 690 with a proof the Lean kernel has checked from start to end, and 682 <= f(8) <= 689 with
a Lean-checked proof from one inequality in the 5-cube, which is accepted on the word of that table or of the
solver runs. The gap is not closed.

## What was checked, and by what

- `prize_work/flagship-erdos-86-f7/lean/Q8Parity.lean`, theorem `Q8Parity.bound_8_690 : Bound 8 690`.
  Lean 4 with Mathlib; 56 declarations; the axiom audit reports the standard axioms only (no `sorry`, no
  `native_decide`). The file stands alone (it imports Mathlib and nothing else).
- `prize_work/flagship-erdos-86-f7/lean/Stars.lean`, theorems `bound_8_691`, `bound_7_306`, `bound_6_133`,
  `bound_5_56`, `bound_4_24` from one general inequality `edge_bound`.
- `prize_work/flagship-erdos-86-f7/lean/Averaging.lean`, theorem `averaging_c4free`: a bound M for Q_n gives
  n * e <= 2 (n + 1) M for Q_(n+1).
- `prize_work/flagship-erdos-86-f7/lean/Q8Bound689.lean`, theorem
  `Q8Bound689.bound_8_689_of_Q5Ineq : Q5Ineq -> Bound 8 689`. It is a conditional statement: 689 follows IF the
  inequality (*) holds in Q5. Lean does not check (*). Statement and limits in the section
  "Thomason-Wagner's table gives 689".
- All four files are checked again by Lean on every change, on GitHub's machines and not on the lab's PC:
  workflow `.github/workflows/lean-check.yml`, result in `prize_work/flagship-erdos-86-f7/lean/check/lean_check.md`
  (first full pass: https://github.com/juw12/open-problems-lab/actions/runs/45603989943; Lean 4.32.0 with Mathlib 4.32.0; no error, no `sorry`, standard axioms only).

The formal statement, so that a reader can judge whether it says what this note says:

    abbrev V (n : Nat) := Fin n -> Bool                       -- vertices of Q_n
    def Adj (x y : V n) : Prop := hamDist x y = 1              -- differ in exactly one coordinate
    def IsEdge (e : Finset (V n)) : Prop := exists x y, Adj x y /\ e = {x, y}
    def HasC4 (H : Finset (Finset (V n))) : Prop :=
      exists a b c d, a != c /\ b != d /\ {a, b} in H /\ {b, c} in H /\ {c, d} in H /\ {d, a} in H
    def Bound (n M : Nat) : Prop :=
      forall H : Finset (Finset (V n)), (forall e in H, IsEdge e) -> not (HasC4 H) -> H.card <= M
    theorem bound_8_690 : Bound 8 690

(ASCII rendering; the file has the exact text.) The kernel checks the proof of this statement. It does not check
that these definitions are the right reading of "C4-free subgraph of the cube"; that is for the reader.

## The proof in words

Let G be C4-free in Q_8 with at least 691 edges; delete edges until it has exactly 691. Let P be one parity
class (128 vertices; every edge has exactly one end in P). A 3-face is a 3-dimensional subcube; there are
C(8,3) * 2^5 = 1792. A star in a 3-face is a vertex of P in the face whose three edges inside the face are all
in G.

1. **One star per 3-face (Chung 1992, Lemmas 4-5).** Two vertices of P in the same 3-face are at distance 2, and two stars there would
   close a 4-cycle. A vertex of degree d is the root of C(d,3) stars in C(d,3) different 3-faces. Hence
   S := sum over v in P of C(d(v),3) <= 1792.
2. **All degrees on P are 5 or 6.** w(d) := C(d,3) - 10 d + 40 takes the values 40, 30, 20, 11, 4, 0, 0, 5, 16
   for d = 0..8: it is 0 exactly at d = 5, 6 and at least 4 elsewhere. The degrees on P add up to 691, so
   S = 10 * 691 - 40 * 128 + sum of w = 1790 + sum of w <= 1792. So every w is 0 and S = 1790.
3. **Exactly two 3-faces carry no star**, since 1792 - 1790 = 2.
4. **Facet parity.** A facet K = {x : x_c = b} is a Q_7 with 560 3-faces. A vertex of P in K has degree
   4, 5 or 6 inside K, and C(4,3) = 4, C(5,3) = 10, C(6,3) = 20 are all even. So the number of stars inside K is
   even, and they sit in different 3-faces of K. 560 is even, so K contains an even number of star-less 3-faces:
   none or both of the two.
5. **Contradiction.** Two different 3-faces are separated by some facet: if their free coordinates differ, take a
   coordinate c free in the second and fixed in the first; if not, their fixed values differ at some c. The facet
   through the first face at c contains the first and not the second. That facet holds exactly one star-less face,
   against step 4.

So no C4-free subgraph of Q_8 has 691 edges. The Lean proof follows the same steps but counts uncovered
(direction, vertex) pairs, 16 of them, instead of naming the two faces.

## Checks against known values

- The same Lean inequality (`Stars.edge_bound`) gives f(4) <= 24 and f(5) <= 56, the exact values, and
  f(6) <= 133, f(7) <= 306 against the true 132 and 304. So the formal statement is neither vacuous nor too strong
  in the cases where the truth is known.
- 690 is above the 682 construction, as it must be.
- Step 4 has no analogue for n = 7 (local degrees 3, 4, 5 give 1, 4, 10 stars, mixed parity). This fits the fact
  that f(7) = 304 needed the much longer argument of arXiv:2611.27284.

## A consequence by arithmetic

With `Averaging.averaging_c4free`, 8 * f(9) <= 18 * f(8) <= 18 * 690 = 12420, so f(9) <= 1552; the same step from
694 gives 1561. (Thomason-Wagner's table gives f(9) <= 1549 directly, see below.) This line is arithmetic on top of two checked files. It is not one Lean statement, because the two
files each carry their own copy of the definitions.

## Thomason-Wagner's table gives 689

A. Thomason and P. Wagner, Bounding the size of square-free subgraphs of the hypercube, Discrete Math. 309
(2009) 1730-1735, read on 2026-11-25. The paper bounds the limit density by 0.62256 and states no bound for
n = 7 or 8. Its Table 1 lists the 33 extreme points of the convex hull of the star vectors (s1, ..., s5),
s_r = sum over v of C(d(v), r), of all C4-free subgraphs of Q5. The table comes from a computer run of about
17 hours. On all 33 points

    12 s1 - 5 s2 + 11 s4 - 20 s5 <= 800                                   (*)

with equality at four of them, among them (112, 144, 80, 16, 0), the 56-edge graph. So (*) holds for every
C4-free subgraph of Q5. Let G be C4-free in Q8. There are C(8,5) * 2^3 = 448 sub-5-cubes, and a star with r
edges lies in C(8-r, 5-r) = 35, 20, 10, 4, 1 of them (r = 1..5). Adding (*) over the sub-5-cubes:

    sum over v of h(d(v)) <= 800 * 448 = 358400,   h(d) = 420 d - 100 C(d,2) + 44 C(d,4) - 20 C(d,5).

For d = 0..8 the values of h(d) - 260 d are 0, 160, 220, 180, 84, 0, 0, 140, 440, never negative. So
520 e <= 358400, e <= 8960/13 = 689.23..., and f(8) <= 689.

- This derivation was made by an AI agent and has not been read by a mathematician. The lab did not find the
  number 689 in the paper or anywhere else, and does not know whether it is written down somewhere.
- **The step from (*) to 689 is checked by Lean.** In `Q8Bound689.lean`:

      def g (d : Nat) : Int := 12 d - 5 C(d,2) + 11 C(d,4) - 20 C(d,5)     -- 0, 12, 19, 21, 29, 45 for d = 0..5
      def Q5Ineq : Prop :=
        forall K : Finset (Finset (V 5)), not (HasC4 K) -> (sum over u : V 5 of g (number of edges of K at u)) <= 800
      theorem bound_8_689_of_Q5Ineq (h : Q5Ineq) : Bound 8 689

  (ASCII rendering; the file has the exact text. `V`, `HasC4`, `Bound` as in the statement above; the sum of
  g over the vertices is the left side of (*).) The proof is the averaging written out above: each sub-5-cube
  of Q8 is a copy of Q5, the 448 copies are added up, and the vertex inequality h(d) >= 260 d is checked by
  the kernel for all 256 sets of directions. No `sorry`, no `native_decide`, standard axioms only
  (https://github.com/juw12/open-problems-lab/actions/runs/45603989943, 17 s). What this buys: nobody has to trust the averaging or the arithmetic any more. What it
  does not buy: (*) is a hypothesis of the theorem, and Lean says nothing about it.
- **(*) itself: two sources, neither of them a proof a person has read.** (a) Table 1, copied from the paper
  and not re-computed. (b) The lab's own check without the table, 2026-11-25: an exact solver (CP-SAT,
  OR-Tools 9.15) was asked for a C4-free subgraph of Q5 with left side >= 801 and answered "infeasible" in
  every case (largest degree 4: 20 s; a vertex of degree 5: 20 min; largest degree <= 3 needs no run, the sum
  is at most 672). Asked for 800 it finds a graph in both cases, so the model is not empty by mistake. A
  second solver (HiGHS) agrees on every case it has finished; two sub-cases were still running when this was
  written. These runs leave no proof file, and the case split uses a symmetry argument nobody has read. Runs:
  https://github.com/juw12/open-problems-lab/actions/runs/45592147713 and .../runs/45596181027; record in
  `prize_work/flagship-erdos-86-f7/sources/q5_inequality_RESULT.md`.
- **A third run that leaves proofs (started 2026-11-25, not finished when this was written).** The same
  question as a formula for a SAT solver (`q5_sat_proof.py`, CaDiCaL 2.1.3). For each "no" the solver writes
  a DRAT proof and drat-trim, a separate small program, checks it; a case counts only if the checker says
  VERIFIED. Then a "no" rests on the encoder (tested against brute force on Q3 and against CP-SAT on Q4) and
  on the checker, not on the solver. The easy cases are verified; the hard ones may not finish in the time
  limit. Run: https://github.com/juw12/open-problems-lab/actions/runs/45604846118; verdict, when the run ends, in
  `prize_work/flagship-erdos-86-f7/sources/q5_sat_runs/q5_sat_verdict.md`. Even a full pass would not be a
  Lean proof of (*): the encoder and the case split would still be trusted.
- (*) is not a fact about one parity class: the largest sum of g(d) = 12 d - 5 C(d,2) + 11 C(d,4) - 20 C(d,5)
  over the 16 vertices of one class is 404, not 400 (`q5_class_bound.py`; the other class then has 382).
- The same table gives f(5) <= 56, f(6) <= 133, f(7) <= 306, f(9) <= 1549 and f(10) <= 3398.
- Scripts: `prize_work/flagship-erdos-86-f7/sources/thomason_wagner_certificate.py` (exact multipliers, checked
  in fractions), `thomason_wagner_star_lp.py` (the linear program), `q5_inequality_cases.py` and `q5_verdict.py`
  (the solver check; `q5_inequality_ilp.py` is the first attempt, which did not finish), `q5_sat_proof.py`
  (the SAT run with checked proofs).

What this note still offers: 690 by a proof of about one page with no computer search, checked by the Lean
kernel.

## What is not known

- **Whether 690 was written down before.** The lab found no upper bound for f(8) in arXiv:2611.27284, arXiv:2603.29127, the OEIS
  entry A245762 or the page of problem 86 (all read on 2026-11-25; the full text of both papers was searched for
  690 to 694, no hit; arXiv:2603.29127 lists "determine ex(Q_8, C_4) exactly, is it 682?" as open). The two
  flag-algebra papers arXiv:1201.0209 and arXiv:1201.3587 give asymptotic densities only.
- **What Chung's 1992 paper has** (F. R. K. Chung, Subgraphs of a hypercube containing no small even cycles,
  J. Graph Theory 16, 273-286; read on 2026-11-25 from the author's page). Her Lemmas 4 and 5 are exact for every
  n and give sum over all v of C(d(v),3) <= 2 C(n,3) 2^(n-3), which is step 1 for both parity classes together.
  So step 1 is hers, and 691 follows from her paper in two lines, although she does not write the number down.
  Her Lemma 11, read as an exact inequality, would give exactly 690 for n = 8. It is not exact: it is one step of
  an argument for large n, its proof uses Goodman's bound in the large-t form, and the printed inequality fails in
  Q_4 for the 14-edge forest made of the edges {0000,0001} and {1110,1111} with all edges at their four ends
  (left side 0, right side 3/4). So 690 is not stated in her paper and does not follow from what is proved
  there. Script: `prize_work/flagship-erdos-86-f7/chung_check.py`.
- **Still not read:** Brass, Harborth and Nienborg (J. Graph Theory 19, 1995). An earlier bound may be there or
  elsewhere. Only the parity step 4 may be new, and it is short enough that it may well be folklore. A reference
  to an earlier source is the most useful reply to this note.
- **The true value.** For 690 edges the same count leaves a slack of 12 instead of 2, and configurations such as
  four star-less faces with the same free coordinates, whose fixed values form an affine plane, pass the parity
  test. Reaching 689 this way needs another idea; the section above gets 689 from the computation of
  Thomason and Wagner. A proof of (*) that Lean or a person can check (from verified SAT proofs, or a short
  argument by hand) would make 689 unconditional; the lab does not have one.
- Nothing here bears on the asymptotic question of problem 86.
