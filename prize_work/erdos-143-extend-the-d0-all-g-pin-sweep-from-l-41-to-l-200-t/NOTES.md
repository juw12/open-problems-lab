# Job lab_erdos_261112_35 (RETRY) — #143: extend d0^all=1/4 g-pin sweep L=41→200

Problem: Erdős #143 / PPL-014 ($500). Dossier slug: `erdos-143-extend-the-d0-all-g-pin-sweep-from-l-41-to-l-200-t`.
Retry note: earlier run (log pool_..._0739.log) died at startup, left zero files. Restarting from job 261111_23 machinery (prize_work/ppl-014/).

## 1. FRONTIER (best published, next to lab's best)
- Published (erdosproblems.com/143, checked 2026-11-12): #143 $500 — countably infinite A⊂(1,∞), |kx−y|≥1 ∀x≠y,k≥1; does A sparse? ∑1/(x log x)<∞? or ∑_{x<n∈A}1/x = o(log n)?
  - o(log n) half: POSITIVE, Koukoulopoulos–Lamzouri–Lichtman (KLL25).
  - 1/(x log x) half: formalised NEGATIVE (Lean 4, non-constructive), github.com/mapice/erdos-143-candidate-disproof @0f4b969.
  - Its quantitative core: periodic-selection lemma per-round layer bound (S.card)/(4L) = d/4 → only needs exhibit density d > 2/3. No published smaller finite-core universal threshold (no published d0^all).
- Lab best: d0^all = min{k/L : g*(L,k)>1/6} = 1/4 for L≤40 at N=5e7 (job 261111_23, referee pass; exact g* L≤16 + analytic worst-class upper bound U'(L,k)=max_r Σ_{j∈S*} W[j][r], S*=k least-weighted classes, L=17..40; tightest L=30,k=7: U'=0.15369943043010403, margin 0.0130 to 1/6).
- This job: extend the U' rule-out to L=41..200 (N=5e6 sweep + full N=5e7 re-check via numpy build_W).

## 3. OPTIONS
1. (CARD) Reuse d0all_harden machinery; sweep L=41..200, k=1..(L-1)//4 (4k<L), at N=5e6 (pure Python build_W, ~0.9 s/L), then full re-check at N=5e7 with numpy-vectorized build_W (gated vs pure Python). Also compute L'(L,k)=max_r Σ(k smallest W[j][r]) = lower bound on g* (min-max ≥ max-min): L'>1/6 at 4k<L would PROVE a genuine crossing below 1/4 (a win). Gain: pin 1/4 to L≤200 or find un-ruled-out L. Odds ~0.45 (card). Cost ~30 min total.
2. (PNT/asymptotics, borrowed) Analytic N→∞ limit of U'(L,k) via PNT-in-AP for Ω-level counts → N-independent pin. Too big for 1 job (rmax truncation, error control); park as follow-up.
3. (LP/game duality, borrowed) g* = min_S max_r F(S,r); zero-sum game value v(L,k)=min_p max_r Σ_j p_j W[j][r] via LP, v≥g*... wait direction: v ≥ g* ≥ L'? No: mixed strategies help the min player → v ≤ g*. If v(L,k) ≤ 1/6 then g* ≤ 1/6 → rules out crossing; v is SHARPER than U' (U' = pure-strategy worst-class). Fallback ONLY if the 5e7 sweep finds U' > 1/6 somewhere.
Chose option 1 (card) first; 3 on standby.

## 4. KILL TEST (≤5 min)
numpy build_W prototype: (a) correctness gate vs pure-Python build_W at N=200k, L∈{4,41,200}, all entries <1e-12; (b) 1-L timing at N=5e7 — if one L > 30 s (or gate fails), numpy path dead → fallback: 5e6-only results + LP on tight L's; per card, report pin at 5e6 with explicit N caveat (partial).

## Reference values (gates)
- G1: level_mass(2e6, omega, L=4, S=(1,2,3)) max-layer == 1.56769507127912 (trusted machinery, per_round_weight_sharpening).
- G2 (5e6): exact g* L=3..16 → d0^all == 1/4, cross rows (4,1),(8,2), bracket_lo 3/13≈0.23077.
- G3 (5e7 code path vs job-23 referee-passed JSON prize_work/ppl-014/d0all_harden.json): U'(30,7)=0.15369943043010403; g*(4,1), g*(8,2) from its exact_rows (full precision).
- THRESH=1/6, RMAX=24, excl={2,3,4} (excl_set(10)=[2,3,4]), wt[n]=1/(n log n) zeroed at excl.

## Verifier design (verify200.py, N=5e7)
- omega_spf (Method 1) vs omega_prime_power (Method 3, nondeg_threshold_verify) ELEMENT-WISE, expect 0 mismatches (job-23 verifier: 0 at 5e7, total run 183 s).
- For L_check = {41, 200, tightest-L from run B, 30}: independent pure-Python direct double-loop W2 (no shared build_W) with Method-3 omega; entry-wise <1e-12 vs numpy W; fresh U' implementation vs run-B stored max_Up <1e-12; recompute g*(4,1), g*(8,2) vs run B + job-23.
- Cross-check run-B JSON per-L max_Up at checked L's.

## 6. CLASSIFY (final)
**extension** — the d0^all = 1/4 pin range extends L≤40 → L≤200 at N=5e7 (160 L's, 0 unruled, 0 proven crossings, worst margin 0.040184 @L=42; N-robust vs the 5e6 sweep). Published frontier (d>2/3, Apicella @0f4b969) unchanged; lab prior job 261111_23. Status **done**: no novel claim (same standard min/max exchange as job 23, now at 5x range) → no needs_review, no Lean statement (finite-N computational bound, no theorem asserted).

## Timing model
- omega_spf(5e7) pure Python ~60–90 s; wt table vectorized ~1 s (numpy) / ~20 s (Python loop); build_W pure-Python 5e7 ~7 s/L (job 23: 403 s total for 38 build_W + sieve + gates); numpy build_W target <2 s/L (prototype gate first).
- numpy arrays at 5e7: om int8 50 MB, wt float64 400 MB, n int32 200 MB. Job 23 ran 650 MB pure-Python fine.
