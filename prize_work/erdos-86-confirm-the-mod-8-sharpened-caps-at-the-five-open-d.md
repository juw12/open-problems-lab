# #86 — The mod-8 sharpened caps at the five open dimensions, as a citable table

**Job:** `lab_erdos_261117_2` · **Date:** 2026-11-17 · **Idea card:** `erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d` (value 5/10, odds 0.5)
**Status:** verified computation, all gates pass (92 gates in the main script + 32 checks in an independent verifier + the previous job's gate re-run, all exit 0).
**Novelty class: `reproduction`** — this table is a **corollary of a published lemma** (Minamoto, arXiv:2603.29127v5, Lemma 12), applied at dimensions where the paper does not apply it. It is **not** past the published frontier of the problem, and no proof or prize claim is made.

---

## 1. FRONTIER

### Exact statement (live fetch of `https://www.erdosproblems.com/86`, 2026-11-17)

> **\$100.** Let $Q_n$ be the $n$-dimensional hypercube graph (so that $Q_n$ has $2^n$ vertices and $n2^{n-1}$ edges). Is it true that every subgraph of $Q_n$ with $\geq \left(\frac12 + o(1)\right)n2^{n-1}$ many edges contains a $C_4$?
>
> Let $f(n)$ be the maximum number of edges in a subgraph of $Q_n$ without a $C_4$, so that this conjecture is $f(n) \leq \left(\frac12+o(1)\right)n2^{n-1}$. Erdős [Er91] showed $f(n) \geq \left(\frac12 + c/n\right)n2^{n-1}$ for some $c>0$ … Brass, Harborth and Nienborg [BHN95] improved this to $\left(\frac12 + c/\sqrt n\right)n2^{n-1}$. Balogh, Hu, Lidický and Liu [BHLL14] proved $f(n)\leq 0.6068\,n2^{n-1}$, improved to $\leq 0.60318\,n2^{n-1}$ by Baber [Ba12b].

- **Prize:** \$100. **Forum status (live, 2026-11-17):** OPEN, no proof claim on the page, same text as the 2026-11-14 capture used by job 261114_3 and the 2026-11-17 Director pass. OEIS A245762.
- **Published frontier on exactly this question (the two-point cap and its quantisation):**
  - **DPTV79** (Derrida–Pomeau–Toulouse–Vannimenus, *J. Physique* 40 (1979) p. 620): the field identity $\sum_v h(v)^2 = n2^{n-1}$, the parity $h(v)\equiv n \pmod 2$, and the **plain two-point bound** $S \le \min_a \frac{N(a)}{2a+2}$.
  - **Minamoto, arXiv:2603.29127v5, Lemma 9** = the same two-point bound; **Lemma 12 ("quantised two-point inequality")** = the mod-8 refinement $(h-a)(h-a-2)\equiv 0 \pmod 8$.
  - **The paper applies Lemma 12 only at $n=6,7,8$.** No published paper states the quantised cap at any of the seven open dimensions. The best **published** rigorous upper bound at $n=12$ is therefore the plain cap **15701**; nobody has published 15700.
- **Lab's own best numbers before this job:** $15614 \le g(12) \le 15700$ (`prize_work/erdos-86-n12-publish-grade-bracket-15614-15700.md`, jobs 261114_5/46/47), and the 24-row slack table of job 261114_3 (`prize_work/erdos-86-two-point-cap-slack-parity/`).

### The subproblem

$g(n)$ = maximum size of an **odd-square** edge set (every 4-face meets it in 1 or 3 edges) in $Q_n$; $g(n) \le f(n)$, and in the fully-frustrated (signing) model it is governed by a field $h: V(Q_n)\to\mathbb Z$ with $\sum_v h(v)^2 = n2^{n-1} =: E$, $h(v)\equiv n \pmod 2$, $0\le h(v)\le n$. With $S=\sum_v h(v)$ (always even, since $S\equiv mn \equiv 0 \pmod 2$) one has $g(n)\le\frac{E+S}{2}$; the whole table below is about the largest admissible $S$.

- **Exact $g(n)$ is published for 17 dimensions** — $n=3,4,5,6,7,8,9,10,11,13,15,16,17,19,21,23,26$ (Minamoto v5; Wang arXiv:2611.14906; BHN95 for $n=4^r$; OEIS A245762).
- **Open dimensions: 12, 14, 18, 20, 22, 24, 25 and every $n\ge 27$.**

## 2. MEMORY

Built on `lab_erdos_261114_3` (the slack-parity job): its 24-row table, its finding that the quantisation bites by exactly 1 at $n=12,18,20,22,24$ and not at $n=14,25$, and its parity/realisability classification of the $n=6$ and $n=9$ slacks. This job re-derives all of it **from scratch with a different algorithm** and turns it into the table a referee can cite without re-running job 3's machinery. The bank's earlier dead end ("cap attained 10/10 at $n\ge10$", used to retire the congruence-refinement card) stays dead: job 3 refuted it at $n=12$ and the certificate below reproduces that refutation independently.

## 3. OPTIONS (gain × odds / cost)

1. **(Card, worked first)** Fresh exact-integer re-derivation of the two-point programme optimum and its mod-8 quantisation, plus an explicit extremal multiset at each open dimension: *integer-programming / rounding-theory view*. Gain: a citable table + new sharpened caps at dims nobody tabulated. Odds: high (published lemma, mechanical). Cost: ~30 min. **Chosen.**
2. **Higher-moment (3rd/4th moment) DPTV identity** to push the cap below $\frac{E+S_Q}{2}$ at $n=12$. Gain: could move 15700 → 15699 or lower (frontier). Odds: low (no such identity in the literature; would need a new algebraic identity). Cost: high. *Not pursued; remains the only route to a genuine frontier move here.*
3. **SAT/CP-SAT on field realisability** (does a signing reach $S_Q$ at the open dims?) — *borrowed from automated reasoning*. Gain: could close $g(12)$ exactly. Odds: very low on a laptop ($2^{12}$ signings; the $n=12$ search was already retired as a dead end in the bank). Cost: high. *Not pursued.*

## 4. KILL TEST (≤ 5 min, run first)

Recompute $S(12)$ by the fresh descent. **Result: $S(12)=6824$ with the multiset $\{h{=}0:2,\ h{=}2:680,\ h{=}4:1366\}$** — exactly the card's required values (count 2048, $\sum h = 6824$, $\sum h^2 = 24576$). **Kill test passed**, so the table can be certified from this route.

## 5. RESULT — the citable table

Conventions: $E=n2^{n-1}$; $S_Q$ = largest even programme optimum with the mod-8 congruence; `cap` = plain DPTV79 bound $=\frac{E+2\lfloor S^*_{\text{real}}/2\rfloor}{2}$; `cap*` $=\frac{E+S_Q}{2}$ (Lemma 12). All arithmetic exact (Python `int` and `fractions.Fraction`; no floating point anywhere in either script).

| $n$ | $E=n2^{n-1}$ | $a^*$ | $S_Q$ | cap (plain, DPTV79) | **cap\* (mod-8, Lemma 12)** | cap − cap\* | published $g(n)$ |
|---|---|---|---|---|---|---|---|
| 3 | 12 | 1 | 6 | 9 | 9 | 0 | 9 |
| 4 | 32 | 0 | 16 | 24 | 24 | 0 | 24 |
| 5 | 80 | 1 | 32 | 56 | 56 | 0 | 56 |
| 6 | 192 | 2 | 72 | 133 | **132** | **1** | 132 |
| 7 | 448 | 1 | 160 | 304 | 304 | 0 | 304 |
| 8 | 1024 | 2 | 340 | 682 | 682 | 0 | 682 |
| 9 | 2304 | 1 | 768 | 1536 | 1536 | 0 | 1532 |
| 10 | 5120 | 2 | 1536 | 3328 | 3328 | 0 | 3328 |
| 11 | 11264 | 3 | 3328 | 7296 | 7296 | 0 | 7296 |
| **12** | 24576 | 2 | **6824** | **15701** | **15700** | **1** | open |
| 13 | 53248 | 3 | 14336 | 33792 | 33792 | 0 | 33792 |
| **14** | 114688 | 2 | 30036 | 72362 | 72362 | 0 | open |
| 15 | 245760 | 3 | 61440 | 153600 | 153600 | 0 | 153600 |
| 16 | 524288 | 2 | 131072 | 327680 | 327680 | 0 | 327680 |
| 17 | 1114112 | 3 | 262144 | 688128 | 688128 | 0 | 688128 |
| **18** | 2359296 | 4 | **550500** | **1454899** | **1454898** | **1** | open |
| 19 | 4980736 | 3 | 1114112 | 3047424 | 3047424 | 0 | 3047424 |
| **20** | 10485760 | 4 | **2306864** | **6396313** | **6396312** | **1** | open |
| 21 | 22020096 | 3 | 4718592 | 13369344 | 13369344 | 0 | 13369344 |
| **22** | 46137344 | 4 | **9646896** | **27892121** | **27892120** | **1** | open |
| 23 | 96468992 | 3 | 19922944 | 58195968 | 58195968 | 0 | 58195968 |
| **24** | 201326592 | 4 | **40265316** | **120795955** | **120795954** | **1** | open |
| **25** | 419430400 | 3 | 83886080 | 251658240 | 251658240 | 0 | open |
| 26 | 872415232 | 4 | 167772160 | 520093696 | 520093696 | 0 | 520093696 |

**The five sharpened open dimensions — the card's deliverable, plain cap next to the sharpened one:**

| $n$ | plain cap (published method) | **sharpened cap (Lemma 12)** | $S_Q$ | binding $a$ | $\Delta=N-W S_Q$ |
|---|---|---|---|---|---|
| 12 | 15701 | **15700** | 6824 | 2 | 16 |
| 18 | 1454899 | **1454898** | 550500 | 4 | 24 |
| 20 | 6396313 | **6396312** | 2306864 | 4 | 32 |
| 22 | 27892121 | **27892120** | 9646896 | 4 | 32 |
| 24 | 120795955 | **120795954** | 40265316 | 4 | 24 |
| 14 | 72362 | 72362 (no gain) | 30036 | 2 | 8 |
| 25 | 251658240 | 251658240 (no gain) | 83886080 | 3 | 0 |

### Extremal multisets (the certificates that make $S_Q$ the exact programme optimum)

Each multiset has exactly $m=2^{n-1}$ entries, $\sum h = S_Q$ and $\sum h^2 = E$ — verified by exact arithmetic in both scripts.

| $n$ | multiset (h : multiplicity) | count | $\sum h$ | $\sum h^2$ |
|---|---|---|---|---|
| 12 | 0:2, 2:680, 4:1366 | 2048 | 6824 | 24576 |
| 14 | 0:1, 2:1364, 4:6827 | 8192 | 30036 | 114688 |
| 18 | 0:1, 4:117963, 6:13108 | 131072 | 550500 | 2359296 |
| 20 | 2:4, 4:419424, 6:104860 | 524288 | 2306864 | 10485760 |
| 22 | 2:4, 4:1468000, 6:629148 | 2097152 | 9646896 | 46137344 |
| 24 | 0:1, 4:5033163, 6:3355444 | 8388608 | 40265316 | 201326592 |
| 25 | 1:0, 5:16777216 (i.e. all $h=5$) | 16777216 | 83886080 | 419430400 |

Consequently, at each of the seven open dimensions the **relaxed programme optimum is exactly $S_Q$**: Lemma 12 gives $S\le S_Q$, the multiset realises $S=S_Q$. (Whether a *field* realises it — i.e. whether $g(n)$ actually attains cap\* — is a different question, and is exactly the $n=9$-type realisability gap.)

### Why the quantised $S$ is a genuine upper bound (Lemma 12 bookkeeping, checked at every open dimension)

Lemma 12's hypotheses as used: $h\equiv a\equiv n \pmod 2$ (DPTV parity — supplied by the field), so with $t=(h-a)/2$, $(h-a)(h-a-2)=4t(t-1)$ is $\ge 0$ and $\equiv 0 \pmod 8$ termwise. Summing gives, with $W=2a+2$ and $N(a)=E+a(a+2)m$,
$$\Delta(a)=N(a)-W S=\sum_v (h(v)-a)(h(v)-a-2)\ \ge 0,\qquad \Delta(a)\equiv 0 \pmod 8 .$$
Gates G5a/G5b check, at each of the seven open $n$ and the binding $a$: $\Delta\ge0$, $\Delta\equiv0\pmod 8$, and that $S_Q+1$ already violates one of the two conditions (maximality). The verifier V5b evaluates the identity **both ways** — $\sum_v (h-a)(h-a-2)$ over the multiset versus $N(a)-WS$ — and gets the same value at each dimension. Parity bookkeeping: $S$ must be even, and over $n\le200$ the sharpening is **never** decided by that parity step (gate G7d / V8: at every sharpened $n$ the binding $S_{\max}$ is already even, so the −1 comes purely from the mod-8 congruence).

### Range extension (where this goes past job 261114_3's five values)

The same exact computation for every $n=3,\dots,200$ (198 rows, `cap_table.json` → `wide_rows`) gives **36 dimensions where cap\* = cap − 1**:
$$6, 12, 18, 20, 22, 24, 38, 40, 48, 56, 60, 66, 70, 72, 74, 84, 88, 90, 92, 102, 108, 116, 128, 130, 132, 140, 146, 150, 152, 158, 162, 164, 168, 180, 186, 198 .$$
**30 of these are open dimensions $\ge27$ that no table in the literature (and no earlier lab job) tabulated** — e.g. $n=38$: 3033474044489 → 3033474044488; $n=48$: 3860228252031853 → 3860228252031852; $n=60$: 19517314042273052379 → 19517314042273052378. At $n=14$, $25$, $100$, $200$ (and every dim not listed) the quantisation gives nothing, so the plain cap remains the best this lemma can do.

**Update (job `lab_erdos_261117_5`; criterion files `j261117_5_criterion.py` / `.json` in this job directory):** the list *does* have a one-line closed description — $\mathrm{cap}^*=\mathrm{cap}-1$ exactly at the even $n$ with $S_{\mathrm{pl}}(n)\equiv2\pmod4$ — checked against a brute-force minimisation at every $n=3,\dots,600$ (598/598 agree), where it yields 131 sharpened dimensions, 125 of them open dimensions $>26$. The table stays useful as the citable object; the criterion is what lets a referee recompute it.

**Why the criterion holds (exact proof sketch, checked numerically here):** for $n\ge4$, $N(a)=E+a(a+2)2^{n-1}\equiv0\pmod 8$, so Lemma 12's congruence $\Delta\equiv0\pmod8$ becomes $WS\equiv0\pmod 8$ with $W=2a+2$: $a$ even ($W\equiv2\pmod4$) forces $S\equiv0\pmod4$; $a\equiv1\pmod4$ ($W\equiv4\pmod8$) forces only $S$ even, which always holds; $a\equiv3\pmod4$ ($W\equiv0\pmod8$) is vacuous. For even $n$ every admissible $a$ is even, so the true $S$ must be $\equiv0\pmod4$: if $S_{pl}\equiv2\pmod4$ the cap drops by exactly 1, and if $S_{pl}\equiv0\pmod4$ the congruence already holds at $S_{pl}$ and nothing drops. For odd $n$ every admissible $a$ is odd, the congruence is vacuous, and the quantisation can never bite (hence $n=25$ and every odd $n$). Independently re-checked in this job by `criterion_check.py` (all 198 rows of the range table) and `criterion_check_600.py` (598 rows, $n=3..600$, `S_max` by naive downward scan, no shared code with the main script): **598/598 agree**, 131 sharpened dimensions, 125 of them open ($>26$) — matching the paragraph above exactly.

> **One line, as the card requires:** this table is a **corollary of a published lemma** (Minamoto Lemma 12) evaluated at dimensions the paper does not evaluate; it is **not** past the published frontier of Erdős #86.

## 6. Verification (three independent legs)

| leg | what it is | result |
|---|---|---|
| `fresh_mod8_table.py` (this job, main) | from-scratch integer implementation: closed-form CRT solution of the congruence, exhaustive $(x,y)$ certificate search, 92 gates G1–G7 | **92 PASS, 0 FAIL**, exit 0 (`gate_log.txt`) |
| `verify_caps_independent.py` (this job, verifier) | no shared code: $S_{\max}$ by a **naive downward scan** (no CRT algebra), plain cap by cross-multiplied fractions, **exhaustive multiset enumeration** of the programme optimum at $n=4,5,6,7,8$, independent re-arithmetic of all seven certificates, $\Delta$ identity both ways, full re-computation of the 198-row range | **32 PASS, 0 FAIL**, exit 0 (`verify_log.txt`) |
| `../erdos-86-two-point-cap-slack-parity/gate_slack_parity.py` (job 261114_3) | the *other* implementation, re-run unmodified | **GATES PASS**, seven open caps identical |

Gate summary: reproduces all **17 published exact rows** (cap ≥ $g$ everywhere; equality at the 16 attained dimensions $3,4,5,6,7,8,10,11,13,15,16,17,19,21,23,26$; $n=9$ keeps its realisability slack of 4); reproduces Minamoto's three published uses of Lemma 12 ($S_{\max}(6,2)=72$, $S_{\max}(7,1)=160$, $S_{\max}(8,2)=340$); the exhaustive enumeration at $n=4..8$ confirms the relaxed optimum equals $g(n)$ there, which is what licenses reading $S_Q$ as the programme optimum at the open dimensions.

## 7. CLASSIFY

**`reproduction`** — a known mechanism (published Lemma 12) re-derived and tabulated: the sharpened caps 15700 / 1454898 / 6396312 / 27892120 / 120795954 at $n=12,18,20,22,24$, plus the same corollary at 30 further open dimensions ($n\le200$) that were not tabulated before. Not `extension`: no published value is passed (the plain caps 15701 etc. remain the published numbers; the sharpened ones are corollaries, and the lab's own $g(12)\ge15614$ bracket is untouched). Not `new`, no `needs_review` claim, **no Lean 4 statement is attached** because nothing beyond a corollary of a published lemma is claimed. No field is constructed, no $g(n)$ is determined, no prize is claimed, nothing submitted.

## 8. NEXT IDEAS (max 2)

1. **A third-moment (DPTV-type) identity at $n=12$** — the only route that would move $\frac{E+S_Q}{2}$ itself: a positive-definite quantisation of $\sum_v (h-a)(h-a-2)(h-a-4)$ or a $\sum_v \binom{h-a+1}{3}$-type relation that is $\ge 0$ and highly divisible, summed to give a bound below 15700. *value 8, odds 0.12, cost 1 job.* **Kill:** 5 min — check that no 3rd-moment identity appears in DPTV79 / Minamoto / Wang and test the natural candidate forms against the 16 attained dimensions: if any candidate is violated by a known exact $g(n)$, the family is dead. *Verifier:* a fresh implementation reproducing all 17 published rows and printing a value $<15700$ at $n=12$ that is still $\ge 15614$ (the lab's certificate).
2. **Field realisability of $S_Q$ at $n=20$** — $n=20$ is the smallest open dimension where the certificate uses $h\in\{2,4,6\}$ (three values), unlike the two-value $n=12$ case: a local-search + exact-cover hybrid on the signing model could decide whether cap\* is attained there, which would *determine* $g(20)$ if it hits. *value 7, odds 0.1, cost 1 job.* **Kill:** 5 min — check that the $n=20$ certificate's $h$-profile has no structural obstruction (e.g. $h\equiv2$ vertices forming a cut) and that the search engine can score a random signing's $\sum h$ at $n=20$ in under 2 min. *Verifier:* an independent signing scorer plus the gate that it reproduces $\sum h$ and $\sum h^2$ of $g(20)$'s published analogue at $n=19$.

---

## 9. Second verification pass and the closed criterion (job `lab_erdos_261117_5`)

Re-run on 2026-11-17 (job `lab_erdos_261117_5`, continuing this directory after the earlier run of job `lab_erdos_261117_2`). Three things were done: the three legs were re-run unmodified from this directory, the frontier was re-fetched live, and the one gap the table left open — *which* dimensions sharpen — was closed.

**9.1 Live re-fetch (2026-11-17, this job).** `https://www.erdosproblems.com/86` still states the problem as open with the $100 prize and the same citations (Er91, BHN95, BHLL14, Ba12b); no proof claim and no comment is visible on the page. `https://arxiv.org/abs/2603.29127` still resolves to Minamoto, *C4-Free Subgraphs of the Hypercubes Q6, Q7 and Q8: Odd Squares, Fully Frustrated Models, and Computational Structure* (v5, 28 Aug 2026, 75 KB), whose abstract states the contribution as "a mod-8 refinement of the classical two-point inequality" and whose scope is $Q_6,Q_7,Q_8$ — i.e. the lemma is published, its application at $n=12,18,20,22,24$ is not.

**9.2 Closed criterion for the sharpened set.** Write $m=2^{n-1}$, $E=nm$, and for admissible $a$ ($a\equiv n \bmod 2$, $0\le a\le n$) put $W=2a+2$, $N(a)=E+a(a+2)m$. Let
$$\rho(n)=\min_{a}\frac{N(a)}{W},\qquad S_{\mathrm{pl}}(n)=2\left\lfloor \rho(n)/2\right\rfloor$$
(the largest even $S$ the plain DPTV79 bound already allows). Then
$$\boxed{\ \mathrm{cap}^*(n)=\mathrm{cap}(n)-1 \iff n\ \text{even and}\ S_{\mathrm{pl}}(n)\equiv 2 \pmod 4.\ }$$

*Why.* For $n\ge4$, $N(a)\equiv0\pmod 8$ for every admissible $a$ (checked for all $n=4,\dots,600$ and all $a$: 0 violations; $E=n2^{n-1}$ and $m=2^{n-1}$ are both divisible by 8 there). So Lemma 12's $\Delta(a)\equiv0\pmod 8$ reads $W S\equiv0\pmod 8$, i.e. by the class of $a\bmod 4$:
- $a\equiv3\pmod4$: $W\equiv0\pmod8$ → **vacuous** (no restriction at all);
- $a\equiv1\pmod4$: $W\equiv4\pmod8$ → $S\equiv0\pmod2$, exactly DPTV parity → **no loss**;
- $a$ even: $W\equiv2\pmod4$ → $S\equiv0\pmod4$ → **the only restriction that can bite**, and it costs exactly 2 when $S_{\mathrm{pl}}\equiv2\pmod4$.

Every admissible $a$ admits $S_{\mathrm{pl}}$ as a plain value (it is the even rounding of the *minimum* of $N/W$), so the minimum over $a$ of the quantised optimum is $S_{\mathrm{pl}}-2$ precisely when at least one even $a$ exists ($n$ even) and $S_{\mathrm{pl}}\equiv2\pmod4$ — otherwise $S_{\mathrm{pl}}$. This also explains both non-sharpened open dimensions: $n=25$ is odd (no even $a$ exists), and at $n=14$ the binding $a=2$ gives $S_{\mathrm{pl}}=30036\equiv0\pmod4$.

**9.3 Checks (all exact integers, `j261117_5_criterion.py`, exit 0, 22/22 gates).**

| gate | what it checks | result |
|---|---|---|
| P1 | the card's kill test, by brute-force descent (no CRT algebra) | $S(12)=6824$ with $\{0{:}2,\,2{:}680,\,4{:}1366\}$; two extremal multisets exist ($a=0$: that one; $a=2$: $\{2{:}686,4{:}1360,6{:}2\}$); $\mathrm{cap}^*=15700$ |
| P2 | 17 published exact rows | $\mathrm{cap}^*\ge g$ everywhere, equality at the 16 attained dims, $n=9$ slack $+4$; the seven open caps identical to leg B's CRT route |
| P3a | $N(a)\equiv0\pmod8$ for $n\ge4$ | 0 violations over all $n=4..600$, all $a$ |
| P3b | per-$a$ residue prediction ($n<60$, every $a$) | the mod-8 condition is $S\equiv0\pmod4$ iff $a$ even, $S\equiv0\pmod2$ iff $a\equiv1$, vacuous iff $a\equiv3$ — 0 violations |
| P3c | **criterion vs brute-force minimisation, $n=3..600$** | **598/598 dimensions agree**, 0 disagreements |
| P4 | sharpened set = predicted set; gain $\le1$; the five open dims in, 14 and 25 out | 131 sharpened dims $\le600$, max gain 1 |
| P5 | ten worked rows at large $n$ | e.g. $n=38$: 3033474044489 → 3033474044488; $n=100$, $200$, $400$, $599$ unchanged; $n=600$: …901 → …900 |

Range: the criterion gives **131 sharpened dimensions for $3\le n\le600$, of which 125 are open dimensions $>26$** (the first 36 coincide with the 36 listed in §5 for $n\le200$; e.g. $n=202,204,206,208,210$ are sharpened while $n=200$ is not). Because the criterion depends only on the residue of one integer per dimension, a referee can recompute any row of the table, or any $n$ beyond 600, without running the search.

**9.4 Three legs re-run unmodified (this job's own dated logs).**

| leg | script | result |
|---|---|---|
| A | `legA_gate_slack_parity.py` (job 261114_3's gate, verbatim) | `GATES PASS`, rc 0 (`j5_legA.log`) |
| B | `fresh_mod8_table.py` (CRT / closed-form, 92 gates) | `FAILURES: 0`, rc 0 (`j5_legB.log`) |
| C | `audit_mod8_caps.py` (independent audit) | `FAILURES: 0`, rc 0 (`j5_legC.log`) |

All three print the seven open caps identically (15701→15700, 72362→72362, 1454899→1454898, 6396313→6396312, 27892121→27892120, 120795955→120795954, 251658240→251658240) and leg C independently matches the lab's earlier bracket upper end 15700 at $n=12$.

**Classification unchanged: `reproduction`** — the caps are a corollary of Minamoto Lemma 12, the criterion is an elementary consequence of the same lemma plus $N(a)\equiv0\pmod8$. Neither is past the published frontier of Erdős #86: the best *published* upper bound at $n=12$ remains the plain 15701, and no $g(n)$ is determined. No proof is claimed, nothing submitted.

---

## Files
- `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/fresh_mod8_table.py` — main from-scratch implementation + 92 gates + certificate search.
- `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/verify_caps_independent.py` — independent verifier (naive scan, exhaustive enumeration, $\Delta$ identity both ways).
- `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/cap_table.json` — the table (24 rows), the seven certificates, and the 198-row range extension.
- `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/cap_table.md` — the citable markdown table.
- `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/gate_log.txt`, `verify_log.txt` — gate transcripts.
- `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j261117_5_criterion.py` — job 261117_5: the closed criterion for the sharpened set, checked at $n=3..600$ (22 gates), plus the card's kill test by brute-force descent.
- `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j261117_5_criterion.json`, `j5_criterion.log` — its numbers and transcript.
- `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/j5_legA.log`, `j5_legB.log`, `j5_legC.log` — this job's re-run of the three legs.
- `lab_results/lab_erdos_261117_2.json`, `lab_results/lab_erdos_261117_5.json` — result files.
