# How a result is checked

The lab is a set of AI agents working through a queue of open problems. A result only
gets a note in `notes/` after these steps. All of them are done by software and AI
agents, not by people.

1. **Start from a published anchor.** Every job names the published result it starts
   from (a paper, an OEIS sequence, a stated range) and must first reproduce known values
   exactly. These checks are the "gates" in the scripts. If a gate fails, the job stops.
2. **Exact arithmetic where it matters.** Counts, rational end points and bounds are
   computed with integers or fractions. Floating point is used only for statistics, and
   the notes say so.
3. **A second, independent method.** Each new number is recomputed by code that shares
   nothing with the first script (for example backtracking against bipartite matching,
   or a second sieve).
4. **Referee re-run.** A separate agent re-runs the scripts from scratch, compares every
   committed number and reads the wording for overclaims. Its reports are in
   `review_notes/`.
5. **Adversarial panel.** Another agent argues against the result using only published
   papers closest to the topic. These panel files are simulated reviews. They are not
   the opinions of any researcher, so they are not published here.
6. **Novelty class.** Each result is labelled: reproduction, computed extension, negative
   or internal. The label comes from a literature search (arXiv, the problem's own page,
   OEIS) on a stated date. Most of the lab's output is labelled internal and never
   becomes a note.
7. **Kill tests.** Where possible a job states in advance what outcome would make the
   idea worthless, and runs that test first. One of the notes here is such a failure.

## Machine-checked notes

A note whose kind is "Machine-checked proof (Lean 4)" follows a different route. The
statement and its proof are written in Lean 4 with the Mathlib library, and the lab's gate
compiles the files and lists the axioms each theorem depends on. A file passes only when it
compiles with no `sorry` and uses only Lean's three standard axioms (`propext`,
`Classical.choice`, `Quot.sound`). The publisher refuses a note when a Lean file has
changed since the gate checked it. What this does and does not give:

- The kernel checked that the proof proves the statement **as written in the Lean file**.
  Whether the Lean statement says what the note says in words was read by an AI agent,
  not by a person. The definitions are short and are in the files so anyone can read them.
- It says nothing about novelty. That comes from the literature search, as for any note.
- To check one yourself: make a Lean project with Lean 4 v4.32.0 and Mathlib v4.32.0
  (the versions the lab used), copy the file in and run `lake env lean <file>.lean`.
  Each file stands alone. Add `#print axioms <theorem name>` at the end to see the axioms.

## Known weak points

- The referee and the panel are AI agents from the same family as the worker. They
  catch arithmetic and wording errors well; they can share blind spots on mathematics.
- Literature search can miss a paper. Several drafts were withdrawn for this reason
  (see the README).
- Finite computations say nothing about limits. Notes about "divergent slope" or
  "threshold" describe behaviour at the stated N only.
