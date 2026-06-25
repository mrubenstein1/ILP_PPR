# PROVENANCE — diagnostic scripts and the production decisions they justify

**Purpose.** This note records *why* the active pipeline (`07`/`08`/`09`) is configured
the way it is. Each decision now baked silently into production was reached through a
diagnostic script archived here. If a committee member or reviewer asks "how do you
know the time cap is harmless / that the slowness wasn't a modelling error / that the
stationary scenario is a valid null?", the answer is the evidence catalogued below.

**Archive contents.** This folder holds the superseded versions of the active scripts
and the diagnostics used to settle the runtime investigation:

- `07__ilp_core.R`, `08__run_models.R`, `09__validation.R` — *pre-revision* versions
  (no solver time cap; no stationary-anchor correction; exact-equality validation
  tolerances). Superseded by the revised versions in the project root.
- `10__diagnose_stationary_gap.R` … `14__run_models_capped.R` — diagnostics (below).
- Diagnostic outputs (CSV/RDS) — the data that actually substantiates the claims.
- Related records (recommended to archive alongside): `P3_runtime_diagnostic_ledger.md`
  (the running diagnostic log), and `test_tie_break.R` (a standalone alternate-tie-break
  experiment — see note under script 12).

**Status convention.** "Backed by" = directly evidenced by an archived artefact.
"Confirmed forward by" = the diagnostic motivated the change; final confirmation comes
from a clean run of the *revised* pipeline (not from the archive).

---

## Decisions now in production, and their evidence

| # | Production decision (where it lives) | Backed by | Key evidence artefact |
|---|---|---|---|
| 1 | Cap each solve at 60 s (`07`, `SOLVER_TIME_LIMIT`) | 11, 13, 14 | `Q4_pick_stability_ladder.csv`; `solver_convergence_capped.csv` |
| 2 | The cap is harmless to the reported result | 13, 14 | `Q4_pick_stability_ladder.csv`; `model_results_capped.csv` |
| 3 | The slow runtime was confined to the stationary null / shared 2020 solve | 14 | `solver_convergence_capped.csv`; `solve_log_capped.csv` |
| 4 | A coefficient / tie-break reformulation was **rejected** as the fix | 12, (test_tie_break.R) | console output recorded in the ledger |
| 5 | Headline: rolling > myopic ~13–15% of value_added, > greedy ~43–46% | 14 | `model_results_capped.csv` |
| 6 | Stationary-anchor fix: ε hazard at 2020 (`07`, `build_scenario_matrices`) | 14 | `model_results_capped.csv` (33% null gap = artefact) |
| 7 | Control runtime via wall-time, not by loosening `MIPGap` | 13, 14 | `Q4_pick_stability_ladder.csv` |
| 8 | Validation tolerances made cap-coherent (`09`, `VA_TOL`; reformulated P3) | 13, 14 | degeneracy characterised in both |

Decisions 6 and 8 are **confirmed forward by** the revised `09` (P3a/P3b) and the
stationary row of the revised `08`: the archive shows *why* the change was needed; the
new run shows the null collapsing.

---

## Per-script record

### 10 — `diagnose_stationary_gap.R`
**What it did.** A single instrumented capped solve on the stationary scenario.
**Finding.** Closed in ~29 s and looked fine — but it used the **rolling** objective,
not the myopic frozen-belief objective, i.e. the wrong target. 
**Why archived (not deleted).** Cautionary: it documents the misfire that the runtime
problem is invisible unless you solve the *myopic* instance. This is the error that
script 14's "different source" tie-break experiment later repeated. Output was console
only; the key numbers are recorded in `P3_runtime_diagnostic_ledger.md`.

### 11 — `diagnose_myopic_vs_rolling_objective.R`
**What it did.** Compared the myopic frozen-belief objective against the rolling
objective on the stationary scenario, per period.
**Finding.** Localised the slowness to the **myopic** policy: the t=1 frozen-baseline
solve carried a large incumbent (~4.6e5) and stalled at the cap at a tiny gap, while
the t≥2 ε-regime solves had a near-zero, maximally degenerate objective. Established the
frozen-belief mechanism behind the degeneracy.
**Backs decision 1.** Output console only → see ledger.

### 12 — `test_tiebreak.R`
**What it did.** Added a symmetry-breaking perturbation to the (correct) instance and
re-solved.
**Finding.** The tie-break did **not** close the optimality gap.
**Backs decision 4** (tie-break/coefficient reformulation is not the fix). 
**Related:** the separately-uploaded `test_tie_break.R` (an "earlier-is-better" ε added
inside `compute_value_vector`) is a second, independent tie-break attempt from a
different source; it tested the **rolling** instance (repeating script 10's target
error) and made runtime *worse*. Archive it with this group as corroborating evidence
that the timing axis is orthogonal to the binding-budget, which-parcel degeneracy.

### 13 — `test_Q4_pick_stability.R`  ★ load-bearing
**What it did.** The "Q4" experiment: solved the myopic t=1 stationary instance across an
effort ladder (5 s / 15 s / 60 s / looser gaps / a long reference), compared the
**implemented (period-1) pick set** across configs, and scored each on the true future;
plus a reproducibility block (pinned vs unpinned threads).
**Finding.** The implemented picks are **not** bit-invariant to effort, but realised
value converges smoothly (3.11% → 1.53% → 0.30% of value_added at 5/15/60 s) and the
60 s and longer solves agree. The contested picks are intrinsic fine structure (not
removable exact ties), and the frozen-belief `MIPGap` is an *optimistic* proxy for
realised decision error.
**Backs decisions 1, 2, 7.** Re-runnable: re-run this if the cap, budget, δ, or dataset
change, to re-certify that the capped decision is stable.
**Outputs:** `Q4_pick_stability_ladder.csv`, `Q4_pick_stability.rds`.

### 14 — `run_models_capped.R`  ★ load-bearing — became the new `08`
**What it did.** Full capped run of all three policies across all five scenarios, with
per-solve convergence logging. The prototype whose logging is now folded into the
revised `08`.
**Finding (the headline).** rolling beats myopic by ~13–15% of value_added and greedy by
~43–46%, consistently across the four climate scenarios; on total J the gaps are ~0.09%
(why the "muted separation" impression arose). The cap binds in only 1–2 of the 9
myopic solves per climate scenario, each at ≤0.05% gap → the realised value_added gap
sits far above the ~0.3–0.5% solver floor. The stationary scenario showed a 33%
value_added gap — diagnosed as the 2020 baseline-hazard anchor artefact, motivating the
fix in decision 6.
**Backs decisions 2, 3, 5, 6.**
**Outputs:** `model_results_capped.csv`, `solver_convergence_capped.csv`,
`solve_log_capped.csv`.

---

## Two-line summary of the resolved question

The myopic solves were slow because the solver reached a near-optimal incumbent in
seconds and then spent the remaining time only *certifying* optimality across a band of
co-valued schedules — search was done, proof was not. Capping that proof is harmless to
the enacted policy (which implements one period and re-optimises), is conservative
(truncation only biases myopic down), and leaves the ~13% rolling-vs-myopic finding
intact; the only scenario that genuinely stalled was the stationary null, whose anchor
has now been corrected.
