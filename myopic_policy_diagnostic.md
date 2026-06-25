# Myopic Policy Runtime — Diagnosis & Resolution

*Formerly `P3_runtime_diagnostic_ledger.md`. Renamed because the problem was never
specific to the P3 validation check — it was a property of the **myopic policy's solves**,
which surface in P3, P4, and `08`. This document is now a **resolved** record: it states
what the problem was, how it was diagnosed, what was decided, and why the decision is
safe. The live-investigation ledger it grew out of (the KNOWN / RULED-OUT / INTUITED /
open-fork format) is preserved in `archive/` for the full blow-by-blow.*

**Status: RESOLVED (2026-06-25).** Production fix is a 60 s per-solve time cap in
`07__ilp_core.R` (`SOLVER_TIME_LIMIT`). The companion correctness issue (the stationary
null) is fixed by the 2020-hazard anchor correction in `build_scenario_matrices`. See the
handoff §8 for the decision record and `archive/PROVENANCE.md` for the evidence map.

Evidence tags used below: `[MEASURED]` (a diagnostic run, numbers cited), `[FROM LOG]`
(real Gurobi console output), `[FROM SOURCE]` (a property of scripts 05/07), `[RESOLVED]`
(a former hypothesis now settled by measurement).

---

## 0. The two problems (both now resolved)

They were always independent; conflating them caused early confusion.

- **Problem A — runtime.** The myopic solves found an excellent answer in seconds, then
  spent minutes *proving* it optimal. **Resolved:** time cap (§3). It was the *policy*, not
  the scenario, and it was *certification*, not search.
- **Problem B — correctness of the stationary null.** P3's bit-exact "myopic == rolling"
  check could fail on real data because the stationary 2020 hazard was the baseline mean,
  not ε `[FROM SOURCE: 05]`. **Resolved:** anchor fix in `07` makes the stationary hazard
  flat, and P3 was reformulated to a solver-independent test (§4).

---

## 1. What was known (measured facts — still the foundation)

### 1a. The instance
- The MIP has **7,569 binary variables** (841 EAUs × 9 periods) and **850 constraints**
  (841 acquire-once + 9 per-period budget knapsacks). `[FROM LOG]` By solver standards
  this is small — size alone never explained the runtime.
- The objective the solver maximizes **is value_added** (the Σ V·y term); `J_baseline` is
  added later in evaluation. `[FROM SOURCE: 07]`

### 1b. The slowness was the *policy*, not the *scenario*
This was the single most important correction to the original framing.
- On the **stationary** scenario, **rolling** closes to OPTIMAL in **28.4 s** (gap
  0.0063%). `[MEASURED: script 11]`
- On the **same** scenario, **myopic** solves hit the cap at tiny gaps. In the real run the
  first myopic solve ran **866 s / ~18M nodes** to 0.01%. `[FROM LOG]`
- Every slow solve observed was a **myopic** solve; the one fast solve on the same scenario
  was rolling. The scenario was not the discriminator; the policy was.

### 1c. The slowness was end-game *proof*, not *search*
- The **LP relaxation is tight**: root 465,194.9 vs integer 464,303.2 — a ~**0.19%**
  starting gap. `[FROM LOG]` An excellent incumbent appears almost immediately; the cost is
  closing the last sliver while the bound creeps down over millions of nodes.
- **Many distinct integer schedules exist at (near-)identical objective value** — the log
  reports several at 464,303, and independent solves returned *different* schedules at
  near-identical value. `[FROM LOG]` `[MEASURED: scripts 10/11]`

### 1d. The mechanism (from source)
- Myopic **freezes** benefit and hazard at the current period and projects them **flat**
  across the horizon, then solves the full remaining-horizon problem — so under its belief
  the hazard is constant, flattening the value coefficients. `[FROM SOURCE: 07]`
- Myopic **implements only the current period** and **discards** the periods-2…9 schedule,
  re-solving next decade. The expensive proof certifies a tail that is thrown away.
  `[FROM SOURCE: 07]`
- Under stationary, the true hazard was the **shared 2020 baseline** at t=0 (> ε), then **ε**
  from 2030 on. `[FROM SOURCE: 05]` So myopic-t1 froze a meaningful hazard (incumbent
  ~4.6e5) while myopic-t≥2 froze ε (objective crushed to [1e-7, 0.33]). `[MEASURED: 11/12]`

---

## 2. The diagnostic fork — and why it turned out not to matter

The live ledger organized everything around one question: is the co-valued band made of
**exact, removable ties manufactured by the ε floor (H1)** or **intrinsic near-flatness
(H2)** — because H1 implied a coefficient-level fix and H2 implied a policy-level one.

**This fork is now moot `[RESOLVED]`.** Two measurements dissolved it:

- **Q4 — pick-stability ladder `[MEASURED: script 13]`.** The *implemented* (period-1)
  myopic picks are **not** bit-invariant to solver effort, but they converge **smoothly**:
  realized value error fell 3.11% → 1.53% → 0.30% of value_added at 5 / 15 / 60 s, and the
  60 s and longer solves agree. Smooth convergence is the signature of **intrinsic
  near-flatness (H2)**, not of exact, permutation-style ties (H1) that a duplicate-collapse
  would remove. It also showed the frozen-belief `MIPGap` is an *optimistic* proxy for
  realized decision error (realized error ran larger than the internal gap).
- **Q3 — climate vs stationary `[MEASURED: script 14]`.** The myopic later-period solves
  *close fast* on the climate scenarios (real, parcel-varying hazard gives the solver value
  to discriminate on); only the stationary null (ε hazard → near-zero, maximally degenerate
  objective) stalls. The hard case is the scenario that matters least.

Because the implemented decision is governed by a converged-enough period-1 solve and the
expensive tail is discarded, the **distinction between H1 and H2 never needed to be
resolved** to fix Problem A. Tie-break / coefficient reformulation was therefore rejected:
`[MEASURED: script 12]` an in-objective tie-break left the solves at the cap and did not
close the gap, and a second, independent tie-break variant (a different source) made
runtime *worse* — corroborating that the contested degeneracy is a binding-budget,
*which-parcel* contest that a timing perturbation cannot break.

---

## 3. Resolution of Problem A — the time cap

**Decision:** cap each solve at **60 s** (`SOLVER_TIME_LIMIT` in `07`); keep `MIPGap` at the
default `1e-4` (the bound is already tight — wall-time, not gap, is the lever); run quiet
(`OutputFlag = 0`).

**Why it is safe — three independent reasons:**
1. **Irrelevance of the discarded tail.** Myopic implements only the current period and
   re-optimizes; the expensive certification spans periods 2–9 that are thrown away, so it
   cannot change the enacted trajectory. `[FROM SOURCE: 07]`
2. **Negligible realized impact.** `[MEASURED: script 14]` The cap binds in only 1–2 of the
   9 myopic solves per climate scenario, each at ≤0.05% gap. The resulting realized
   value_added uncertainty is <~0.5%, against a **~13%** rolling-vs-myopic effect.
3. **Conservative direction.** Truncation only ever biases the myopic incumbent *down*, so
   it can only *widen*, never manufacture, the rolling advantage. A small measured gap is
   thus strong evidence the true gap is at least that large.

**Reproducibility.** A capped, multi-threaded solve is not bit-reproducible run-to-run
(parallel tie-breaking). For reported numbers, run `08` with `REPRODUCIBLE = TRUE`
(single-thread, fixed seed) or confirm stability across two multi-threaded runs; the
~0.3%-scale wobble is far below the effect, so the result is reportable to ~2 sig figs
either way.

---

## 4. Resolution of Problem B — the stationary null

**Cause `[FROM SOURCE: 05/07]`.** `build_scenario_matrices` stitched the shared 2020
baseline row onto every scenario, so the stationary scenario's 2020 hazard was the
baseline mean (> ε), not ε. Myopic froze that non-ε hazard at the anchor and mis-targeted
the one decision (2020→2030) that carries essentially all of the stationary value_added —
producing a spurious **~33%** myopic gap. `[MEASURED: script 14]`

**Fix (in `07`).** For the stationary scenario only, set the 2020 hazard to the scenario's
own ε floor: `lam[,1] <- lam[,2]`. Applied at assembly time, this preserves the
shared-baseline convention for the climate scenarios and touches neither script `05` nor
the data panel.

**Validation reformulation (in `09`).** A capped, degenerate null cannot satisfy bit-exact
schedule identity even when correct. P3 now tests the **solver-independent** content of the
null: **(P3a)** stationary benefit and (non-terminal) hazard are flat; **(P3b)** the myopic
*frozen* value vector equals the rolling *true* value vector. The realized myopic-vs-rolling
gap is reported for information, not asserted. P2/P4 tolerances were set to value_added
scale (coherent with the cap); P1/P5 were unaffected.

**One open check.** The fix flattens the *hazard*. P3a additionally checks that *benefit*
is flat — the other precondition of the null. If P3a's benefit check fails, the stationary
abundance trajectory is not flat, which is an upstream (01–04) data-construction matter,
not a solver one. Confirm on the run.

---

## 5. What this produced (the payoff)

The investigation's real output was not just a runtime fix but the **headline result**,
made trustworthy by the convergence evidence above: across the four climate scenarios,
rolling beats myopic by **~13–15% of value_added** and greedy by **~43–46%**, with myopic
also systematically under-acquiring (~55 vs ~65 parcels) and under-spending (~10.5B vs
~13.4B). `[MEASURED: script 14]` The "nearly indistinguishable" expectation was an artifact
of reading total J (~0.09% gap) and the near-zero stationary null; on value_added the
separation is large and well clear of the solver-precision floor.

---

## 6. One-line summary

The myopic solves were slow because the solver reached a near-optimal incumbent in seconds
and then only *certified* it across a band of co-valued schedules spanning periods the
policy discards — search was done, proof was not. Capping that proof at 60 s is harmless to
the enacted policy, conservative, and leaves the ~13% effect intact; the only genuinely
slow scenario was the stationary null, whose anchor is now corrected and whose validation
test is now solver-independent. The "exact-ties vs near-flatness" fork that organized the
original investigation was rendered moot by the pick-stability and climate-vs-stationary
measurements.
