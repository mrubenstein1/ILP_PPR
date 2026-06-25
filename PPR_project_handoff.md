# PPR Land Acquisition Under Climate Non-Stationarity — Project Handoff

*Context document for a future collaborator (human or AI) joining this project. It states
the research question, data, modeling approach, current state, the decisions that are
already settled (so they are not reopened), and the known limitations.*

**Last updated 2026-06-25**, after the first complete R/Gurobi run on real data. The
solver-runtime investigation that dominated the prior session is now **resolved**; the
headline comparison has been produced. See §4 for the current state and §8 for the
resolution record. Sections §1–§3 (question, data, locked design) are unchanged.

---

## 0. Orientation (read this first)

A graduate thesis in conservation ecology / spatial optimization (author: Madeleine
Rubenstein). It builds three land-acquisition decision models for the U.S. Prairie
Pothole Region (PPR) and compares them to answer one question: **does accounting for
projected future ecological conditions (non-stationarity) lead to better conservation
acquisition decisions than assuming present conditions persist?**

A few decisions are **locked and should not be reopened without explicit instruction**
(each was reached deliberately, several after substantial validation):

1. The objective is **landscape-total**, not a protected-portfolio objective. The
   portfolio framing was explicitly rejected as ecologically unrealistic.
2. The conversion-risk survival clock is **global from t = 0** (not forward-replanned).
   The forward alternative was tested and rejected because it breaks the optimality
   ordering (myopic can beat rolling).
3. Conversion risk is treated as **modular and provisional** — it enters the models
   through exactly one data column (`trans_prob`) and may be replaced by a new dataset
   or metric later. Model code must not bake in conversion assumptions elsewhere.

**Headline finding (new this session).** Foresight matters, and substantially. Across the
four climate scenarios the foresighted (rolling) policy beats the myopic policy by
**~13–15% of value_added** and the greedy status-quo heuristic by **~43–46%**. This
*overturns* the earlier expectation (former §4.3) that the models would be "nearly
indistinguishable." The separation is masked on total landscape J (where the gap is
~0.09%, because only ~60 of 841 parcels are ever acquired) and is visible on
**value_added = J − J_baseline**, the welfare acquisition actually creates. Read §4.3 for
the numbers and the mechanism.

---

## 1. Project and Research Question

### 1.1 The question
Conservation agencies acquire land (or perpetual easements) to protect waterfowl
breeding habitat. Classic acquisition prioritization assumes the landscape is
stationary — that today's habitat value and today's threats describe the future. Under
climate change this is false: habitat suitability, duck abundance, and conversion
pressure all shift over the century. The thesis asks whether a decision-maker who
plans against *projected* trajectories outperforms one who (a) re-optimizes each period
but assumes the present persists, or (b) follows a simple status-quo heuristic.

### 1.2 The three decision models
All three choose, period by period, which EAUs (spatial units, §2.1) to acquire under a
budget. They differ only in the information and criterion used:

- **Greedy (status quo).** Each period, rank available parcels by current benefit-to-cost
  ratio (`scaled_abundance / cost`) and buy down the list until the budget is exhausted.
  No foresight, no conversion-risk weighting, no future periods considered.
- **Myopic ILP (fabricated future).** Each period, freeze both benefit and hazard at the
  *current* period's values and assume they persist for the rest of the horizon; solve
  the optimal acquisition program under that (false) stationary belief; implement only
  the current period's purchases; advance one period and observe the true next state.
- **Rolling-horizon ILP (true future).** Same re-solve structure, but uses the *true*
  projected benefit and hazard trajectories. Implements only the current period, then
  re-solves. This is the foresighted benchmark — the best achievable under perfect
  foresight (see §3.4, P2).

### 1.3 Relationship to prior work
This replicates the *spirit* (not the exact structure) of a prior Markov Decision
Process toy-problem study by the same author (Rubenstein et al., *Conservation Biology*,
in review). That study found the models perform similarly under stationarity, but the
foresighted (optimal) model outperforms under non-stationarity, with the advantage
growing as volatility increases (reported gaps up to ~7.34% over myopic and ~14.42%
over greedy on the toy problem). The present work asks whether that signal survives at
realistic landscape scale with empirically derived inputs. **It does — and on the
decision-relevant metric (value_added) it is larger than the toy precedent:** ~13–15%
over myopic and ~43–46% over greedy (§4.3). **Important contrast:** the toy problem had a
small landscape with a large acquisition fraction, so differences showed up directly in
the total objective; the real model has ~841 parcels with a tiny acquisition fraction, so
the same signal is invisible on total J and must be read on value_added (§4.3).

---

## 2. Data: Sources and Structure

The spatial framework and data layers are built from real sources via a six-script
pipeline (`01__…R` through `06__…R`). The conversion-risk layer is provisional.

### 2.1 Spatial and temporal frame
- **Spatial unit: the EAU (Equal Area Unit).** Script 01 grids the USFWS Wetland
  Management District (WMD) shapefile into equal-area cells (~282 km² each; original
  gridding script credited to Heini Kujala). After dropping 4 WMDs, **841 EAUs across 20
  active WMDs** remain. Each EAU is the atomic acquire/don't-acquire unit.
- **Temporal frame.** Decadal steps **2020–2100 = 9 periods** (period index t = (year −
  2020)/10, so t ∈ {0,…,8}). 2020 is an observed anchor; 2030–2100 are projections.

### 2.2 Scenario structure
Futures are indexed by an emissions pathway (`rcp`) crossed with a climate-moisture
realization (`gcm`). Factor levels actually used in the data:
- `rcp ∈ {"baseline", "45", "85", "stationary"}`
- `gcm ∈ {"baseline", "wet", "dry", "stationary"}`

The five analysis scenarios are the four climate futures **(45,wet), (45,dry),
(85,wet), (85,dry)** plus a **stationary** scenario **(stationary,stationary)** that
holds 2020 habitat constant from 2030 onward (the null against which non-stationarity
is judged). Every scenario shares a single **2020 baseline row** (`rcp="baseline",
gcm="baseline"`) that is stitched on as period 0.

> **Note (new this session):** the shared-baseline stitch is correct for the *climate*
> scenarios but was wrong for the stationary null — the baseline 2020 hazard is not ε, so
> the stationary scenario was not actually stationary at the anchor year. This is now
> corrected at scenario-assembly time for the stationary scenario only (§8.3); the
> shared-baseline convention is otherwise preserved.

### 2.3 The `data_panel` object
The pipeline's master output (`input_data/data_panel.rds` / `.csv`) is **long format,
one row per EAU × year × rcp × gcm**: per EAU, one 2020 baseline row plus, for each
future decade, five rows (the four climate futures + stationary).

| Column            | Meaning                                              | Source script |
|-------------------|------------------------------------------------------|---------------|
| `eau_id`          | EAU identifier                                       | 01            |
| `wmd_id`          | parent Wetland Management District                   | 01/03         |
| `year`            | 2020, 2030, …, 2100                                  | 03            |
| `rcp`, `gcm`      | scenario coordinates (levels above)                  | 03            |
| `prop_suitable`   | proportion of EAU that is suitable habitat           | 02            |
| `abs_abundance`   | absolute projected breeding duck pairs               | 04            |
| `scaled_abundance`| **benefit b** — EAU-level abundance used by models   | 04            |
| `trans_prob`      | **hazard λ** — per-period conversion probability     | 05            |
| `cost`            | acquisition cost (USD), scenario-invariant           | 06            |

### 2.4 Per-layer provenance
- **`prop_suitable` (script 02).** From FOREsce land-cover projections under the PPR
  GCAM Reference scenario, for RCP 4.5 and 8.5; suitable-habitat fraction per EAU. (Note
  in the script: slight area misalignment between EAUs and the FOREsce land-cover grid.)
- **`abs_abundance` / `scaled_abundance` (script 04).** Duck abundance anchored to a 2020
  observational value per WMD and projected/interpolated to decadal steps, restructured
  by gcm/rcp/timestep. `scaled_abundance` is the benefit coefficient the models use.
- **`trans_prob` (script 05).** A **net-loss proxy** for conversion risk:
  `trans_prob = max(ε, (prop_suitable_t − prop_suitable_{t+1}) / prop_suitable_t)`. Special
  cases: 0 where there is no habitat or at the 2100 terminal period; ε (a small
  background-risk floor) under the stationary scenario; the 2020 baseline value is the
  mean of the RCP 4.5 and 8.5 2020→2030 transitions. ε defaults to one order of
  magnitude below the smallest observed positive loss rate and is a sensitivity knob.
  **This layer is provisional** (§3.7).
- **`cost` (script 06).** From the PLACES fair-market-value vacant-land raster (Nolte
  2020, *PNAS* 117(47):29577; Dryad doi:10.5061/dryad.np5hqbzq9), 2017 USD. Mean FMV per
  hectare per EAU × EAU area × (1.02)^(year − 2017). Scenario-invariant; varies only by
  EAU and year (2%/yr inflation).

### 2.5 Pipeline scripts (data construction, already complete)
`01` establish EAUs from the WMD shapefile; `02` land-cover/suitable-habitat
characterization; `03` assemble the data panel (incl. wet/dry duplication and the
stationary scenario); `04` duck-abundance benefit; `05` conversion-risk hazard; `06`
acquisition cost. These are written, run, and verified upstream; the modeling scripts
consume `data_panel`.

---

## 3. Modeling Approach and Key Assumptions

### 3.1 Objective — landscape-total (LOCKED)
Maximize the **time-discounted total breeding duck pairs on the whole landscape** over
2020–2100. An acquired parcel contributes its (climate-driven) abundance `b[i,t]` with
**certainty** for every period at or after acquisition — acquisition removes conversion
risk but the abundance still follows its climate trajectory. An unacquired-but-not-yet-
converted parcel contributes `b[i,t] · S[i,t]` in expectation, where `S` is its survival
probability. Total welfare decomposes as

&nbsp;&nbsp;&nbsp;&nbsp;`J(policy) = J_baseline + Σ_i V[i, τ_i]`

where `J_baseline` (the do-nothing landscape total) is a constant and `V[i,τ]` is the
incremental value of acquiring parcel i at period τ (§3.3). Maximizing the acquisition
term maximizes `J`.

**Why not a protected-portfolio objective.** An alternative formulation counts only
ducks on *protected* land (unprotected-but-extant habitat is worth nothing). The author
**explicitly rejected this** as ecologically unrealistic and uninteresting: ducks on the
landscape are valuable whether or not the parcel is in the acquisition portfolio. Do not
reintroduce the portfolio objective without explicit instruction.

### 3.2 Decision program (ILP)
Binary `y[i,τ] = 1` iff parcel i is acquired at period τ.

&nbsp;&nbsp;&nbsp;&nbsp;maximize&nbsp;&nbsp;`Σ_{i,τ} V[i,τ] · y[i,τ]`
&nbsp;&nbsp;&nbsp;&nbsp;s.t.&nbsp;&nbsp;`Σ_τ y[i,τ] ≤ 1` (each parcel acquired at most once)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`Σ_i cost[i,τ] · y[i,τ] ≤ B_τ` (per-period budget, no rollover)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`y ∈ {0,1}`

### 3.3 Acquisition value and survival
- **Survival (global from t = 0):** `S[i,t] = Π_{t'=0}^{t-1} (1 − λ[i,t'])`, with `S[i,0]=1`.
  The terminal-period hazard never enters S (no subsequent period; `trans_prob` is 0 there
  by construction).
- **Acquisition value (expected loss prevented):**
  `V[i,τ] = Σ_{t=τ}^{T} δ^t · b[i,t] · (1 − S[i,t])`.
  This is the welfare gained by removing conversion risk from parcel i starting at τ.
- **Discounting** is absolute present value to 2020: `δ^t` with **δ = 0.95 per decade**
  (so 2100 is worth 0.95^8 ≈ 0.66). δ is a free parameter for sensitivity analysis. Note
  that within any single per-period re-solve, relative vs. absolute discounting does not
  change the chosen acquisitions (uniform scaling preserves the argmax), so absolute
  discounting is used for consistency with evaluation.

### 3.4 The survival-clock decision (LOCKED — and subtle)
A central correctness question was whether, in each per-period re-solve, the survival
clock should **reset to the current period** (forward re-planning) or stay **global from
t = 0**. The resolution is **global for both the myopic and rolling models**, for two
reasons established both theoretically and by validation:

1. **Global survival makes the rolling policy the provable welfare optimum.** With global
   `V`, the rolling re-solve loop reproduces the single full-horizon optimum exactly
   (property P2), so rolling is the correct upper-benchmark and is never beaten.
2. **Global survival makes the stationary null exact.** Under a *genuinely* stationary
   scenario the myopic frozen-future belief is *correct*, so myopic and rolling pose the
   identical program every period (property P3). NB: this requires the stationary inputs
   to be flat across the whole horizon, including the 2020 anchor — see §8.3 for the fix
   that ensures this on the real data, and the reformulated P3 in §4.2.

The forward-replanning alternative (`reset_at = t_now`) was implemented and tested. It
**causes order violations** — myopic occasionally beats rolling — which is disqualifying,
since rolling is defined as the best-with-foresight benchmark. (Mechanically,
forward-conditional re-planning optimizes a slightly wrong objective, off by a factor of
1/S at the current period.) Hence it was rejected. Do not switch to forward re-planning.

### 3.5 The three policies, precisely
- **Greedy:** per period, sort available parcels by `b[i,t]/cost[i,t]`, acquire greedily
  within `B_t`. (Benefit-cost ratio — a *different* criterion from the welfare objective,
  so greedy can differ from rolling even under stationarity.)
- **Myopic:** per period `t_now`, set `b[i,t] := b[i,t_now]` and `λ[i,t] := λ[i,t_now]` for
  all remaining t; compute `V` under global survival with that frozen belief; solve over
  parcels still available and periods `t_now…T`; implement only the `t_now` purchases;
  advance and observe the true next period.
- **Rolling:** per period, use the true `b` and `λ` trajectories (global `V`); solve over
  available parcels and `t_now…T`; implement only `t_now`. Equivalent to one full-horizon
  solve under deterministic foresight (kept as a loop for fidelity and future stochastic
  extension).

### 3.6 "Option A" — conversion in the objective only
All unacquired parcels are treated as deterministically available at each re-solve;
conversion risk enters **only** through the `(1 − S)` weighting in the objective (and the
`S` weighting in evaluation), **not** as a constraint that removes converted parcels from
the choice set. This is justified because the **per-period** conversion probability is
small, so few parcels actually convert out of availability within a decade — even though
the *cumulative*, discounted loss-prevented value (value_added) is substantial enough for
the policies to separate clearly (§4.3). If a high-conversion dataset is adopted, revisit
whether converted parcels should be removed from availability.

### 3.7 Conversion risk is modular (LOCKED stance)
Conversion enters the models in exactly two places — `compute_value_vector` (the `1 − S`
term) and `evaluate_policy` (the `S` term) — both fed by `trans_prob`. **A change to the
conversion metric or dataset is a script-05 change only; no model code changes.** A
methods memo (`methods_memo_conversion_risk.md`, produced earlier) documents the debate
over the net-loss proxy vs. a vote-counting alternative (the two were found empirically
equivalent — median churn ratio 1.0 — so the existing net-loss method was kept), its
drawbacks, and a broader finding: in this system the dominant axis of non-stationarity is
**waterfowl redistribution** (abundance shifting across the landscape), not habitat loss
per se. This is now a candidate explanation for *why* the models separate as much as they
do despite low per-period conversion (§4.3, §6 item 5).

---

## 4. Current State

### 4.1 Modeling scripts (written and RUN on real data with Gurobi)
Three R scripts continue the pipeline numbering:
- **`07__ilp_core.R`** — shared engine: parameter block (δ, budget, scenario coding, and
  a **Solver control** block — see §8.2), `build_scenario_matrices` (with the stationary
  anchor correction, §8.3), `compute_survival_matrix`, `compute_value_vector`,
  `make_budget`, `solve_acquisition_ilp` (Gurobi wrapper: cost-scaling, a per-solve time
  cap, and an optional solve-log hook), the three policies + `run_full_horizon`,
  `evaluate_policy`, and synthetic-instance + brute-force helpers for validation.
- **`08__run_models.R`** — driver: runs all three models on all five scenarios, scores
  each against its scenario's true future, computes `J_baseline`, and writes results +
  per-period trajectories + a per-run solver-convergence summary (RDS + CSV). Has a
  `REPRODUCIBLE` toggle (single-thread, fixed seed) for bit-reproducible reported numbers.
- **`09__validation.R`** — correctness suite (run first, §4.2), plus a documented stub for
  the deferred MDP comparison.

### 4.2 Validation status (real data, Gurobi, capped solver)
The suite is run with the production time cap active, so its tolerances are now coherent
with how the policy is actually solved (see §8.2). The five properties:

- **P1** ILP optimum equals brute-force optimum (tiny synthetic instance; exact).
- **P2** rolling ≈ full-horizon single solve (real scenario; checked up to solver
  tolerance, since the single full-horizon solve is the hardest instance and may stop at
  the cap while the rolling loop converges).
- **P3** *reformulated.* The stationary null is now tested by its **solver-independent**
  content rather than bit-exact schedule identity (which a capped, degenerate null solve
  cannot satisfy even when correct): **(P3a)** the stationary benefit and (non-terminal)
  hazard trajectories are flat; **(P3b)** the myopic *frozen* value vector equals the
  rolling *true* value vector. The realized myopic-vs-rolling gap is reported for
  information only. *Watch P3a's benefit check:* a FAIL there would mean the stationary
  abundance trajectory is not flat — an upstream (01–04) data-construction matter, not a
  solver one (§6 item 6).
- **P4** rolling ≥ myopic on every scenario, checked up to solver tolerance (a margin of
  ~2% of value_added — above capped-solver noise, below any genuine order violation).
- **P5** LP relaxation ≥ ILP (holds for free under the cap: the LP bound dominates any
  feasible incumbent).

Solver is **Gurobi via the native R `gurobi` package**; CBC blows up beyond ~40 parcels,
which is why Gurobi is required at N = 841.

### 4.3 Key empirical finding (the headline)
Under landscape-total with the current data, **foresight provides a large, consistent
advantage on value_added** across all four climate scenarios:

| metric (climate scenarios) | greedy | myopic | rolling |
|---|---|---|---|
| value_added gap vs rolling | ~43–46% | **~13–15%** | 0 (reference) |
| gap on total landscape J    | ~0.3%  | ~0.09% | 0 |
| parcels acquired (of 841)   | ~60–64 | ~55–57 | ~65–67 |
| total spend                 | ~13.0B | ~10.5B | ~13.4B |

Two things to internalize:

1. **Read value_added, not total J.** Because only ~60 of 841 parcels are acquired,
   `J_baseline` dominates `J`, so the gap on total J is ~0.09% and looks negligible. The
   decision-relevant quantity is `value_added = J − J_baseline` (the welfare acquisition
   creates), where the ordering is large and clear. The driver reports both.
2. **The mechanism is the cost of delayed *and* under-deployed protection.** A myopic
   manager re-plans each decade and self-corrects, but it protects later than the
   foresighted manager and — notably — **systematically under-acquires and under-spends**
   (~55 parcels / ~10.5B vs rolling's ~65 / ~13.4B). Under its frozen belief it does not
   "see" enough value to deploy the full budget on parcels whose worth depends on a future
   it cannot anticipate. Greedy diverges more still, because its benefit-cost criterion
   differs from the welfare objective.

**This supersedes the prior "nearly indistinguishable" expectation.** That impression came
from looking at total-J gaps and at near-null diagnostics; on the climate scenarios with
value_added, the separation is robust and well clear of any solver-precision floor (§8.2).
A plausible driver of the magnitude is waterfowl **redistribution** (§3.7): even with low
per-period conversion, anticipating *where abundance moves* is worth a lot — a candidate
decomposition for the thesis (§6 item 5).

### 4.4 The stationary null
With the §8.3 anchor fix, the stationary scenario is a genuine null. Pre-fix it showed a
spurious ~33% myopic value_added gap caused entirely by the non-ε 2020 baseline hazard;
post-fix, value_added under stationarity is near zero (there is nothing to foresee) and
myopic and rolling coincide up to solver tie-breaking. **Confirm this on your run** via
P3a/P3b and the stationary row of `08`'s output.

---

## 5. Assumptions and Limitations (consolidated)

- **Conversion metric is provisional.** Net-loss proxy from habitat erosion; low in
  per-period magnitude; may be replaced. Methods memo notes redistribution (not loss) is
  the dominant non-stationarity signal — now a candidate explanation for the *size* of the
  observed separation, and a candidate direction for a better risk metric.
- **Deterministic expected-value evaluation.** Survival `S` is used as an expected weight;
  there is no Monte Carlo over realized conversion events. "Option A" treats unacquired
  parcels as always available (justified by small per-period λ, §3.6).
- **No spatial interactions.** Each EAU is independent — no adjacency, connectivity, or
  contiguity constraints, and no spatial spillover in benefit.
- **Acquisition is permanent and one-shot.** Acquire-once, no divestment; per-period
  budget with no rollover.
- **Greedy uses benefit-cost ratio, not the welfare objective** — an intentional
  status-quo baseline, not a "welfare-optimal under wrong beliefs" model the way myopic is.
- **Reported numbers carry a small solver-precision band.** Each solve is capped at 60 s
  (§8.2); on the climate scenarios this leaves <~0.5% of value_added of uncertainty on the
  myopic numbers, far below the ~13% effect, and it biases myopic *down* (so the comparison
  is conservative). For exact reproducibility, run `08` with `REPRODUCIBLE = TRUE`
  (single-thread).
- **Exact MDP benchmark is intractable.** A joint-landscape MDP has ~3^841 states (§ stub
  in `09`).
- **Stationary abundance flatness is asserted, not yet confirmed on this machine.** The
  anchor fix flattens the stationary *hazard*; P3a additionally checks that *benefit* is
  flat. If it is not, the null is not exact and the cause is upstream (01–04).

---

## 6. Open Items / Next Steps

1. **Confirm the run.** `09__validation.R` (Gurobi) → P1–P5, with attention to P3a
   (stationary benefit flatness) → `08__run_models.R` → inspect `model_results.csv`, the
   trajectories, and `solver_convergence.csv`.
2. **Lock reportable numbers.** Run `08` with `REPRODUCIBLE = TRUE` for the figures that
   go in the thesis; report value_added gaps to ~2 significant figures (≈13% myopic,
   ≈45% greedy) with the solver band noted.
3. **Sensitivity analyses.** δ (discounting), `BUDGET_EAUS_PER_PERIOD` (default 5, range
   2–10), and ε (conversion floor). Parameter blocks are at the top of `07` and `05`.
4. **Conversion data / metric.** Possibly adopt a redistribution-aware or higher-conversion
   risk layer; absorbed through `trans_prob` alone.
5. **Decompose the foresight signal (thesis enrichment).** Myopic freezes *both* benefit
   and hazard. Re-run freezing only `b` vs only `λ` to attribute the ~13% to anticipating
   abundance redistribution vs habitat loss. If redistribution dominates, that is a sharper
   and more novel framing than "conversion risk."
6. **Stationary abundance flatness (if P3a fails).** If the stationary benefit trajectory
   is not flat, decide upstream (01–04) whether the stationary scenario should hold 2020
   abundance constant; this is a data-construction decision, not a model change.
7. **MDP comparison (deferred).** Interface stub documented in `09`.

---

## 7. Thesis framing

Communicate the result precisely: model choice matters in proportion to how much
foresight changes the *timing and amount* of protection, and that effect is large here
(~13% myopic, ~45% greedy on value_added) even though it is invisible on the landscape
total. The value of foresight is the value of protecting at-risk, value-shifting parcels
*earlier and in greater number* than a present-assuming manager will. The honest caveats
are the provisional conversion metric and the small solver-precision band on the capped
solves; neither threatens the finding.

---

## 8. Resolution Record — First Full R/Gurobi Run (2026-06-24 → 2026-06-25)

*This section replaces the prior "session addendum." It records what was diagnosed and
decided. The blow-by-blow diagnostic trail lives in `myopic_policy_diagnostic.md` (the
former `P3_runtime_diagnostic_ledger.md`) and `archive/PROVENANCE.md`.*

**Environment.** Gurobi 13.0.2 (academic license), R `gurobi` package, Apple M2 Pro /
16 GB / macOS 15.7.7. `input_data/data_panel.rds` consumed successfully.

### 8.1 The runtime problem, and its resolution
**Symptom (prior session).** The myopic solves on the stationary scenario reached an
excellent incumbent in seconds, then took many minutes (the first real solve ran 866 s,
~18M nodes) *proving* optimality.

**Diagnosis (this session, measured).** The slowness was the **myopic policy, not the
scenario** (rolling closes the same instance in ~28 s); it was the **optimality-proof
phase, not search** (the LP relaxation is tight, ~0.19%; the incumbent saturates early
while only the bound descends); and it was driven by a large band of **co-valued
schedules** produced by the flattened frozen-belief objective. A pick-stability experiment
(the "Q4" ladder) showed the implemented period-1 picks converge smoothly with effort
(realized error 3.11% → 1.53% → 0.30% of value_added at 5/15/60 s) rather than being
exact, removable ties; and the climate scenarios' later-period solves close quickly. The
earlier "exact ties vs intrinsic near-flatness" fork that organized the diagnostic ledger
was thereby rendered **moot**: tie-break/coefficient reformulation was rejected (it did not
close the gap and, in one variant, made runtime worse), and the resolution is policy-level.

**Decision (in production).** Cap each solve at **60 s** (`SOLVER_TIME_LIMIT` in `07`).
This is justified because (i) each policy implements only the current period and
re-optimizes, so the expensive certification of the discarded tail never affects the
enacted trajectory; (ii) the cap binds in only 1–2 of the 9 myopic solves per climate
scenario, each at ≤0.05% gap, leaving <~0.5% realized uncertainty on value_added against a
~13% effect; and (iii) truncation only biases myopic *down*, so the rolling-vs-myopic
comparison is conservative under the cap. `MIPGap` was **not** loosened (the bound is
already tight; the lever is wall-time, not gap).

### 8.2 Solver control block (`07`) — KEEP
The Gurobi call now reads its parameters from a documented **Solver control** block:
`SOLVER_TIME_LIMIT = 60`, `SOLVER_MIP_GAP = 1e-4`, `SOLVER_THREADS = 0L` (set 1L for
reproducibility), `SOLVER_SEED = 1L`, `SOLVER_OUTPUT = 0L` (quiet). The earlier
**cost-scaling** fix is retained: inside `solve_acquisition_ilp`, costs and the budget RHS
are divided by `cost_scale = 1e6` before Gurobi (identity-preserving; only improves
conditioning; reporting stays in raw USD). A small, inert **solve-log hook** records each
solve's status/gap/runtime when a driver creates a `.SOLVE_LOG` (used by `08`'s
convergence summary); it is otherwise a no-op. *Do not re-add `NumericFocus`/`ScaleFlag`*
(tried last session, reverted — they slowed every solve for no benefit once data was
scaled).

### 8.3 Stationary null correction (`07`, `build_scenario_matrices`) — KEEP
Prior session identified that the stationary scenario inherited the shared 2020 baseline
hazard (the mean of floored RCP45/85 2020→2030 transitions, > ε) instead of ε, making it
non-stationary at the anchor and producing a spurious ~33% myopic gap. Fix: for the
stationary scenario only, set the 2020 hazard to the scenario's own ε floor
(`lam[,1] <- lam[,2]`). This is applied at matrix-assembly time and **resolves the §8.5
design tension cleanly** — it does not touch script `05` or the data panel, and it
preserves the shared-baseline convention for the climate scenarios (which correctly use the
baseline 2020 hazard). It flattens the *hazard*; the *benefit*-flatness precondition is now
certified separately by P3a.

### 8.4 Validation reformulation (`09`) — KEEP
With the cap active, the original `1e-6`/bit-exact tolerances were incompatible with how the
model is solved. P3 was reformulated to the solver-independent null test (§4.2); P2 and P4
were given tolerances coherent with the cap (≈ value_added scale); P1 and P5 were unchanged.
The realized stationary myopic-vs-rolling gap is now reported as information, not asserted.

### 8.5 What changed vs the prior plan
The prior session's recommended next lever (`MIPFocus = 3` / `Cuts = 2` to attack the
bound) was **not** adopted: the pick-stability and convergence measurements showed a time
cap is provably harmless and sufficient, so the simpler, better-justified fix won. The
prior worry that loosening precision could "contaminate the muted effect" was also
re-examined and found to conflate the solver's internal gap with the reported value_added;
the cap controls wall-time, not gap, and the headline effect is an order of magnitude
larger than the residual either way.
