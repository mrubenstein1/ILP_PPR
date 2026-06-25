# PPR Land Acquisition Under Climate Non-Stationarity — Project Handoff

*Prepared as a context document for a future AI agent joining this project. It states
the research question, data, modeling approach, current state, and the decisions that
are already settled (so they are not reopened) along with the known limitations.*

---

## 0. Orientation (read this first)

A graduate thesis in conservation ecology / spatial optimization (author: Madeleine
Rubenstein). It builds three land-acquisition decision models for the U.S. Prairie
Pothole Region (PPR) and compares them to answer one question: **does accounting for
projected future ecological conditions (non-stationarity) lead to better conservation
acquisition decisions than assuming present conditions persist?**

A few decisions are **locked and should not be reopened without explicit instruction**
(they were each reached deliberately, some after substantial validation):

1. The objective is **landscape-total**, not a protected-portfolio objective. The
   portfolio framing was explicitly rejected as ecologically unrealistic.
2. The conversion-risk survival clock is **global from t = 0** (not forward-replanned).
   The forward alternative was tested and rejected because it breaks the optimality
   ordering (myopic can beat rolling).
3. Conversion risk is treated as **modular and provisional** — it enters the models
   through exactly one data column (`trans_prob`) and may be replaced by a new dataset
   or metric later. Model code must not bake in conversion assumptions elsewhere.

A consequence worth internalizing early (Section 4.3): under the *current* conversion
data the three models separate only modestly, and this is expected, not a bug.

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
realistic landscape scale with empirically derived inputs. **Important contrast:** the
toy problem had a small landscape with a large acquisition fraction, so differences
showed up directly in the total objective; the real model has ~841 parcels with a tiny
acquisition fraction, which changes how (and where) the model differences are visible
(§4.3).

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
2. **Global survival makes the stationary null exact.** Under a stationary scenario the
   myopic frozen-future belief is *correct*, so myopic and rolling produce identical
   schedules and identical J (property P3, difference = 0).

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
the choice set. This simplification is justified because realized conversion is near zero
in the current data. If a high-conversion dataset is adopted, revisit whether converted
parcels should be removed from availability.

### 3.7 Conversion risk is modular (LOCKED stance)
Conversion enters the models in exactly two places — `compute_value_vector` (the `1 − S`
term) and `evaluate_policy` (the `S` term) — both fed by `trans_prob`. **A change to the
conversion metric or dataset is a script-05 change only; no model code changes.** A
methods memo (`methods_memo_conversion_risk.md`, produced earlier) documents the debate
over the net-loss proxy vs. a vote-counting alternative (the two were found empirically
equivalent — median churn ratio 1.0 — so the existing net-loss method was kept), its
drawbacks, and a broader finding: in this system the dominant axis of non-stationarity is
**waterfowl redistribution** (abundance shifting across the landscape), not habitat loss
per se. This is why conversion magnitude is low and provisional.

---

## 4. Current State

### 4.1 Modeling scripts (written; not yet run in R/Gurobi by the author)
Three R scripts continue the pipeline numbering and have been delivered:
- **`07__ilp_core.R`** — shared engine: parameter block (δ, budget, scenario coding),
  `build_scenario_matrices`, `compute_survival_matrix`, `compute_value_vector`,
  `make_budget`, `solve_acquisition_ilp` (Gurobi wrapper, isolated "SOLVER CALL" block),
  `run_greedy` / `run_myopic_ilp` / `run_rolling_horizon` / `run_full_horizon`,
  `evaluate_policy`, and synthetic-instance + brute-force helpers for validation.
- **`08__run_models.R`** — driver: runs all three models on all five scenarios, scores
  each against its scenario's true future, computes `J_baseline`, and writes a results
  table + per-period landscape trajectories (RDS + CSV).
- **`09__validation.R`** — correctness suite (run first, see §4.2), plus a documented
  stub for the deferred MDP comparison.

### 4.2 Validation status
The full modeling logic was validated in a **separate Python/PuLP solver harness**
(R and Gurobi were unavailable in the build environment), including a line-for-line
Python mirror of the R indexing to catch off-by-one/translation bugs. Five properties
pass:

- **P1** ILP optimum equals brute-force optimum (small instance).
- **P2** rolling re-solve loop equals the full-horizon single solve.
- **P3** stationary null: myopic equals rolling exactly (identical schedules, J diff 0).
- **P4** rolling ≥ myopic on every tested instance (no order violations).
- **P5** LP relaxation ≥ ILP (small integrality gap, ~0.04–0.5%).

**The author must run `09__validation.R` first with a live Gurobi academic license to
confirm these on the real data, then run `08__run_models.R`.** Solver is **Gurobi via the
native R `gurobi` package**; CBC (the Python stand-in) blows up beyond ~40 parcels, which
is why Gurobi is required at N = 841.

### 4.3 Key empirical finding (set expectations accordingly)
Under landscape-total with the **current low-conversion data, the myopic and rolling
models are nearly indistinguishable**, and this is a structural feature, not an error.
When little habitat is ever lost, an unacquired parcel still contributes ≈ `b` (because
`S ≈ 1`), so the acquisition schedule barely moves the landscape total. The myopic-vs-
rolling gap is real but modest and **grows with (a) conversion magnitude, (b) abundance
"surprise" / back-loading, and (c) budget tightness**. The mechanism of foresight value
is the **cost of delayed protection**: a myopic manager can self-correct by re-planning
(it acquires a parcel once its value becomes apparent), but it protects later than the
foresighted manager and loses some abundance to conversion in the interim.

Two reporting consequences, both handled in `08`:
- Because only ~5 of 841 parcels are acquired per period, `J_baseline` dominates `J`, so
  gaps measured on **total J look tiny**. The driver therefore also reports
  **`value_added = J − J_baseline`** (the welfare acquisition actually creates) and the gap
  on that quantity, where the model ordering is visible.
- Greedy diverges more than myopic does, because its benefit-cost criterion differs from
  the welfare objective; greedy can differ from rolling even under stationarity.

If a higher-conversion dataset or a redistribution-aware risk metric is adopted, the
separation should emerge with no code change (only `trans_prob` moves).

---

## 5. Assumptions and Limitations (consolidated)

- **Conversion metric is provisional.** Net-loss proxy from habitat erosion; low in
  magnitude; may be replaced. Methods memo notes redistribution (not loss) is the
  dominant non-stationarity signal — a candidate reason the current model separation is
  muted, and a candidate direction for a better risk metric.
- **Deterministic expected-value evaluation.** Survival `S` is used as an expected
  weight; there is no Monte Carlo over realized conversion events. "Option A" treats
  unacquired parcels as always available.
- **No spatial interactions.** Each EAU is independent — no adjacency, connectivity, or
  contiguity constraints, and no spatial spillover in benefit. Acquisition value is
  purely per-parcel.
- **Acquisition is permanent and one-shot.** Acquire-once, no divestment; per-period
  budget with no rollover.
- **Greedy uses benefit-cost ratio, not the welfare objective** — an intentional
  status-quo baseline, but it means greedy is not a "welfare-optimal under wrong beliefs"
  model the way myopic is.
- **Scale/structure differs from the toy precedent.** Huge landscape + tiny acquisition
  fraction means total-J differences are small even when per-decision differences exist;
  interpret via `value_added`.
- **Exact MDP benchmark is intractable.** A joint-landscape MDP has ~3^841 states.
- **Not yet executed end-to-end in R/Gurobi.** Logic is Python-validated; the R/Gurobi run
  on real data is the immediate next step and the first place real numbers appear.
- **Verify scenario strings.** The `SCENARIOS` block in `07` hard-codes `rcp`/`gcm` levels
  ("45","85","wet","dry","stationary"); confirm they match `data_panel` before trusting
  output.

---

## 6. Open Items / Next Steps

1. **Run the models.** `09__validation.R` (with Gurobi) → confirm P1–P5 on real data →
   `08__run_models.R` → inspect `model_results.csv` and trajectories.
2. **Sensitivity analyses.** δ (discounting), `BUDGET_EAUS_PER_PERIOD` (default 5, range
   2–10), and ε (conversion floor). The parameter blocks are at the top of `07` and `05`.
3. **Conversion data / metric.** Possibly adopt a higher-conversion or redistribution-
   aware risk layer; the model is built to absorb this through `trans_prob` alone. This is
   the most likely lever to move the headline result.
4. **MDP comparison (deferred).** Decide whether to pursue a like-for-like MDP benchmark;
   it would require either a toy-sized landscape or a per-EAU/approximate value-function
   decomposition. Interface stub is documented in `09`.
5. **Thesis framing.** Communicate the muted-separation finding honestly: model choice
   matters in proportion to how much the landscape is actually losing, and the value of
   foresight here is the value of protecting at-risk, value-shifting parcels *earlier*.

---

## Appendix A — Locked decisions (quick reference)
1. Objective: landscape-total. (Portfolio rejected.)
2. Survival clock: global from t = 0, both models. (Forward re-planning rejected — order
   violations.)
3. Acquisition value `V[i,τ] = Σ_{t≥τ} δ^t b (1−S)` = expected loss prevented.
4. Discount δ = 0.95/decade, absolute present value.
5. Budget `B_t = BUDGET_EAUS_PER_PERIOD × median(cost_t)`; default 5, sensitivity 2–10.
6. ILP: acquire-once + per-period budget (no rollover) + binary.
7. Option A: conversion enters objective only, not the constraint set.
8. Conversion proxy: existing net-loss method (script 05); kept modular.
9. Solver: Gurobi via R `gurobi` package (academic license).
10. Models evaluated against each scenario's *true* future (myopic is scored on reality
    despite choosing under a fabricated future).

## Appendix B — Key formulas
- Survival: `S[i,t] = Π_{t'=0}^{t-1}(1 − λ[i,t'])`, `S[i,0] = 1`.
- Acquisition value: `V[i,τ] = Σ_{t=τ}^{T} δ^t · b[i,t] · (1 − S[i,t])`.
- Welfare decomposition: `J = J_baseline + Σ_i V[i,τ_i]`.
- Evaluation: `J = Σ_i Σ_t δ^t · b[i,t] · w[i,t]`, `w = 1` if protected by t, else `S[i,t]`.
- Hazard: `trans_prob = max(ε, (ps_t − ps_{t+1})/ps_t)`; 0 at terminal/no-habitat; ε under
  stationary; 2020 = mean of RCP45/85 2020→2030.
- Cost: `cost[i,t] = mean_fmv_per_ha[i] × area_ha × 1.02^(year − 2017)`.

## Appendix C — File map
- Data pipeline (complete): `01__create_EAUs.R`, `02__EAU_prop_suitable.R`,
  `03__create_data_panel.R`, `04__benefit_data.R`, `05__risk_of_loss.R`, `06__cost_data.R`.
- Models (written, pending R/Gurobi run): `07__ilp_core.R`, `08__run_models.R`,
  `09__validation.R`.
- Master data object: `input_data/data_panel.rds` / `.csv`.
- Model outputs (created by `08`): `output_data/model_results.{rds,csv}`,
  `output_data/model_trajectories.csv`.
- Supporting docs: `methods_memo_conversion_risk.md` (conversion-risk rationale);
  `Land_Acquisition_as_a_time.docx` (project write-up); `README.md`.
- Prior MDP toy-problem scripts (precedent, reference only): `PPR_toyproblem.R`,
  `PPR_nonStationary_example.r`, `greedy_solver.r`, `mdp_finite_horizon_nonStationary.r`,
  `mdp_myopic_forward_look_policy.r`, `significance_tests.R`, `volatility.R`.

---

## 8. Session Addendum — First R/Gurobi Run Attempt (2026-06-24)

*Appended after a working session that attempted the first live run of `07`/`08`/`09`
on the real data with Gurobi. This records what was changed, what passed, and the open
questions — especially one structural finding about the stationary null (P3) that the
earlier Python/PuLP validation never exercised. The §1–6 body above is unchanged; treat
this section as the current front line.*

**Environment.** Gurobi 13.0.2 (academic license, active), R `gurobi` package working,
Apple M2 Pro / 16 GB / macOS 15.7.7. `input_data/data_panel.rds` present and consumed
successfully by the models.

### 8.1 Confirmed this session (can be trusted)
- **Scenario strings verified.** The `SCENARIOS` block in `07` matches `data_panel`:
  `rcp ∈ {"45","85","stationary"}`, `gcm ∈ {"wet","dry","stationary"}`, plus the shared
  `baseline/baseline` period-0 row. **The §5/§6 "verify scenario strings" caveat can be
  retired.**
- **P1 PASS** (synthetic, instant).
- **P2 PASS** (real scenario `rcp85_dry`; fast once the solver missteps below were undone).

### 8.2 Code change made to `07__ilp_core.R` — cost scaling (KEEP)
The first live run threw two Gurobi warnings: *large matrix coefficients* (`[1e+00,
1e+10]`) and *large RHS* (`[1e+00, 3e+09]`), because `cost` is raw USD (billions for some
EAUs). Inside `solve_acquisition_ilp`, costs and the budget RHS are now divided by
`cost_scale <- 1e6` **before** being passed to Gurobi; the objective vector `V` is **not**
scaled. This is identity-preserving (dividing both sides of the budget constraint by a
constant leaves the feasible region and the optimum unchanged) and only improves
conditioning — confirmed: matrix range dropped to `[1e+00, 1e+04]`, RHS to `[1e+00,
3e+03]`, warnings gone. **Reporting is unaffected:** `data_panel$cost`, `make_budget`,
and `evaluate_policy`'s `total_spend` all still operate in raw USD; the 1e6 factor lives
and dies inside the solver wrapper. A short dead-code duplication (an unscaled
`A`/`sense`/`rhs` build) was deleted; there is now a single scaled construction.

### 8.3 Solver-parameter missteps and the corrected understanding
- **Tried and REVERTED:** `NumericFocus = 2` and `ScaleFlag = 3`. These were added on the
  (mistaken) theory that the warnings were causing the slowdown. They forced slower
  careful arithmetic on *every* solve and made the previously-fast **P2 hang**. Once the
  data was scaled there was no numerical *failure* to fix, only warnings, so they bought
  nothing. Removing them restored P2. **Do not re-add them.**
- **Key correction on `MIPGap`.** The real slowdown is **not numerical** — it is a MIP
  optimality-gap-closing problem on the **stationary** scenario: Gurobi finds a
  near-optimal incumbent within seconds, then explores millions of nodes proving
  optimality because (a) the per-period budget knapsack creates many near-tied solutions
  and (b) the objective range stays wide (`[1e-07, 2e+04]`, untouched by cost scaling
  since those are the `V` values), giving a weak LP bound that descends very slowly.
  **Gurobi's default `MIPGap` is `1e-4`. The grind observed IS the `1e-4` grind** —
  setting `MIPGap = 1e-4` explicitly is a no-op and will NOT speed anything up. Loosening
  to `1e-3` was considered and **rejected**: a 0.1% solver gap is comparable to the muted
  myopic-vs-rolling effect size (§4.3) and could contaminate the headline result in `08`.
- **Recommended next lever (keeps the tight `1e-4` bound):** attack the *bound*, not the
  gap — `MIPFocus = 3` (directs effort at the best bound, which is the bottleneck) and
  optionally `Cuts = 2`. These are correctness-neutral (affect runtime only). They may
  only partly help a fundamentally degenerate instance; the fallback is to accept the
  runtime (only the stationary scenario grinds, and only its first full solve is
  expensive).
- **Housekeeping:** `OutputFlag` is currently `1` for debugging. Set it back to `0` in the
  production run, especially in `08` (rolling/myopic call the solver ~9× per scenario ×
  5 scenarios — verbose logs otherwise).

### 8.4 Structural finding on P3 (IMPORTANT — unresolved)
P3 requires myopic == rolling **exactly** on the stationary null. The premise (§3.4) is
that under stationarity the myopic frozen-future belief is *correct*. This holds in the
**synthetic** instances the Python harness used (`make_synthetic_instance` builds the
stationary hazard constant across **all** periods including t = 0), which is why P3
validated there.

On the **real** data the premise breaks at exactly one period. The stationary scenario's
period-0 (2020) row is the **shared 2020 baseline** (`rcp="baseline"`), whose `trans_prob`
is the mean of the floored RCP 4.5/8.5 2020→2030 transitions — **not ε**. From 2030 onward
the stationary hazard is ε (confirmed in script `05`). So at the **2020 decision**, myopic
freezes that baseline hazard and projects it forward, while rolling sees the hazard drop
to ε from 2030 on. Myopic's belief is correct from 2030 onward but **wrong at 2020**.

**Consequence:** myopic and rolling can make different 2020 purchases, which would trip
P3's bit-exact schedule-identity assertion **for a real structural reason, independent of
any solver gap**. This is a synthetic-vs-real-data mismatch, not necessarily a model bug —
and it is potentially interesting for the thesis: the stationary "null" is not perfectly
stationary at the anchor year because of the shared-baseline stitch (§2.2).

P2 and P5 do **not** share this exposure and need no tolerance changes: rolling's first
solve is byte-identical to the full-horizon solve (so they stay deterministically equal at
any gap), and the gap only makes P5's `LP ≥ ILP` bound easier to satisfy.

### 8.5 Decision tree for P3 (deferred to the author)
Do **not** pre-emptively loosen P3's tolerance; the correct fix depends on the failure
mode observed. After re-running P3, read its two `report(...)` lines:
- **`J_myopic ≈ J_rolling` to ~`1e-4`, schedules differ** → the 2020 baseline-vs-ε
  asymmetry and/or knapsack near-ties. Author's choice: **(a)** relax P3 to "equal within
  tolerance" and document the anchor-year wrinkle, or **(b)** make the stationary scenario
  use ε at 2020 as well — a **script-`05`** change that collides with the locked "every
  scenario shares one 2020 baseline row" convention (§2.2), so a genuine design tension,
  not a free fix.
- **`J` differs by more than ~`1e-4`** → more substantive divergence; understand it before
  changing any tolerance.

### 8.6 Where the project is being left / immediate next step
1. Add `MIPFocus = 3` to the Gurobi params in `solve_acquisition_ilp` (keep
   `OutputFlag = 1` for now), leave `MIPGap` at default, and **re-run P3 alone**. Report
   the two `report(...)` lines (the myopic/rolling `J` values and the schedule-identity
   PASS/FAIL). That single observation determines whether the remaining work is a
   tolerance change, a data decision, or nothing.
2. Then complete **P4** and confirm **P5** on real data.
3. Once P3–P5 are settled, reconcile the body: §4.1/§4.2 ("not yet run in R/Gurobi") and
   §6 item 1 should be updated to reflect the live run; the §5/§6 scenario-string caveat
   is already retired (§8.1).
4. `08__run_models.R` has **not** been run yet. Set `OutputFlag = 0` there first.

**Open questions carried forward:** (i) Does `MIPFocus = 3`/`Cuts` make the stationary
gap close in acceptable time at `1e-4`, or must the runtime simply be accepted? (ii) Is
the 2020 baseline-hazard-under-stationary an artifact to absorb in P3's check, or a data
decision (ε at 2020) that touches the shared-baseline convention? (iii) Do P4/P5 hold on
real data under the scaled, default-gap solver?
