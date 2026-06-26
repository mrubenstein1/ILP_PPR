# PPR Land Acquisition Under Climate Non-Stationarity — Project Handoff

*Context document for a future collaborator (human or AI) joining this project. It states
the research question, data, modeling approach, current state, the decisions that are
already settled (so they are not reopened), and the known limitations.*

**Last updated 2026-06-26** (P3a closed — the stationary null is now fully certified, §4.4;
**`08` now persists the per-parcel schedule and reporting is consolidated into a single
`10_results_figs.R` that produces tables + figures + maps**, §4.1/§6 items 3 & 7), after (a)
the first complete R/Gurobi run on
real data and (b) a follow-on **budget-deployment ("spend-down") session**. The solver-runtime
investigation is **resolved** (§8); the headline comparison has been produced and now
holds under a *realistic full-budget-deployment* rule (§9). See §4 for the current state,
§8 for the runtime resolution record, and **§9 for the spend-down decision record and its
central finding — that the foresight gap is a *targeting* effect, not a *deployment*
one.** Sections §1–§3 (question, data, locked design) are unchanged except for the new
§3.8 (spend-down).

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
4. Each policy **deploys its full per-period budget** (the realistic "spend-down" rule,
   §3.8, `SPEND_DOWN_MODE = "spend"` in `07`, default on). This was added after the first
   run; it does **not** change the headline (§9), and it replaces a MIPGap-sensitive
   under-spend artifact with an explicit behavioural rule. Set `SPEND_DOWN_MODE = "off"`
   to reproduce the pre-spend-down numbers.

**Headline finding.** Foresight matters, and substantially. Across the four climate
scenarios the foresighted (rolling) policy beats the myopic policy by **~13–15% of
value_added** and the greedy status-quo heuristic by **~43–46%**. This *overturns* the
earlier expectation (former §4.3) that the models would be "nearly indistinguishable." The
separation is masked on total landscape J (where the gap is ~0.09%, because only ~60 of
841 parcels are ever acquired) and is visible on **value_added = J − J_baseline**, the
welfare acquisition actually creates. Read §4.3 for the numbers and §9 for the mechanism.

**Mechanism finding (spend-down session).** Forcing myopic to deploy its entire budget
(it previously left ~22% idle) lifts its acquisitions to ~full-budget levels but leaves
the gap essentially unchanged (it even widens slightly in 3 of 4 climate scenarios). So
**the cost of myopia is mis-*targeting*, not under-*deployment*.** Under-spending in the
first run was a *symptom* of the blind belief, not its cause; a manager that cannot see
which currently-safe parcels will become at-risk gains nothing from spending more, because
it spends the freed budget blindly. This is a *stronger* result than the original framing
and answers the obvious reviewer objection ("isn't the gap just an unspent-budget
artifact?"). See §9.

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

> **Risk gates the redistribution signal (clarified, spend-down session).** Acquisition
> value is `V = Σ δ^t · b · (1 − S)`: you can only "save" a parcel to the extent it is
> *both* abundant *and* at conversion risk. A parcel gaining abundance via redistribution
> but facing no conversion risk has `(1 − S) ≈ 0` and contributes ~nothing to acquisition
> value. So even though redistribution is the larger *physical* non-stationarity signal,
> the *decision-relevant* signal is **risk** (where loss falls, on parcels that would have
> been valuable). The §6-item-5 freeze-b-only vs freeze-λ-only decomposition quantifies
> this split; the structural prediction is that λ gates it.

### 3.8 Budget deployment — "spend-down" (decided, default on)
Each ILP solve optimises acquisition value `V` as the **primary** objective exactly as in
§3.2, then — among value-co-optimal plans — runs a **strictly lower-priority secondary**
objective that deploys the *implemented* period's budget as fully as possible. This is the
realistic rule (PPR managers cannot roll budget over, so they do not leave money idle) and
it removes a MIPGap-sensitive artifact: under value-only solving, myopic left ~22% of its
budget unspent because, under its frozen belief, the marginal parcels looked near-valueless
and the solver was indifferent to them.

- **Implementation (`07`):** Gurobi native multi-objective (`model$multiobj`), priority-2 =
  value (`V`) with `SPEND_DOWN_RELTOL = 1e-4` allowed degradation, priority-1 = the
  secondary. `SPEND_DOWN_MODE ∈ {"off","spend","count"}` (default `"spend"`): `"spend"`
  maximises dollars deployed; `"count"` maximises parcels acquired ("more lottery tickets").
  The secondary targets only `implement_periods` (the single enacted period in the
  rolling/myopic loop) — deliberately, so it deploys the budget actually about to be spent
  and avoids a perverse "defer to a pricier later period" tie under cost inflation.
- **Scope:** affects only the ILP policies — materially for myopic, negligibly for rolling
  (which already nearly exhausts the budget). **Greedy is untouched** (it deploys its budget
  by construction, buying down its ranked list until nothing affordable remains).
- **Important subtlety (do not over-trust monotonicity).** Forcing deployment was expected
  to *weakly raise* myopic's true J (pure addition of non-negative-value parcels) and thus
  only narrow the gap. In practice **it slightly *widens* the gap in 3 of 4 climate
  scenarios.** The reason: `SPEND_DOWN_RELTOL` bounds degradation of the *frozen* objective
  (the belief), not the *true* value used to score the result, and frozen and true value
  diverge precisely on the currently-safe / future-at-risk parcels. Within the frozen-tie
  band the cost-driven tiebreak can swap a truly-valuable parcel for a truly-marginal one
  with no measurable change in the frozen objective. (Part of the wobble is also
  multi-threaded run-to-run non-determinism, ~0.3% scale.) If a strictly monotone "never
  hurts" variant is ever needed, use **fix-and-fill** (hard-fix the value-optimal core as
  `y = 1`, then maximise spend over the remaining budget); it would hold the gap flat but
  still would not *close* it (§9).

---

## 4. Current State

### 4.1 Modeling scripts (written and RUN on real data with Gurobi)
Three R scripts continue the pipeline numbering, plus two reporting/extraction scripts:
- **`07__ilp_core.R`** — shared engine: parameter block (δ, budget, scenario coding, a
  **Solver control** block — see §8.2 — and the **spend-down** block, §3.8),
  `build_scenario_matrices` (with the stationary anchor correction, §8.3),
  `compute_survival_matrix`, `compute_value_vector`, `make_budget`, `solve_acquisition_ilp`
  (Gurobi wrapper: cost-scaling, a per-solve time cap, the optional value-first/spend-second
  multi-objective tiebreak, and an optional solve-log hook), the three policies +
  `run_full_horizon` (all now forward `spend_down`/`implement_periods`), `evaluate_policy`,
  and synthetic-instance + brute-force helpers for validation.
- **`08__run_models.R`** — driver: runs all three models on all five scenarios, scores
  each against its scenario's true future, computes `J_baseline`, and writes results +
  per-period trajectories + a per-run solver-convergence summary (RDS + CSV). Records
  `spend_down_mode` in the saved params and prints it at run start. Has a `REPRODUCIBLE`
  toggle (single-thread, fixed seed) for bit-reproducible reported numbers. **Now also
  persists the per-parcel acquisition schedule (2026-06-26):** it joins each policy's
  `acquired` vector to EAU coordinates (from `input_data/eau_wmd_lookup.rds`) and writes
  `output_data/acquisition_schedule_spatial.csv`, also embedding `schedule` in the saved
  `.rds`. This resolved the former "08 discards schedules" gap (§6 item 7) and retired the
  separate extract step. The persisted schedule reflects the run's mode, so run with
  `REPRODUCIBLE <- TRUE` for a stable map footprint.
- **`09__validation.R`** — correctness suite (run first, §4.2), plus a documented stub for
  the deferred MDP comparison. Formulation-certifying tests (P1/P2/P5) run with
  `spend_down = "off"` (they concern the primary objective, which the tiebreak never
  changes); P4 runs in the active mode (production ordering).
- **`10_results_figs.R`** (reporting — single consolidated script) — pure reader of `08`'s
  outputs (no solver, no `07`). Reads `model_results.csv`, `model_trajectories.csv`, and
  `acquisition_schedule_spatial.csv` and writes to `output_figs/`: two results tables + six
  figures (`ggplot2`/`flextable`: value-added gaps; absolute Δpairs; wet/dry
  stakes-vs-premium; J-vs-value_added; cumulative gap over time; landscape decline) **and**
  three maps (footprint by decade; rolling-vs-myopic foresight difference; wet-vs-dry).
  Restyleable config block on top; `DRAW_MAPS` toggles the map section. The maps need only
  `ggplot2`; `terra`/`sf` are used solely for the optional WMD outline overlay (it dissolves
  `input_data/wmd_raster_equal_area.tif` to district polygons) and are skipped gracefully if
  absent — so tables + figures still build without the spatial stack.
- **Retired (2026-06-26):** `10_extract_schedules.R` / `11_extract_schedules.R` (the schedule
  now comes from `08`) and the standalone `11_maps.R` (folded into `10_results_figs.R`).
  *Naming note:* the archived diagnostics (formerly `10`–`14`, §6) live in `archive/` and are
  separate from this reporting script.

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
  information only. **P3a now passes on the production run (2026-06-26):** the stationary
  benefit and non-terminal hazard trajectories are flat, so the null is exact and this loop
  is closed. (Had P3a's benefit check failed, it would have meant the stationary abundance
  trajectory is not flat — an upstream (01–04) data-construction matter, not a solver one;
  see §6 item 8.)
- **P4** rolling ≥ myopic on every scenario, checked up to solver tolerance (a margin of
  ~2% of value_added — above capped-solver noise, below any genuine order violation). Run
  in the active `SPEND_DOWN_MODE` (production ordering); spend-down lifts myopic toward but
  never past the foresight optimum, so the ordering holds with margin.
- **P5** LP relaxation ≥ ILP (holds for free under the cap: the LP bound dominates any
  feasible incumbent).

Solver is **Gurobi via the native R `gurobi` package**; CBC blows up beyond ~40 parcels,
which is why Gurobi is required at N = 841.

### 4.3 Key empirical finding (the headline)
Under landscape-total with the current data, **foresight provides a large, consistent
advantage on value_added** across all four climate scenarios. Numbers below are the
**production (spend-down on)** run:

| metric (climate scenarios) | greedy | myopic | rolling |
|---|---|---|---|
| value_added gap vs rolling | ~43–46% | **~13–15%** | 0 (reference) |
| gap on total landscape J    | ~0.24–0.31% | ~0.08–0.10% | 0 |
| parcels acquired (of 841)   | ~60–64 | ~61–65 | ~65–67 |
| total spend                 | ~13.0B | ~13.3–13.4B | ~13.4B |

Per-scenario value_added gaps (spend-down run): myopic 13.8 / 13.4 / 14.9 / 15.3 %, greedy
45.1 / 45.1 / 45.5 / 43.4 % for (45wet, 45dry, 85wet, 85dry). Stationary value_added ≈ 9
duck-pairs out of 37.8M (a true null; its %-gaps are noise on a near-zero base).

> **Spend-down vs value-only.** Pre-spend-down, myopic acquired ~55–57 parcels and spent
> ~10.5B (≈78% of budget); spend-down lifts it to ~61–65 / ~13.3–13.4B (≈full budget) but
> the **gap is unchanged** (§9). The headline ~13–15% holds either way; the production run
> uses spend-down because it is the realistic rule and forecloses the "unspent-budget
> artifact" objection.

Two things to internalize:

1. **Read value_added, not total J.** Because only ~60 of 841 parcels are acquired,
   `J_baseline` dominates `J`, so the gap on total J is ~0.09% and looks negligible. The
   decision-relevant quantity is `value_added = J − J_baseline` (the welfare acquisition
   creates), where the ordering is large and clear. The driver reports both.
2. **The mechanism is mis-*targeting*, not under-*deployment*.** A myopic manager re-plans
   each decade and self-corrects, but it cannot identify the parcels that look safe under
   today's hazard yet will become valuable losses later — so it protects the wrong set.
   Forcing it to deploy its full budget (spend-down, §3.8/§9) does **not** recover the gap:
   it spends the freed budget blindly, and the gap is unchanged-to-slightly-wider. The
   original under-spend was a *symptom* of the blind belief, not the cause. Greedy diverges
   more still, because its benefit-cost criterion differs from the welfare objective
   entirely (and it spends the full budget too — same lesson, larger).

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
myopic and rolling coincide up to solver tie-breaking. **Confirmed on the production run
(2026-06-26):** P3a and P3b both pass, and the stationary row of `08`'s output shows
value_added ≈ 0 (≈9 duck-pairs out of 37.8M) with myopic ≈ rolling up to solver
tie-breaking. The null is exact.

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
  budget with no rollover (the no-rollover assumption is what motivates the spend-down
  rule, §3.8).
- **Greedy uses benefit-cost ratio, not the welfare objective** — an intentional
  status-quo baseline, not a "welfare-optimal under wrong beliefs" model the way myopic is.
- **Budget-deployment rule is a modelling choice (spend-down, §3.8).** Production uses
  `SPEND_DOWN_MODE = "spend"`. It changes *which/how many* parcels myopic buys but **not the
  headline gap** (§9). Two caveats for a future agent: (i) the "forced deployment can only
  help" intuition is *false* here — the `1e-4` reltol bounds the belief's objective, not
  true value, so the gap can drift slightly the wrong way; (ii) `"spend"` vs `"count"` give
  different schedules and a fix-and-fill variant would behave differently again — none close
  the gap, but the maps and any per-parcel claims depend on the choice, so state it.
- **The cap's conservative direction holds only value-side, not deployment-side.** Under
  value-only solving, truncation biased myopic *down* (conservative). Under spend-down that
  one-directional guarantee is weaker (see above); rely on the *magnitude* of the effect
  (~13–15% ≫ ~0.3% wobble), not on a strict monotonicity argument.
- **Reported numbers carry a small solver-precision band.** Each solve is capped at 60 s
  (§8.2); on the climate scenarios this leaves <~0.5% of value_added of uncertainty on the
  myopic numbers, far below the ~13% effect. Multi-threaded runs are not bit-reproducible
  (~0.3% wobble); for exact reproducibility — and a stable map footprint, since the maps
  now read the schedule `08` persists — run `08` with single-thread + fixed seed
  (`REPRODUCIBLE <- TRUE`).
- **Exact MDP benchmark is intractable.** A joint-landscape MDP has ~3^841 states (§ stub
  in `09`).
- **Stationary abundance flatness — confirmed (2026-06-26).** The anchor fix flattens the
  stationary *hazard*; P3a additionally certifies that *benefit* is flat, and it now
  **passes**, so the stationary null is exact. (Had it failed, the null would not have been
  exact and the cause would have been upstream in 01–04.)

---

## 6. Open Items / Next Steps

1. **Confirm the run.** `09__validation.R` (Gurobi) → P1–P5. *(P3a, stationary benefit
   flatness, now passes — §4.4; this item is satisfied on the production run.)* →
   `08__run_models.R` → inspect `model_results.csv`, the trajectories, and
   `solver_convergence.csv`.
2. **Lock reportable numbers.** Run `08` (spend-down on) single-thread/seeded for the
   figures that go in the thesis; report value_added gaps to ~2 significant figures (≈13–15%
   myopic, ≈43–46% greedy) with the solver band noted.
3. **Figures + maps (done, 2026-06-26 — consolidated).** `08` now persists
   `acquisition_schedule_spatial.csv` (per scenario × model × EAU: acquired period/year +
   EAU coordinates, ESRI:102039 USGS Conus Albers, on a regular ~18.5 km lattice). The single
   `10_results_figs.R` reads `08`'s three CSVs and produces the two tables, six figures, and
   three maps: (a) acquisition-footprint per model × scenario, shaded by acquisition decade;
   (b) rolling-vs-myopic "what foresight adds/avoids" difference (both-/rolling-only/
   myopic-only/neither); (c) wet-vs-dry footprint comparison. WMD outlines are overlaid from
   `input_data/wmd_raster_equal_area.tif` (dissolved to district polygons via `terra`/`sf`;
   optional). *Remaining polish (optional):* a state basemap, and the final scenario subset
   for the thesis plates.
4. **Decompose the foresight signal (thesis enrichment, HIGH VALUE).** Myopic freezes
   *both* `b` and `λ`. Re-run with selective freezing — (frozen b, true λ) and (true b,
   frozen λ) — to attribute the gap to anticipating *risk* vs *abundance redistribution*.
   Structural prediction (§3.8 box): **λ gates it** — "sees risk, blind to abundance" should
   recover most of rolling's value, "sees abundance, blind to risk" little, because
   `V = b·(1−S)` pays only for prevented loss. This is the cleanest next analysis and pairs
   naturally with the targeting-vs-deployment finding (§9). NB cost is *not* a channel:
   myopic already uses the true cost trajectory, so cost foresight contributes nothing to
   the rolling-vs-myopic gap.
5. **Conversion-risk magnitude sensitivity (separate branch).** Planned: scale the hazard
   up (the PPR's true gross conversion likely exceeds the FOREsce net-loss proxy). Do **not**
   multiply `trans_prob` directly (it is a per-decade probability in [0,1]; ×10 overflows and
   corrupts the survival recursion). Scale the *hazard rate*: `λ' = 1 − exp(−k·(−ln(1−λ)))`
   (≈ k·λ for small λ, saturates gracefully). Watch the **Option A** assumption (§3.6): at
   high λ, "unacquired ⇒ always available" may break and converted parcels may need removal
   from the choice set. Expect larger gaps on *both* metrics and *faster* solves (higher,
   more-differentiated λ un-flattens the degenerate objective).
6. **Other sensitivity analyses.** δ (discounting), `BUDGET_EAUS_PER_PERIOD` (default 5,
   range 2–10), ε (conversion floor), and `SPEND_DOWN_MODE` ("spend" vs "count").
7. **Persist schedules in `08` — DONE (2026-06-26).** `08` now joins each policy's
   `acquired` vector to EAU coordinates and writes `acquisition_schedule_spatial.csv` (and
   embeds `schedule` in `model_results.rds`), so maps read a saved artifact and never
   re-solve. The interim extract script is retired.
8. **Stationary abundance flatness — RESOLVED (2026-06-26).** P3a passes: the stationary
   benefit trajectory is flat, so the null is exact and no upstream (01–04) change is needed.
   (Kept for history: had P3a failed, the fix would have been an upstream decision about
   whether the stationary scenario holds 2020 abundance constant — a data-construction
   choice, not a model change.)
9. **MDP comparison (deferred).** Interface stub documented in `09`.

---

## 7. Thesis framing

Communicate the result precisely: model choice matters in proportion to how well
foresight identifies *which* parcels to protect, and that effect is large here (~13–15%
myopic, ~43–46% greedy on value_added) even though it is invisible on the landscape total.
The value of foresight is the value of protecting the parcels that look safe under today's
conditions but will become valuable losses later — a **targeting** advantage. The
spend-down result sharpens this: forcing the myopic manager to spend its whole budget does
not recover the gap (§9), so the story is not "myopia under-invests" but "myopia invests in
the wrong places." The honest caveats are the provisional conversion metric and the small
solver-precision band on the capped solves; neither threatens the finding.

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
certified separately by P3a — **which passes on the production run (2026-06-26)**, so the
stationary null is fully exact.

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

---

## 9. Spend-down Decision Record — Budget-Deployment Session (2026-06-25)

*This section records the budget-deployment ("spend-down") work that followed the first
full run. It is a design + finding record, parallel to §8. The central result revised the
mechanism narrative in §0/§4.3/§7.*

### 9.1 Motivation
In the first (value-only) run, myopic spent ~78% of its budget and acquired ~55–57 of the
parcels rolling did (~65–67). Two reasons to force full deployment: (i) **realism** — PPR
managers cannot roll budget over, so leaving ~22% idle is not a real behaviour; (ii)
**rigour** — the under-spend was partly a *MIPGap artifact* (the marginal parcels had
frozen value below the solver's stopping threshold, so the solver was indifferent and left
them unbought), which a reviewer could dismiss as a tolerance choice. Spend-down replaces
that artifact with an explicit, defensible rule.

### 9.2 What was implemented (`07`/`08`/`09`)
See §3.8 for the mechanism. In short: a Gurobi lexicographic multi-objective — value first,
then (among value-co-optimal plans) maximise deployment of the *implemented* period's
budget — exposed as `SPEND_DOWN_MODE ∈ {"off","spend","count"}` (default `"spend"`) with
`SPEND_DOWN_RELTOL = 1e-4`. `08` records the mode; `09` runs formulation tests (P1/P2/P5)
with `"off"` and the ordering test (P4) in the active mode. Greedy is unaffected.

### 9.3 The central finding — targeting, not deployment
With spend-down on, myopic deploys ≈full budget (~13.3–13.4B, ~61–65 parcels), but **the
value_added gap is unchanged** — ~13–15%, and it *widens slightly* in 3 of 4 climate
scenarios (myopic value_added moved −1.7% to +0.1%; rolling ≈ unchanged):

| scenario | myopic gap, value-only → spend-down | myopic spend → |
|---|---|---|
| RCP4.5 wet | 13.1% → 13.8% | 10.5B → 13.4B |
| RCP4.5 dry | 13.5% → 13.4% | 10.5B → 13.3B |
| RCP8.5 wet | 13.5% → 14.9% | 10.5B → 13.4B |
| RCP8.5 dry | 15.1% → 15.3% | 10.5B → 13.4B |

**Interpretation.** The cost of myopia is **mis-targeting**, not under-deployment. Myopic
has no signal for which currently-safe parcels will become valuable losses, so spending the
freed budget buys a roughly value-neutral basket and recovers ~none of the gap. Under-spend
was a *symptom* of the blind belief, not its cause. This is a *stronger* and more defensible
result than the original "under-acquires/under-spends" framing, and it forecloses the
"unspent-budget artifact" objection.

### 9.4 Why it can *worsen* slightly — a correction to record
An earlier verbal claim that forced deployment "can only weakly raise myopic's true J / only
narrow the gap" was **too strong**. `SPEND_DOWN_RELTOL` bounds degradation of the *frozen*
(belief) objective, not the *true* value used to score the result; frozen and true value
diverge exactly on the contested currently-safe parcels, so within the frozen-tie band the
cost-driven tiebreak can swap a truly-valuable parcel for a truly-marginal one at no
measurable frozen-objective cost. Hence true J (and the gap) can drift slightly the wrong
way. The clean "never hurts" variant is **fix-and-fill** (hard-fix the value-optimal core,
then fill remaining budget with positive-value parcels); it would hold the gap flat, but it
still would **not close** it — so the substantive conclusion is unchanged, and fix-and-fill
is optional polish, not a correction needed for the result.

### 9.5 Tooling produced
- `10_results_figs.R` — single `ggplot2`/`flextable` reporting script: two tables + six
  figures + three maps, all from `08`'s saved CSVs → `output_figs/` (no solver). *(Updated
  2026-06-26: this consolidates the former `10_figures.R` and the retired extract/map
  scripts; the schedule it maps now comes directly from `08` — §4.1, §6 items 3 & 7.)*

### 9.6 Decision
**Adopt spend-down (`"spend"`) as the production rule.** It is realistic, it leaves the
headline intact, and it strengthens the framing. Keep `"off"` available to reproduce the
first-run numbers. Consider running the §6-item-4 freeze decomposition next — it is the
natural complement to this finding (deployment ruled out ⇒ isolate which *information*,
risk vs abundance, drives the targeting advantage).
