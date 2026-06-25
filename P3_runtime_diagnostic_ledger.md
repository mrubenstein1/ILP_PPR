# P3 Runtime — Diagnostic State Ledger

**Purpose.** One place that records what we *know* about the P3 slow-runtime problem versus what we *intuit*, so the next decision is made on measured fact rather than inference. Every claim below carries a tag for how it is known:

- `[MEASURED]` — observed directly in a diagnostic run during this investigation, with the numbers cited.
- `[FROM LOG]` — read from the real P3 Gurobi console output.
- `[FROM SOURCE]` — a property of the code itself (scripts 05/07), not an inference.
- `[REPORTED]` — established in the prior §8 session / P-suite, taken as given here.
- `[HYPOTHESIS — UNMEASURED]` — plausible, consistent with evidence, but **not yet tested**. Do not treat as fact.

Date of state: 2026-06-24.

---

## 0. The two problems, kept separate

These are independent and must not be conflated; past confusion came from treating them as one.

- **Problem A — runtime.** The myopic solves on the stationary scenario find an excellent answer in seconds, then take many minutes (or time out) *proving* it optimal. This is the subject of this ledger.
- **Problem B — correctness.** P3's bit-exact "myopic == rolling" schedule check can fail on real stationary data for a legitimate structural reason (the 2020 baseline hazard is not ε). `[FROM SOURCE: 05]` Fixing A does nothing for B and vice versa. The agreed fix for B (split P3 into a strict synthetic invariant + a tolerance-based real-data check) stands regardless of how A resolves.

---

## 1. KNOWN — observed facts

### 1a. The instance

| # | Fact | How known |
|---|------|-----------|
| 1 | The MIP has **7,569 binary variables** (841 EAUs × 9 periods) and **850 constraints** (841 acquire-once + 9 per-period budget knapsacks); 15,138 nonzeros. | `[FROM LOG]` header: "850 rows, 7569 columns, 15138 nonzeros" |
| 2 | The objective the solver maximizes **is `value_added`** (the acquisition term Σ V·y). `J_baseline` is added later in evaluation, not inside the solve. | `[FROM SOURCE: 07]` |
| 3 | By solver standards this is a **small** model. Size alone does not explain the runtime. | inference from #1, uncontested |

### 1b. Where the slowness lives — it is the *policy*, not the *scenario*

This is the single most important correction to the original framing (handoff §8 attributed it to the stationary scenario).

| # | Fact | How known |
|---|------|-----------|
| 4 | On the **stationary** scenario, **rolling** (true V) closes to OPTIMAL in **28.4 s** (incumbent 100,873, gap 0.0063%). | `[MEASURED: script 11]` |
| 5 | On the **same** scenario, **myopic-t1** and **myopic-t2** both **hit the 60 s cap** (TIME_LIMIT) at tiny gaps (0.10% and 0.04%). | `[MEASURED: script 11]` |
| 6 | In the real P3 run, the first solve (myopic, t=1) ran **866 s** to a 0.01% gap, exploring **~18 million nodes**. | `[FROM LOG]` |
| 7 | P1 (synthetic invariants) and P2 (full-horizon, on rcp85_dry) run quickly and cleanly. **Neither runs the myopic loop** — P3 is the only test that exercises myopic. | `[REPORTED]` + `[FROM SOURCE: 09]` |

> Net: every slow solve observed is a **myopic** solve. The one fast solve on the same scenario is rolling. The scenario is not the discriminator; the policy is.

### 1c. The *nature* of the slowness — it is end-game proof, not search

| # | Fact | How known |
|---|------|-----------|
| 8 | The **LP relaxation is tight**: root relaxation 465,194.885 vs integer optimum 464,303.15 — a starting gap of **~0.19%**. | `[FROM LOG]` |
| 9 | An excellent incumbent appears almost immediately; the cost is spent closing the **last sliver** (≈0.1% → 0.01%), with the bound creeping down hundredths of a percent per minute over millions of nodes. | `[FROM LOG]` |
| 10 | **Multiple distinct integer schedules exist at (near-)identical objective value** — the log reports several solutions all at 464,303, and independent diagnostic solves returned *different* schedules at near-identical value. | `[FROM LOG]` + `[MEASURED: scripts 10/11]` |

### 1d. What the myopic code actually does (mechanism, from source)

| # | Fact | How known |
|---|------|-----------|
| 11 | Myopic **freezes** benefit and hazard at the current period and projects them **flat across all remaining periods**, then solves the full remaining-horizon problem. So under myopic's belief the hazard is **constant over the horizon**. | `[FROM SOURCE: 07]` |
| 12 | Myopic **implements only the current period's purchases** each iteration and **discards** the periods-2…9 schedule, re-solving from scratch next decade. The expensive optimality proof therefore spans periods that are **thrown away**. | `[FROM SOURCE: 07]` |
| 13 | Under the stationary scenario the true hazard is the **shared 2020 baseline** at t=0 (period 1), then **ε** from 2030 (period 2) on. The baseline row is the mean of floored RCP45/85 2020→2030 transitions, generally **> ε**. | `[FROM SOURCE: 05]` |
| 14 | Consequence of #11+#13, measured: myopic-t1 freezes the baseline hazard (objective up to ~2e4, incumbent ~4.64e5); myopic-t2 onward freezes ε, crushing the objective into **[1e-7, 0.33]** (incumbent ~7.9). | `[MEASURED: scripts 11/12]` |
| 15 | The current script-07 Gurobi param block is **clean**: only `OutputFlag` is set (no `NumericFocus`, `ScaleFlag`, `MIPGap`, `MIPFocus`, `Cuts`); `cost_scale = 1e6` is present. | `[FROM SOURCE: 07]` |

---

## 2. RULED OUT — tested and rejected

Each of these was proposed, then falsified by a measurement. They are dead ends.

| Hypothesis | Why it's dead | How known |
|------------|---------------|-----------|
| Stale `NumericFocus`/`ScaleFlag` were forcing slow arithmetic | Param block confirmed clean; still slow. | `[FROM SOURCE: 07]` + `[FROM LOG]` |
| Weak LP relaxation / slow-descending bound | Relaxation is tight (0.19% start, #8). | `[FROM LOG]` |
| A swarm of near-zero-value variables bloats the model (→ pruning) | `share < 1e-4` is ~0.4–0.5% in **all** objectives, fast and slow alike — not a discriminator. | `[MEASURED: script 11]` |
| The degeneracy is about **timing** (buy-now vs buy-later) (→ earlier-is-better tie-break) | Tie-break at κ=1e-6 and κ=1e-8 left **all three solves at the time limit**. | `[MEASURED: script 12]` |
| "Just wait it out" | Bound stalls; one of *many* solves alone is 866 s; P3 runs the full myopic **and** rolling loops. Not viable. | `[FROM LOG]` |
| Certified-bound / loosened gap as a standalone runtime fix | The gap is *already* tight; the cost is the proof itself. Stopping early helps only if the gap is loosened, which risks contaminating the muted myopic-vs-rolling effect. (Note: a time limit may still be justifiable on *different* grounds — see Open Question Q3.) | reasoning from #2, #8, #9 |

**Also established (not a dead end, but a resolved contradiction in the record):** an early single-solve read (script 10) closed in 29 s and suggested "no problem." That solve used the **rolling** objective, not myopic — it was the wrong target, which is why it disagreed with the real P3 run. The contradiction is fully explained by #4 vs #5. `[MEASURED: scripts 10 vs 11]`

**Confirmed safe but ineffective:** the tie-break did not distort the answer — the true (un-perturbed) value of every tie-break schedule stayed inside the reference interval [4.6407e5, 4.6454e5] (κ=1e-6 reproduced the reference incumbent exactly; κ=1e-8 found 4.6426e5). It simply didn't close the gap. `[MEASURED: script 12]`

---

## 3. INTUITED — not yet measured

Everything in this section is **inference**, consistent with the evidence but **unconfirmed**. The repeated failures above all came from acting on items like these as if they were facts. The items are grouped by what *kind* of claim they make — a distinction the earlier draft blurred:

- **3A — Diagnoses of the cause** (H1, H2): competing explanations of *why* the solver stalls.
- **3B — Mitigations** (H3): not explanations at all, but ways to *live with* the stall.
- **3C — Consequence for Problem B** (H6): not a cause, but the bridge showing that the A-side cause and the B-side correctness fix are coupled.

### 3A. Diagnoses of the cause

H1 and H2 both attribute the stall to a band of co-valued schedules the solver must rule out. They differ on **where those ties come from and whether they are exact** — and that difference is decisive, because they respond to **opposite fixes**.

- **H1 — "Manufactured, exact ties" (imposed symmetry).** `[HYPOTHESIS — UNMEASURED]` The ε floor (`[FROM SOURCE: 05]`) clamps every parcel at or below ε to the *identical* value ε. Identical hazard → identical (or proportional) V → parcels that are **perfectly interchangeable**, forcing the solver to branch through permutations of them to prove no swap helps. These ties are an **artifact of the floor**: remove or jitter it and they vanish. **If true, the ties are exact and breakable** — collapse the duplicates or add symmetry handling and the proof closes with *no* effect on the answer. **Untested:** the count of *distinct* (vs duplicated) hazard / cost / V values.

- **H2 — "Natural near-flatness" (intrinsic).** `[HYPOTHESIS — UNMEASURED]` Even with no floor, the *natural* values are clustered so tightly that schedules are **near-equal, though not identical** — the landscape genuinely does not care much which of these similar parcels you buy. **If true, the ties are not cleanly breakable**: any perturbation large enough to separate them is large enough to start choosing for reasons unrelated to conservation value (exactly the bias risk that sank the tie-break). "The optimal myopic schedule" is then close to ill-defined, no solver setting fixes it, and the resolution belongs in the **policy/test definition**, not the solver. **Untested:** how many integer solutions lie within ~0.1% of the optimum (solution-pool count).

  > **They are not mutually exclusive.** The floor can *convert* mild natural near-flatness (H2) into exact ties (H1) by clamping genuinely-different parcels onto one value. So the real question is a mix: *how much* of the co-valued band is exact-and-removable vs natural-and-intrinsic. Q2 measures the exact half; Q3 measures the near half. This is the central fork (Section 4).

### 3B. Mitigations (not diagnoses)

- **H3 — "The used decision is already settled."** `[HYPOTHESIS — UNMEASURED]` **This makes no claim about the cause.** It observes that myopic *implements only the current period* (#12) and re-solves next decade, so the long proof is spent certifying periods 2–9 that are **discarded**. *If* the period-1 pick set is invariant to further solver effort, a time limit is **provably harmless to the policy's behavior and to its realized, scored-on-reality value** — regardless of whether the cause is H1 or H2. H3 is therefore a candidate **response** that sits alongside "justified time limit," not a competing explanation. **Untested (and the test is a *safety* check, not a diagnosis):** is the period-1 pick set stable across time limits? (Q4)

### 3C. Consequence for Problem B (the bridge)

- **H6 — "Uniqueness couples A and B."** `[HYPOTHESIS — UNMEASURED]` **Categorically distinct from both the causes (3A) and the mitigation (3B):** it is a claim about the *correctness test*, not the runtime. P3 asserts myopic and rolling pick **bit-identical schedules**, which only makes sense if "the optimal schedule" is a single, well-defined object. If H2 holds (a large equivalence class of co-optimal schedules), myopic has no canonical schedule — the solver returns an arbitrary member — so a bit-identity check tests a **non-unique object** and can "fail" merely because two solves picked different members of the same tie. Thus resolving the A-side fork (H1 vs H2) **changes what the correct B-side test is**: under H1 (exact, removable) you can canonicalize and salvage a schedule check; under H2 (intrinsic) you abandon schedule-identity for value-equality — *independent* of the 2020-anchor reason already on record. H6 is the link, not a diagnosis.

---

## 4. The one fork that determines the fix

Everything reduces to: **is the co-valued band exact-and-removable (H1) or natural-and-intrinsic (H2)** — and in what mix?

- To the extent **H1**, the fix is structural/data-side: collapse the manufactured duplicate coefficients (or jitter the floor / add symmetry-aware handling) and the proof closes, with no effect on the answer.
- To the extent **H2**, no solver fix exists; the answer is to accept that the myopic schedule is non-unique and decide *at the policy/test level* — a justified time limit (pending the H3 safety check) for Problem A, and the P3 split for Problem B (per H6).

We currently **cannot say the mix**, and that is the honest state. Past attempts implicitly assumed H1-style breakability without checking — the tie-break result (it stayed safe but did not close) is mild evidence that at least *part* of the band is H2, but that is not yet measured directly.

---

## 5. OPEN QUESTIONS → the measurements that resolve them

None of these is a fix; each is a **read-only** measurement chosen to discriminate among the hypotheses, written as a decision rule ("if X → conclude \_\_; if Y → conclude \_\_") so the section reads as a procedure, not a to-do list. Q4 and Q5 are the two former "hypotheses" H4/H5, restated as the experiments they always were. Recommended order reflects information-per-effort.

**Q1 — Count the *exact* duplication in the coefficients.**
*Measurement:* number of distinct values in `sc$lam`, `sc$cost`, and the myopic V, and the size of the largest identical-value groups; in particular, how many parcels sit *exactly* at ε.
*Mechanism it tests:* H1 says the floor manufactures exact ties. *Decision:* if a large share of parcels collapse onto a handful of identical hazard/V values → **H1 is real and large**, the ties are exact, and collapsing duplicates is a legitimate, answer-preserving fix. If values are nearly all distinct → **H1 is small**, and the band must be H2 (near-flatness), which is *not* cleanly breakable.

**Q2 — Count the schedules within ~0.1% of the optimum (solution pool).**
*Measurement:* Gurobi `PoolSearchMode` + `PoolGap≈1e-3`; report how many distinct integer schedules populate the near-optimal band.
*Mechanism it tests:* H2 (intrinsic near-flatness) and the input to H6. *Decision:* a handful → the optimum is effectively unique → tuning/canonicalization can work and a schedule-based P3 is salvageable. Thousands → "the optimal myopic schedule" is **not a well-defined object** → no perturbation fixes it, the runtime resolution is policy-level (time limit), and **Problem B must drop bit-identity for value-equality** (this is H6 turning from hypothesis into finding).

> Q1 and Q2 are the **fork-resolving pair**. Q1 measures the *exact-and-removable* half of the band; Q2 measures the *near-and-intrinsic* half. Together they give the mix in Section 4. Run these two first.

**Q3 — Run myopic on a climate scenario (rcp85_dry).** *(was H4)*
*Measurement:* myopic-t1 and myopic-t2 runtime/gap on rcp85_dry, vs stationary.
*Mechanism it tests:* this is a **controlled experiment**, not just a scope check — it holds the policy fixed and varies only the hazard's spread, the exact variable H1/H2 implicate. Two diagnostic facts come out of it:
- *t=1* freezes the **shared 2020 baseline**, which is *identical across all scenarios* (`[FROM SOURCE: 05]`). So myopic-t1 should be **just as slow on rcp85_dry as on stationary**. If it isn't, the cause is not the baseline row and our mechanism is wrong.
- *t≥2* freezes the scenario's own hazard — **ε under stationary, but real spread-out values under climate.** *Decision:* if climate myopic-t2 **closes fast** while stationary myopic-t2 stalls → direct evidence that **low hazard variance is the driver** (supports H1/H2 and pins the mechanism to the ε regime), *and* the project bottom line is reassured (the stall is contained to the null, not the whole results table). If climate myopic-t2 **also stalls** → variance is *not* the driver and a mechanism is still missing.

**Q4 — Is the period-1 pick set invariant to solver effort?** *(safety check for the H3 mitigation, not a diagnosis)*
*Measurement:* solve myopic-t1 at several TimeLimits (e.g. 5 s, 30 s, full); compare the *implemented* (period-1) pick sets.
*Decision:* identical across limits → a time limit is **provably harmless** to what myopic actually does, and Problem A can be closed with a justified cap regardless of the H1/H2 mix. Pick set wobbles → the implemented decision genuinely needs solver effort, and a time limit is *not* safe.

**Q5 — Does each period's budget bind?** *(was H5)*
*Measurement:* budget used vs available, per period, at the incumbent.
*Mechanism it tests:* whether there is a hard combinatorial problem *at all*, or only symmetry. *Decision:* budgets mostly **slack** → there is no real packing decision; the only thing left is choosing among interchangeable parcels → the hardness **is** symmetry (H1/H2 are the whole story). Budgets tightly **binding** → a genuine knapsack with many similar value-per-dollar parcels, which is hard *independent* of any floor artifact → symmetry is only *part* of the story and duplicate-collapsing alone will not be enough.

**Leading expectation (explicitly an intuition, not a finding):** Q1 + Q3 together likely explain it — the shared baseline plus the ε floor producing low-variance, many-identical coefficients. But that exact instinct has misled this investigation more than once, so it stays here until the numbers land.

---

## 6. One-line summary

We know **with certainty**: the problem is the *myopic* policy (not the scenario), the relaxation is *tight*, the cost is *proving optimality* across *many co-valued schedules* over periods myopic *discards*, and pruning / tie-break / waiting / stale-params are all *ruled out*. We do **not** yet know the **mix** between two causes — **exact ties manufactured by the ε floor (H1, removable)** and **natural near-flatness (H2, intrinsic)** — and that mix decides whether Problem A's fix is a coefficient-level reformulation or a policy-level time limit. Two further items are *not* causes and were re-labelled as such: the time-limit idea (H3) is a **mitigation** whose safety is checked by Q4, and schedule non-uniqueness (H6) is the **bridge** by which the H1/H2 answer also determines Problem B's correct test. Q1 + Q2 resolve the fork; run them first.
