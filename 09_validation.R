# ══════════════════════════════════════════════════════════════════════════════
# 09__validation.R  —  Correctness checks for the acquisition models
# ══════════════════════════════════════════════════════════════════════════════
#
# Run this FIRST (with Gurobi active) to certify the formulation on your machine.
# It mirrors the property tests that were validated in a separate solver harness:
#
#   P1  ILP == brute-force optimum            (tiny synthetic instance; exact)
#   P2  rolling ~= full-horizon single solve  (real scenario; up to solver tolerance)
#   P3  stationary null: myopic belief correct (solver-independent: flat inputs +
#                                               frozen value == true value)
#   P4  ordering: rolling >= myopic           (all scenarios; up to solver tolerance)
#   P5  LP relaxation >= ILP                   (real scenario; valid bound)
#
# Each check appends a PASS/FAIL to `checks`; the script stops if any fails.
#
# ── SPEND-DOWN AND VALIDATION ──────────────────────────────────────────────────
# The formulation-certifying tests run with spend_down = "off" so they certify the
# pure value optimum, independent of the budget-deployment tiebreak: P1 (ILP == brute
# force), P2 (rolling loop == full-horizon solve), and P5 (LP >= ILP) all concern the
# primary objective, which the spend tiebreak never changes. The behavioural tests run
# in the ACTIVE SPEND_DOWN_MODE: P4 (rolling >= myopic) checks the ordering of the
# policy as deployed, and P3's informational realised gap reflects production. P3a/P3b
# are solver-independent (properties of the input and value vectors) and unaffected by
# either setting.
#
# ── TOLERANCES UNDER THE PRODUCTION TIME CAP ──────────────────────────────────
# 07 now caps each solve (SOLVER_TIME_LIMIT) so the production policy is computed
# the same way it is reported. Two consequences for validation:
#   (1) A capped solve returns a near-optimal incumbent, not a certified optimum, so
#       equalities that hold only under exact solves (P2, P4) are checked up to a
#       tolerance set to exceed capped-solver noise (~sub-0.5% of value_added) while
#       still catching genuine, large (survival-clock-bug-scale) errors.
#   (2) The stationary null is maximally degenerate (epsilon hazard => near-zero,
#       near-flat objective), so two capped solves pick different CO-OPTIMAL
#       schedules. Asserting bit-identical schedules would fail on a correct model.
#       P3 is therefore reformulated to test the SOLVER-INDEPENDENT content of the
#       null — that the stationary inputs are flat and the myopic frozen value vector
#       equals the rolling true value vector — and reports the realised myopic-vs-
#       rolling gap for information only. (P1/P5 are unaffected: P1's instance is tiny
#       and solves exactly; P5's LP bound dominates any feasible incumbent, capped or
#       not.) To validate at full convergence instead, raise SOLVER_TIME_LIMIT before
#       sourcing — but note that uncapping the stationary myopic solve reinstates the
#       multi-minute runtime the cap was introduced to bound.
#
# WHY a tiny synthetic instance for P1.  Brute force enumerates (n_t+1)^n_eau
# assignments, so it is only tractable for a handful of EAUs. The real panel
# (n_eau = 879) is far beyond enumeration — P1 certifies the FORMULATION on a small
# instance; P2/P3/P5 then exercise the same solver on the real data.
# ══════════════════════════════════════════════════════════════════════════════

source("07_ilp_core.R")
if (!exists(".SETUP_DONE")) source("00_setup.R")

TOL    <- 1e-6     # exact-arithmetic tolerance (P1, and the solver-independent P3 checks)
VA_TOL <- 0.02    # solver-tolerance checks (P2, P4): fraction of value_added that a
                  # capped solve may move J by; exceeds capped-solver noise (<~0.5%)
                  # yet flags any genuine order/equivalence violation (which would be
                  # many % of value_added). value_added = J - J_baseline.
checks <- list()
report <- function(name, ok, detail = "") {
  checks[[name]] <<- isTRUE(ok)
  cat(sprintf("  [%s] %-46s %s\n", ifelse(isTRUE(ok), "PASS", "FAIL"), name, detail))
}

data_panel <- readRDS(file.path(DIR_DERIVED, "data_panel.rds"))


# ── P1. ILP == brute-force optimum (tiny synthetic instance) ──────────────────####
cat("\nP1  ILP vs brute force (synthetic)\n")
{
  inst <- make_synthetic_instance(n_eau = 6, n_t = 5, seed = 11, severity = 1.5)
  V <- t(vapply(seq_len(nrow(inst$b)),
                function(i) compute_value_vector(inst$b[i, ], inst$lam[i, ]),
                numeric(ncol(inst$b))))
  bud  <- make_budget(inst$cost, n_eaus = 2)
  ilp  <- solve_acquisition_ilp(V, inst$cost, bud,
                                avail = seq_len(nrow(V)), periods = seq_len(ncol(V)),
                                spend_down = "off")
  bf   <- brute_force_optimum(V, inst$cost, bud)
  report("ILP objective equals brute-force optimum",
         abs(ilp$objval - bf) < 1e-4,
         sprintf("(ILP = %.6f, brute = %.6f)", ilp$objval, bf))
}


# ── P2. rolling == full-horizon single solve (real scenario) ──────────────────####
cat("\nP2  rolling-loop vs full-horizon solve (real scenario rcp85_dry)\n")
{
  sc  <- build_scenario_matrices(data_panel, SCENARIOS[["rcp85_dry"]])
  bud <- make_budget(sc$cost)
  roll <- run_rolling_horizon(sc$b, sc$lam, sc$cost, bud, spend_down = "off")
  full <- run_full_horizon(sc$b, sc$lam, sc$cost, bud, spend_down = "off")
  Jr <- evaluate_policy(roll, sc$b, sc$lam, sc$cost)$J
  Jf <- evaluate_policy(full, sc$b, sc$lam, sc$cost)$J
  Jb <- evaluate_policy(rep(NA_integer_, nrow(sc$b)), sc$b, sc$lam, sc$cost)$J
  va  <- abs(Jf - Jb)
  # Under the cap, the single full-horizon solve is the hardest instance and may stop
  # at the time limit, while the rolling loop decomposes into easier per-period solves
  # that converge — so rolling can even edge slightly ABOVE the capped full solve.
  # Equivalence is therefore checked up to VA_TOL of value_added, not exactly.
  report("rolling J ~= full-horizon J (up to solver tolerance)",
         abs(Jr - Jf) <= VA_TOL * max(va, 1),
         sprintf("(rolling = %.2f, full = %.2f, |gap| = %.3f%% of value_added)",
                 Jr, Jf, 100 * abs(Jr - Jf) / max(va, 1)))
}


# ── P3. stationary null: under stationarity the myopic belief is exactly correct ──####
cat("\nP3  stationary null (solver-independent: flat inputs + frozen value == true value)\n")
{
  sc  <- build_scenario_matrices(data_panel, SCENARIOS[["stationary"]])
  bud <- make_budget(sc$cost)
  n_t <- ncol(sc$b)
  nonterm <- seq_len(n_t - 1)   # terminal hazard is 0 by construction, excluded

  # The content of the null is that the world does not change, so the myopic frozen
  # belief (current values persist) coincides with reality and the two policies pose
  # identical programs every period. We test that directly and solver-independently,
  # rather than asserting bit-identical schedules from a degenerate capped solver.

  # P3a — the stationary scenario IS stationary: flat benefit and flat (non-terminal)
  # hazard. The 2020-anchor fix in build_scenario_matrices guarantees flat hazard; a
  # FAIL on benefit would mean the stationary abundance trajectory is not flat (an
  # upstream data-construction issue in 01-04), not a solver problem.
  b_flat   <- all(abs(sc$b - sc$b[, 1])                <= TOL * pmax(1, abs(sc$b[, 1])))
  lam_flat <- all(abs(sc$lam[, nonterm] - sc$lam[, 1]) <= 1e-12)
  report("stationary: benefit trajectory is flat",          b_flat)
  report("stationary: (non-terminal) hazard is flat",       lam_flat)

  # P3b — given stationarity, the myopic FROZEN value vector equals the rolling TRUE
  # value vector at every period (exact, solver-free). This is the precise statement
  # that myopic and rolling optimise the identical objective each period.
  V_true <- t(vapply(seq_len(nrow(sc$b)),
                     function(i) compute_value_vector(sc$b[i, ], sc$lam[i, ]),
                     numeric(n_t)))
  V_froz <- t(vapply(seq_len(nrow(sc$b)),
                     function(i) compute_value_vector(rep(sc$b[i, 1L], n_t),
                                                      rep(sc$lam[i, 1L], n_t)),
                     numeric(n_t)))
  dV <- max(abs(V_true - V_froz))
  report("stationary: myopic frozen value == rolling true value", dV <= TOL * max(1, max(abs(V_true))),
         sprintf("(max |dV| = %.2e)", dV))

  # Informational ONLY (not gated): realised myopic vs rolling under the production
  # cap. The stationary null is maximally degenerate, so the capped solver may return
  # different co-optimal schedules; the realised gap should be a negligible fraction
  # of the climate-scenario effect. Reported for transparency.
  myo  <- run_myopic_ilp(sc$b, sc$lam, sc$cost, bud)
  roll <- run_rolling_horizon(sc$b, sc$lam, sc$cost, bud)
  Jm <- evaluate_policy(myo,  sc$b, sc$lam, sc$cost)$J
  Jr <- evaluate_policy(roll, sc$b, sc$lam, sc$cost)$J
  Jb <- evaluate_policy(rep(NA_integer_, nrow(sc$b)), sc$b, sc$lam, sc$cost)$J
  va_r <- Jr - Jb
  cat(sprintf("    [info] realised: myopic J = %.2f  rolling J = %.2f  Δ = %.2f (%s of rolling value_added)\n",
              Jm, Jr, Jr - Jm,
              if (abs(va_r) < 1e-9) "n/a" else sprintf("%.3f%%", 100 * (Jr - Jm) / va_r)))
}


# ── P4. ordering: rolling >= myopic on every scenario (no order violations) ───####
cat("\nP4  ordering rolling >= myopic across all scenarios\n")
{
  # Run in the ACTIVE SPEND_DOWN_MODE (the production configuration): the foresight
  # ordering must hold for the policy as actually deployed. Spend-down lifts myopic
  # toward rolling (narrows the gap) but cannot make myopic exceed the foresight
  # optimum, so rolling >= myopic should still hold with margin.
  ok_all <- TRUE
  for (sc_name in names(SCENARIOS)) {
    sc  <- build_scenario_matrices(data_panel, SCENARIOS[[sc_name]])
    bud <- make_budget(sc$cost)
    Jm <- evaluate_policy(run_myopic_ilp(sc$b, sc$lam, sc$cost, bud),
                          sc$b, sc$lam, sc$cost)$J
    Jr <- evaluate_policy(run_rolling_horizon(sc$b, sc$lam, sc$cost, bud),
                          sc$b, sc$lam, sc$cost)$J
    Jb <- evaluate_policy(rep(NA_integer_, nrow(sc$b)), sc$b, sc$lam, sc$cost)$J
    va_r <- Jr - Jb
    # Rolling is the foresight optimum, so rolling >= myopic should hold up to solver
    # suboptimality. Flag a violation only if myopic exceeds rolling by more than
    # VA_TOL of value_added — beyond capped-solver noise, indicating a real order
    # violation (e.g. a survival-clock error letting myopic genuinely beat rolling).
    margin <- Jr - Jm
    viol <- margin < -VA_TOL * max(abs(va_r), 1)
    if (viol) ok_all <- FALSE
    cat(sprintf("    %-12s rolling = %.2f  myopic = %.2f  margin = %+.3f%% of value_added  %s\n",
                sc_name, Jr, Jm,
                if (abs(va_r) < 1e-9) 0 else 100 * margin / abs(va_r),
                ifelse(viol, "<-- VIOLATION", "ok")))
  }
  report("rolling >= myopic on every scenario (up to solver tolerance)", ok_all)
}


# ── P5. LP relaxation is a valid upper bound on the ILP (real scenario) ───────####
cat("\nP5  LP relaxation >= ILP (real scenario rcp45_wet)\n")
{
  sc  <- build_scenario_matrices(data_panel, SCENARIOS[["rcp45_wet"]])
  bud <- make_budget(sc$cost)
  V <- t(vapply(seq_len(nrow(sc$b)),
                function(i) compute_value_vector(sc$b[i, ], sc$lam[i, ]),
                numeric(ncol(sc$b))))
  ip <- solve_acquisition_ilp(V, sc$cost, bud, seq_len(nrow(V)), seq_len(ncol(V)),
                              relax = FALSE, spend_down = "off")
  lp <- solve_acquisition_ilp(V, sc$cost, bud, seq_len(nrow(V)), seq_len(ncol(V)),
                              relax = TRUE)
  report("LP relaxation >= ILP optimum",
         lp$objval >= ip$objval - 1e-6 * max(1, abs(ip$objval)),
         sprintf("(LP = %.4f, ILP = %.4f, gap = %.4f%%)",
                 lp$objval, ip$objval, 100 * (lp$objval - ip$objval) / ip$objval))
}


# ── Summary ───────────────────────────────────────────────────────────────────####
cat("\n══ VALIDATION SUMMARY ══\n")
n_pass <- sum(unlist(checks)); n_tot <- length(checks)
cat(sprintf("  %d / %d checks passed\n", n_pass, n_tot))
if (n_pass < n_tot) {
  cat("  FAILED:", paste(names(checks)[!unlist(checks)], collapse = "; "), "\n")
  stop("Validation failed — do not trust 08__run_models.R output until resolved.")
} else {
  cat("  All correctness checks passed.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# STUB — MDP comparison (DEFERRED)
# ══════════════════════════════════════════════════════════════════════════════
# An exact finite-horizon MDP over the joint landscape state is intractable here:
# the state is the protected/unprotected (and converted) status of every EAU, so
# the state space is ~3^n_eau — already astronomical beyond n_eau ~ 12-15, and
# hopeless at 841. A like-for-like MDP benchmark would require either (a) a drastic
# toy-sized landscape, or (b) per-EAU decomposition with an approximate value
# function. Whether to pursue this is an open decision. The hook below documents
# the intended interface so it can be filled in later without disturbing 07/08.
#
# run_mdp_optimal <- function(b_mat, lam_mat, cost_mat, budget, delta = DELTA) {
#   stop("MDP comparison not implemented — see STUB note in 09__validation.R.")
# }
# ══════════════════════════════════════════════════════════════════════════════
