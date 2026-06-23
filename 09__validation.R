# ══════════════════════════════════════════════════════════════════════════════
# 09__validation.R  —  Correctness checks for the acquisition models
# ══════════════════════════════════════════════════════════════════════════════
#
# Run this FIRST (with Gurobi active) to certify the formulation on your machine.
# It mirrors the property tests that were validated in a separate solver harness:
#
#   P1  ILP == brute-force optimum            (tiny synthetic instance)
#   P2  rolling == full-horizon single solve  (real scenario; deterministic foresight)
#   P3  stationary null: myopic == rolling    (stationary scenario; exact)
#   P4  ordering: rolling >= myopic           (all scenarios; no order violations)
#   P5  LP relaxation >= ILP                   (real scenario; valid bound)
#
# Each check appends a PASS/FAIL to `checks`; the script stops if any fails.
#
# WHY a tiny synthetic instance for P1.  Brute force enumerates (n_t+1)^n_eau
# assignments, so it is only tractable for a handful of EAUs. The real panel
# (n_eau = 841) is far beyond enumeration — P1 certifies the FORMULATION on a small
# instance; P2/P3/P5 then exercise the same solver on the real data.
# ══════════════════════════════════════════════════════════════════════════════

source("07__ilp_core.R")

TOL    <- 1e-6
checks <- list()
report <- function(name, ok, detail = "") {
  checks[[name]] <<- isTRUE(ok)
  cat(sprintf("  [%s] %-46s %s\n", ifelse(isTRUE(ok), "PASS", "FAIL"), name, detail))
}

data_panel <- readRDS("input_data/data_panel.rds")


# ── P1. ILP == brute-force optimum (tiny synthetic instance) ──────────────────####
cat("\nP1  ILP vs brute force (synthetic)\n")
{
  inst <- make_synthetic_instance(n_eau = 6, n_t = 5, seed = 11, severity = 1.5)
  V <- t(vapply(seq_len(nrow(inst$b)),
                function(i) compute_value_vector(inst$b[i, ], inst$lam[i, ]),
                numeric(ncol(inst$b))))
  bud  <- make_budget(inst$cost, n_eaus = 2)
  ilp  <- solve_acquisition_ilp(V, inst$cost, bud,
                                avail = seq_len(nrow(V)), periods = seq_len(ncol(V)))
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
  roll <- run_rolling_horizon(sc$b, sc$lam, sc$cost, bud)
  full <- run_full_horizon(sc$b, sc$lam, sc$cost, bud)
  Jr <- evaluate_policy(roll, sc$b, sc$lam, sc$cost)$J
  Jf <- evaluate_policy(full, sc$b, sc$lam, sc$cost)$J
  report("rolling J equals full-horizon J",
         abs(Jr - Jf) < 1e-4 * abs(Jf),
         sprintf("(rolling = %.4f, full = %.4f)", Jr, Jf))
}


# ── P3. stationary null: myopic == rolling, exactly (stationary scenario) ─────####
cat("\nP3  stationary null (myopic must equal rolling on the stationary scenario)\n")
{
  sc  <- build_scenario_matrices(data_panel, SCENARIOS[["stationary"]])
  bud <- make_budget(sc$cost)
  myo  <- run_myopic_ilp(sc$b, sc$lam, sc$cost, bud)
  roll <- run_rolling_horizon(sc$b, sc$lam, sc$cost, bud)
  Jm <- evaluate_policy(myo,  sc$b, sc$lam, sc$cost)$J
  Jr <- evaluate_policy(roll, sc$b, sc$lam, sc$cost)$J
  same_schedule <- all((is.na(myo) & is.na(roll)) | (myo == roll), na.rm = FALSE)
  report("stationary: myopic J equals rolling J",
         abs(Jm - Jr) < 1e-6 * max(1, abs(Jr)),
         sprintf("(myopic = %.6f, rolling = %.6f)", Jm, Jr))
  report("stationary: myopic and rolling pick identical schedules", same_schedule)
}


# ── P4. ordering: rolling >= myopic on every scenario (no order violations) ───####
cat("\nP4  ordering rolling >= myopic across all scenarios\n")
{
  ok_all <- TRUE
  for (sc_name in names(SCENARIOS)) {
    sc  <- build_scenario_matrices(data_panel, SCENARIOS[[sc_name]])
    bud <- make_budget(sc$cost)
    Jm <- evaluate_policy(run_myopic_ilp(sc$b, sc$lam, sc$cost, bud),
                          sc$b, sc$lam, sc$cost)$J
    Jr <- evaluate_policy(run_rolling_horizon(sc$b, sc$lam, sc$cost, bud),
                          sc$b, sc$lam, sc$cost)$J
    viol <- Jr < Jm - 1e-6 * max(1, abs(Jr))
    if (viol) ok_all <- FALSE
    cat(sprintf("    %-12s rolling = %.2f  myopic = %.2f  %s\n",
                sc_name, Jr, Jm, ifelse(viol, "<-- VIOLATION", "ok")))
  }
  report("rolling >= myopic on every scenario", ok_all)
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
                              relax = FALSE)
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
