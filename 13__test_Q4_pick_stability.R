# ══════════════════════════════════════════════════════════════════════════════
# 13__test_Q4_pick_stability.R
#   Q4: Is the IMPLEMENTED (period-1) myopic pick set invariant to solver effort?
# ══════════════════════════════════════════════════════════════════════════════
#
# WHY THIS IS THE DECISION-DOMINATING MEASUREMENT
# ------------------------------------------------------------------------------
# Myopic implements ONLY the current period's purchases and re-solves next decade
# (07: run_myopic_ilp). The expensive optimality proof on the t=1 stationary solve
# is spent ordering a periods-2..9 schedule that is thrown away on the next line.
# So the question that actually matters for Problem A is NOT "exact ties vs natural
# near-flatness" (the H1/H2 fork) but simply:
#
#     Does the set of parcels myopic ACQUIRES AT t=1 change as we give the solver
#     more effort (longer TimeLimit) or a looser MIPGap?
#
# If that set is invariant, then a TimeLimit (or a looser gap) is provably harmless
# to what the policy DOES and to its realized, scored-on-reality value -- regardless
# of the H1/H2 mix -- and Problem A closes with a one-line solver-param change.
# If it wobbles, we need to know by how much, whether the swapped parcels are ties,
# and whether the wobble touches the realized number.
#
# WHAT THIS SCRIPT DOES
#   Experiment 1 (effort ladder, DETERMINISTIC: Threads=1 + fixed Seed)
#     Solves the exact myopic t=1 stationary instance under a ladder of stopping
#     rules {5s,15s,60s @ default gap; 1e-3 and 1e-2 gaps @ ample time; a longer
#     best-effort REFERENCE}. For each it reports status/runtime/gap, the period-1
#     pick set, set-difference vs the reference, and the TRUE-future value of the
#     period-1 implementation. Pinning Threads/Seed makes differences attributable
#     to the stopping rule, not to thread races (a confound that has misled before).
#
#   Experiment 2 (reproducibility at a fixed cap)
#     Repeats one capped solve N times, once pinned (Threads=1, fixed Seed -> should
#     be identical) and once UNPINNED (all threads, varying Seed -> production-style).
#     Answers: is a capped solve reproducible run-to-run? If not, capping injects
#     policy irreproducibility and a loosened gap is the safer lever.
#
# NON-INVASIVE: solve_acquisition_ilp is overridden LOCALLY with a faithful superset
# of the 07 version (identical model assembly + 1e6 cost scaling; accepts extra
# Gurobi params; returns objbound/mipgap/runtime/status). Script 07 is untouched.
#
# EXACTNESS NOTE: with FREEZE_COL = 1, avail = all parcels and periods = 1..n_t, so
# this is byte-identical to run_myopic_ilp's FIRST iteration. (For FREEZE_COL > 1 the
# real loop would shrink `avail` by the parcels already bought at t<FREEZE_COL; this
# script keeps avail = all, which probes the same solve SHAPE but is not the exact
# loop state. Keep FREEZE_COL = 1 for the decision-dominating, exact case.)
# ══════════════════════════════════════════════════════════════════════════════

source("07__ilp_core.R")

# ══ PARAMETERS ═══════════════════════════════════════════════════════════════════
FREEZE_COL <- 1L     # period myopic freezes AND implements (1 = 2020). Exact at 1.
THREADS    <- 1L     # Exp-1 ladder: single-thread => deterministic given the cap
SEED       <- 1L     # Exp-1 ladder: fixed seed   => removes run-to-run jitter
REPRO_N    <- 3L     # Exp-2: number of repeats
REPRO_TL   <- 30     # Exp-2: the fixed cap (s) to repeat

# Effort ladder. The LAST entry is the REFERENCE (most effort); every other config
# is compared against it. Crank REFERENCE's TimeLimit if you want a tighter anchor.
configs <- list(
  list(label = "cap_5s_gap1e-4",  TimeLimit =   5, MIPGap = 1e-4),
  list(label = "cap_15s_gap1e-4", TimeLimit =  15, MIPGap = 1e-4),
  list(label = "cap_60s_gap1e-4", TimeLimit =  60, MIPGap = 1e-4),
  list(label = "gap_1e-3_ample",  TimeLimit = 120, MIPGap = 1e-3),
  list(label = "gap_1e-2_ample",  TimeLimit = 120, MIPGap = 1e-2),
  list(label = "REFERENCE_300s",  TimeLimit = 300, MIPGap = 1e-4)
)
REF_LABEL <- "REFERENCE_300s"
# ══════════════════════════════════════════════════════════════════════════════


# ── Instrumented solver: faithful superset of the 07 version ────────────────────
solve_acquisition_ilp <- function(V_mat, cost_mat, budget, avail, periods,
                                   relax = FALSE, extra_params = list()) {
  if (length(avail) == 0 || length(periods) == 0)
    return(list(picks = data.frame(eau_idx = integer(0), period = integer(0)),
                objval = 0, x = numeric(0), grid = NULL,
                objbound = NA_real_, mipgap = NA_real_, runtime = 0, status = "EMPTY"))

  grid  <- expand.grid(eau_idx = avail, period = periods)
  nvar  <- nrow(grid)
  obj   <- V_mat[cbind(grid$eau_idx, grid$period)]
  costs <- cost_mat[cbind(grid$eau_idx, grid$period)]

  n_av <- length(avail); n_pe <- length(periods)
  row_once <- match(grid$eau_idx, avail)
  row_budg <- n_av + match(grid$period, periods)

  if (!requireNamespace("gurobi", quietly = TRUE))
    stop("The 'gurobi' R package is required (activate an academic licence).")

  cost_scale <- 1e6
  costs_s    <- costs / cost_scale

  A <- Matrix::sparseMatrix(
    i = c(row_once, row_budg), j = c(seq_len(nvar), seq_len(nvar)),
    x = c(rep(1, nvar), costs_s), dims = c(n_av + n_pe, nvar)
  )
  sense <- rep("<=", n_av + n_pe)
  rhs   <- c(rep(1, n_av), budget[periods] / cost_scale)

  model <- list(A = A, obj = obj, sense = sense, rhs = rhs,
                lb = rep(0, nvar), ub = rep(1, nvar),
                vtype = if (relax) "C" else "B", modelsense = "max")

  params <- modifyList(list(OutputFlag = 0), extra_params)
  res    <- gurobi::gurobi(model, params = params)

  gx <- function(v, d = NA_real_) if (is.null(v)) d else v
  x <- res$x
  chosen <- if (relax) integer(0) else which(x > 0.5)
  list(picks = grid[chosen, , drop = FALSE], objval = res$objval, x = x, grid = grid,
       objbound = gx(res$objbound), mipgap = gx(res$mipgap),
       runtime  = gx(res$runtime, 0), status = gx(res$status, "UNKNOWN"))
}

# small set helper
symdiff <- function(a, b) union(setdiff(a, b), setdiff(b, a))


# ── Build the exact myopic t=1 stationary instance ──────────────────────────────
data_panel <- readRDS("input_data/data_panel.rds")
sc    <- build_scenario_matrices(data_panel, SCENARIOS[["stationary"]])
bud   <- make_budget(sc$cost)
n_eau <- nrow(sc$b); n_t <- ncol(sc$b)

# Frozen-belief objective EXACTLY as run_myopic_ilp builds it for tnow = FREEZE_COL:
# freeze benefit & hazard at the current period, flat across the whole horizon.
Vf <- matrix(0, n_eau, n_t)
for (i in seq_len(n_eau)) {
  b_frozen   <- rep(sc$b[i,   FREEZE_COL], n_t)
  lam_frozen <- rep(sc$lam[i, FREEZE_COL], n_t)
  Vf[i, ] <- compute_value_vector(b_frozen, lam_frozen)
}
avail   <- seq_len(n_eau)        # at t=1 every parcel is available
periods <- FREEZE_COL:n_t        # solve over the remaining horizon

# Do-nothing landscape total on the TRUE stationary future (denominator for value_added)
J_baseline <- evaluate_policy(rep(NA_integer_, n_eau), sc$b, sc$lam, sc$cost)$J


# ── Solve one config, return the implemented (period-1) set + its true-future value ─
solve_cfg <- function(cfg) {
  ep <- list(OutputFlag = 0, Threads = THREADS, Seed = SEED,
             TimeLimit = cfg$TimeLimit, MIPGap = cfg$MIPGap)
  sol <- solve_acquisition_ilp(Vf, sc$cost, bud, avail, periods, extra_params = ep)

  # IMPLEMENTED set = picks at the freeze/implement period only (the rest is discarded)
  p1 <- sort(sol$picks$eau_idx[sol$picks$period == FREEZE_COL])

  # Score "protect exactly these from FREEZE_COL onward, nothing else" on the TRUE future
  acq <- rep(NA_integer_, n_eau); if (length(p1)) acq[p1] <- FREEZE_COL
  trueJ <- evaluate_policy(acq, sc$b, sc$lam, sc$cost)$J

  list(label = cfg$label, status = sol$status, runtime = sol$runtime,
       gap_pct = 100 * sol$mipgap, incumbent = sol$objval,
       n_pick1 = length(p1), pick1 = p1, trueJ = trueJ)
}


# ══ EXPERIMENT 1 — effort ladder (deterministic) ════════════════════════════════
cat("\n══════════════════════════════════════════════════════════\n")
cat("  Q4 EXPERIMENT 1: does solver effort change the t=1 picks?\n")
cat(sprintf("  (stationary, freeze period %d; Threads=%d, Seed=%d)\n",
            FREEZE_COL, THREADS, SEED))
cat("══════════════════════════════════════════════════════════\n")

runs <- lapply(configs, function(cfg) {
  cat(sprintf("  solving %-18s (TimeLimit=%gs, MIPGap=%g) ...\n",
              cfg$label, cfg$TimeLimit, cfg$MIPGap))
  solve_cfg(cfg)
})
names(runs) <- vapply(configs, `[[`, "", "label")

if (!REF_LABEL %in% names(runs)) stop("REF_LABEL not found among configs.")
ref <- runs[[REF_LABEL]]
ref_VA <- ref$trueJ - J_baseline   # reference value_added on the true future

# ── Comparison table vs the reference ───────────────────────────────────────────
cmp <- do.call(rbind, lapply(runs, function(r) {
  inter   <- length(intersect(r$pick1, ref$pick1))
  uni     <- length(union(r$pick1, ref$pick1))
  data.frame(
    config       = r$label,
    status       = r$status,
    runtime_s    = round(r$runtime, 1),
    gap_pct      = signif(r$gap_pct, 3),
    incumbent    = signif(r$incumbent, 7),
    n_pick1      = r$n_pick1,
    n_added      = length(setdiff(r$pick1, ref$pick1)),   # in this config, not in ref
    n_dropped    = length(setdiff(ref$pick1, r$pick1)),   # in ref, not in this config
    jaccard      = if (uni == 0) 1 else round(inter / uni, 4),
    trueJ        = round(r$trueJ, 2),
    d_trueJ      = round(r$trueJ - ref$trueJ, 4),
    d_trueVA_pct = if (abs(ref_VA) < 1e-12) NA_real_
                   else round(100 * (r$trueJ - ref$trueJ) / ref_VA, 4),
    stringsAsFactors = FALSE
  )
}))
rownames(cmp) <- NULL

cat("\n── Period-1 pick set vs REFERENCE (", REF_LABEL, ") ──\n", sep = "")
cat("   n_added/n_dropped vs reference; jaccard=1 & both 0 => identical set.\n")
cat("   d_trueVA_pct = change in realized value_added from those pick differences.\n\n")
print(cmp, row.names = FALSE)

# ── Tie structure of any swapped parcels (added or dropped anywhere in the ladder) ─
diff_idx <- sort(unique(unlist(lapply(runs, function(r) symdiff(r$pick1, ref$pick1)))))
if (length(diff_idx) > 0) {
  tie_tbl <- data.frame(
    eau_id        = sc$eau_ids[diff_idx],
    V_freeze      = signif(Vf[cbind(diff_idx, rep(FREEZE_COL, length(diff_idx)))], 6),
    cost_freeze   = signif(sc$cost[cbind(diff_idx, rep(FREEZE_COL, length(diff_idx)))], 6),
    in_reference  = diff_idx %in% ref$pick1,
    stringsAsFactors = FALSE
  )
  tie_tbl <- tie_tbl[order(-tie_tbl$V_freeze), ]
  cat("\n── Parcels that move across the ladder (tie diagnosis) ──\n")
  cat("   If swapped parcels share (near-)identical V_freeze, the differences are\n")
  cat("   ties => harmless. If V_freeze differs materially, effort is buying real value.\n\n")
  print(utils::head(tie_tbl, 25), row.names = FALSE)
  if (nrow(tie_tbl) > 25) cat(sprintf("   ... (%d more)\n", nrow(tie_tbl) - 25))
} else {
  cat("\n── No parcel ever moves: every config's t=1 set == reference. ──\n")
}


# ══ EXPERIMENT 2 — reproducibility at a fixed cap ═══════════════════════════════
cat("\n\n══════════════════════════════════════════════════════════\n")
cat(sprintf("  Q4 EXPERIMENT 2: is a %gs-capped solve reproducible?\n", REPRO_TL))
cat("══════════════════════════════════════════════════════════\n")

repro_sets <- function(threads, seeds, label) {
  cat(sprintf("  %s: Threads=%d, %d repeats ...\n", label, threads, length(seeds)))
  lapply(seeds, function(s) {
    ep <- list(OutputFlag = 0, Threads = threads, Seed = s,
               TimeLimit = REPRO_TL, MIPGap = 1e-4)
    sol <- solve_acquisition_ilp(Vf, sc$cost, bud, avail, periods, extra_params = ep)
    sort(sol$picks$eau_idx[sol$picks$period == FREEZE_COL])
  })
}

summarise_sets <- function(sets) {
  base   <- sets[[1]]
  all_id <- all(vapply(sets, function(s) setequal(s, base), TRUE))
  maxsym <- max(vapply(sets, function(s) length(symdiff(s, base)), 0L))
  c(all_identical = all_id, max_symdiff_vs_run1 = maxsym, n_pick1_run1 = length(base))
}

det_sets  <- repro_sets(1L, rep(SEED, REPRO_N),  "PINNED   (Threads=1, fixed Seed)")
ndet_sets <- repro_sets(0L, seq_len(REPRO_N),    "UNPINNED (all threads, varying Seed)")

det_s  <- summarise_sets(det_sets)
ndet_s <- summarise_sets(ndet_sets)

cat("\n  PINNED   (Threads=1, fixed Seed)  -> should be identical across repeats:\n")
cat(sprintf("    all_identical = %s | max symdiff vs run1 = %d | n_pick1 = %d\n",
            det_s["all_identical"] == 1, det_s["max_symdiff_vs_run1"], det_s["n_pick1_run1"]))
cat("  UNPINNED (all threads, varying Seed) -> production-style run-to-run wobble:\n")
cat(sprintf("    all_identical = %s | max symdiff vs run1 = %d | n_pick1 = %d\n",
            ndet_s["all_identical"] == 1, ndet_s["max_symdiff_vs_run1"], ndet_s["n_pick1_run1"]))


# ══ VERDICT ═════════════════════════════════════════════════════════════════════
cat("\n\n══════════════════════════════════════════════════════════\n")
cat("  VERDICT\n")
cat("══════════════════════════════════════════════════════════\n")

ladder_stable <- all(cmp$n_added == 0 & cmp$n_dropped == 0)
max_dVA <- max(abs(cmp$d_trueVA_pct), na.rm = TRUE)

if (ladder_stable && ndet_s["all_identical"] == 1) {
  cat("  STABLE. The implemented t=1 pick set is invariant to solver effort across\n")
  cat("  the whole ladder AND reproducible under production threading. A TimeLimit\n")
  cat("  (or a looser MIPGap) is provably harmless to what myopic does at t=1 and to\n")
  cat("  its realized value. Problem A can be closed with a solver-param change; the\n")
  cat("  H1/H2 fork need not be resolved.\n")
} else if (ladder_stable) {
  cat("  EFFORT-STABLE BUT THREAD-SENSITIVE. More effort doesn't change the picks,\n")
  cat("  but unpinned runs disagree (max symdiff vs run1 = ", ndet_s["max_symdiff_vs_run1"],
      ").\n", sep = "")
  cat("  => Capping is safe ONLY if Threads/Seed are pinned in production, OR prefer a\n")
  cat("     looser gap. Check the tie table: thread-swapped parcels should be ties.\n")
} else {
  cat(sprintf("  WOBBLE. The t=1 pick set changes with effort (max realized value_added\n"))
  cat(sprintf("  shift = %.4f%% of reference value_added). The picks are NOT effort-\n", max_dVA))
  cat("  invariant: inspect the tie table. If the moving parcels are genuine ties\n")
  cat("  (equal V_freeze) and d_trueVA_pct is negligible, a cap is still defensible;\n")
  cat("  if not, effort is buying real value and a cap/loose gap is NOT safe -- the\n")
  cat("  H1/H2 question then genuinely matters.\n")
}

cat(sprintf("\n  (reference value_added on true future = %.4f; J_baseline = %.2f)\n",
            ref_VA, J_baseline))


# ── Save for the record ─────────────────────────────────────────────────────────
out_dir <- "output_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(cmp, file.path(out_dir, "Q4_pick_stability_ladder.csv"), row.names = FALSE)
saveRDS(list(runs = runs, reference = REF_LABEL, J_baseline = J_baseline,
             det_sets = det_sets, ndet_sets = ndet_sets,
             freeze_col = FREEZE_COL),
        file.path(out_dir, "Q4_pick_stability.rds"))
cat(sprintf("\nSaved:\n  %s\n  %s\n",
            file.path(out_dir, "Q4_pick_stability_ladder.csv"),
            file.path(out_dir, "Q4_pick_stability.rds")))
