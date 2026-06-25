# ══════════════════════════════════════════════════════════════════════════════
# 12__test_tiebreak.R
#   Does an earlier-is-better tie-break close the degenerate myopic solves WITHOUT
#   sacrificing true conservation value?
# ══════════════════════════════════════════════════════════════════════════════
#
# The myopic frozen-hazard objective is near-symmetric in TIMING (buy-now vs
# buy-later ties to within rounding), so Gurobi finds a great incumbent instantly
# but cannot prove optimality. A tiny "earlier is better" perturbation,
#       V_pert[i,t] = V[i,t] - rho * (t - 1),     rho = kappa * max(V),
# makes earlier strictly preferred, collapsing the timing ties to a unique optimum.
#
# Sizing: rho scales with each solve's own magnitude, so the MAX possible total
# distortion is ~1.5e-5 of value_added whether value_added is 4.6e5 (myopic t=1)
# or 7.9 (myopic t=2). Two things are checked per solve:
#   (1) EFFECTIVE? status -> OPTIMAL, runtime drops sharply.
#   (2) SAFE?      the TRUE (un-perturbed) value of the chosen schedule still lands
#                  inside the reference interval [un-perturbed incumbent, bound].
#
# If both hold, the tie-break is the right fix for Problem A. It must be applied
# IDENTICALLY to myopic and rolling in production; this script only proves it out.
# Script 07 is untouched (solver overridden locally as before).
# ══════════════════════════════════════════════════════════════════════════════

source("07__ilp_core.R")

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
  if (!requireNamespace("gurobi", quietly = TRUE)) stop("gurobi package required.")
  cost_scale <- 1e6
  A <- Matrix::sparseMatrix(
    i = c(row_once, row_budg), j = c(seq_len(nvar), seq_len(nvar)),
    x = c(rep(1, nvar), costs / cost_scale), dims = c(n_av + n_pe, nvar))
  model <- list(A = A, obj = obj, sense = rep("<=", n_av + n_pe),
                rhs = c(rep(1, n_av), budget[periods] / cost_scale),
                lb = rep(0, nvar), ub = rep(1, nvar),
                vtype = if (relax) "C" else "B", modelsense = "max")
  res <- gurobi::gurobi(model, params = modifyList(list(OutputFlag = 0), extra_params))
  gx <- function(v, d = NA_real_) if (is.null(v)) d else v
  chosen <- if (relax) integer(0) else which(res$x > 0.5)
  list(picks = grid[chosen, , drop = FALSE], objval = res$objval, x = res$x, grid = grid,
       objbound = gx(res$objbound), mipgap = gx(res$mipgap),
       runtime = gx(res$runtime, 0), status = gx(res$status, "UNKNOWN"))
}

# ── Setup ──────────────────────────────────────────────────────────────────────
data_panel <- readRDS("input_data/data_panel.rds")
sc  <- build_scenario_matrices(data_panel, SCENARIOS[["stationary"]])
bud <- make_budget(sc$cost)
n_eau <- nrow(sc$b); n_t <- ncol(sc$b)

myopic_V <- function(tnow) {
  Vf <- matrix(0, n_eau, n_t)
  for (i in seq_len(n_eau))
    Vf[i, ] <- compute_value_vector(rep(sc$b[i, tnow], n_t), rep(sc$lam[i, tnow], n_t))
  Vf
}
V1 <- myopic_V(1)   # freezes baseline hazard  -> the 866s solve
V2 <- myopic_V(2)   # freezes epsilon hazard    -> the stalled solve

# ── Reference: un-perturbed 60s solve (incumbent + bound to bracket the optimum) ─
reference <- function(V, avail, periods, tl = 60) {
  s <- solve_acquisition_ilp(V, sc$cost, bud, avail, periods,
                             extra_params = list(OutputFlag = 0, TimeLimit = tl))
  c(incumbent = s$objval, bound = s$objbound)
}

# ── Tie-break solve, then score the chosen schedule on the TRUE objective ───────
tiebreak <- function(label, V, avail, periods, kappa, tl = 180) {
  rho    <- kappa * max(V[is.finite(V)])
  V_pert <- sweep(V, 2, rho * (seq_len(n_t) - 1), `-`)   # subtract rho*(period-1)
  s <- solve_acquisition_ilp(V_pert, sc$cost, bud, avail, periods,
                             extra_params = list(OutputFlag = 0, TimeLimit = tl))
  p <- s$picks
  true_va <- if (nrow(p) > 0) sum(V[cbind(p$eau_idx, p$period)]) else 0  # TRUE value
  data.frame(objective = label, kappa = kappa, rho = signif(rho, 3),
             status = s$status, runtime_s = round(s$runtime, 1),
             pert_gap_pct = signif(100 * s$mipgap, 3),
             true_value_added = true_va, n_picks = nrow(p),
             stringsAsFactors = FALSE)
}

all_eau <- seq_len(n_eau); per1 <- seq_len(n_t); per2 <- 2:n_t

cat("\n── Reference intervals (un-perturbed, 60s cap) ────────────────────────────\n")
ref1 <- reference(V1, all_eau, per1); ref2 <- reference(V2, all_eau, per2)
cat(sprintf("  myopic t=1 : true optimum in [%.4f , %.4f]\n", ref1["incumbent"], ref1["bound"]))
cat(sprintf("  myopic t=2 : true optimum in [%.6f , %.6f]\n", ref2["incumbent"], ref2["bound"]))

cat("\n── Tie-break solves ───────────────────────────────────────────────────────\n")
res <- rbind(
  tiebreak("myopic t=1", V1, all_eau, per1, kappa = 1e-6),
  tiebreak("myopic t=1", V1, all_eau, per1, kappa = 1e-8),
  tiebreak("myopic t=2", V2, all_eau, per2, kappa = 1e-6)
)
print(res, row.names = FALSE)

cat("\nRead it like this:\n")
cat("  EFFECTIVE? status OPTIMAL with a small runtime_s = the grind is gone.\n")
cat("  SAFE?      true_value_added must sit INSIDE the matching reference interval\n")
cat("             above. If it does, the nudge cost no real conservation value;\n")
cat("             it only chose a canonical (earliest) representative among ties.\n")
cat("  kappa 1e-6 vs 1e-8: if both close and both stay in-interval, the result is\n")
cat("             robust to the nudge size (good). If 1e-6 drifts but 1e-8 doesn't,\n")
cat("             that tells us how small rho must be.\n")

out_dir <- "output_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(res, file.path(out_dir, "tiebreak_test.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(out_dir, "tiebreak_test.csv")))
