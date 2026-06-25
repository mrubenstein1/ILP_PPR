# ══════════════════════════════════════════════════════════════════════════════
# 11__diagnose_myopic_vs_rolling_objective.R
#   Why does P3 grind when the script-10 single-solve read closed in 29s?
# ══════════════════════════════════════════════════════════════════════════════
#
# HYPOTHESIS being tested:
#   Script 10 solved the ROLLING first-iteration objective (the TRUE stationary V).
#   That is well-conditioned and fast. But P3 runs MYOPIC first, and myopic FREEZES
#   benefit and hazard at the current period and projects them flat across the whole
#   horizon. Under the stationary scenario that frozen objective is far more
#   degenerate, and IS the slow solve(s) in the P3 console output.
#
# This script builds the stationary scenario and, WITHOUT running the full loops,
# constructs the three objective vectors that the policies actually hand Gurobi:
#   (A) rolling, t=1        : true V                    -> what script 10 solved
#   (B) myopic, t=1         : freeze baseline hazard    -> log block 1 (the 866s solve)
#   (C) myopic, t=2         : freeze epsilon hazard      -> log block 2 (the stall)
# For each it prints the objective SHAPE (range + how many near-zero entries) and
# does a TIME-CAPPED solve so you can see which ones close and which grind — none
# can hang, because every solve has a 60s TimeLimit.
#
# Nothing here touches script 07; the solver is overridden locally (faithful
# superset that returns objbound/mipgap/runtime/status and accepts extra params).
# ══════════════════════════════════════════════════════════════════════════════

source("07__ilp_core.R")

TIME_LIMIT <- 60   # seconds per solve; enough to reveal "closes" vs "stalls"

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
    stop("The 'gurobi' R package is required.")

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

  params <- modifyList(list(OutputFlag = 0), extra_params)   # quiet by default
  res    <- gurobi::gurobi(model, params = params)

  gx <- function(v, d = NA_real_) if (is.null(v)) d else v
  x <- res$x
  chosen <- if (relax) integer(0) else which(x > 0.5)
  list(picks = grid[chosen, , drop = FALSE], objval = res$objval, x = x, grid = grid,
       objbound = gx(res$objbound), mipgap = gx(res$mipgap),
       runtime  = gx(res$runtime, 0), status = gx(res$status, "UNKNOWN"))
}

# ── Build the stationary scenario ──────────────────────────────────────────────
data_panel <- readRDS("input_data/data_panel.rds")
sc  <- build_scenario_matrices(data_panel, SCENARIOS[["stationary"]])
bud <- make_budget(sc$cost)
n_eau <- nrow(sc$b); n_t <- ncol(sc$b)

# ── Construct the three objective vectors exactly as the policies do ───────────
# (A) Rolling, t=1: true V over the real stationary trajectories
V_roll <- t(vapply(seq_len(n_eau),
                   function(i) compute_value_vector(sc$b[i, ], sc$lam[i, ]),
                   numeric(n_t)))

# Myopic freezes benefit & hazard at column `tnow`, flat across the horizon
myopic_V <- function(tnow) {
  Vf <- matrix(0, n_eau, n_t)
  for (i in seq_len(n_eau)) {
    b_frozen   <- rep(sc$b[i,   tnow], n_t)
    lam_frozen <- rep(sc$lam[i, tnow], n_t)
    Vf[i, ] <- compute_value_vector(b_frozen, lam_frozen)
  }
  Vf
}
V_myo1 <- myopic_V(1)   # freezes the 2020 BASELINE hazard
V_myo2 <- myopic_V(2)   # freezes the 2030 hazard (epsilon, under stationary)

# ── Helper: describe an objective vector's shape ───────────────────────────────
describe_obj <- function(label, V_mat, avail, periods) {
  vals <- V_mat[cbind(rep(avail, length(periods)),
                      rep(periods, each = length(avail)))]
  vals <- vals[is.finite(vals)]
  cat(sprintf("\n%s\n", label))
  cat(sprintf("  n vars            : %d\n", length(vals)))
  cat(sprintf("  objective range   : [%.3e , %.3e]\n", min(vals), max(vals)))
  cat(sprintf("  median |value|    : %.3e\n", median(abs(vals))))
  cat(sprintf("  share < 1e-4      : %.1f%%   (near-zero swarm)\n",
              100 * mean(abs(vals) < 1e-4)))
  cat(sprintf("  share < 1e-6      : %.1f%%\n", 100 * mean(abs(vals) < 1e-6)))
}

# ── Helper: time-capped solve, report close-vs-stall ───────────────────────────
capped_solve <- function(label, V_mat, avail, periods) {
  t0  <- Sys.time()
  sol <- solve_acquisition_ilp(V_mat, sc$cost, bud, avail, periods,
                               relax = FALSE,
                               extra_params = list(OutputFlag = 0,
                                                   TimeLimit = TIME_LIMIT))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  data.frame(
    objective    = label,
    incumbent    = sol$objval,
    bound        = sol$objbound,
    rel_gap_pct  = 100 * sol$mipgap,
    runtime_s    = round(sol$runtime, 1),
    wall_s       = round(wall, 1),
    status       = sol$status,           # OPTIMAL = closed; TIME_LIMIT = stalled
    stringsAsFactors = FALSE
  )
}

all_eau  <- seq_len(n_eau)
per_full <- seq_len(n_t)      # periods 1..9 (t = 0..8)
per_t2   <- 2:n_t             # periods 2..9 (myopic's second iteration)

cat("\n══════════════════════════════════════════════════════════\n")
cat("  OBJECTIVE SHAPE (what each policy hands Gurobi)\n")
cat("══════════════════════════════════════════════════════════")
describe_obj("(A) ROLLING t=1  (true V)            -> script-10 read",
             V_roll, all_eau, per_full)
describe_obj("(B) MYOPIC  t=1  (freeze baseline)   -> P3 log block 1",
             V_myo1, all_eau, per_full)
describe_obj("(C) MYOPIC  t=2  (freeze epsilon)    -> P3 log block 2",
             V_myo2, all_eau, per_t2)

cat("\n\n══════════════════════════════════════════════════════════\n")
cat(sprintf("  TIME-CAPPED SOLVES (%ds limit each; none can hang)\n", TIME_LIMIT))
cat("══════════════════════════════════════════════════════════\n")
res <- rbind(
  capped_solve("(A) rolling t=1",  V_roll, all_eau, per_full),
  capped_solve("(B) myopic  t=1",  V_myo1, all_eau, per_full),
  capped_solve("(C) myopic  t=2",  V_myo2, all_eau, per_t2)
)
print(res, row.names = FALSE)

cat("\nReading this:\n")
cat("  status OPTIMAL  + small runtime_s   -> closes fine (not the problem).\n")
cat("  status TIME_LIMIT (hit the cap)     -> this objective is the grind.\n")
cat("  A tiny rel_gap_pct with TIME_LIMIT  -> great incumbent, can't PROVE it:\n")
cat("                                         the near-symmetry / degeneracy.\n")
cat("  Compare the 'share < 1e-4' column above: the bigger the near-zero swarm,\n")
cat("  the more the solve stalls -> that is exactly what pruning targets.\n")

out_dir <- "output_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(res, file.path(out_dir, "myopic_vs_rolling_objective.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(out_dir, "myopic_vs_rolling_objective.csv")))
