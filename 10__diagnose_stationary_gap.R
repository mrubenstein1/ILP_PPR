# ══════════════════════════════════════════════════════════════════════════════
# 10__diagnose_stationary_gap.R  —  Read the certified gap on the stationary solve
# ══════════════════════════════════════════════════════════════════════════════
#
# PURPOSE (one-time diagnostic; NOT part of the production pipeline)
# ------------------------------------------------------------------------------
# The stationary-scenario solve finds an excellent incumbent in seconds, then the
# best BOUND stalls — Gurobi cannot prove optimality because the instance is
# near-symmetric. Before changing anything, this answers one question:
#
#     How far from optimal is the answer we ALREADY have, and is that gap small
#     enough to simply stop and report it?  (Option A2 — "certified bound")
#
# It re-solves the exact expensive instance (stationary, full horizon, all 841
# parcels) but with a short TimeLimit, capturing the incumbent and the best bound.
# Because the bound is stalled, a ~2-minute read returns essentially the same gap
# a 45-minute run would — you are not waiting 45 minutes again.
#
# It runs the solve TWICE — default params, then MIPFocus=3 / Cuts=2 — so the same
# script also tells you, in one shot, whether Option A1 moves the bound at all.
#
# NON-INVASIVE: solve_acquisition_ilp is overridden LOCALLY here with a faithful
# superset of the 07 version (identical model assembly, identical 1e6 cost scaling)
# that merely (a) accepts extra Gurobi params and (b) returns objbound / mipgap /
# runtime / status. Script 07 is left untouched. If you adopt A2 for production,
# fold those same three additions into the real solve_acquisition_ilp in 07.
# ══════════════════════════════════════════════════════════════════════════════

source("07__ilp_core.R")

# ── Tunable: how long to let each read run. The bound stalls within seconds, so
#    120s is ample. Lower it for an even faster read. ────────────────────────────
TIME_LIMIT <- 120

# ── Instrumented solver: faithful superset of the 07 version ────────────────────
#    The model assembly is byte-for-byte the same as 07. Only the lines marked
#    "# <- added / changed" differ, and none of them alter the optimum.
solve_acquisition_ilp <- function(V_mat, cost_mat, budget, avail, periods,
                                   relax = FALSE, extra_params = list()) {   # <- added arg
  if (length(avail) == 0 || length(periods) == 0)
    return(list(picks = data.frame(eau_idx = integer(0), period = integer(0)),
                objval = 0, x = numeric(0), grid = NULL,
                objbound = NA_real_, mipgap = NA_real_,                      # <- added
                runtime = 0, status = "EMPTY"))                             # <- added

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
    i    = c(row_once, row_budg),
    j    = c(seq_len(nvar), seq_len(nvar)),
    x    = c(rep(1, nvar), costs_s),
    dims = c(n_av + n_pe, nvar)
  )
  sense <- rep("<=", n_av + n_pe)
  rhs   <- c(rep(1, n_av), budget[periods] / cost_scale)

  model <- list(
    A = A, obj = obj, sense = sense, rhs = rhs,
    lb = rep(0, nvar), ub = rep(1, nvar),
    vtype = if (relax) "C" else "B", modelsense = "max"
  )

  params <- modifyList(list(OutputFlag = 1), extra_params)                  # <- added
  res    <- gurobi::gurobi(model, params = params)                          # <- changed

  gx <- function(v, d = NA_real_) if (is.null(v)) d else v                  # <- added
  x <- res$x
  chosen <- if (relax) integer(0) else which(x > 0.5)
  picks  <- grid[chosen, , drop = FALSE]
  list(picks = picks, objval = res$objval, x = x, grid = grid,
       objbound = gx(res$objbound), mipgap = gx(res$mipgap),                # <- added
       runtime  = gx(res$runtime, 0), status = gx(res$status, "UNKNOWN"))   # <- added
}

# ── Build the exact expensive instance: stationary, full horizon, all parcels ──
#    This reproduces the rolling policy's FIRST solve (t = 0, all 9 periods),
#    which is the solve that stalls.
data_panel <- readRDS("input_data/data_panel.rds")
sc  <- build_scenario_matrices(data_panel, SCENARIOS[["stationary"]])
bud <- make_budget(sc$cost)
V   <- t(vapply(seq_len(nrow(sc$b)),
                function(i) compute_value_vector(sc$b[i, ], sc$lam[i, ]),
                numeric(ncol(sc$b))))

# ── Read the certified gap under one config ────────────────────────────────────
read_gap <- function(label, extra_params) {
  cat("\n────────────────────────────────────────────────────────\n")
  cat("  ", label, "\n", sep = "")
  cat("────────────────────────────────────────────────────────\n")
  sol <- solve_acquisition_ilp(
    V, sc$cost, bud,
    avail = seq_len(nrow(V)), periods = seq_len(ncol(V)),
    relax = FALSE, extra_params = extra_params
  )
  inc <- sol$objval     # incumbent  = best value_added found so far
  bnd <- sol$objbound   # best bound = proven UPPER bound on value_added (max sense)
  data.frame(
    config      = label,
    incumbent   = inc,
    bound       = bnd,
    abs_gap     = bnd - inc,            # value_added units (bound is an upper bound)
    rel_gap_pct = 100 * sol$mipgap,     # Gurobi's MIPGap, as a percent
    runtime_s   = round(sol$runtime, 1),
    status      = sol$status,
    stringsAsFactors = FALSE
  )
}

r1 <- read_gap("Config 1: default Gurobi params",
               list(OutputFlag = 1, TimeLimit = TIME_LIMIT))
r2 <- read_gap("Config 2: MIPFocus = 3, Cuts = 2   (Option A1)",
               list(OutputFlag = 1, TimeLimit = TIME_LIMIT, MIPFocus = 3, Cuts = 2))

summary_tbl <- rbind(r1, r2)

# ── Report ─────────────────────────────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════════════════\n")
cat("  CERTIFIED-GAP SUMMARY   (stationary, full-horizon solve)\n")
cat("══════════════════════════════════════════════════════════\n")
print(summary_tbl, row.names = FALSE)

cat("\nHow to read this:\n")
cat("  incumbent  : the value_added the answer you ALREADY have achieves.\n")
cat("  bound      : best provable ceiling; the true optimum lies in [inc, bound].\n")
cat("  abs_gap    : bound - incumbent = the MOST extra value_added the true\n")
cat("               optimum could have. This IS your reporting uncertainty band.\n")
cat("  rel_gap_pct: abs_gap as a percent of incumbent (Gurobi's MIPGap).\n")

cat("\nDecision guide (Option A2 — can we just stop and report?):\n")
cat("  - rel_gap_pct already <= 0.01%  -> it effectively converged: stop, done.\n")
cat("  - else judge abs_gap against the rolling-vs-myopic value_added differences\n")
cat("    you intend to report. abs_gap much smaller than those -> A2 suffices:\n")
cat("    stop and report 'value_added = incumbent (+0 / +abs_gap)'.\n")
cat("  - abs_gap comparable to the effect size -> A2 not enough; escalate to\n")
cat("    pruning, then the tie-break.\n")

cat("\nDid Option A1 help?  Compare 'bound' across the two configs:\n")
cat("  - bound drops toward incumbent (rel_gap_pct shrinks) -> MIPFocus/Cuts works.\n")
cat("  - bound essentially unchanged                        -> A1 is not the lever.\n")
cat("  (If a config solved within the limit, its status reads OPTIMAL and\n")
cat("   runtime_s < TIME_LIMIT — that alone tells you it closed.)\n")

# ── Save for the record ────────────────────────────────────────────────────────
out_dir <- "output_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(summary_tbl, file.path(out_dir, "stationary_gap_diagnostic.csv"),
          row.names = FALSE)
cat("\nSaved: ", file.path(out_dir, "stationary_gap_diagnostic.csv"), "\n", sep = "")
