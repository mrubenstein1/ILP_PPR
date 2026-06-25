# ══════════════════════════════════════════════════════════════════════════
# test_tie_break.R — standalone test of the symmetry-breaking eps fix
# ══════════════════════════════════════════════════════════════════════════
# Run this from your project root (where input_data/ and 07__ilp_core.R live).
# It does NOT modify 07__ilp_core.R. It locally re-defines compute_value_vector
# with the eps argument, builds the exact stalling instance (stationary, full
# horizon, all 841 parcels), and compares solve behavior with eps = 0 vs the
# chosen eps, both capped at a short TimeLimit so you get a fast read either way.
# ══════════════════════════════════════════════════════════════════════════

source("07__ilp_core.R")

# ---- locally override with the eps-aware version (07 itself is untouched) ----
compute_value_vector <- function(b_vec, lam_vec, delta = DELTA, eps = 0) {
  n_t <- length(b_vec)
  S <- numeric(n_t); S[1] <- 1
  if (n_t >= 2) for (t in 2:n_t) S[t] <- S[t - 1] * (1 - lam_vec[t - 1])
  disc    <- delta^(0:(n_t - 1))
  contrib <- disc * b_vec * (1 - S)
  V <- rev(cumsum(rev(contrib)))
  V - eps * (seq_len(n_t) - 1)
}

choose_tie_break_eps <- function(V_mat, n_t, safety = 1e-6) {
  vals <- abs(V_mat[V_mat != 0])
  if (length(vals) == 0) return(0)
  (safety * median(vals)) / n_t
}

TIME_LIMIT <- 120  # short capped read, same logic as 10__diagnose_stationary_gap.R

# ---- build the exact stalling instance ----
data_panel <- readRDS("input_data/data_panel.rds")
sc  <- build_scenario_matrices(data_panel, SCENARIOS[["stationary"]])
bud <- make_budget(sc$cost)
n_t <- ncol(sc$b)

V0 <- t(vapply(seq_len(nrow(sc$b)),
               function(i) compute_value_vector(sc$b[i, ], sc$lam[i, ], eps = 0),
               numeric(n_t)))
eps <- choose_tie_break_eps(V0, n_t)
cat("Chosen eps:", eps,
    "| max cumulative perturbation (period 9):", eps * (n_t - 1),
    "| median |V|:", median(abs(V0[V0 != 0])), "\n")

Veps <- t(vapply(seq_len(nrow(sc$b)),
                 function(i) compute_value_vector(sc$b[i, ], sc$lam[i, ], eps = eps),
                 numeric(n_t)))

# ---- instrumented solve (same as 10__diagnose_stationary_gap.R) ----
solve_instrumented <- function(V_mat, label) {
  cat("\n──", label, "──\n")
  grid  <- expand.grid(eau_idx = seq_len(nrow(V_mat)), period = seq_len(n_t))
  nvar  <- nrow(grid)
  obj   <- V_mat[cbind(grid$eau_idx, grid$period)]
  costs <- sc$cost[cbind(grid$eau_idx, grid$period)]
  n_av  <- nrow(V_mat)
  row_once <- match(grid$eau_idx, seq_len(n_av))
  row_budg <- n_av + match(grid$period, seq_len(n_t))
  cost_scale <- 1e6
  costs_s <- costs / cost_scale

  A <- Matrix::sparseMatrix(
    i = c(row_once, row_budg), j = c(seq_len(nvar), seq_len(nvar)),
    x = c(rep(1, nvar), costs_s), dims = c(n_av + n_t, nvar)
  )
  model <- list(
    A = A, obj = obj, sense = rep("<=", n_av + n_t),
    rhs = c(rep(1, n_av), bud / cost_scale),
    lb = rep(0, nvar), ub = rep(1, nvar), vtype = "B", modelsense = "max"
  )
  res <- gurobi::gurobi(model, params = list(OutputFlag = 1, TimeLimit = TIME_LIMIT))
  data.frame(label = label, incumbent = res$objval, bound = res$objbound,
             rel_gap_pct = 100 * res$mipgap, runtime_s = round(res$runtime, 1),
             status = res$status)
}

r_baseline <- solve_instrumented(V0,   "eps = 0 (current behavior)")
r_tiebreak <- solve_instrumented(Veps, paste0("eps = ", signif(eps, 4), " (tie-break)"))

cat("\n══ COMPARISON ══\n")
print(rbind(r_baseline, r_tiebreak), row.names = FALSE)
cat("\nIf rel_gap_pct drops sharply (e.g. to ~0%) with the tie-break and runtime_s\n")
cat("< TimeLimit (status OPTIMAL), the symmetry was indeed the bottleneck and this\n")
cat("fix is safe to fold into 07__ilp_core.R.\n")
