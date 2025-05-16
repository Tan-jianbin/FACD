#' Canonical Cross-Covariance Analysis (CCVA)
#'
#' Performs canonical cross-covariance analysis between two high-dimensional longitudinal datasets.
#'
#' @param grid_num Integer. Number of grid points for time discretization (default: 201).
#' @param dat_X List. Contains components `n`, `p`, and `dat` for dataset X.
#' @param dat_Y List. Contains components `n`, `p`, and `dat` for dataset Y.
#' @param mean_zero Logical. Center data to zero mean before analysis (default: FALSE).
#' @param print_tuning Logical. Print tuning information during sparse fitting (default: TRUE).
#' @param length_tun Integer. Number of candidate tuning parameters for sparse CCA (default: 30).
#' @param max_var Numeric vector of length 2. Maximum variance parameters for (X, Y) groups (default: c(min(n,p_X), min(n,p_Y))).
#' @param Sparse Logical. If TRUE, perform sparse CCA with variable selection (default: TRUE).
#' @param comp Integer. Number of canonical components to compute.
#' @param FVE Numeric. Fraction of variance explained threshold for basis estimation (default: 0.95).
#' @return A list containing:
#' \item{Mean}{Estimated means for each feature.}
#' \item{Time_mark}{Mapping from subject times to grid indices.}
#' \item{Optimal_basis}{Learned basis functions and component numbers for X and Y.}
#' \item{var_X, var_Y}{Selected variables for each dataset.}
#' \item{load_X, load_Y}{Longitudinal/functional loading for each dataset.}
#' \item{Z_X, Z_Y}{Subject-specific canonical variables.}
#' @export
CCVA <- function(grid_num    = 201,
                 dat_X, dat_Y,
                 mean_zero    = FALSE,
                 print_tuning = TRUE,
                 length_tun   = 30,
                 max_var,
                 Sparse       = TRUE,
                 comp,
                 FVE          = 0.95) {

  ## 1. SETUP ---------------------------------------------------------------
  n   <- dat_X$n
  p_X <- dat_X$p
  p_Y <- dat_Y$p
  time_grid <- seq(0, 1, length.out = grid_num)

  ## 2. BASIS & GROUP INDICES ----------------------------------------------
  Optimal_basis <- Optimal_basis_est(dat_X, dat_Y,
                                     grid_num,
                                     max_comp_num = 5,
                                     mean_zero, FVE)
  group_x <- rep(seq_len(p_X),
                 each = Optimal_basis$Optimal_bais$X$comp_num)
  group_y <- rep(seq_len(p_Y),
                 each = Optimal_basis$Optimal_bais$Y$comp_num)

  ## 3. PROJECT DATA ONTO BASIS --------------------------------------------
  dat_x_g <- t(sapply(seq_len(n), function(i) {
    unlist(lapply(seq_len(p_X), function(j) {
      phi_mat <- Optimal_basis$Optimal_bais$X$phi[
        Optimal_basis$Time_mark$X[[i]],
        1:Optimal_basis$Optimal_bais$X$comp_num
      ]
      vals <- dat_X$dat[[j]][[i]]$val -
        Optimal_basis$Mean$X[[j]]$mu[
          Optimal_basis$Time_mark$X[[i]]
        ]
      crossprod(phi_mat, vals) / length(dat_X$dat[[j]][[i]]$time)
    }))
  }))
  dat_y_g <- t(sapply(seq_len(n), function(i) {
    unlist(lapply(seq_len(p_Y), function(j) {
      phi_mat <- Optimal_basis$Optimal_bais$Y$phi[
        Optimal_basis$Time_mark$Y[[i]],
        1:Optimal_basis$Optimal_bais$Y$comp_num
      ]
      vals <- dat_Y$dat[[j]][[i]]$val -
        Optimal_basis$Mean$Y[[j]]$mu[
          Optimal_basis$Time_mark$Y[[i]]
        ]
      crossprod(phi_mat, vals) / length(dat_Y$dat[[j]][[i]]$time)
    }))
  }))

  ## 4. CROSS‐COVARIANCE & INITIAL SVD --------------------------------------
  Optimal_basis$G_matrix <- sapply(seq_len(n), function(i) {
    tcrossprod(dat_x_g[i,], dat_y_g[i,])
  }, simplify = "array")
  Cov <- apply(Optimal_basis$G_matrix, c(1,2), mean)
  svd_fit <- svd(Cov)
  fit_cca <- list(
    X_coef = svd_fit$u,
    Y_coef = svd_fit$v,
    cov    = svd_fit$d
  )
  tol_cov <- sqrt(sum(svd_fit$d^2))
  rm(svd_fit); gc()

  ## 5. SPARSE CCA (OPTIONAL) -----------------------------------------------
  if (Sparse) {
    cv_num  <- 5
    cv_mark <- lapply(seq_len(cv_num),
                      function(k) floor(seq(k, n, cv_num)))

    Initial_cca <- list(
      max_sin = fit_cca$cov[1],
      X_coef  = fit_cca$X_coef[,1],
      Y_coef  = fit_cca$Y_coef[,1],
      Cov     = Cov,
      CV      = lapply(seq_len(cv_num), function(k) {
        train_idx <- setdiff(seq_len(n), cv_mark[[k]])
        test_idx  <- cv_mark[[k]]
        list(
          train = rowMeans(Optimal_basis$G_matrix[,,train_idx], dims=2),
          test  = rowMeans(Optimal_basis$G_matrix[,,test_idx ], dims=2)
        )
      }),
      cv_num = cv_num
    )
    Optimal_basis$G_matrix <- NULL; gc()

    fit_sp <- SLCCA_tr(group_x, group_y,
                       print_tuning, length_tun, max_var,
                       Initial_cca, n, p_X, p_Y)
    gc()

    var_X <- list(fit_sp$var_X)
    var_Y <- list(fit_sp$var_Y)
    fit_cca$X_coef[,1] <- fit_sp$fit_opt$fit_cca$X_coef
    fit_cca$Y_coef[,1] <- fit_sp$fit_opt$fit_cca$Y_coef
    fit_cca$cov[1]     <- crossprod(fit_cca$X_coef[,1], Cov %*% fit_cca$Y_coef[,1])

    ## Deflation for multiple components
    if (comp > 1) {
      for (r in 2:comp) {
        fit_mat <- fit_cca$cov[r-1] *
          tcrossprod(fit_cca$X_coef[,r-1],
                     fit_cca$Y_coef[,r-1])
        Cov <- Cov - fit_mat
        svd_fit <- svd(Cov)
        Initial_cca <- list(
          max_sin = svd_fit$d[1],
          X_coef  = svd_fit$u[,1],
          Y_coef  = svd_fit$v[,1],
          Cov     = Cov,
          CV      = lapply(Initial_cca$CV, function(cv) {
            list(
              train = cv$train - fit_mat,
              test  = cv$test  - fit_mat
            )
          }),
          cv_num = cv_num
        ); gc()

        fit_sp <- SLCCA_tr(group_x, group_y,
                           print_tuning, length_tun, max_var,
                           Initial_cca, n, p_X, p_Y); gc()

        var_X[[r]] <- fit_sp$var_X
        var_Y[[r]] <- fit_sp$var_Y
        fit_cca$X_coef[,r] <- fit_sp$fit_opt$fit_cca$X_coef
        fit_cca$Y_coef[,r] <- fit_sp$fit_opt$fit_cca$Y_coef
        fit_cca$cov[r]     <- crossprod(fit_cca$X_coef[,r],
                                        Cov %*% fit_cca$Y_coef[,r])
      }
    }
  } else {
    var_X <- replicate(comp, seq_len(p_X), simplify = FALSE)
    var_Y <- replicate(comp, seq_len(p_Y), simplify = FALSE)
  }
  rm(Cov)

  ## 6. COMPUTE LOADING SURFACES & SCORES -------------------------------
  load_X <- lapply(seq_len(comp), function(k) {
    matrix(fit_cca$X_coef[,k], p_X, byrow=TRUE) %*%
      t(Optimal_basis$Optimal_bais$X$phi[,1:Optimal_basis$Optimal_bais$X$comp_num])
  })
  load_Y <- lapply(seq_len(comp), function(k) {
    matrix(fit_cca$Y_coef[,k], p_Y, byrow=TRUE) %*%
      t(Optimal_basis$Optimal_bais$Y$phi[,1:Optimal_basis$Optimal_bais$Y$comp_num])
  })

  Z_X <- lapply(seq_len(comp), function(k) {
    sapply(seq_len(n), function(i) {
      sum(sapply(seq_len(p_X), function(j) {
        mean(
          (dat_X$dat[[j]][[i]]$val -
             Optimal_basis$Mean$X[[j]]$mu[
               Optimal_basis$Time_mark$X[[i]]
             ]) *
            load_X[[k]][j, Optimal_basis$Time_mark$X[[i]]]
        )
      }))
    })
  })
  Z_Y <- lapply(seq_len(comp), function(k) {
    sapply(seq_len(n), function(i) {
      sum(sapply(seq_len(p_Y), function(j) {
        mean(
          (dat_Y$dat[[j]][[i]]$val -
             Optimal_basis$Mean$Y[[j]]$mu[
               Optimal_basis$Time_mark$Y[[i]]
             ]) *
            load_Y[[k]][j, Optimal_basis$Time_mark$Y[[i]]]
        )
      }))
    })
  })

  ## 7. RETURN RESULTS ----------------------------------------------------
  list(
    Mean          = Optimal_basis$Mean,
    Time_mark     = Optimal_basis$Time_mark,
    Optimal_basis = Optimal_basis$Optimal_bais,
    var_X         = var_X,
    var_Y         = var_Y,
    load_X        = load_X,
    load_Y        = load_Y,
    Z_X           = Z_X,
    Z_Y           = Z_Y
  )
}

#' Estimate Data-Adaptive Basis Functions
#'
#' Computes mean functions, time mappings, and optimal basis functions for two longitudinal datasets.
#'
#' @param dat_X List with elements:
#'   \describe{
#'     \item{n}{Number of subjects}
#'     \item{p}{Number of features}
#'     \item{dat}{Nested list of length \code{p}, each entry a list of length \code{n} with \code{$time} and \code{$val}}}
#' @param dat_Y List with the same structure as \code{dat_X} for dataset Y.
#' @param grid_num Integer. Number of grid points for discretizing the time domain.
#' @param max_comp_num Integer. Maximum number of basis components to compute (default: 5).
#' @param mean_zero Logical. If \code{TRUE}, forces mean curves to zero (default: \code{FALSE}).
#' @param FVE Numeric. Fraction of variance explained threshold for selecting components.
#'
#' @return A list with:
#'   \describe{
#'     \item{Mean}{List with \code{X} and \code{Y}, each a list of mean curves and their workGrid.}
#'     \item{Time_mark}{List with \code{X} and \code{Y}, each mapping observed times to grid indices.}
#'     \item{Optimal_basis}{List with \code{X} and \code{Y}, each containing \code{phi}, \code{workGrid}, \code{eigenvalue}, and \code{comp_num}.}}
#'
#' @importFrom stats spline
#' @importFrom fdapace GetMeanCurve
#' @export
Optimal_basis_est <- function(dat_X,
                              dat_Y,
                              grid_num,
                              max_comp_num = 5,
                              mean_zero     = FALSE,
                              FVE) {
  # Setup
  n        <- dat_X$n
  p_X      <- dat_X$p
  p_Y      <- dat_Y$p
  fda_X    <- dat_X$dat
  fda_Y    <- dat_Y$dat
  work_grid <- seq(0, 1, length.out = grid_num)

  # Helper: estimate mean curve for a single feature
  estimate_mean <- function(feature_list) {
    times_all <- unique(unlist(lapply(feature_list, `[[`, "time")))
    if (length(times_all) <= 5) {
      vals_mat <- sapply(feature_list, function(obs) {
        sapply(times_all, function(t) {
          idx <- which(obs$time == t)
          if (length(idx)) obs$val[idx] else NA
        })
      })
      mu_vals <- rowMeans(vals_mat, na.rm = TRUE)
      mu      <- spline(x = times_all, y = mu_vals, xout = work_grid)$y
    } else {
      Ly <- lapply(feature_list, `[[`, "val")
      Lt <- lapply(feature_list, `[[`, "time")
      mu <- GetMeanCurve(
        Ly    = do.call(cbind, Ly),
        Lt    = Lt,
        optns = list(
          error          = TRUE,
          nRegGrid       = grid_num,
          methodBwMu     = "GCV",
          methodMuCovEst = "smooth"
        )
      )$mu
    }
    list(mu = mu, workGrid = work_grid)
  }

  # Compute mean functions
  if (mean_zero) {
    fit_mean_X <- replicate(p_X,
                            list(mu = rep(0, grid_num), workGrid = work_grid),
                            simplify = FALSE)
    fit_mean_Y <- replicate(p_Y,
                            list(mu = rep(0, grid_num), workGrid = work_grid),
                            simplify = FALSE)
  } else {
    fit_mean_X <- lapply(fda_X, estimate_mean)
    fit_mean_Y <- lapply(fda_Y, estimate_mean)
  }

  # Map observed times to grid indices
  map_time <- function(fda_obj, fit_mean) {
    lapply(seq_len(n), function(i) {
      sapply(fda_obj[[1]][[i]]$time, function(t) which.min(abs(t - fit_mean[[1]]$workGrid)))
    })
  }
  time_mark_X <- map_time(fda_X, fit_mean_X)
  time_mark_Y <- map_time(fda_Y, fit_mean_Y)

  # Build subject‐wise covariance arrays
  build_G <- function(fda_obj, fit_mean, time_mark) {
    lapply(seq_len(n), function(i) {
      Tlen <- length(time_mark[[i]])
      Gsub <- array(0, dim = c(Tlen, Tlen, p_X, p_Y))
      for (j1 in seq_len(p_X)) {
        for (j2 in seq_len(p_Y)) {
          xv <- fda_obj[[j1]][[i]]$val - fit_mean[[j1]]$mu[time_mark[[i]]]
          yv <- fda_obj[[j2]][[i]]$val - fit_mean[[j2]]$mu[time_mark[[i]]]
          Gsub[,,j1,j2] <- tcrossprod(xv, yv)
        }
      }
      apply(Gsub, c(1,2), mean)
    })
  }
  G_X_list <- build_G(fda_X, fit_mean_X, time_mark_X)
  G_Y_list <- build_G(fda_Y, fit_mean_Y, time_mark_Y)

  # Estimate basis via FSVD_r()
  basis_X <- FSVD_r(M        = G_X_list,
                    time     = work_grid,
                    time_mark= time_mark_X,
                    R_max    = max_comp_num,
                    FVE      = FVE)
  basis_Y <- FSVD_r(M        = G_Y_list,
                    time     = work_grid,
                    time_mark= time_mark_Y,
                    R_max    = max_comp_num,
                    FVE      = FVE)

  # Return results
  list(
    Mean         = list(X = fit_mean_X, Y = fit_mean_Y),
    Time_mark    = list(X = time_mark_X, Y = time_mark_Y),
    Optimal_basis= list(X = basis_X, Y = basis_Y)
  )
}

#' Sparse Canonical Correlation Analysis
#'
#' Performs grid search and cross-validation to tune sparse CCA penalties.
#'
#' @param group_x Integer vector of group labels for X coefficients.
#' @param group_y Integer vector of group labels for Y coefficients.
#' @param print_tuning Logical; if TRUE, prints CV error matrix (default: FALSE).
#' @param length_tun Integer; number of tuning values per margin (default: 30).
#' @param max_var Numeric(2); maximum nonzero groups allowed for (X, Y).
#' @param Initial_cca List; initial CCA state containing Cov, Y_coef, X_coef, and CV splits.
#' @param n Integer; number of subjects.
#' @param p_X Integer; number of X groups.
#' @param p_Y Integer; number of Y groups.
#' @return A list with elements:
#'   \item{fit_opt}{Best CV-chosen CCA fit (from cv_slcca).}
#'   \item{error}{CV error of the chosen fit.}
#'   \item{CV_error}{Matrix of CV errors over the tuning grid.}
#'   \item{var_X}{Indices of selected X groups.}
#'   \item{var_Y}{Indices of selected Y groups.}
#' @export
#'
#' @note Ensure the following helper functions are implemented or imported:
#'   * \code{cca_AX(lambda, abs, Initial_cca, n, group_x, group_y)}
#'   * \code{cv_slcca(lambda, Initial_cca, n, group_x, group_y, max_var)}
#' @importFrom stats sd
#' @export
SLCCA_tr <- function(group_x, group_y, print_tuning = F, length_tun = 30, max_var, Initial_cca, n, p_X, p_Y){

  # Use grid search to determine the initial tuning parameters
  vec <- c(Initial_cca$Cov %*% Initial_cca$Y_coef)
  max_t <- max(sapply(unique(group_x), function(j) sqrt(sum(vec[which(group_x == j)] ^ 2))))

  fit <- cca_AX(lambda = c(max_t, max_t), abs = 0.001, Initial_cca, n, group_x, group_y)
  var_num_X <- length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0))
  var_num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0))

  while(var_num_X + var_num_Y != 0){
    max_t <- max_t * 5
    fit <- cca_AX(lambda = c(max_t, max_t), abs = 0.001, Initial_cca, n, group_x, group_y)
    var_num_X <- length(which(sapply(unique(group_x), function(j) sum(fit$fit_cca$X_coef[group_x == j] ^ 2)) != 0))
    var_num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit$fit_cca$Y_coef[group_y == j] ^ 2)) != 0))
  }

  lambda_list <- c(max_t / 2, max_t / 10, max_t / 100, max_t / 1000, max_t / 10000)
  var_num_X <- matrix(NA, length(lambda_list), length(lambda_list))
  var_num_Y <- matrix(NA, length(lambda_list), length(lambda_list))

  for(a in 1:length(lambda_list)){
    for(b in 1:length(lambda_list)){
      fit <- cca_AX(lambda = c(lambda_list[a], lambda_list[b]), abs = 0.001, Initial_cca, n, group_x, group_y)
      var_num_X[b,a] <- length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0))
      var_num_Y[b,a] <- length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0))
    }
  }

  A <- var_num_X >= max_var[1]
  if(length(which(colMeans(A) > 0)) == 0){
    min_X <- max_t / 100000
    fit <- cca_AX(lambda = c(min_X, max_t), abs = 0.001, Initial_cca, n, group_x, group_y)
    while(length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0)) <= max_var[1]){
      min_X <- min_X / 10
      fit <- cca_AX(lambda = c(min_X, max_t), abs = 0.001, Initial_cca, n, group_x, group_y)
    }
  }else{
    min_X <- lambda_list[min(which(colMeans(A) > 0))]
  }

  A <- var_num_Y >= max_var[2]
  if(length(which(rowMeans(A) > 0)) == 0){
    min_Y <- max_t / 100000
    fit <- cca_AX(lambda = c(min_X, min_Y), abs = 0.001, Initial_cca, n, group_x, group_y)
    while(length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0)) <= max_var[2]){
      min_Y <- min_Y / 10
      fit <- cca_AX(lambda = c(min_X, min_Y), abs = 0.001, Initial_cca, n, group_x, group_y)
    }
  }else{
    min_Y <- lambda_list[min(which(rowMeans(A) > 0))]
  }

  lambda_list_X <- exp(seq(log(max_t), log(min_X), length.out = length_tun + 1))[-1]
  lambda_list_Y <- matrix(NA, length_tun, length_tun)
  CV_error <- matrix(NA, length(lambda_list_X), length(lambda_list_X))

  for(a in 1:length(lambda_list_X)){
    a_X <- c(Initial_cca$Cov %*% Initial_cca$Y_coef)
    a_X <- c(sapply(1:p_X, function(j){
      norm <- sqrt(sum(a_X[group_x == j] ^ 2))
      if(norm == 0){
        rep(0, sum(group_x == j))
      }else{
        max((1 - lambda_list_X[a] / norm), 0) * a_X[group_x == j]
      }
    }))
    x_norm <- c(sqrt(crossprod(a_X)))
    if(x_norm != 0){
      a_X <- a_X / x_norm
    }

    a_Y <- c(t(a_X) %*% Initial_cca$Cov)
    max_Y <- max(c(sapply(1:p_Y, function(j){
      sqrt(sum(a_Y[group_y == j] ^ 2))
    })))
    if(max_Y == 0){
      max_Y <- max_t
    }

    lambda_list_Y[,a] <- exp(seq(log(max_Y), log(min_Y), length.out = length_tun + 1))[-1]
    for(b in 1:length(lambda_list_Y[,a])){
      fit <- cv_slcca(lambda = c(lambda_list_X[a], lambda_list_Y[b,a]), Initial_cca, n, group_x, group_y, max_var)
      CV_error[b,a] <- fit$error
      if(print_tuning == TRUE) print(CV_error)
    }
  }

  CV_error[CV_error == -10^10] <- NA
  CV_error[CV_error == -10^11] <- NA
  sd_CV <- 0.1 * sd(CV_error[is.finite(CV_error)])

  mark <- (CV_error >= max(CV_error, na.rm = TRUE) - sd_CV)
  row_mark <- min(which(apply(mark, 2, function(x) mean(x, na.rm = TRUE)) != 0))
  lambda_1 <- c(lambda_list_X[row_mark], lambda_list_Y[min(which(mark[,row_mark] == TRUE)), row_mark])
  col_mark <- min(which(apply(mark, 1, function(x) mean(x, na.rm = TRUE)) != 0))
  lambda_2 <- c(lambda_list_X[min(which(mark[col_mark,] == TRUE))], lambda_list_Y[col_mark, min(which(mark[col_mark,] == TRUE))])

  fit_opt_1 <- cv_slcca(lambda = lambda_1, Initial_cca, n, group_x, group_y, max_var)
  fit_opt_2 <- cv_slcca(lambda = lambda_2, Initial_cca, n, group_x, group_y, max_var)
  fit_opt <- if(sum(sapply(unique(group_x), function(j) sum(fit_opt_1$fit_cca$X_coef[group_x == j] ^ 2)) != 0) +
                sum(sapply(unique(group_y), function(j) sum(fit_opt_1$fit_cca$Y_coef[group_y == j] ^ 2)) != 0) <=
                sum(sapply(unique(group_x), function(j) sum(fit_opt_2$fit_cca$X_coef[group_x == j] ^ 2)) != 0) +
                sum(sapply(unique(group_y), function(j) sum(fit_opt_2$fit_cca$Y_coef[group_y == j] ^ 2)) != 0)) {
    fit_opt_1
  } else {
    fit_opt_2
  }

  list(
    fit_opt  = fit_opt,
    error    = fit_opt$error,
    CV_error = CV_error,
    var_X    = which(sapply(unique(group_x), function(j) sum(fit_opt$fit_cca$X_coef[group_x == j]^2)) != 0),
    var_Y    = which(sapply(unique(group_y), function(j) sum(fit_opt$fit_cca$Y_coef[group_y == j]^2)) != 0)
  )
}

#' SVD Basis Extraction
#'
#' Computes SVD basis functions for cross-covariance arrays.
#'
#' @param M List of subject-level covariance matrices (time×time).
#' @param time Numeric vector of grid time points.
#' @param time_mark List mapping each subject’s observed times to grid indices.
#' @param R_max Integer. Maximum number of basis components.
#' @param FVE Numeric. Fraction of variance explained threshold.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{phi}{Matrix of basis functions (columns).}
#'     \item{workGrid}{Time grid vector.}
#'     \item{eigenvalue}{Eigenvalues corresponding to basis functions.}
#'     \item{comp_num}{Selected number of components.}}
#'
#' @note Requires implementation of `fit_basis(M, time_mark, time)` that returns:
#'   \describe{
#'     \item{tol_nor}{Total norm for FVE calculation}
#'     \item{eigendocomposition}{List with `values` and `vectors`}}
#' @export
FSVD_r <- function(M, time, time_mark, R_max, FVE){
  fit <- fit_basis(M, time_mark, time)
  tol_norm <- fit$tol_nor
  op_basis <- list()
  op_basis$phi <- fit$eigendocomposition$vectors[,1:R_max] * sqrt(length(time))
  op_basis$workGrid <- time
  op_basis$eigenvalue <- fit$eigendocomposition$values[1:R_max] / length(time)
  op_basis$comp_num <- max(sum(cumsum(op_basis$eigenvalue) / tol_norm <= FVE) + 1, 2)
  return(op_basis)
}

#' Compute Covariance Basis via Eigen Decomposition
#'
#' Builds the covariance surface and computes its eigendecomposition,
#' returning key components for functional basis extraction.
#'
#' @param M List of subject-level covariance matrices (time×time arrays).
#' @param time_mark List mapping each subject’s observed times to grid indices.
#' @param time Numeric vector of grid time points.
#'
#' @return A list containing:
#'   \describe{
#'     \item{eigenvalue}{Leading eigenvalue of the covariance surface, normalized by grid size.}
#'     \item{eigenfunction}{Leading eigenfunction scaled by \code{sqrt(nrow(cov))}.}
#'     \item{eigendocomposition}{Full result of \code{eigen()}, with \code{$values} and \code{$vectors}.}
#'     \item{tol_nor}{Total sum of eigenvalues divided by grid size (for FVE calculation).}}
#'
#' @note Requires implementation of \code{fit_cov_surface()} to assemble the covariance matrix.
#' @examples
#' # fit_cov <- fit_basis(M_list, time_mark_list, time_grid)
#'
#' @export
fit_basis <- function(M, time_mark, time) {
  # Build the full covariance surface from subject-level matrices
  cov_mat <- fit_cov_surface(M, time_mark, time_mark, time)
  # Perform symmetric eigen decomposition
  fit <- eigen(cov_mat, symmetric = TRUE)
  # Leading eigenfunction scaled by grid length
  primary_phi <- fit$vectors[, 1] * sqrt(nrow(cov_mat))

  # Compute normalized leading eigenvalue
  leading_ev <- as.numeric(t(primary_phi) %*% cov_mat %*% primary_phi / (nrow(cov_mat)^2))

  # Return components for basis extraction
  list(
    eigenvalue        = leading_ev,
    eigenfunction     = primary_phi,
    eigendocomposition= fit,
    tol_nor           = sum(fit$values) / nrow(cov_mat)
  )
}


#' Build and Smooth Covariance Surface
#'
#' Assembles the covariance surface matrix from subject-level cross-covariance arrays
#' and fits a smooth surface using either thin-plate splines or a GAM.
#'
#' @param M List of length \code{n}, each element a time×time covariance matrix for a subject.
#' @param time_mark_X List of length \code{n}, mapping X observation indices to grid indices.
#' @param time_mark_Y List of length \code{n}, mapping Y observation indices to grid indices.
#' @param time Numeric vector of grid time points.
#'
#' @return A symmetric covariance matrix of dimension \code{length(time)}×\code{length(time)}.
#'
#' @importFrom fields    Tps predictSurface
#' @importFrom mgcv     gam te
#' @importFrom stats    gaussian predict
#' @export
fit_cov_surface <- function(M, time_mark_X, time_mark_Y, time) {
  n <- length(M)
  matrix_sample <- data.frame()

  # Collect observed (x,y,z) with weights
  for (i in seq_len(n)) {
    x_idx <- time_mark_X[[i]]
    y_idx <- time_mark_Y[[i]]
    for (k in seq_along(x_idx)) {
      x_val <- time[x_idx[k]]
      y_vals <- time[y_idx]
      z_vals <- M[[i]][k, ]
      df_i <- data.frame(
        x = rep(x_val, length(y_vals)),
        y = y_vals,
        z = z_vals,
        w = 1 / (length(x_idx) * length(y_idx))
      )
      matrix_sample <- rbind(matrix_sample, df_i)
    }
  }

  uni_x <- sort(unique(matrix_sample$x))
  num_k <- min(length(uni_x), 15)

  # Create prediction grid
  grid_x <- rep(time, length(time))
  grid_y <- rep(time, each = length(time))

  # Fit thin-plate spline if grid small, else GAM
  if (length(uni_x) * length(uni_x) <= 25) {
    fit_tps <- Tps(x = cbind(matrix_sample$x, matrix_sample$y),
                   Y = matrix_sample$z,
                   df = nrow(matrix_sample) - 2,
                   method = "GCV.model")
    Cov_mat <- predictSurface(fit_tps,
                              grid.list = list(x = time, y = time),
                              extrap = TRUE)$z
  } else {
    fit_gam <- gam(z ~ te(x, y, bs = "cr", k = c(num_k, num_k)),
                   data = matrix_sample,
                   weights = matrix_sample$w,
                   family = gaussian(),
                   method = "GCV.Cp")
    pred <- predict(fit_gam, newdata = data.frame(x = grid_x, y = grid_y))
    Cov_mat <- matrix(pred, nrow = length(time))
  }

  return(Cov_mat)
}

#' Penalized CCA Core Solver
#'
#' Performs alternating updates for penalized canonical correlation analysis with group lasso penalties.
#'
#' @param lambda Numeric vector of length 2; penalties for X and Y groups.
#' @param abs Numeric; convergence threshold for relative change in coefficients.
#' @param Initial_cca List; initial CCA state containing `Cov`, `X_coef`, and `Y_coef`.
#' @param n Integer; number of subjects (unused in this function, kept for interface consistency).
#' @param group_x Integer vector; group labels for X coefficients.
#' @param group_y Integer vector; group labels for Y coefficients.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{X_coef}{Updated X canonical coefficients (numeric vector).}
#'     \item{Y_coef}{Updated Y canonical coefficients (numeric vector).}}
#'
#' @export
cca_AX <- function(lambda, abs, Initial_cca, n, group_x, group_y) {
  a_X <- Initial_cca$X_coef
  a_Y <- Initial_cca$Y_coef
  error <- 1
  ite <- 0
  num_X <- length(unique(group_x))
  num_Y <- length(unique(group_y))

  while ((error > abs) && (ite < 200)) {
    ite <- ite + 1
    a_X_t <- a_X
    a_Y_t <- a_Y

    # Update X coefficients
    a_X <- as.numeric(Initial_cca$Cov %*% a_Y)
    a_X <- unlist(lapply(seq_len(num_X), function(j) {
      idx <- which(group_x == j)
      norm_j <- sqrt(sum(a_X[idx]^2))
      if (norm_j == 0) {
        rep(0, length(idx))
      } else {
        shrink <- max(1 - lambda[1] / norm_j, 0)
        shrink * a_X[idx]
      }
    }))
    x_norm <- sqrt(sum(a_X^2))
    if (x_norm != 0) a_X <- a_X / x_norm

    # Update Y coefficients
    a_Y <- as.numeric(t(a_X) %*% Initial_cca$Cov)
    a_Y <- unlist(lapply(seq_len(num_Y), function(j) {
      idx <- which(group_y == j)
      norm_j <- sqrt(sum(a_Y[idx]^2))
      if (norm_j == 0) {
        rep(0, length(idx))
      } else {
        shrink <- max(1 - lambda[2] / norm_j, 0)
        shrink * a_Y[idx]
      }
    }))
    y_norm <- sqrt(sum(a_Y^2))
    if (y_norm != 0) a_Y <- a_Y / y_norm

    # Compute convergence error
    err_x <- if (x_norm != 0) sum(abs(a_X_t - a_X)) / sum(abs(a_X_t)) else 0
    err_y <- if (y_norm != 0) sum(abs(a_Y_t - a_Y)) / sum(abs(a_Y_t)) else 0
    error <- err_x + err_y
  }

  list(X_coef = a_X, Y_coef = a_Y)
}

#' Cross-Validation Wrapper for Sparse CCA
#'
#' Performs K-fold cross-validation for penalized canonical correlation analysis
#' using the core solver `cca_AX`, and computes mean absolute correlation on held-out folds.
#'
#' @param lambda Numeric vector of length 2; penalties for X and Y groups.
#' @param Initial_cca List; initial CCA state containing `X_coef`, `Y_coef`, and `Cov`.
#' @param n Integer; number of subjects.
#' @param zx Numeric matrix; X-projected training data (subjects × features).
#' @param zy Numeric matrix; Y-projected training data (subjects × features).
#' @param group_x Integer vector; group labels for X coefficients.
#' @param group_y Integer vector; group labels for Y coefficients.
#' @param max_var Numeric(2); maximum nonzero groups allowed for (X, Y).
#' @param cv_num Integer; number of folds for cross-validation (default: min(5, n)).
#'
#' @return A list with:
#'   \describe{
#'     \item{error}{Mean absolute correlation across validation folds.}
#'     \item{fit_cca}{Result of `cca_AX` on full data.}
#'     \item{lambda}{The penalty parameters used.}
#'   }
#'
#' @importFrom stats cor
#' @export
cv_slcca <- function(lambda, Initial_cca, n,
                     zx, zy,
                     group_x, group_y,
                     max_var,
                     cv_num = min(5, n)) {
  cv_mark <- lapply(seq_len(cv_num), function(k) floor(seq(k, n, cv_num)))

  # Fit on full data
  fit_cca <- cca_AX(lambda = lambda,
                    abs    = 0.001,
                    Initial_cca,
                    n, group_x, group_y)

  # Count selected groups
  num_X <- sum(sapply(unique(group_x), function(j) sum(fit_cca$X_coef[group_x == j]^2) != 0))
  num_Y <- sum(sapply(unique(group_y), function(j) sum(fit_cca$Y_coef[group_y == j]^2) != 0))

  # If outside allowed range, return zero
  if (num_X > max_var[1] || num_X <= 0 ||
      num_Y > max_var[2] || num_Y <= 0) {
    error <- 0
  } else {
    # CV loop
    errors <- sapply(seq_len(cv_num), function(k) {
      test_idx  <- cv_mark[[k]]
      train_idx <- setdiff(seq_len(n), test_idx)

      cv_state <- list(
        X_coef = fit_cca$X_coef,
        Y_coef = fit_cca$Y_coef,
        Cov    = crossprod(zx[train_idx, ], zy[train_idx, ]) / length(train_idx)
      )
      cv_fit <- cca_AX(lambda = lambda,
                       abs    = 0.001,
                       Initial_cca = cv_state,
                       n, group_x, group_y)

      X_val <- zx[test_idx, , drop = FALSE] %*% cv_fit$X_coef
      Y_val <- zy[test_idx, , drop = FALSE] %*% cv_fit$Y_coef
      cr <- cor(X_val, Y_val)
      if (!is.finite(cr)) 0 else abs(cr)
    })
    error <- mean(errors)
  }

  list(error   = error,
       fit_cca = fit_cca,
       lambda  = lambda)
}
