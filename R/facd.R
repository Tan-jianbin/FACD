#' Functional-Aggregated Cross-covariance Decomposition (FACD)
#'
#' Performs canonical cross-covariance analysis between two high-dimensional longitudinal datasets.
#'
#' @param grid_num Integer. Number of grid points for time discretization (default: 201).
#' @param dat_X List. Contains components `n`, `p`, and `dat` for dataset X.
#'   All observed times must be scaled to [0, 1]; see Details for preprocessing.
#' @param dat_Y List. Contains components `n`, `p`, and `dat` for dataset Y,
#'   with the same preprocessing requirements as \code{dat_X}.
#' @param mean_zero Logical. If TRUE, assumes that the mean curves are zero;
#'   otherwise estimates and subtracts them (default: FALSE). This does not
#'   standardize feature variances.
#' @param print_tuning Logical. Print tuning information during sparse fitting (default: TRUE).
#' @param length_tun Integer. Number of candidate tuning parameters for sparse CCA (default: 30).
#' @param max_var Numeric vector of length 2. Maximum variance parameters for (X, Y) groups.
#' @param Sparse Logical. If TRUE, perform sparse CCA with variable selection (default: TRUE).
#' @param comp Integer. Number of canonical components to compute.
#' @param FVE Numeric. Fraction of variance explained threshold for basis estimation (default: 0.95).
#' @param max_comp_num Integer. Maximum number of basis components (default: 5).
#' @return A list containing:
#' \item{Mean}{Estimated means for each feature.}
#' \item{Time_mark}{Mapping from subject times to grid indices.}
#' \item{Optimal_basis}{Learned basis functions and component numbers for X and Y.}
#' \item{var_X, var_Y}{Selected variables for each dataset.}
#' \item{load_X, load_Y}{Longitudinal/functional loading for each dataset.}
#' \item{Z_X, Z_Y}{Subject-specific canonical variables.}
#' \item{eta}{Estimated component cross-covariances.}
#' \item{tol_cov}{Frobenius norm of the initial cross-covariance matrix.}
#' 
#' @details
#' Before calling \code{FACD()}, scale every observed \code{$time} and the
#' corresponding \code{time_grid} to [0, 1], for example using
#' \code{(time - time_min) / (time_max - time_min)} with \code{time_max > time_min}.
#' Use one time mapping across all subjects and features within each dataset,
#' rather than rescaling each subject separately. When X and Y share a physical
#' time axis, use the same units and reference interval to preserve alignment.
#'
#' Feature-wise standardization is recommended, especially when features have
#' different units or scales. For each feature, pool observations across subjects
#' and time points, then subtract its mean and divide by its standard deviation
#' (z-score scaling). Remove features with no finite observations or zero variance
#' before standardizing. Do not standardize each subject's trajectory separately.
#' \code{FACD()} does not automatically rescale time or standardize feature variances.
#'
#' \code{\link{wrap_dat}} automatically rescales time to [0, 1] and centers each
#' feature, dividing by the root mean square (RMS) of its uncentered values
#' rather than its standard deviation. For z-score preprocessing, standardize the raw feature
#' columns before wrapping. Each call to \code{wrap_dat()} uses its own observed
#' time range; check that these ranges match when X and Y share a time axis.
#' The bundled datasets are already prepared for the example.
#'
#' @example example.R
#' @export
FACD <- function(grid_num = 201, dat_X, dat_Y, mean_zero = F, print_tuning = T, length_tun = 30, 
                 max_var, Sparse = T, comp, FVE = 0.95, max_comp_num = 5){
  
  ## First-step estimation
  time_grid <- seq(0, 1, length.out = grid_num)
  n <- dat_X$n
  p_X <- dat_X$p
  p_Y <- dat_Y$p
  
  Optimal_basis <- Optimal_basis_est(dat_X, dat_Y, grid_num, max_comp_num = max_comp_num, mean_zero, FVE)
  
  group_x <- rep(1:p_X, each = Optimal_basis$Optimal_bais$X$comp_num)
  group_y <- rep(1:p_Y, each = Optimal_basis$Optimal_bais$Y$comp_num)
  
  ## SVD on cross-covariance
  dat_x_g <- t(sapply(1:n, function(i){
    c(sapply(1:p_X, function(j){
      t(Optimal_basis$Optimal_bais$X$phi[Optimal_basis$Time_mark$X[[i]],1:Optimal_basis$Optimal_bais$X$comp_num]) %*% c(dat_X$dat[[j]][[i]]$val - Optimal_basis$Mean$X[[j]]$mu[Optimal_basis$Time_mark$X[[i]]]) / length(dat_X$dat[[j]][[i]]$time)
    }))
  }))
  
  dat_y_g <- t(sapply(1:n, function(i){
    c(sapply(1:p_Y, function(j){
      t(Optimal_basis$Optimal_bais$Y$phi[Optimal_basis$Time_mark$Y[[i]],1:Optimal_basis$Optimal_bais$Y$comp_num]) %*% c(dat_Y$dat[[j]][[i]]$val - Optimal_basis$Mean$Y[[j]]$mu[Optimal_basis$Time_mark$Y[[i]]]) / length(dat_Y$dat[[j]][[i]]$time)
    }))
  }))
  
  Optimal_basis$G_matrix <- sapply(1:n, function(i){
    tcrossprod(dat_x_g[i,], dat_y_g[i,])
  }, simplify = 'array')
  
  Cov <- apply(Optimal_basis$G_matrix, c(1, 2), mean)
  fit <- svd(Cov)
  fit_cca <- list()
  fit_cca$X_coef <- fit$u
  fit_cca$Y_coef <- fit$v
  fit_cca$cov <- fit$d
  tol_cov <- sqrt(sum(fit$d ^ 2))
  rm(fit)
  gc()
  
  if(Sparse == T){
    cv_num <- 5
    cv_mark <- lapply(1:cv_num, function(k){
      floor(seq(k, n, cv_num))
    })
    
    Initial_cca <- list(max_sin = fit_cca$cov[1], X_coef = fit_cca$X_coef[,1], Y_coef = fit_cca$Y_coef[,1],
                        Cov = Cov,
                        CV = lapply(1:cv_num, function(k){
                          list(train = rowMeans(Optimal_basis$G_matrix[,,seq_len(dim(Optimal_basis$G_matrix)[3])[-cv_mark[[k]]]], dims = 2),
                               test = rowMeans(Optimal_basis$G_matrix[,,seq_len(dim(Optimal_basis$G_matrix)[3])[cv_mark[[k]]]], dims = 2))
                        }), cv_num = cv_num
    )
    Optimal_basis$G_matrix <- NULL
    gc()
    
    fit <- SLCCA_tr(group_x = group_x, group_y = group_y, print_tuning, length_tun, max_var, Initial_cca, n, p_X, p_Y)
    gc()
    
    var_X <- list()
    var_Y <- list()
    
    var_X[[1]] <- fit$var_X
    var_Y[[1]] <- fit$var_Y
    
    fit_cca$X_coef[,1] <- fit$fit_opt$fit_cca$X_coef
    
    fit_cca$Y_coef[,1] <- fit$fit_opt$fit_cca$Y_coef
    fit_cca$cov[1] <- t(fit_cca$X_coef[,1]) %*% Cov %*% fit_cca$Y_coef[,1]
    
    if(comp > 1){
      for(r in 2:comp){
        
        u_previous <- fit_cca$X_coef[, r - 1]
        v_previous <- fit_cca$Y_coef[, r - 1]
        
        eta_previous <- drop(
          crossprod(
            u_previous,
            Cov %*% v_previous
          )
        )
        
        fit_matrix <- eta_previous *
          tcrossprod(u_previous, v_previous)
        
        Cov <- Cov - fit_matrix
        
        CV_previous <- Initial_cca$CV
        
        CV_residual <- lapply(seq_len(cv_num), function(k){
          
          eta_train <- drop(
            crossprod(
              u_previous,
              CV_previous[[k]]$train %*% v_previous
            )
          )
          
          eta_test <- drop(
            crossprod(
              u_previous,
              CV_previous[[k]]$test %*% v_previous
            )
          )
          
          list(
            train = CV_previous[[k]]$train -
              eta_train * tcrossprod(u_previous, v_previous),
            
            test = CV_previous[[k]]$test -
              eta_test * tcrossprod(u_previous, v_previous)
          )
        })
        
        svd_fit <- svd(Cov)
        
        Initial_cca <- list(
          max_sin = svd_fit$d[1],
          X_coef = svd_fit$u[, 1],
          Y_coef = svd_fit$v[, 1],
          Cov = Cov,
          CV = CV_residual,
          cv_num = cv_num
        )
        
        fit <- SLCCA_tr(
          group_x = group_x,
          group_y = group_y,
          print_tuning,
          length_tun,
          max_var,
          Initial_cca,
          n,
          p_X,
          p_Y
        )
        
        gc()
        
        var_X[[r]] <- fit$var_X
        var_Y[[r]] <- fit$var_Y
        
        fit_cca$X_coef[, r] <- fit$fit_opt$fit_cca$X_coef
        
        fit_cca$Y_coef[, r] <- fit$fit_opt$fit_cca$Y_coef
        
        fit_cca$cov[r] <- drop(
          crossprod(
            fit_cca$X_coef[, r],
            Cov %*% fit_cca$Y_coef[, r]
          )
        )
      }
    }
    rm(Initial_cca)
    
  }else{
    var_X <- lapply(1:comp, function(k) 1:p_X)
    var_Y <- lapply(1:comp, function(k) 1:p_Y)
  }
  rm(Cov)
  
  load_X <- lapply(1:comp, function(k){
    matrix(fit_cca$X_coef[,k], p_X, byrow = T) %*% t(Optimal_basis$Optimal_bais$X$phi[,1:Optimal_basis$Optimal_bais$X$comp_num])
  })
  load_Y <- lapply(1:comp, function(k){
    matrix(fit_cca$Y_coef[,k], p_Y, byrow = T) %*% t(Optimal_basis$Optimal_bais$Y$phi[,1:Optimal_basis$Optimal_bais$Y$comp_num])
  })
  
  if(comp > 1){
    for(r in 2:comp){
      for(s in 1:(r - 1)){
        load_X[[r]] <- load_X[[r]] -
          sum(load_X[[r]] * load_X[[s]]) /
          sum(load_X[[s]]^2) * load_X[[s]]
        
        load_Y[[r]] <- load_Y[[r]] -
          sum(load_Y[[r]] * load_Y[[s]]) /
          sum(load_Y[[s]]^2) * load_Y[[s]]
      }
      
      load_X[[r]] <- load_X[[r]] / sqrt(sum(rowMeans(load_X[[r]]^2)))
      load_Y[[r]] <- load_Y[[r]] / sqrt(sum(rowMeans(load_Y[[r]]^2)))
    }
  }
  
  Z_X <- vector("list", comp)
  
  for(k in 1:comp){
    Z_X[[k]] <- sapply(1:n, function(i){
      time_i <- Optimal_basis$Time_mark$X[[i]]
      score <- 0
      
      for(j in 1:p_X){
        residual <- dat_X$dat[[j]][[i]]$val - Optimal_basis$Mean$X[[j]]$mu[time_i]
        
        if(k > 1){
          for(r in 1:(k - 1)){
            residual <- residual - Z_X[[r]][i] * load_X[[r]][j, time_i]
          }
        }
        
        score <- score + mean(residual * load_X[[k]][j, time_i])
      }
      
      score
    })
  }
  
  Z_Y <- vector("list", comp)
  
  for(k in 1:comp){
    Z_Y[[k]] <- sapply(1:n, function(i){
      time_i <- Optimal_basis$Time_mark$Y[[i]]
      score <- 0
      
      for(j in 1:p_Y){
        residual <- dat_Y$dat[[j]][[i]]$val - Optimal_basis$Mean$Y[[j]]$mu[time_i]
        
        if(k > 1){
          for(r in 1:(k - 1)){
            residual <- residual - Z_Y[[r]][i] * load_Y[[r]][j, time_i]
          }
        }
        
        score <- score + mean(residual * load_Y[[k]][j, time_i])
      }
      
      score
    })
  }
  
  return(list(Mean = Optimal_basis$Mean,
              Time_mark = Optimal_basis$Time_mark,
              Optimal_basis = Optimal_basis$Optimal_bais,
              var_X = var_X,
              var_Y = var_Y,
              load_X = load_X,
              load_Y = load_Y,
              Z_X = Z_X,
              Z_Y = Z_Y,
              eta = fit_cca$cov[1:comp],
              tol_cov = tol_cov
  ))
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
Optimal_basis_est <- function(dat_X, dat_Y, grid_num, max_comp_num = 5, mean_zero = F, FVE){
  
  n <- dat_X$n
  p_X <- dat_X$p
  p_Y <- dat_Y$p
  
  fda_object_X <- dat_X$dat
  fda_object_Y <- dat_Y$dat
  
  # Estimation of mean functions
  if(mean_zero == T){
    Fit_mean_X <- lapply(1:p_X, function(i){
      list(mu = rep(0, grid_num), workGrid = seq(0, 1, length.out = grid_num))
    })
    
    Fit_mean_Y <- lapply(1:p_Y, function(i){
      list(mu = rep(0, grid_num), workGrid = seq(0, 1, length.out = grid_num))
    })
  }else{
    
    Lt <- unique(unlist((lapply(1:n, function(j) fda_object_X[[1]][[j]]$time))))
    
    if(length(Lt) <= 5){
      Fit_mean_X <- lapply(1:p_X, function(i, fda_object, grid_num, n){
        Lt <- unique(unlist((lapply(1:n, function(j) fda_object[[1]][[j]]$time))))
        Ly <- sapply(1:n, function(j) sapply(1:length(Lt), function(k){
          mark <- which(fda_object[[i]][[j]]$time == Lt[k])
          if(length(mark) > 0){
            fda_object[[i]][[j]]$val[mark]
          }else{
            NA
          }
        }))
        
        dis_mu <- sapply(1:length(Lt), function(t) mean(Ly[t,], na.rm = T))
        mu <- spline(x = Lt, y = dis_mu, xout = seq(0, 1, length.out = grid_num))$y
        
        return(list(mu = mu, workGrid = seq(0, 1, length.out = grid_num)))
      }, fda_object = fda_object_X, grid_num = grid_num, n = n)
      
    }else{
      
      Fit_mean_X <- lapply(1:p_X, function(i, fda_object, grid_num, n){
        Ly <- sapply(1:n, function(j) fda_object[[i]][[j]]$val)
        Lt <- sapply(1:n, function(j) fda_object[[i]][[j]]$time)
        
        mu <- GetMeanCurve(Ly, Lt, optns = list(error = T, nRegGrid = grid_num, 
                                                methodBwMu = "GCV",
                                                methodMuCovEst = "smooth"
        ))$mu
        
        return(list(mu = mu, workGrid = seq(0, 1, length.out = grid_num)))
      }, fda_object = fda_object_X, grid_num = grid_num, n = n)
    }
    
    Lt <- unique(unlist((lapply(1:n, function(j) fda_object_Y[[1]][[j]]$time))))
    
    if(length(Lt) <= 5){
      Fit_mean_Y <- lapply(1:p_Y, function(i, fda_object, grid_num, n){
        Lt <- unique(unlist((lapply(1:n, function(j) fda_object[[i]][[j]]$time))))
        Ly <- sapply(1:n, function(j) sapply(1:length(Lt), function(k){
          mark <- which(fda_object[[i]][[j]]$time == Lt[k])
          if(length(mark) > 0){
            fda_object[[i]][[j]]$val[mark]
          }else{
            NA
          }
        }))
        
        dis_mu <- sapply(1:length(Lt), function(t) mean(Ly[t,], na.rm = T))
        mu <- spline(x = Lt, y = dis_mu, xout = seq(0, 1, length.out = grid_num))$y
        
        return(list(mu = mu, workGrid = seq(0, 1, length.out = grid_num)))
      }, fda_object = fda_object_Y, grid_num = grid_num, n = n)
      
    }else{
      Fit_mean_Y <- lapply(1:p_Y, function(i, fda_object, grid_num, n){
        Ly <- sapply(1:n, function(j) fda_object[[i]][[j]]$val)
        Lt <- sapply(1:n, function(j) fda_object[[i]][[j]]$time)
        
        mu <- GetMeanCurve(Ly, Lt, optns = list(error = T, nRegGrid = grid_num, 
                                                methodBwMu = "GCV",
                                                methodMuCovEst = "smooth"
        ))$mu
        
        return(list(mu = mu, workGrid = seq(0, 1, length.out = grid_num)))
      }, fda_object = fda_object_Y, grid_num = grid_num, n = n)
    }
  }
  
  time_mark_X <- lapply(1:n, function(i){
    sapply(1:length(fda_object_X[[1]][[i]]$time), function(h){
      which.min(abs(fda_object_X[[1]][[i]]$time[h] - Fit_mean_X[[1]]$workGrid))
    })
  }) 
  time_mark_Y <- lapply(1:n, function(i){
    sapply(1:length(fda_object_Y[[1]][[i]]$time), function(h){
      which.min(abs(fda_object_Y[[1]][[i]]$time[h] - Fit_mean_Y[[1]]$workGrid))
    })
  }) 
  
  # Optimal basis estimation
  G_X <- list()
  G_Y <- list()
  
  time <- seq(0, 1, length.out = grid_num)
  
  for(i in 1:n){
    G <- sapply(1:p_X, function(j_1){
      sapply(1:p_Y, function(j_2){
        tcrossprod(fda_object_X[[j_1]][[i]]$val - Fit_mean_X[[j_1]]$mu[time_mark_X[[i]]],
                   fda_object_Y[[j_2]][[i]]$val - Fit_mean_Y[[j_2]]$mu[time_mark_Y[[i]]])
      }, simplify = "array")
    }, simplify = "array")
    
    G_X[[i]] <- sapply(1:length(time_mark_X[[i]]), function(k){
      sapply(1:length(time_mark_X[[i]]), function(g){
        mean(G[k,,,] * G[g,,,])
      })
    }, simplify = "array")
    
    G_Y[[i]] <- sapply(1:length(time_mark_Y[[i]]), function(k){
      sapply(1:length(time_mark_Y[[i]]), function(g){
        mean(G[,k,,] * G[,g,,])
      })
    }, simplify = "array")
    
    rm(G)
    gc()
  }
  
  ## X data
  op_basis_X <- FSVD_r(M = G_X, time, time_mark = time_mark_X, R_max = max_comp_num, FVE)
  
  ## Y data
  op_basis_Y <- FSVD_r(M = G_Y, time, time_mark = time_mark_Y, R_max = max_comp_num, FVE)
  
  return(list(
    Mean = list(X = Fit_mean_X, Y = Fit_mean_Y),
    Time_mark = list(X = time_mark_X, Y = time_mark_Y),
    Optimal_bais = list(X = op_basis_X, Y = op_basis_Y))
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
#'   * \code{cca_AX(lambda, tol, Initial_cca, n, group_x, group_y)}
#'   * \code{cv_slcca(lambda, Initial_cca, n, group_x, group_y, max_var)}
#' @importFrom stats sd
#' @export
SLCCA_tr <- function(group_x, group_y, print_tuning = F, length_tun = 30, max_var, Initial_cca, n, p_X, p_Y){
  
  # Use grid search to determine the initial tuning parameters
  vec <- c(Initial_cca$Cov %*% Initial_cca$Y_coef) 
  max_t <- max(sapply(unique(group_x), function(j) sqrt(sum(vec[which(group_x == j)] ^ 2)))) 
  
  fit <- cca_AX(lambda = c(max_t, max_t), tol = 0.001, Initial_cca, n, group_x, group_y)
  var_num_X <- length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0))
  var_num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0))
  
  while(var_num_X + var_num_Y != 0){
    max_t <- max_t * 5
    fit <- cca_AX(lambda = c(max_t, max_t), tol = 0.001, Initial_cca, n, group_x, group_y)
    var_num_X <- length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0))
    var_num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0))
  }
  
  lambda_list <- c(max_t / 2, max_t / 10, max_t / 100, max_t / 1000, max_t / 10000)
  var_num_X <- matrix(NA, length(lambda_list), length(lambda_list))
  var_num_Y <- matrix(NA, length(lambda_list), length(lambda_list))
  
  for(a in 1:length(lambda_list)){
    for(b in 1:length(lambda_list)){
      fit <- cca_AX(lambda = c(lambda_list[a], lambda_list[b]), tol = 0.001, Initial_cca, n, group_x, group_y)
      var_num_X[b,a] <- length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0))
      var_num_Y[b,a] <- length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0))
    }
  }
  
  A <- var_num_X >= max_var[1]
  if(length(which(colMeans(A) > 0)) == 0){
    min_X <- max_t / 100000
    fit <- cca_AX(lambda = c(min_X, max_t), tol = 0.001, Initial_cca, n, group_x, group_y)
    while(length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0)) <= max_var[1]){
      min_X <- min_X / 10
      fit <- cca_AX(lambda = c(min_X, max_t), tol = 0.001, Initial_cca, n, group_x, group_y)
    }
  }else{
    min_X <- lambda_list[min(which(colMeans(A) > 0))]
  }
  
  A <- var_num_Y >= max_var[2]
  if(length(which(rowMeans(A) > 0)) == 0){
    min_Y <- max_t / 100000
    fit <- cca_AX(lambda = c(min_X, min_Y), tol = 0.001, Initial_cca, n, group_x, group_y)
    while(length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0)) <= max_var[2]){
      min_Y <- min_Y / 10
      fit <- cca_AX(lambda = c(min_X, min_Y), tol = 0.001, Initial_cca, n, group_x, group_y)
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
      if(b > 1){
        if(CV_error[b-1,a] != -10 ^ (11)){
          fit <- cv_slcca(lambda = c(lambda_list_X[a], lambda_list_Y[b,a]), Initial_cca, n, group_x, group_y, max_var)
          CV_error[b,a] <- fit$error
          if(CV_error[b,a] == max(CV_error, na.rm = T)){
            lambda <- c(lambda_list_X[a], lambda_list_Y[b,a])
          }
        }else{
          CV_error[b,a] <- -10 ^ (11)
          if(CV_error[b,a] == max(CV_error, na.rm = T)){
            lambda <- c(lambda_list_X[a], lambda_list_Y[b,a])
          }
        }
        
      }else{
        fit <- cv_slcca(lambda = c(lambda_list_X[a], lambda_list_Y[b,a]), Initial_cca, n, group_x, group_y, max_var)
        CV_error[b,a] <- fit$error
        if(CV_error[b,a] == max(CV_error, na.rm = T)){
          lambda <- c(lambda_list_X[a], lambda_list_Y[b,a])
        }
      }
      
      if(print_tuning == T){
        print(CV_error)
      }
      # gc()
    }
  }
  
  CV_error[CV_error == - 10 ^ 10] <- NA
  CV_error[CV_error == - 10 ^ 11] <- NA
  sd_CV <- 0.1 * sd(CV_error[is.finite(CV_error)])
  
  mark <- (CV_error >= max(CV_error, na.rm = T) - sd_CV)
  
  row_mark <- min(which(apply(mark, c(2), function(x) mean(x, na.rm = T)) != 0))
  lambda_1 <- c(lambda_list_X[row_mark], lambda_list_Y[min(which(mark[,row_mark] == T)), row_mark])
  
  col_mark <- min(which(apply(mark, c(1), function(x) mean(x, na.rm = T)) != 0))
  lambda_2 <- c(lambda_list_X[min(which(mark[col_mark,] == T))], lambda_list_Y[col_mark, min(which(mark[col_mark,] == T))])
  
  fit_opt_1 <- cv_slcca(lambda = lambda_1, Initial_cca, n, group_x, group_y, max_var)
  fit_opt_2 <- cv_slcca(lambda = lambda_2, Initial_cca, n, group_x, group_y, max_var)
  
  if(sum(sapply(unique(group_x), function(j) sum(fit_opt_1$fit_cca$X_coef[group_x == j] ^ 2)) != 0) + 
     sum(sapply(unique(group_y), function(j) sum(fit_opt_1$fit_cca$Y_coef[group_y == j] ^ 2)) != 0) <=
     sum(sapply(unique(group_x), function(j) sum(fit_opt_2$fit_cca$X_coef[group_x == j] ^ 2)) != 0) + 
     sum(sapply(unique(group_y), function(j) sum(fit_opt_2$fit_cca$Y_coef[group_y == j] ^ 2)) != 0)
  ){
    fit_opt <- fit_opt_1
  }else{
    fit_opt <- fit_opt_2
  }
  cor <- fit_opt$error
  
  return(list(fit_opt = fit_opt,
              error = cor,
              CV_error = CV_error,
              var_X = which(sapply(unique(group_x), function(j) sum(fit_opt$fit_cca$X_coef[group_x == j] ^ 2)) != 0),
              var_Y = which(sapply(unique(group_y), function(j) sum(fit_opt$fit_cca$Y_coef[group_y == j] ^ 2)) != 0)
  ))
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
  
  op_basis$comp_num <- min(max(sum(cumsum(op_basis$eigenvalue) / tol_norm <= FVE) + 1, 2), R_max)
  
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
fit_basis <- function(M, time_mark, time){
  
  fit_cov <- fit_cov_surface(M, time_mark, time_mark, time)
  fit <- base::eigen(fit_cov, symmetric = T)
  eigenfunction <- fit$vectors[,1] * sqrt(nrow(fit_cov))
  
  return(
    list(
      eigenvalue = c(t(eigenfunction) %*% fit_cov %*% eigenfunction / nrow(fit_cov) ^ 2),
      eigenfunction = eigenfunction,
      eigendocomposition = fit,
      tol_nor = sum(fit[[1]]) / nrow(fit_cov)
    )
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
fit_cov_surface <- function(M, time_mark_X, time_mark_Y, time){
  
  n <- length(M)
  matrix_sample <- vector()
  for(i in 1:n){
    A <- sapply(1:length(time_mark_X[[i]]), function(k){
      cbind(rep(time[time_mark_X[[i]]][k], length(time_mark_Y[[i]])),
            time[time_mark_Y[[i]]],
            M[[i]][k,]
      )
    }, simplify = "array")
    
    matrix_sample <- rbind(matrix_sample, data.frame(
      x = c(A[,1,]),
      y = c(A[,2,]),
      z = c(A[,3,]),
      w = sapply(1:length(A[,1,]), function(k) 1 / (length(time_mark_X[[i]]) * length(time_mark_Y[[i]])))
    ))
  }
  
  uni_time_x <- sort(unique(matrix_sample$x))
  uni_time_y <- sort(unique(matrix_sample$y))
  num_k <- min(c(length(uni_time_x), length(uni_time_y), 15))
  
  x <- rep(time, length(time))
  y <- c(sapply(time, function(k) rep(k, length(time))))
  
  if(length(uni_time_x) * length(uni_time_x) <= 25){
    dat_tran <- sapply(1:length(uni_time_x), function(k){
      sapply(1:length(uni_time_y), function(j){
        obs <- matrix_sample$z[intersect(which(matrix_sample$x == uni_time_x[k]), which(matrix_sample$y == uni_time_y[j]))]
        return(c(uni_time_x[k], uni_time_y[j], mean(obs), length(obs)))
      }, simplify = "array")
    }, simplify = "array")
    
    matrix_sample <- data.frame(
      x = c(dat_tran[1,,]),
      y = c(dat_tran[2,,]),
      z = c(dat_tran[3,,])
    )
    rm(dat_tran)
    
    fit <- Tps(x = cbind(matrix_sample$x, matrix_sample$y), Y = matrix_sample$z, df = nrow(matrix_sample) - 2, method = "GCV.model")
    
    Cov_fit <- c(predictSurface(fit, grid.list = list(x = time, y = time),
                                extrap = T)$z)
    
  }else{
    
    fit <- gam(z ~ te(x, y, bs = "cr", k = c(num_k, num_k)), data = matrix_sample, weights = matrix_sample$w, family = "gaussian", method="GCV.Cp")
    Cov_fit <- predict(fit, newdata = list(x = x, y = y))
    
  }
  rm(fit)
  
  return(matrix(Cov_fit, nrow = length(time)))
}

#' Penalized CCA Core Solver
#'
#' Performs alternating updates for penalized canonical correlation analysis with group lasso penalties.
#'
#' @param lambda Numeric vector of length 2; penalties for X and Y groups.
#' @param tol Numeric; convergence threshold for relative change in coefficients.
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
cca_AX <- function(lambda, tol, Initial_cca, n, group_x, group_y){
  
  a_X <- Initial_cca$X_coef
  a_Y <- Initial_cca$Y_coef
  
  error <- 1
  ite <- 0
  
  num_X <- length(unique(group_x))
  num_Y <- length(unique(group_y))
  
  while((error > tol) & (ite < 200)){
    ite <- ite + 1
    
    a_X_t <- a_X
    a_Y_t <- a_Y
    
    a_X <- c(Initial_cca$Cov %*% a_Y) 
    a_X <- c(sapply(1:num_X, function(j){
      norm <- sqrt(sum(a_X[group_x == j] ^ 2))
      if(norm == 0){
        rep(0, sum(group_x == j))
      }else{
        max((1 - lambda[1] / norm), 0) * a_X[group_x == j]
      }
    })) 
    x_norm <- c(sqrt(crossprod(a_X)))
    if(x_norm != 0){
      a_X <- a_X / x_norm
    }
    
    a_Y <- c(t(a_X) %*% Initial_cca$Cov) 
    a_Y <- c(sapply(1:num_Y, function(j){
      norm <- sqrt(sum(a_Y[group_y == j] ^ 2))
      if(norm == 0){
        rep(0, sum(group_y == j))
      }else{
        max((1 - lambda[2] / norm), 0) * a_Y[group_y == j]
      }
    })) 
    y_norm <- c(sqrt(crossprod(a_Y)))
    if(y_norm != 0){
      a_Y <- a_Y / y_norm
    }
    
    if(x_norm != 0){
      error_x <- (sum(abs(a_X_t - a_X)) / sum(abs(a_X_t))) 
    }else{
      error_x <- 0
    }
    
    if(y_norm != 0){
      error_y <- (sum(abs(a_Y_t - a_Y)) / sum(abs(a_Y_t)))
    }else{
      error_y <- 0
    }
    
    error <- error_x + error_y
  }
  
  
  return(list(X_coef = as.numeric(a_X),
              Y_coef = as.numeric(a_Y)
  ))
}

#' Cross-Validation Wrapper for Sparse CCA
#'
#' Performs K-fold cross-validation for penalized canonical correlation analysis
#' using the core solver `cca_AX`, and computes mean absolute correlation on held-out folds.
#'
#' @param lambda Numeric vector of length 2; penalties for X and Y groups.
#' @param Initial_cca List; initial CCA state containing `X_coef`, `Y_coef`, and `Cov`.
#' @param n Integer; number of subjects.
#' @param group_x Integer vector; group labels for X coefficients.
#' @param group_y Integer vector; group labels for Y coefficients.
#' @param max_var Numeric(2); maximum nonzero groups allowed for (X, Y).
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
cv_slcca <- function(lambda, Initial_cca, n, group_x, group_y, max_var){
  
  fit_cca <- cca_AX(lambda, tol = 0.001, Initial_cca, n, group_x, group_y)
  num_X <- length(which(sapply(unique(group_x), function(j) sum(fit_cca$X_coef[group_x == j] ^ 2)) != 0))
  num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit_cca$Y_coef[group_y == j] ^ 2)) != 0))
  
  if((num_X > max_var[1]) | (num_Y > max_var[2])){
    error <- - 10 ^ 11
  }else if((num_X == 0) | (num_Y == 0)){
    error <- - 10 ^ 10
  }else{
    error <- mean(sapply(1:Initial_cca$cv_num, function(k){
      cca_cv_fit <- cca_AX(lambda, tol = 0.001, 
                           Initial_cca = list(X_coef = fit_cca$X_coef, Y_coef = fit_cca$Y_coef, 
                                              Cov = Initial_cca$CV[[k]]$train), n, 
                           group_x, group_y)
      
      cor <- c(t(cca_cv_fit$X_coef) %*% Initial_cca$CV[[k]]$test %*% cca_cv_fit$Y_coef)
      
      return(cor)
    }))
  }
  
  return(list(error = error, 
              fit_cca = fit_cca,
              lambda = lambda
  ))
}

#' Plot FACD Loading Surfaces for X and Y Features
#'
#' Creates a heatmap visualization of the top canonical loading surfaces for X and Y
#' from a FACD analysis, highlighting the most important features over specified time points.
#'
#' @param fit A list returned by \code{FACD()}, containing loading surfaces and basis information.
#' @param dat_X A list with component \code{feature_label}: character vector of X feature names.
#' @param dat_Y A list with component \code{feature_label}: character vector of Y feature names.
#' @param hight_feature Integer; maximum number of top features to display (default: 5).
#' @param k Integer; the component index to visualize (e.g., 1 for the first canonical component).
#' @param name_X Character; label for the X data facet (e.g., "X").
#' @param name_Y Character; label for the Y data facet (e.g., "Y").
#' @param time_illu Numeric vector; specific grid time points to illustrate.
#' @param time_label Character vector; axis labels corresponding to \code{time_illu}.
#' @param title_name Character; plot title.
#'
#' @return A list with:
#'   \describe{
#'     \item{plot}{A \code{ggplot} heatmap object.}
#'     \item{feature_X}{Character vector of selected X feature names.}
#'     \item{feature_Y}{Character vector of selected Y feature names.}}
#'
#' @import ggplot2 ggrepel
#' @export
plot_FACD <- function(fit, dat_X, dat_Y, hight_feature = 5, k, name_X, name_Y, time_illu, time_label, title_name){
  # Columns evaluated by ggplot2 inside dat_plot.
  time <- label <- val <- NULL
  
  ## Longitudinal canonical vector
  Canonical_X <- fit$load_X[[k]]
  Canonical_X <- Canonical_X * sign(sum(Canonical_X))
  Canonical_Y <- fit$load_Y[[k]]
  Canonical_Y <- Canonical_Y * sign(sum(Canonical_Y))
  
  ### X
  order_X <- order(rowSums(Canonical_X ^ 2), decreasing = T)
  nonzeo_X <- sum(sort(rowSums(Canonical_X ^ 2), decreasing = T) > 0)
  feature_X <- dat_X$feature_label[order_X[1:nonzeo_X]]
  
  mark <- sapply(1:length(time_illu), function(h) which.min(abs(fit$Optimal_basis$X$workGrid - time_illu[h]))) 
  Lcv_X <- matrix(sapply(1:min(hight_feature, nonzeo_X), function(l){
    Canonical_X[order_X[l],mark]
  }), ncol = min(hight_feature, nonzeo_X))
  
  imp_fea <- sapply(1:dim(Lcv_X)[2], function(k) dat_X$feature_label[order_X[k]])
  imp_fea <- factor(c(matrix(imp_fea, nrow = nrow(Lcv_X), ncol = length(imp_fea), byrow = T)),
                    levels = rev(imp_fea))
  dat_plot_X <- data.frame(time = rep(time_illu, ncol(Lcv_X)),
                           val  = c(Lcv_X), 
                           label = imp_fea,
                           data_type = rep(name_X, length(imp_fea))
  )
  
  ### Y
  order_Y <- order(rowSums(Canonical_Y ^ 2), decreasing = T)
  nonzeo_Y <- sum(sort(rowSums(Canonical_Y ^ 2), decreasing = T) > 0)
  feature_Y <- dat_Y$feature_label[order_Y[1:nonzeo_Y]]
  
  Lcv_Y <- matrix(sapply(1:min(hight_feature, nonzeo_Y), function(l){
    Canonical_Y[order_Y[l],mark]
  }), ncol = min(hight_feature, nonzeo_Y))
  
  imp_fea <- sapply(1:dim(Lcv_Y)[2], function(k) dat_Y$feature_label[order_Y[k]])
  imp_fea <- factor(c(matrix(imp_fea, nrow = dim(Lcv_Y)[1], ncol = length(imp_fea), byrow = T)),
                    levels = rev(imp_fea))
  dat_plot_Y <- data.frame(time = rep(time_illu, ncol(Lcv_Y)),
                           val  = c(Lcv_Y), 
                           label = imp_fea,
                           data_type = rep(name_Y, length(imp_fea))
  )
  
  dat_plot <- rbind(dat_plot_X, dat_plot_Y)
  
  pl <- ggplot(dat_plot) + 
    geom_tile(aes(x = time, y = label, fill = val)) +
    scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "orange", midpoint = 0) +
    # geom_line(aes(x = time, y = val, group = label), color = "gray", alpha = 0.3) +
    # geom_line(aes(x = time, y =  high_fea, group = label), color = "black") +
    # ggrepel::geom_text_repel(aes(x = high_x, y = high_y, label = label)) +
    facet_wrap(~data_type, nrow = 1, scales = "free_y") +
    # geom_hline(yintercept = 0, color = "blue", linetype = 2) +
    theme_bw(base_family = "Times") +
    scale_x_continuous(breaks = time_illu, labels = time_label) +
    xlab("") + ylab("") + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "right",
          text = element_text(size = 12),
          panel.border = element_blank(),
          # text = element_text(family = "STHeiti"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0)) +
    labs(title = title_name, fill = "")
  
  return(list(plot = pl,
              feature_X = feature_X,
              feature_Y = feature_Y))
}

#' Wrap Raw Data into FACD-Compatible Format
#'
#' Converts a wide-format data matrix and accompanying sample metadata into the
#' nested list structure required by FACD, scaling values and normalizing time.
#'
#' @param dat A data frame with column names as feature labels.
#' @param dat_count A data frame containing at least columns \code{SubjectID} and \code{timepoint},
#'   mapping each row of \code{dat} to a specific subject and time index.
#' @param time_label A finite numeric vector used to map \code{timepoint} indices
#'   to actual times. Observed times must contain at least two distinct values.
#'
#' @return A list with components:
#'   \describe{
#'     \item{dat}{Nested list of length \code{p}, each a list of length \code{n} with elements \code{time} and \code{val}.}
#'     \item{n}{Number of subjects.}
#'     \item{p}{Number of features.}
#'     \item{sub_label}{Vector of unique subject identifiers.}
#'     \item{feature_label}{Vector of feature names (column names of \code{dat}).}
#'     \item{time_grid}{Sorted, scaled unique time points across all features.}
#'     \item{time_min}{Minimum raw time before scaling.}
#'     \item{time_max}{Maximum raw time before scaling.}}
#'
#' @details
#' This function rescales observed times to [0, 1] using the minimum and maximum
#' across all subjects and features in this dataset. The original bounds are
#' returned as \code{time_min} and \code{time_max}. For datasets sharing a physical
#' time axis, check that separate calls use the same units and time bounds.
#'
#' After averaging duplicate subject/timepoint measurements, each feature is
#' centered by its pooled mean and divided by \code{sqrt(mean(values^2))}, the
#' RMS of its uncentered values. This is not unit-variance z-score scaling.
#' To apply the recommended z-score preprocessing, standardize the raw feature
#' columns across all subjects and time points before calling \code{wrap_dat()}.
#' Remove features with no finite observations or zero variance beforehand.
#'
#' @export
wrap_dat <- function(dat, dat_count, time_label) {
  feature_label <- colnames(dat)
  sub_label     <- unique(dat_count$SubjectID)
  
  n <- length(sub_label)
  p <- length(feature_label)
  
  # Build nested list of time/value pairs
  fda_object <- lapply(seq_len(p), function(i) {
    lapply(seq_len(n), function(j) {
      idx <- which(dat_count$SubjectID == sub_label[j])
      times_raw <- dat_count$timepoint[idx]
      vals_raw  <- dat[idx, feature_label[i]]
      
      # Aggregate by timepoint
      uniq_t <- sort(unique(times_raw))
      vals   <- sapply(uniq_t, function(t) mean(vals_raw[times_raw == t]))
      times  <- time_label[uniq_t]
      
      data.frame(time = times, val = vals)[is.finite(vals), , drop = FALSE]
    })
  })
  
  # Scale and normalize time to [0,1]
  all_times <- sort(unique(unlist(lapply(fda_object, function(feat) {
    unlist(lapply(feat, `[[`, "time"))
  })))
  )
  min_t <- min(all_times)
  max_t <- max(all_times)
  time_grid <- (all_times - min_t) / (max_t - min_t)
  
  # Center & scale values per feature
  for (i in seq_len(p)) {
    vals_all <- unlist(lapply(fda_object[[i]], `[[`, "val"))
    mu   <- mean(vals_all)
    sigma <- sqrt(mean(vals_all^2))
    for (j in seq_len(n)) {
      if (nrow(fda_object[[i]][[j]]) > 0) {
        fda_object[[i]][[j]]$time <-
          (fda_object[[i]][[j]]$time - min_t) / (max_t - min_t)
        fda_object[[i]][[j]]$val  <-
          (fda_object[[i]][[j]]$val  - mu) / sigma
      }
    }
  }
  
  list(
    dat           = fda_object,
    n             = n,
    p             = p,
    sub_label     = sub_label,
    feature_label = feature_label,
    time_grid     = time_grid,
    time_min      = min_t,
    time_max      = max_t
  )
}