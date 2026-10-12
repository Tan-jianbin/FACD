test_that("FACD fits and plots bundled data without global configuration", {
  e <- new.env()
  data("dat_X", "dat_Y", package = "FACD", envir = e)
  subset_features <- function(dat, p) {
    dat$dat <- dat$dat[seq_len(p)]
    dat$feature_label <- dat$feature_label[seq_len(p)]
    dat$p <- p
    dat
  }
  x <- subset_features(e$dat_X, 6L)
  y <- subset_features(e$dat_Y, 4L)

  for (sparse in c(FALSE, TRUE)) {
    fit <- FACD(dat_X = x, dat_Y = y, grid_num = 21,
                max_var = c(3, 2), length_tun = 5,
                print_tuning = FALSE, Sparse = sparse, comp = 2)
    expect_length(fit$Z_X, 2)
    expect_length(fit$Z_Y, 2)
    expect_length(fit$eta, 2)
    expect_true(all(is.finite(unlist(fit[c("Z_X", "Z_Y", "eta", "tol_cov")]))))
    for (k in 1:2) {
      expect_equal(dim(fit$load_X[[k]]), c(6L, 21L))
      expect_equal(dim(fit$load_Y[[k]]), c(4L, 21L))
    }
    expect_equal(sum(fit$load_X[[1]] * fit$load_X[[2]]), 0, tolerance = 1e-8)
    expect_equal(sum(fit$load_Y[[1]] * fit$load_Y[[2]]), 0, tolerance = 1e-8)
    plotted <- plot_FACD(fit, x, y, hight_feature = 2, k = 1,
                         name_X = "X", name_Y = "Y", time_illu = seq(0, 1, .25),
                         time_label = as.character(0:4), title_name = "FACD")
    expect_s3_class(plotted$plot, "ggplot")
    expect_true(nrow(ggplot2::ggplot_build(plotted$plot)$data[[1]]) > 0)
  }
})

test_that("cross-validation counts Y groups using Y's dimension", {
  # Five Y coefficients versus three X coefficients: an X mask would produce NA.
  initial <- list(Cov = matrix(c(4, 0, 0, 0, 0,
                                 0, 3, 0, 0, 0,
                                 0, 0, 2, 0, 0), 3, 5, byrow = TRUE),
                  X_coef = c(1, 0, 0), Y_coef = c(1, 0, 0, 0, 0))
  initial$cv_num <- 2
  initial$CV <- rep(list(list(train = initial$Cov, test = initial$Cov)), 2)
  fit <- cv_slcca(c(0, 0), initial, n = 10,
                  group_x = 1:3, group_y = 1:5, max_var = c(3, 5))
  expect_true(is.finite(fit$error))
  expect_equal(fit$error, 4)
})
