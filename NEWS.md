# FACD 0.2.0

- Rename the package and public analysis/plotting functions from CCVA to FACD.
  Use `library(FACD)`, `FACD()`, and `plot_FACD()`.
- Update multi-component deflation, loading orthogonalization, and residual scores;
  return component cross-covariances (`eta`) and total covariance norm (`tol_cov`).
- Expose `max_comp_num = 5` in `FACD()` so basis estimation has a defined limit.
- Add regression coverage for cross-validation with unequal X/Y dimensions.
- Regenerate function/data documentation and update package installation tests.
