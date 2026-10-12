# FACD: Functional-Aggregated Cross-covariance Decomposition for High-Dimensional Longitudinal Data

**FACD** is an R package implementing the methodology from the paper **"Associating High-Dimensional Longitudinal Datasets through an Efficient Cross-Covariance Decomposition"** by *Jianbin Tan* and *Pixu Shi*, available on [arXiv](https://arxiv.org/abs/2601.13405).

## Description: 
 FACD encodes associations between high-dimensional longitudinal datasets by providing a data-adaptive basis to capture longitudinal signals and incorporating variable selection to facilitate the discovery of key associations.

## Installation

```r
install.packages("remotes")
remotes::install_github("Tan-jianbin/FACD")
```

## Data preprocessing

Before calling `FACD()`:

- **Scale time to `[0, 1]`.** Transform all observed `$time` values and
  `time_grid`, for example with `(time - time_min) / (time_max - time_min)`.
  The time range must be nonzero. Use one mapping across all subjects and
  features within a dataset. If X and Y share a physical time axis, use the
  same units and reference interval for both.
- **Standardize features before analysis**, especially when their units or
  scales differ. A recommended choice is z-score scaling: subtract each
  feature's mean and divide by its standard deviation, pooling observations
  across subjects and time points. Do not scale each subject separately.
  Remove features with no finite observations or zero variance first.
- **Match subjects across datasets** and keep them in the same order.

`FACD()` does not automatically rescale time or standardize feature variances.
`mean_zero = FALSE` estimates and subtracts mean curves; it does not replace
feature standardization.

`wrap_dat()` automatically rescales time to `[0, 1]` and centers feature values,
but divides by the RMS of the uncentered values rather than their standard
deviation. For z-score preprocessing, apply `scale()` to the raw feature columns
before wrapping, as illustrated below. Each call uses its own observed time
bounds, so check `time_min` and `time_max` when X and Y share a physical time axis.
The bundled example datasets are already prepared and need no further scaling.

## Example

The full worked example is stored in **example.R** and is reproduced below:

```r
#-------------------------------------------------------------
# Example code
#-------------------------------------------------------------

# Load required packages
library(FACD)
library(ggplot2)
library(ggrepel)

# 1. Load data --------------------------------------------------
# The dat_X and dat_Y objects should have the following structure:
#   n            : sample size
#   p            : number of features (dimension)
#   dat          : nested list of length p, each a list of length n with elements $time and $val
#   sub_label    : vector of subject identifiers
#   feature_label: vector of feature names
#   time_grid    : transformed unique time points on [0,1]
#   time_min     : minimum original time
#   time_max     : maximum original time
# These objects can be created from raw data via wrap_dat().
#
# Before FACD(): all observed $time values and time_grid must be in [0, 1].
# Use one time mapping across subjects/features, not one per subject.
# Feature-wise z-score standardization across subjects and times is recommended.
# wrap_dat() rescales time automatically and centers/scales values by their RMS;
# its value scaling is not the same as dividing by the standard deviation.
# The bundled dat_X and dat_Y are already prepared; do not preprocess them again.
#
# For your own raw data (rows = observations; columns = numeric features):
# Remove all-missing/non-finite and zero-variance features before scaling.
# raw_X <- as.data.frame(scale(raw_X))
# raw_Y <- as.data.frame(scale(raw_Y))
# dat_X <- wrap_dat(raw_X, count_X, time_label_X)
# dat_Y <- wrap_dat(raw_Y, count_Y, time_label_Y)
# count_X/count_Y contain SubjectID and integer timepoint indices into the
# corresponding numeric time_label vector (e.g., c(0, 2, 15, 30, 60)).
# Each dataset must have at least two distinct observed times. If X and Y share
# a physical time axis, check that wrap_dat() uses the same time bounds for both.
# Subjects must also appear in the same order in dat_X and dat_Y.

data(dat_X)
data(dat_Y)

# 2. FACD estimation ------------------------------------------------------

# Parameters
grid_num   <- 201                        # number of grid points in [0,1]
max_var    <- floor(c(min(dat_X$p/2, dat_X$n),
                      min(dat_Y$p/2, dat_Y$n)))  # max selected features per dataset

# Run FACD (first canonical component only)
Result <- FACD(
  grid_num     = grid_num,
  dat_X        = dat_X,
  dat_Y        = dat_Y,
  max_var      = max_var,
  length_tun   = 30,
  print_tuning = FALSE,
  comp         = 1,
  FVE          = 0.95
)

# 3. Correlation of canonical variables ------------------------------------
X_cca <- Result$Z_X[[1]]
Y_cca <- Result$Z_Y[[1]]
n <- dat_X$n

# Prepare data frame for plotting
cor <- abs(crossprod(X_cca, Y_cca) / sqrt(crossprod(X_cca, X_cca) * crossprod(Y_cca, Y_cca)))
dat_plot <- data.frame(
  X = X_cca, Y = Y_cca,
  cor = c(cor, rep(NA, n - 1)), var1 = rep("Metabolite", n), var2 = rep("Lipid", n),
  text_x = c(min(X_cca), rep(NA, n - 1)), text_y = c(max(Y_cca) + 0.2, rep(NA, n - 1))
)

ggplot(dat_plot) + 
  geom_point(aes(x = X, y = Y)) +  
    ggrepel::geom_text_repel(
      aes(
        x     = text_x,
        y     = text_y,
        label = paste0("Correlation ", round(cor, 3))
      ),
      size         = 4,
      max.overlaps = 100
    ) +
  theme_bw(base_family = "Times") +
  xlab("Metabolite") + ylab("Lipid") + 
  theme(legend.position = "none",
        text = element_text(size = 15),
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))

# 4. Heatmap of top loading surfaces --------------------------------------

# Define time points and labels
time_illu <- 0:4/4
time_lbl  <- c("Baseline", "2 min", "15 min", "30 min", "60 min")
h_feat    <- 10  # number of top features to display

# Generate heatmap plot and retrieve important features
res_plot <- plot_FACD(
  fit           = Result,
  dat_X         = dat_X,
  dat_Y         = dat_Y,
  hight_feature = h_feat,
  k             = 1,
  name_X        = "Metabolite",
  name_Y        = "Lipid",
  time_illu     = time_illu,
  time_label    = time_lbl,
  title_name    = "Longitudinal loadings"
)

# Display heatmap
print(res_plot$plot)

# Important features
cat("Top Metabolite features:\n", paste(res_plot$feature_X, collapse = ", "), "\n")
cat("Top Lipid features:\n",      paste(res_plot$feature_Y, collapse = ", "), "\n")
```
