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
#   dat          : nested list of length p, each a list of length n with 
# elements $time and $val
#   sub_label    : vector of subject identifiers
#   feature_label: vector of feature names
#   time_grid    : transformed unique time points on [0,1]
#   time_min     : minimum original time
#   time_max     : maximum original time
# These objects can be created from raw data via wrap_dat().

data(dat_X)
data(dat_Y)

# 2. FACD estimation ------------------------------------------------------

# Parameters
grid_num   <- 201                        # number of grid points in [0,1]
max_var    <- floor(c(min(dat_X$p/2, dat_X$n),
                      min(dat_Y$p/2, dat_Y$n)))  # max selected features 

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
cor <- abs(crossprod(X_cca, Y_cca) / 
             sqrt(crossprod(X_cca, X_cca) * crossprod(Y_cca, Y_cca)))
dat_plot <- data.frame(
  X = X_cca, Y = Y_cca,
  cor = c(cor, rep(NA, n - 1)), var1 = rep("Metabolite", n), 
  var2 = rep("Lipid", n),
  text_x = c(min(X_cca), rep(NA, n - 1)), 
  text_y = c(max(Y_cca) + 0.2, rep(NA, n - 1))
)

ggplot(dat_plot) + 
  geom_point(aes(x = X, y = Y)) +  
    geom_text_repel(
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
cat("Top Meta features:\n", paste(res_plot$feature_X, collapse = ", "), "\n")
cat("Top Lipid features:\n", paste(res_plot$feature_Y, collapse = ", "), "\n")
