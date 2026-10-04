# ZIBB_PGLMM_Squamata.R
#
# Squamata Zero-Inflated Beta-Binomial Phylogenetic GLMM Pipeline
# Model neoplasia and malignancy COUNTS across Squamata while accounting for the below:
#
#   1. Different numbers of records/necropsies among species
#   2. Excess zero tumor counts
#   3. Overdispersion beyond a normal binomial model
#   4. Shared evolutionary history among species
#   5. Life-history traits:
#        - maximum longevity
#        - gestation length
#        - adult body mass
#        - litter size
#
# RESPONSE DEFINITIONS & VARS
# --------------------
# NeoplasiaCases  = NeoplasiaWithDenominators
# MalignancyCases = Malignant
# Trials          = RecordsWithDenominators
#
# MODEL FAMILY
# ------------
# True Zero Inflated Beta Binomial:
#
# cases | trials(Trials) ~ predictors +
#     (1 | gr(Species, cov = A_model))
#
# zi ~ log_trials_s
#
# IMPORTANT
# ---------
# We deliberately DO NOT use prevalence as the response!!
#
# Instead of:
#
#     NeoplasiaPrevalence ~ predictor
#
# we use:
#
#     NeoplasiaCases | trials(Trials) ~ predictor
#
# This allows the model to distinguish:
#
#     2 tumors / 20 records
#
# from
#
#     20 tumors / 200 records
#
# even though both have prevalence = 0.10 !!
# Hence the whole point of ZIBB
#

# ==================
# 0. PACKAGES
# ==================

# Install these once if needed:
#
# install.packages(c(
#   "tidyverse",
#   "ape",
#   "brms",
#   "posterior",
#   "bayesplot"
# ))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(brms)
  library(posterior)
  library(bayesplot)
})


# Use up to 4 CPU cores for Stan
N_CORES <- min(4, parallel::detectCores())

options(mc.cores = N_CORES)

set.seed(1234)

# 1. IMPORTANT INPUTS AND DIRECTORIES ! 

# Input files
squamata_file <- "squamata_data.csv"
full_data_file <- "min20-2022.05.16.csv"
tree_file <- "min20Fixed516.nwk"


# Output directory
out_dir <- "ZIBB_PGLMM_Squamata"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# Usually we would like to see a minimum number of species required before fitting a model 
#
# For instance the old multivariable modeling reportedly has only ~14 species
# A ZIBB PGLMM with several predictors + phylogeny + zero inflation 
# is a bit...challenging at this size but oh well
# 
# 20 is a conservative default though
#
MIN_SPECIES <- 20


# Zero inflation model
#
# "log_trials":
#     asks whether sampling effort predicts whether a species ends up in the excess zero component
#
# "intercept":
#     estimates one overall excess zero probability
#
ZI_MODE <- "log_trials"

# Alternative:
# ZI_MODE <- "intercept"


# Stan settings
N_CHAINS <- 4
N_ITER   <- 4000
N_WARMUP <- 2000

ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15



# 2. CHECK IF TRUE ZIBB IS AVAILABLE! 


if (!exists(
  "zero_inflated_beta_binomial",
  envir = asNamespace("brms"),
  inherits = FALSE
)) {
  
  stop(
    "\nYour installed version of brms does not contain ",
    "zero_inflated_beta_binomial().\n",
    "Update brms before running this pipeline.\n\n",
    "This script intentionally does NOT automatically replace the model ",
    "with a regular zero inflated binomial because the requested model ",
    "is a true ZIBB.\n"
  )
}



# 3. LOAD THE DATA 

# VERY IMPORTANT 
# We would like to really prefere the Squamata only CSV BUT
# If for whatever reason someone doesnt have it or it isn't present we should
# always fall back to the complete AZA dataset from Zach and Walker --> Cancer Prevalence across Vertebrates
# and filter Orders == Squamata.

if (file.exists(squamata_file)) {
  
  message("Reading: ", squamata_file)
  
  raw <- read.csv(
    squamata_file,
    check.names = TRUE,
    stringsAsFactors = FALSE
  )
  
} else if (file.exists(full_data_file)) {
  
  message(
    squamata_file,
    " not found. Using complete dataset instead: ",
    full_data_file
  )
  
  raw <- read.csv(
    full_data_file,
    check.names = TRUE,
    stringsAsFactors = FALSE
  )
  
  if (!"Orders" %in% names(raw)) {
    stop("Cannot find the Orders column needed to identify Squamata.")
  }
  
  raw <- raw %>%
    filter(grepl("Squamata", Orders, ignore.case = TRUE))
  
} else {
  
  stop(
    "Neither squamata_data.csv nor min20-2022.05.16.csv could be found."
  )
}


message("Rows loaded: ", nrow(raw))



# 4. VERIFY OUR REQUIRED COLUMNS


required_columns <- c(
  "Species",
  "RecordsWithDenominators",
  "NeoplasiaWithDenominators",
  "Malignant",
  "adult_weight.g.",
  "max_longevity.months.",
  "Gestation.months.",
  "litter_size"
)


missing_columns <- setdiff(
  required_columns,
  names(raw)
)


if (length(missing_columns) > 0) {
  
  cat("\nColumns available in the CSV:\n\n")
  print(names(raw))
  
  stop(
    "\n\nMissing required columns:\n",
    paste(missing_columns, collapse = ", "),
    "\n\nCheck the CSV column names before continuing."
  )
}


# Family is useful for inspection/plots but not required for the model ! 
if (!"Family" %in% names(raw)) {
  raw$Family <- NA_character_
}



# 5. SANITIZE RAW DATA


# There are the old existing scripts that treat negative numerical values as missing
#
# Do this ONLY to numeric columns ! 
#
raw <- raw %>%
  mutate(
    across(
      where(is.numeric),
      ~ replace(.x, .x < 0, NA_real_)
    )
  )


# Clean species names
raw$Species <- trimws(raw$Species)
raw$Species <- gsub(" ", "_", raw$Species)


# Check duplicate species
duplicate_species <- raw$Species[
  duplicated(raw$Species)
]


if (length(duplicate_species) > 0) {
  
  cat("\nDuplicated species detected:\n")
  print(unique(duplicate_species))
  
  stop(
    "\nResolve duplicated species before fitting the ZIBB models."
  )
}



# 6. BUILD THE ZIBB ANALYSIS DATASET


dat <- raw %>%
  transmute(
    
    Species = Species,
    Family = Family,
    
    # --------------------------------------
    # COUNTS
    # --------------------------------------
    
    NeoplasiaCases = NeoplasiaWithDenominators,
    
    MalignancyCases = Malignant,
    
    # denominator / sampling effort
    Trials = RecordsWithDenominators,
    
    
    # --------------------------------------
    # LIFE HISTORY TRAITS
    # --------------------------------------
    
    adult_weight_g = adult_weight.g.,
    
    max_longevity_months = max_longevity.months.,
    
    gestation_months = Gestation.months.,
    
    litter_size = litter_size
  )



# 7. AUDIT COUNTS AND DENOMINATORS

# VERY IMPORTANT 
# ZIBB requires counts !!!!
#
# Trials and cases should therefore be integer valued !!!

is_integerish <- function(x) {
  
  all(
    is.na(x) |
      abs(x - round(x)) < 1e-8
  )
}


if (!is_integerish(dat$Trials)) {
  stop("RecordsWithDenominators contains non-integer values.")
}

if (!is_integerish(dat$NeoplasiaCases)) {
  stop("NeoplasiaWithDenominators contains non-integer values.")
}

if (!is_integerish(dat$MalignancyCases)) {
  stop("Malignant contains non-integer values.")
}


# Impossible count checks:
#
# cases cannot exceed trials ! 

bad_neoplasia <- dat %>%
  filter(
    !is.na(NeoplasiaCases),
    !is.na(Trials),
    NeoplasiaCases > Trials
  )


bad_malignancy <- dat %>%
  filter(
    !is.na(MalignancyCases),
    !is.na(Trials),
    MalignancyCases > Trials
  )


if (nrow(bad_neoplasia) > 0) {
  
  print(bad_neoplasia)
  
  stop(
    "At least one species has NeoplasiaCases > Trials. ",
    "Do not fit the model until this is resolved."
  )
}


if (nrow(bad_malignancy) > 0) {
  
  print(bad_malignancy)
  
  stop(
    "At least one species has MalignancyCases > Trials. ",
    "Do not fit the model until this is resolved."
  )
}



# 8. LOAD AND SANITIZE PHYLOGENY


if (!file.exists(tree_file)) {
  stop("Cannot find phylogenetic tree: ", tree_file)
}


tree <- ape::read.tree(tree_file)

tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- gsub(" ", "_", tree$tip.label)


# Match Squamata species to phylogeny

species_overlap <- intersect(
  dat$Species,
  tree$tip.label
)


message(
  "Species in Squamata data: ",
  length(unique(dat$Species))
)

message(
  "Species found in phylogeny: ",
  length(species_overlap)
)


not_in_tree <- setdiff(
  dat$Species,
  tree$tip.label
)


if (length(not_in_tree) > 0) {
  
  write.csv(
    data.frame(Species = not_in_tree),
    file.path(
      out_dir,
      "species_missing_from_tree.csv"
    ),
    row.names = FALSE
  )
  
  message(
    length(not_in_tree),
    " species were not found in the tree. ",
    "See species_missing_from_tree.csv"
  )
}


# Keep only species represented in phylogeny
dat <- dat %>%
  filter(Species %in% species_overlap)



# 9. TRANSFORM LIFE-HISTORY VARS


# The existing Squamata analyses frequently use log-transformed life history
# variables. Here we are going to log10-transform positive traits and then standardize them !
#
# So standardization means:
#
#     mean = 0
#     SD   = 1
#
# Thus a coefficient corresponds to roughly a 1-SD change in that predictor.


z_scale <- function(x) {
  
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  
  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  
  (x - m) / s
}


dat <- dat %>%
  mutate(
    
    # --------------------------------------
    # LOG TRANSFORMS
    # --------------------------------------
    
    log_mass = ifelse(
      adult_weight_g > 0,
      log10(adult_weight_g),
      NA_real_
    ),
    
    log_longevity = ifelse(
      max_longevity_months > 0,
      log10(max_longevity_months),
      NA_real_
    ),
    
    log_gestation = ifelse(
      gestation_months > 0,
      log10(gestation_months),
      NA_real_
    ),
    
    log_litter_size = ifelse(
      litter_size > 0,
      log10(litter_size),
      NA_real_
    ),
    
    # --------------------------------------
    # STANDARDIZED PREDICTORS
    # --------------------------------------
    
    mass_s = z_scale(log_mass),
    
    longevity_s = z_scale(log_longevity),
    
    gestation_s = z_scale(log_gestation),
    
    litter_size_s = z_scale(log_litter_size),
    
    
    # --------------------------------------
    # ZERO INFLATION PREDICTOR !!!!
    # --------------------------------------
    
    log_trials = ifelse(
      Trials > 0,
      log(Trials),
      NA_real_
    ),
    
    log_trials_s = z_scale(log_trials)
  )



# 10. COVERAGE REPORT


coverage <- tibble(
  
  Variable = c(
    "Trials",
    "NeoplasiaCases",
    "MalignancyCases",
    "Adult body mass",
    "Maximum longevity",
    "Gestation",
    "Litter size"
  ),
  
  N_nonmissing = c(
    sum(!is.na(dat$Trials) & dat$Trials > 0),
    sum(!is.na(dat$NeoplasiaCases)),
    sum(!is.na(dat$MalignancyCases)),
    sum(!is.na(dat$mass_s)),
    sum(!is.na(dat$longevity_s)),
    sum(!is.na(dat$gestation_s)),
    sum(!is.na(dat$litter_size_s))
  )
)


print(coverage)


write.csv(
  coverage,
  file.path(
    out_dir,
    "data_coverage.csv"
  ),
  row.names = FALSE
)


# Basic tumor-count audit

tumor_audit <- tibble(
  
  Outcome = c(
    "Neoplasia",
    "Malignancy"
  ),
  
  N_with_data = c(
    
    sum(
      !is.na(dat$NeoplasiaCases) &
        !is.na(dat$Trials) &
        dat$Trials > 0
    ),
    
    sum(
      !is.na(dat$MalignancyCases) &
        !is.na(dat$Trials) &
        dat$Trials > 0
    )
  ),
  
  N_zero_cases = c(
    
    sum(
      dat$NeoplasiaCases == 0,
      na.rm = TRUE
    ),
    
    sum(
      dat$MalignancyCases == 0,
      na.rm = TRUE
    )
  )
)


print(tumor_audit)


write.csv(
  tumor_audit,
  file.path(
    out_dir,
    "tumor_count_audit.csv"
  ),
  row.names = FALSE
)


# 11. PRIORS


# Predictors have been standardized, so now Normal(0,1) is a reasonable
# regularizing prior for regression coefficients
#
# phi controls beta-binomial dispersion

make_priors <- function(zi_mode) {
  
  p <- c(
    
    # Mean model
    brms::set_prior(
      "normal(0, 1.5)",
      class = "Intercept"
    ),
    
    brms::set_prior(
      "normal(0, 1)",
      class = "b"
    ),
    
    # Phylogenetic random effect
    brms::set_prior(
      "normal(0, 1)",
      class = "sd"
    ),
    
    # Beta binomial dispersion
    brms::set_prior(
      "exponential(1)",
      class = "phi"
    ),
    
    # Zero inflation intercept
    brms::set_prior(
      "normal(0, 1.5)",
      class = "Intercept",
      dpar = "zi"
    )
  )
  
  # IMPORTANT !!!!
  # Only add a zi regression prior if the zi model has a predictor !!!
  if (zi_mode == "log_trials") {
    
    p <- c(
      p,
      
      brms::set_prior(
        "normal(0, 1)",
        class = "b",
        dpar = "zi"
      )
    )
  }
  
  
  p
}



# 12. EFFECT PLOT FUNCTION


# This creates our:
#
#   observed prevalence points
#   posterior median prediction
#   95% Bayesian credible interval
#
# IMPORTANT!!!

# posterior_epred() returns predicted tumor COUNTS
#
# We divide predictions by Trials so that the plot is truly showing
# predicted PREVALENCE rather than predicted counts ! ! ! 


predictor_labels <- c(
  
  longevity_s =
    "Standardized log10 Maximum Longevity",
  
  gestation_s =
    "Standardized log10 Gestation Length",
  
  mass_s =
    "Standardized log10 Adult Body Mass",
  
  litter_size_s =
    "Standardized log10 Litter Size"
)


export_effect_plot <- function(
    fit,
    df_model,
    focal_predictor,
    all_predictors,
    response_label,
    model_label,
    out_dir) {
  
  
  if (!focal_predictor %in% names(df_model)) {
    return(invisible(NULL))
  }
  
  
  # Prediction range = range actually observed in the model
  x_seq <- seq(
    min(df_model[[focal_predictor]], na.rm = TRUE),
    max(df_model[[focal_predictor]], na.rm = TRUE),
    length.out = 150
  )
  
  
  # Median sampling effort for visualization
  median_trials <- round(
    median(
      df_model$Trials,
      na.rm = TRUE
    )
  )
  
  
  # Start prediction dataframe
  newdat <- tibble(
    
    Trials = rep(
      median_trials,
      length(x_seq)
    ),
    
    # zero = average standardized sampling effort
    log_trials_s = rep(
      0,
      length(x_seq)
    ),
    
    Species = factor(
      rep(
        levels(df_model$Species)[1],
        length(x_seq)
      ),
      levels = levels(df_model$Species)
    )
  )
  
  
  # Hold all standardized predictors at their mean (= 0)
  for (pv in all_predictors) {
    
    newdat[[pv]] <- 0
  }
  
  
  # Vary focal predictor
  newdat[[focal_predictor]] <- x_seq
  
  
  # Expected posterior counts
  epred_counts <- brms::posterior_epred(
    
    fit,
    
    newdata = newdat,
    
    re_formula = NA,
    
    allow_new_levels = TRUE
  )
  
  
  # Convert expected counts -> expected prevalence
  epred_prev <- sweep(
    
    epred_counts,
    
    MARGIN = 2,
    
    STATS = newdat$Trials,
    
    FUN = "/"
  )
  
  
  prediction_df <- tibble(
    
    x = x_seq,
    
    Estimate = apply(
      epred_prev,
      2,
      median
    ),
    
    Lower = apply(
      epred_prev,
      2,
      quantile,
      probs = 0.025
    ),
    
    Upper = apply(
      epred_prev,
      2,
      quantile,
      probs = 0.975
    )
  )
  
  
  observed_df <- df_model %>%
    mutate(
      ObservedPrevalence =
        cases / Trials
    )
  
  
  x_label <- predictor_labels[[focal_predictor]]
  
  if (is.null(x_label)) {
    x_label <- focal_predictor
  }
  
  
  g <- ggplot() +
    
    # Observed prevalence
    geom_point(
      data = observed_df,
      aes(
        x = .data[[focal_predictor]],
        y = ObservedPrevalence,
        size = Trials
      ),
      color = "grey40",
      alpha = 0.55
    ) +
    
    # 95% credible interval
    geom_ribbon(
      data = prediction_df,
      aes(
        x = x,
        ymin = Lower,
        ymax = Upper
      ),
      fill = "grey80",
      alpha = 0.8
    ) +
    
    # Posterior median
    geom_line(
      data = prediction_df,
      aes(
        x = x,
        y = Estimate
      ),
      color = "black",
      linewidth = 1.2
    ) +
    
    scale_y_continuous(
      limits = c(0, NA),
      labels = scales::label_percent(
        accuracy = 1
      )
    ) +
    
    scale_size_continuous(
      name = "Records",
      range = c(1.5, 6)
    ) +
    
    labs(
      
      title = paste0(
        response_label,
        " — ",
        model_label
      ),
      
      subtitle =
        "Zero-Inflated Beta-Binomial phylogenetic GLMM",
      
      x = x_label,
      
      y = paste0(
        "Observed / predicted ",
        tolower(response_label),
        " prevalence"
      )
    ) +
    
    theme_classic(base_size = 13) +
    
    theme(
      legend.position = "right",
      plot.title =
        element_text(face = "bold")
    )
  
  
  filename <- paste0(
    "plot_",
    response_label,
    "_",
    model_label,
    "_",
    focal_predictor,
    ".png"
  )
  
  
  ggsave(
    filename = file.path(
      out_dir,
      filename
    ),
    plot = g,
    width = 8,
    height = 5.5,
    dpi = 300,
    bg = "white"
  )
  
  
  invisible(g)
}


# 13. ZERO-FREQUENCY POSTERIOR CHECK


# This checks whether the fitted model reproduces the observed proportion
# of species with zero cases 


export_zero_check <- function(
    fit,
    observed_cases,
    response_label,
    model_label,
    out_dir) {
  
  
  y_rep <- brms::posterior_predict(
    fit,
    ndraws = 500
  )
  
  
  predicted_zero_fraction <- rowMeans(
    y_rep == 0
  )
  
  
  observed_zero_fraction <- mean(
    observed_cases == 0
  )
  
  
  check_df <- tibble(
    PredictedZeroFraction =
      predicted_zero_fraction
  )
  
  
  g <- ggplot(
    check_df,
    aes(x = PredictedZeroFraction)
  ) +
    
    geom_density(
      fill = "#79A7D3",
      alpha = 0.55
    ) +
    
    geom_vline(
      xintercept =
        observed_zero_fraction,
      color = "red",
      linewidth = 1.2,
      linetype = 2
    ) +
    
    labs(
      
      title = paste0(
        response_label,
        " — zero-count posterior check"
      ),
      
      subtitle =
        "Red line = observed proportion of species with zero cases",
      
      x =
        "Proportion of species with zero tumor cases",
      
      y =
        "Posterior density"
    ) +
    
    theme_classic(base_size = 13)
  
  
  ggsave(
    
    file.path(
      out_dir,
      paste0(
        "zero_check_",
        response_label,
        "_",
        model_label,
        ".png"
      )
    ),
    
    g,
    
    width = 7,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}



# 14. MAIN ZIBB-PGLMM FITTING FUNCTION


fit_zibb_model <- function(
    response_col,
    response_label,
    predictor_vars,
    model_label,
    dat,
    full_tree,
    out_dir,
    zi_mode = "log_trials",
    min_species = 20) {
  
  
  cat(
    "\n\n============================================================\n"
  )
  
  cat(
    "MODEL:",
    model_label,
    "\n"
  )
  
  cat(
    "OUTCOME:",
    response_label,
    "\n"
  )
  
  cat(
    "============================================================\n"
  )
  
  

  # Prepare data for THIS MODEL ONLY
  #
  # This is very important ! ! ! 
  #
  # A longevity-only model should NOT lose a species just because
  # that species lacks gestation data.

  
  needed <- c(
    "Species",
    "Family",
    "Trials",
    "log_trials_s",
    response_col,
    predictor_vars
  )
  
  
  df_model <- dat %>%
    select(
      all_of(needed)
    )
  
  
  # Rename outcome to generic "cases"
  names(df_model)[
    names(df_model) == response_col
  ] <- "cases"
  
  
  # Drop missing values ONLY for variables required by this model
  df_model <- df_model %>%
    drop_na(
      all_of(
        c(
          "Species",
          "Trials",
          "log_trials_s",
          "cases",
          predictor_vars
        )
      )
    ) %>%
    
    filter(
      Trials > 0,
      cases >= 0,
      cases <= Trials
    )
  
  
  # --------------------------------------
  # Match to tree
  # --------------------------------------
  
  species_use <- intersect(
    full_tree$tip.label,
    df_model$Species
  )
  
  
  tree_model <- ape::keep.tip(
    full_tree,
    species_use
  )
  
  
  df_model <- df_model %>%
    filter(
      Species %in%
        tree_model$tip.label
    )
  
  
  # Reorder dataset EXACTLY to tree order
  df_model <- df_model[
    match(
      tree_model$tip.label,
      df_model$Species
    ),
  ]
  
  
  if (anyNA(df_model$Species)) {
    stop("Species/tree ordering failed.")
  }
  
  
  # Factor levels must match tree
  df_model$Species <- factor(
    df_model$Species,
    levels =
      tree_model$tip.label
  )
  
  
  N_species <- nrow(
    df_model
  )
  
  
  zero_prop <- mean(
    df_model$cases == 0
  )
  
  
  median_trials <- median(
    df_model$Trials
  )
  
  
  cat(
    "Species in model:",
    N_species,
    "\n"
  )
  
  cat(
    "Zero-case species:",
    sum(df_model$cases == 0),
    "/",
    N_species,
    "(",
    round(zero_prop * 100, 1),
    "%)\n"
  )
  
  cat(
    "Median records per species:",
    median_trials,
    "\n"
  )
  
  
  # Don't fit very small complex ZIBB models automatically
  if (N_species < min_species) {
    
    warning(
      response_label,
      " / ",
      model_label,
      " skipped because N = ",
      N_species,
      " < MIN_SPECIES = ",
      min_species,
      "."
    )
    
    
    return(
      list(
        fitted = FALSE,
        N_species = N_species,
        reason = "Too few species"
      )
    )
  }
  
  
  # Save exact data used for reproducibility
  write.csv(
    
    df_model,
    
    file.path(
      out_dir,
      paste0(
        "data_used_",
        response_label,
        "_",
        model_label,
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # --------------------------------------
  # Build phylogenetic correlation matrix ! 
  # --------------------------------------
  
  A_model <- ape::vcv(
    tree_model,
    corr = TRUE
  )
  
  
  A_model <- A_model[
    tree_model$tip.label,
    tree_model$tip.label
  ]
  
  
  stopifnot(
    identical(
      rownames(A_model),
      tree_model$tip.label
    )
  )
  
  
  stopifnot(
    identical(
      colnames(A_model),
      tree_model$tip.label
    )
  )
  
  
  # --------------------------------------
  # Mean formula
  # --------------------------------------
  
  predictor_string <- paste(
    predictor_vars,
    collapse = " + "
  )
  
  
  mean_formula <- as.formula(
    
    paste0(
      
      "cases | trials(Trials) ~ ",
      
      predictor_string,
      
      " + (1 | gr(Species, cov = A_model))"
    )
  )
  
  
  # --------------------------------------
  # Zero-inflation formula
  # --------------------------------------
  
  if (zi_mode == "log_trials") {
    
    zi_formula <-
      as.formula(
        "zi ~ log_trials_s"
      )
    
  } else {
    
    zi_formula <-
      as.formula(
        "zi ~ 1"
      )
  }
  
  
  model_formula <- brms::bf(
    mean_formula,
    zi_formula
  )
  
  # --------------------------------------
  # Priors
  # --------------------------------------
  
  priors <- make_priors(
    zi_mode
  )
  
  # --------------------------------------
  # FIT TRUE ZIBB
  # --------------------------------------
  
  fit <- brms::brm(
    
    formula = model_formula,
    
    family =
      brms::zero_inflated_beta_binomial(),
    
    data = df_model,
    
    data2 = list(
      A_model = A_model
    ),
    
    prior = priors,
    
    chains = N_CHAINS,
    
    cores = N_CORES,
    
    iter = N_ITER,
    
    warmup = N_WARMUP,
    
    seed = 1234,
    
    control = list(
      adapt_delta =
        ADAPT_DELTA,
      max_treedepth =
        MAX_TREEDEPTH
    ),
    
    save_pars =
      brms::save_pars(
        all = TRUE
      ),
    
    refresh = 500
  )
  
  # --------------------------------------
  # Save fitted model
  # --------------------------------------
  
  model_file <- file.path(
    
    out_dir,
    
    paste0(
      "fit_",
      response_label,
      "_",
      model_label,
      ".rds"
    )
  )
  
  
  saveRDS(
    fit,
    model_file
  )
  
  # --------------------------------------
  # Save summary ! 
  # --------------------------------------
  
  summary_file <- file.path(
    
    out_dir,
    
    paste0(
      "summary_",
      response_label,
      "_",
      model_label,
      ".txt"
    )
  )
  
  
  writeLines(
    
    capture.output(
      print(
        summary(fit),
        digits = 3
      )
    ),
    
    summary_file
  )
  
  # --------------------------------------
  # MCMC Diagnostics
  # --------------------------------------
  
  diag_df <- posterior::summarise_draws(
    
    posterior::as_draws_array(
      fit
    ),
    
    "mean",
    "sd",
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
  
  
  write.csv(
    
    diag_df,
    
    file.path(
      out_dir,
      paste0(
        "diagnostics_",
        response_label,
        "_",
        model_label,
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  max_rhat <- max(
    diag_df$rhat,
    na.rm = TRUE
  )
  
  # --------------------------------------
  # Divergence count
  # --------------------------------------
  
  nuts <- tryCatch(
    
    brms::nuts_params(
      fit
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  
  if (!is.null(nuts)) {
    
    N_divergent <- sum(
      nuts$Parameter ==
        "divergent__" &
        nuts$Value == 1
    )
    
  } else {
    
    N_divergent <- NA_integer_
  }
  
  
  cat(
    "Maximum Rhat:",
    round(max_rhat, 4),
    "\n"
  )
  
  cat(
    "Divergent transitions:",
    N_divergent,
    "\n"
  )
  
  
  # --------------------------------------
  # Posterior draws
  # --------------------------------------
  
  post <-
    posterior::as_draws_df(
      fit
    )
  
  
  # --------------------------------------
  # Beta-binomial dispersion parameter phi
  # --------------------------------------
  
  if ("phi" %in% names(post)) {
    
    phi_q <- quantile(
      post$phi,
      probs = c(
        0.025,
        0.50,
        0.975
      )
    )
    
  } else {
    
    phi_q <- c(
      NA,
      NA,
      NA
    )
  }
  
  
  # --------------------------------------
  # Extract predictor effects
  # --------------------------------------
  
  result_rows <- list()
  
  
  for (pv in predictor_vars) {
    
    
    parameter_name <-
      paste0(
        "b_",
        pv
      )
    
    
    if (!parameter_name %in%
        names(post)) {
      
      next
    }
    
    
    beta <-
      post[[parameter_name]]
    
    
    OR <-
      exp(beta)
    
    
    beta_q <- quantile(
      beta,
      probs = c(
        0.025,
        0.50,
        0.975
      )
    )
    
    
    OR_q <- quantile(
      OR,
      probs = c(
        0.025,
        0.50,
        0.975
      )
    )
    
    
    result_rows[[pv]] <-
      tibble(
        
        Response =
          response_label,
        
        Model =
          model_label,
        
        Predictor =
          pv,
        
        N_species =
          N_species,
        
        Zero_fraction =
          zero_prop,
        
        Median_trials =
          median_trials,
        
        Beta_lower95 =
          beta_q[1],
        
        Beta_median =
          beta_q[2],
        
        Beta_upper95 =
          beta_q[3],
        
        OddsRatio_lower95 =
          OR_q[1],
        
        OddsRatio_median =
          OR_q[2],
        
        OddsRatio_upper95 =
          OR_q[3],
        
        Posterior_P_positive =
          mean(beta > 0),
        
        Posterior_P_negative =
          mean(beta < 0),
        
        CrI_excludes_zero =
          beta_q[1] > 0 |
          beta_q[3] < 0,
        
        Phi_lower95 =
          phi_q[1],
        
        Phi_median =
          phi_q[2],
        
        Phi_upper95 =
          phi_q[3],
        
        Max_Rhat =
          max_rhat,
        
        Divergences =
          N_divergent
      )
  }
  
  
  result_df <-
    bind_rows(
      result_rows
    )
  
  
  # --------------------------------------
  # Zero-inflation sampling effort effect
  # --------------------------------------
  
  zi_parameter <-
    "b_zi_log_trials_s"
  
  
  if (
    zi_mode == "log_trials" &&
    zi_parameter %in%
    names(post)
  ) {
    
    
    zi_beta <-
      post[[zi_parameter]]
    
    
    zi_q <-
      quantile(
        zi_beta,
        probs = c(
          0.025,
          0.50,
          0.975
        )
      )
    
    
    cat(
      "\nZero-inflation sampling-effort coefficient:\n"
    )
    
    print(
      zi_q
    )
  }
  
  
  # --------------------------------------
  # Effect plots
  # --------------------------------------
  
  for (pv in predictor_vars) {
    
    export_effect_plot(
      
      fit = fit,
      
      df_model =
        df_model,
      
      focal_predictor =
        pv,
      
      all_predictors =
        predictor_vars,
      
      response_label =
        response_label,
      
      model_label =
        model_label,
      
      out_dir =
        out_dir
    )
  }
  
  
  # --------------------------------------
  # Zero-frequency posterior predictive check
  # --------------------------------------
  
  export_zero_check(
    
    fit = fit,
    
    observed_cases =
      df_model$cases,
    
    response_label =
      response_label,
    
    model_label =
      model_label,
    
    out_dir =
      out_dir
  )
  
  
  # --------------------------------------
  # Standard posterior predictive check
  # --------------------------------------
  
  pp_plot <- brms::pp_check(
    
    fit,
    
    type = "stat",
    
    stat = "mean",
    
    ndraws = 100
  )
  
  
  ggsave(
    
    file.path(
      out_dir,
      paste0(
        "ppcheck_",
        response_label,
        "_",
        model_label,
        ".png"
      )
    ),
    
    pp_plot,
    
    width = 7,
    height = 5,
    
    dpi = 300,
    
    bg = "white"
  )
  
  
  # ---------
  # Return
  # ---------
  
  list(
    
    fitted =
      TRUE,
    
    fit =
      fit,
    
    results =
      result_df,
    
    N_species =
      N_species,
    
    zero_fraction =
      zero_prop,
    
    divergences =
      N_divergent
  )
}


# ===================
# 15. DEFINE MODELS
# ===================

# We are STARTING with univariate models because they preserve more species
#
# Then we examine biologically relevant multivariate combinations
#
# We are NOT doing automatic stepwise selection here ! 


model_specs <- list(
  
  
  # ------------
  # UNIVARIATE
  # ------------
  
  list(
    label =
      "longevity_only",
    predictors =
      c("longevity_s")
  ),
  
  
  list(
    label =
      "gestation_only",
    predictors =
      c("gestation_s")
  ),
  
  
  list(
    label =
      "mass_only",
    predictors =
      c("mass_s")
  ),
  
  
  list(
    label =
      "litter_size_only",
    predictors =
      c("litter_size_s")
  ),
  
  
  # ---------------
  # MULTIVARIATE
  # ---------------
  
  list(
    label =
      "longevity_gestation",
    predictors =
      c(
        "longevity_s",
        "gestation_s"
      )
  ),
  
  
  list(
    label =
      "longevity_gestation_mass",
    predictors =
      c(
        "longevity_s",
        "gestation_s",
        "mass_s"
      )
  ),
  
  
  list(
    label =
      "full_life_history",
    predictors =
      c(
        "longevity_s",
        "gestation_s",
        "mass_s",
        "litter_size_s"
      )
  )
)


# ======================
# 16. DEFINE OUR RESPONSES
# ======================

responses <- list(
  
  Neoplasia = "NeoplasiaCases",
  
  Malignancy = "MalignancyCases"
)


# =================
# 17. RUN MODELS
# =================

results_list <- list()


for (
  response_label
  in names(responses)
) {
  
  
  response_col <-
    responses[[response_label]]
  
  
  for (
    spec
    in model_specs
  ) {
    
    
    model_label <-
      spec$label
    
    
    predictor_vars <-
      spec$predictors
    
    
    result <-
      fit_zibb_model(
        
        response_col =
          response_col,
        
        response_label =
          response_label,
        
        predictor_vars =
          predictor_vars,
        
        model_label =
          model_label,
        
        dat =
          dat,
        
        full_tree =
          tree,
        
        out_dir =
          out_dir,
        
        zi_mode =
          ZI_MODE,
        
        min_species =
          MIN_SPECIES
      )
    
    
    if (
      isTRUE(
        result$fitted
      )
    ) {
      
      
      key <-
        paste(
          response_label,
          model_label,
          sep = "_"
        )
      
      
      results_list[[key]] <-
        result$results
    }
  }
}


# ==========================
# 18. COMBINE MODEL RESULTS
# ==========================

if (
  length(results_list) > 0
) {
  
  
  results_all <-
    bind_rows(
      results_list
    )
  
  
  print(
    results_all,
    n = Inf
  )
  
  
  write.csv(
    
    results_all,
    
    file.path(
      out_dir,
      "ZIBB_PGLMM_Squamata_summary.csv"
    ),
    
    row.names = FALSE
  )
  
  
} else {
  
  
  warning(
    "No models were fitted. ",
    "Check missing data and MIN_SPECIES."
  )
}


# ====================================
# 19. IMPORTANT MODEL COMPARISON NOTE
# ====================================

cat(
  "\n\nIMPORTANT MODEL COMPARISON NOTE:\n"
)

cat(
  "Do NOT directly compare LOO/AIC values between models that contain different species.\n"
)

cat(
  "For example, a longevity-only model may use far more species than a\n"
)

cat(
  "longevity + gestation + mass model.\n"
)

cat(
  "Formal model comparison should use the SAME species in every compared model.\n\n"
)


# ==============================
# 20. SAVE SESSION INFORMATION
# ==============================

writeLines(
  
  capture.output(
    sessionInfo()
  ),
  
  file.path(
    out_dir,
    "sessionInfo.txt"
  )
)


cat(
  "\n============================================================\n"
)

cat(
  "ZIBB-PGLMM PIPELINE COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "Results saved to:\n"
)

cat(
  normalizePath(
    out_dir
  ),
  "\n"
)