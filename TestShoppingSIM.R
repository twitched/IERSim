# Load the MASS package for mvrnorm
# install.packages("MASS")
library(MASS)
library(tidyverse)

# --------------------------------------------------------------------------
# 1. Define Model Parameters and Inputs
# --------------------------------------------------------------------------

# Number of observations to simulate
N_observations <- 500 # You can change this

# --- Latent Variable (Factor) Parameters (Eta: η) ---
factor_names <- c("SSR", "RI", "WOM", "SAT", "DJ", "IJ", "PJ")
K_factors <- length(factor_names)

# Latent factor covariance matrix (Psi: Ψ) - Provided by user
Psi_matrix <- matrix(c(
  1.020100, 0.478538, 0.650440, 0.556510, 0.892032, 0.705485, 0.555399,
  0.478538, 1.060900, 0.454848, 0.970775, 0.710700, 0.614807, 0.530244,
  0.650440, 0.454848, 0.846400, 0.453560, 0.672888, 0.665988, 0.462852,
  0.556510, 0.970775, 0.453560, 2.102500, 1.020510, 0.957580, 1.000935,
  0.892032, 0.710700, 0.672888, 1.020510, 1.904400, 0.998982, 0.888030,
  0.705485, 0.614807, 0.665988, 0.957580, 0.998982, 1.612900, 0.817245,
  0.555399, 0.530244, 0.462852, 1.000935, 0.888030, 0.817245, 1.368900
), nrow = K_factors, byrow = TRUE, dimnames = list(factor_names, factor_names))

# Latent factor mean vector (Beta: β) - Provided by user
Beta_vector <- c(SSR = 5.49, RI = 3.25, WOM = 4.03, SAT = 4.27, DJ = 3.44, IJ = 3.89, PJ = 3.93)

# Target Reliabilities for each factor's measurement (interpreted as McDonald's Omega)
target_reliabilities <- c(
  SSR = 0.88, RI = 0.88, WOM = 0.88,
  SAT = 0.90, DJ = 0.89, IJ = 0.89, PJ = 0.87
)

# --- Observed Variable (Indicator) Parameters (Y, Epsilon: ε) ---
indicator_definition <- list(
  SSR = c("ssr1", "ssr2", "ssr3"),
  RI  = c("ri1", "ri2", "ri3"),
  WOM = c("wom1", "wom2", "wom3"),
  SAT = c("sat1", "sat2", "sat3"),
  DJ  = c("dj1", "dj2", "dj3", "dj4"),
  IJ  = c("ij1", "ij2", "ij3", "ij4"),
  PJ  = c("pj1", "pj2", "pj3", "pj4")
)
indicator_names <- unlist(indicator_definition)
I_indicators <- length(indicator_names)

# Factor Loading Matrix (Lambda: Λ)
# I_indicators x K_factors matrix (24x7)
Lambda_matrix <- matrix(0, nrow = I_indicators, ncol = K_factors,
                        dimnames = list(indicator_names, factor_names))

# ** CRITICAL ASSUMPTION: Define factor loadings here **
# We assume the first indicator for each factor has loading 1.0, others 0.7.
# You can change 0.7 to another value or specify a unique pattern.
assumed_secondary_loading <- 0.7

for (factor in factor_names) {
  inds <- indicator_definition[[factor]]
  Lambda_matrix[inds[1], factor] <- 1.0 # First indicator loading fixed to 1
  if (length(inds) > 1) {
    Lambda_matrix[inds[2:length(inds)], factor] <- assumed_secondary_loading
  }
}

# --- Calculate Residual Variances (diagonal of Theta: Θ) based on reliabilities ---
residual_variances <- numeric(I_indicators)
names(residual_variances) <- indicator_names

for (factor_idx in 1:K_factors) {
  factor_name <- factor_names[factor_idx]
  factor_var <- Psi_matrix[factor_name, factor_name]
  rel_k <- target_reliabilities[factor_name]
  
  # Get indicators and their loadings for the current factor
  current_indicators <- indicator_definition[[factor_name]]
  current_loadings <- Lambda_matrix[current_indicators, factor_name]
  
  sum_lambda_k <- sum(current_loadings)
  L_k <- (sum_lambda_k^2) * factor_var
  
  # Sum of residual variances for indicators of this factor
  E_k <- (L_k * (1 - rel_k)) / rel_k
  
  if (E_k < 0) {
    stop(paste("Calculated sum of residual variances for factor", factor_name, "is negative.",
               "This might indicate inconsistent reliability targets or loading assumptions."))
  }
  
  # Allocate the sum E_k equally among the indicators of this factor
  num_inds_for_factor <- length(current_indicators)
  individual_residual_variance <- E_k / num_inds_for_factor
  
  if (individual_residual_variance < 0) {
    stop(paste("Calculated individual residual variance for indicators of factor", factor_name, "is negative."))
  }
  
  residual_variances[current_indicators] <- individual_residual_variance
}

Theta_matrix <- diag(residual_variances)
dimnames(Theta_matrix) <- list(indicator_names, indicator_names)

# Mean vector for residuals (typically all zeros)
epsilon_means <- rep(0, I_indicators)

# NEW: Define target mean for observed indicators
target_observed_mean <- 4.0 # Example: aiming for a mean of 4 for all indicators

# NEW: Calculate indicator intercepts (tau_vector)
# tau_vector will be an I_indicators x 1 vector
tau_vector <- numeric(I_indicators)
names(tau_vector) <- indicator_names

for (factor_k_name in factor_names) {
  # Get the mean of the current latent factor
  beta_k <- Beta_vector[factor_k_name]
  
  # Get indicators for this factor
  indicators_for_factor_k <- indicator_definition[[factor_k_name]]
  
  for (indicator_i_name in indicators_for_factor_k) {
    # Get the loading of indicator i on factor k
    lambda_ik <- Lambda_matrix[indicator_i_name, factor_k_name]
    
    # Calculate tau_i
    tau_vector[indicator_i_name] <- target_observed_mean - (lambda_ik * beta_k)
  }
}

print("Calculated Indicator Intercepts (Tau Vector):")
print(tau_vector)

# --------------------------------------------------------------------------
# 2. Simulate Latent Factor Scores (η)
# --------------------------------------------------------------------------
set.seed(123) # for reproducibility

# Eta (η) will be an N_observations x K_factors matrix
eta_scores <- mvrnorm(n = N_observations, mu = Beta_vector, Sigma = Psi_matrix)

# --------------------------------------------------------------------------
# 3. Simulate Residuals (ε)
# --------------------------------------------------------------------------
# Epsilon (ε) will be an N_observations x I_indicators matrix
epsilon_scores <- mvrnorm(n = N_observations, mu = epsilon_means, Sigma = Theta_matrix)

# --------------------------------------------------------------------------
# 4. Calculate Observed Indicator Scores (Y)
# --------------------------------------------------------------------------
# Y = 1*τ' + η*Λ' + ε
# We need to add the tau_vector to each row of (eta_scores %*% t(Lambda_matrix) + epsilon_scores)
# Or, more directly: Y_ij = tau_i + (eta_j %*% Lambda_i_T) + epsilon_ij

# Create an N_observations x I_indicators matrix of intercepts
# Each row is tau_vector, repeated N_observations times
tau_matrix_repeated <- matrix(tau_vector, nrow = N_observations, ncol = I_indicators, byrow = TRUE)

Y_observed_scores <- tau_matrix_repeated + eta_scores %*% t(Lambda_matrix) + epsilon_scores

simulated_data_matrix_algebra <- as.data.frame(Y_observed_scores)
colnames(simulated_data_matrix_algebra) <- indicator_names

# --------------------------------------------------------------------------
# 5. Inspect the Simulated Data (Optional)
# --------------------------------------------------------------------------
print("First few rows of simulated data:")
head(simulated_data_matrix_algebra)

print("Calculated individual residual variances (Theta diagonal):")
print(diag(Theta_matrix))


print("Sample Means of Observed Variables:")
colMeans(simulated_data_matrix_algebra)

# (Assuming the previous script has been run and the following objects exist:
# simulated_data_matrix_algebra, Psi_matrix, Beta_vector, factor_names)

# --------------------------------------------------------------------------
# 6. Define and Fit the Model using lavaan
# --------------------------------------------------------------------------
# Install and load lavaan if you haven't already
# install.packages("lavaan")
library(lavaan)

# Define the lavaan model for fitting.
# For identification, lavaan will automatically fix the loading of the first
# indicator for each latent variable to 1.0. All other parameters
# (other loadings, residual variances, latent variances, latent covariances,
# and latent means) will be estimated from the data.
lavaan_model_syntax <- "
 # Measurement model
  SSR =~ ssr1 + ssr2 + ssr3
  RI  =~ ri1 + ri2 + ri3
  WOM =~ wom1 + wom2 + wom3
  SAT =~ sat1 + sat2 + sat3
  DJ  =~ dj1 + dj2 + dj3 + dj4
  IJ  =~ ij1 + ij2 + ij3 + ij4
  PJ  =~ pj1 + pj2 + pj3 + pj4

  # Variances
  SSR ~~ SSR
  RI ~~ RI
  WOM ~~ WOM
  SAT ~~ SAT
  DJ ~~ DJ
  IJ ~~ IJ
  PJ ~~ PJ
  
  # Covariances
  SSR ~~ RI  
  SSR ~~ WOM 
  SSR ~~ SAT 
  SSR ~~ DJ  
  SSR ~~ IJ  
  SSR ~~ PJ  
  RI ~~ WOM 
  RI ~~ SAT 
  RI ~~ DJ  
  RI ~~ IJ  
  RI ~~ PJ  
  WOM ~~ SAT 
  WOM ~~ DJ  
  WOM ~~ IJ  
  WOM ~~ PJ  
  SAT ~~ DJ  
  SAT ~~ IJ  
  SAT ~~ PJ  
  DJ ~~ IJ  
  DJ ~~ PJ  
  IJ ~~ PJ
  
  #intercepts
  SSR ~ 1
  RI  ~ 1
  WOM ~ 1
  SAT ~ 1
  DJ  ~ 1
  IJ  ~ 1
  PJ  ~ 1
  
  # Means
  ssr1 ~ ssr_int * 1
  ssr2 ~ ssr_int * 1
  ssr3 ~ ssr_int * 1 
  ri1 ~ ri_int * 1
  ri2 ~ ri_int * 1
  ri3 ~ ri_int * 1
  wom1 ~ wom_int * 1
  wom2 ~ wom_int * 1
  wom3 ~ wom_int * 1
  sat1 ~ sat_int * 1
  sat2 ~ sat_int * 1
  sat3 ~ sat_int * 1
  dj1 ~ dj_int * 1
  dj2 ~ dj_int * 1
  dj3 ~ dj_int * 1
  dj4 ~ dj_int * 1
  ij1 ~ ij_int * 1
  ij2 ~ ij_int * 1
  ij3 ~ ij_int * 1
  ij4 ~ ij_int * 1
  pj1 ~ pj_int * 1
  pj2 ~ pj_int * 1
  pj3 ~ pj_int * 1
  pj4 ~ pj_int * 1
"

# Fit the CFA model to the simulated data
# lavaan's cfa() function will automatically include a mean structure.
fit_lavaan <- cfa(lavaan_model_syntax, data = simulated_data_matrix_algebra)

# Check model summary (optional, but good for seeing overall fit and parameter estimates)
summary(fit_lavaan, fit.measures = TRUE, standardized = TRUE)

#Function to plot histograms of observed variables
histograms <- function(data) {
  data <- data |> 
    pivot_longer(cols = everything())
  
  #ensure measure names retain their order
  data$name <- as.character(data$name)
  data$name <- factor(data$name, levels=unique(data$name))
  
  #plot
  ggplot(data, aes(value)) +
    geom_histogram(binwidth = .5) +
    scale_x_continuous(breaks = seq(1, 7, 1), limits = c(0.5, 7.5)) +
    facet_wrap(vars(name))
}

# Plot histograms of observed variables
histograms(simulated_data_matrix_algebra)

# --------------------------------------------------------------------------
# 7. Compare Given (Population) vs. Lavaan Fitted Latent Parameters
# --------------------------------------------------------------------------

# --- Latent Variable Covariance Matrix ---
# Given (Population) Latent Covariance Matrix (Psi_matrix from data generation)
print("Given (Population) Latent Variable Covariance Matrix (Psi):")
print(Psi_matrix)

# Lavaan Fitted Latent Covariance Matrix
fitted_cov_lv <- lavInspect(fit_lavaan, "est")$psi # More direct way for CFA
# Or, if it was a general SEM model: fitted_cov_lv <- lavInspect(fit_lavaan, "cov.lv")
# For CFA, psi contains the symmetric variance-covariance matrix of the factors
dimnames(fitted_cov_lv) <- list(factor_names, factor_names) # Ensure names for comparison
print("Lavaan Fitted Latent Variable Covariance Matrix (Estimated Psi):")
print(fitted_cov_lv)

print("Difference (Population - Fitted) for Latent Covariance Matrix:")
print(Psi_matrix - fitted_cov_lv)


# --- Latent Variable Means ---
# Given (Population) Latent Mean Vector (Beta_vector from data generation)
print("Given (Population) Latent Variable Mean Vector (Beta):")
print(Beta_vector)

# Lavaan Fitted Latent Mean Vector
fitted_mean_lv <- lavInspect(fit_lavaan, "est")$alpha # More direct way for CFA
#fitted_mean_lv <- lavInspect(fit_lavaan, "mean.lv")
# For CFA, alpha contains the latent means
names(fitted_mean_lv) <- factor_names # Ensure names for comparison
print("Lavaan Fitted Latent Variable Mean Vector (Estimated Beta):")
print(fitted_mean_lv)

print("Difference (Population - Fitted) for Latent Means:")
print(Beta_vector - fitted_mean_lv)

# --- Additionally, compare observed variable covariance matrix and means ---
# Population observed covariance matrix (model-implied based on generation parameters)
# Lambda_matrix, Psi_matrix, Theta_matrix should be from the data generation script
if (exists("Lambda_matrix") && exists("Theta_matrix")) {
  population_observed_cov <- Lambda_matrix %*% Psi_matrix %*% t(Lambda_matrix) + Theta_matrix
  print("Population Model-Implied Observed Variable Covariance Matrix (first 6x6):")
  print(round(population_observed_cov[1:6, 1:6], 4))
}

# Lavaan fitted observed covariance matrix (model-implied by the fitted model)
fitted_observed_cov <- fitted.values(fit_lavaan)$cov
if (!is.null(fitted_observed_cov)) {
  print("Lavaan Model-Implied (Fitted) Observed Variable Covariance Matrix (first 6x6):")
  print(round(fitted_observed_cov[1:6, 1:6], 4))
} else {
  print("Could not retrieve fitted observed covariance matrix from lavaan object.")
}

# Population observed means (model-implied based on generation parameters)
if (exists("Lambda_matrix")) {
  population_observed_means <- Lambda_matrix %*% Beta_vector
  # Note: this formula assumes indicator intercepts are zero.
  # More accurately, the mean of y_j = Lambda * beta_j.
  # The mean of simulated_data_matrix_algebra should be close to this if the epsilon means are zero.
  # Let's directly use the means from the data generation process if available or sample means
  # as a proxy for "population" means from which the sample was drawn.
  print("Population Model-Implied Observed Variable Means (first 6):")
  # Assuming indicator intercepts are effectively zero in generation, E(y) = Lambda %*% Beta
  print(round(population_observed_means[1:6,1], 4))
}

# Lavaan fitted observed means (model-implied by the fitted model)
fitted_observed_means <- fitted.values(fit_lavaan)$mean
if (!is.null(fitted_observed_means)) {
  print("Lavaan Model-Implied (Fitted) Observed Variable Means (first 6):")
  print(round(fitted_observed_means[1:6], 4))
} else {
  print("Could not retrieve fitted observed means from lavaan object.")
}

# Sample observed means from the generated data
print("Sample Observed Variable Means from Generated Data (first 6):")
print(round(colMeans(simulated_data_matrix_algebra)[1:6], 4))
