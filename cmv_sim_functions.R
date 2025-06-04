# Functions for CMV simulations with Lavaan
library(lavaan)
library(purrr)
library(dplyr)


# given a CFA lavaan model, add a CMV marker variabled
# that loads on all measures returns
# lavaan parameter table suitable for use with `cfa`

add_cmv_marker <- function(model, marker_name = "CMV", marker_measures = 7, marker_loading = .9) {
  if (inherits(model, "string")) { # is probably model syntax string
    model <- lavaanify(model) # change to parameter table
  }
  model <- model[, -which(names(model) %in% c("id", "plabel"))] #remove the ids so we can add them back in order
  measures <- model[model$op == "=~", ] # get just the measures
  cmv_measures <- measures[, ] #make a copy
  cmv_measures$lhs <- marker_name #add the marker name to the left hand side
  new_model <- rbind(measures,                #put the free parameters back together
                     cmv_measures)
  new_model$free <- seq_len(nrow(new_model)) #number them
  new_model <- rbind(new_model, model[model$op != "=~", ]) #add in the rest of the parameters

  new_model <- rbind(new_model, seq_len(marker_measures) |>
    purrr::map_dfr(\(x) {
      rbind(
        data.frame(
          lhs = marker_name,
          # add in measures for CMV
          op = "=~",
          rhs = paste0(marker_name, x),
          user = 0,
          block = 1,
          group = 1,
          free = 0,
          exo = 0,
          ustart = marker_loading,
          label = ""
        )
      )
    }))

  # add in a variance term for the CMV measures
  new_model <- rbind(new_model, 1:marker_measures |>
    purrr::map_dfr(\(x) {
      rbind(
        data.frame(
          lhs = paste0(marker_name, x),
          op = "~~",
          rhs = paste0(marker_name, x),
          user = 0,
          block = 1,
          group = 1,
          free = 0,
          exo = 0,
          ustart = 0,
          label = ""
        )
      )
    }))

  # add in a variance term for the CMV variable
  new_model <-
    rbind(
      new_model,
      data.frame(
        lhs = marker_name,
        op = "~~",
        rhs = marker_name,
        user = 0,
        block = 1,
        group = 1,
        free = 0,
        exo = 0,
        ustart = 1,
        label = ""
      )
    )
  new_model$id <- seq_len(nrow(new_model)) #add the id back in
  new_model <-
    new_model[, c("id", setdiff(names(new_model), "id"))] #make it the first column
  new_model$plabel <-
    paste0(".p", new_model$id, ".") #add the plabel back in
  return(data.frame(new_model))
}

get_test_model <- function() {
  return("
     # measurement model
  X =~ x1 + x2 + x3 + x4 + x5
  Y =~ y1 + y2 + y3 + y4 + x5
  Z =~ z1 + z2 + z3 + z4 + z5
  
  # covariances
  X ~~ Y
  Y ~~ Z
  X ~~ Z")
}

test_add_cmv_marker <- function() {
  get_test_model() |> lavaanify() |> add_cmv_marker() |> print()
}

# given a lavaan model in parameter table form, fix the loadings on the factor to loading
# todo: vectorize
fix_loadings <- function(model, factor, loading) {
  model[model$op == "=~" & model$lhs == factor, "ustart"] <- loading
  return(model)
}

test_fix_loadings <- function() {
  get_test_model() |>
    lavaanify() |>
    fix_loadings("X", 0.7) |>
    fix_loadings("Y", 0.8) |>
    fix_loadings("Z", 0.9) |>
    print()
}

# fix all of the CMV marker variable loadings (or all the other loadings if invert=TRUE) to loading
# lavaan model should be in parameter table format
fix_cmv_loadings <- function(model, loading, invert = FALSE) {
  if (invert) {
    model[model$op == "=~" & model$lhs != "CMV", "ustart"] <- loading
    return(model)
  } else {
    return(fix_loadings(model, "CMV", loading))
  }
}

# fix all covariances in the lavaan model in parameter table format to covariance
fix_covariances <- function(model, covariance) {
  model[model$op == "~~" & model$lhs != model$rhs, "ustart"] <- covariance
  return(model)
}

test_fix_covariances <- function() {
  get_test_model() |> lavaanify() |> fix_covariances(0) |> print()
}

# ChatGPT made this function to convert from a parameter table to a lavaan model syntax string
convert_to_lavaan_syntax <- function(param_table) {
  # Extract measurement model paths
  measurement_model <- param_table |>
    filter(op == "=~") |>
    group_by(lhs) |>
    summarize(rhs = paste(rhs, collapse = " + ")) |>
    rowwise() |>
    mutate(measurement = paste(lhs, "=~", rhs)) |>
    pull(measurement)

  # Extract variances and covariances
  var_cov <- param_table |>
    filter(op == "~~") |>
    rowwise() |>
    mutate(var_cov = paste(lhs, "~~", rhs)) |>
    pull(var_cov)

  # Combine the paths
  model_syntax <- c(measurement_model, var_cov)

  # Convert to a single string
  model_string <- paste(model_syntax, collapse = "\n  ")

  return(model_string)
}

test_convert_to_lavaan_syntax <- function() {
  test_cfa_model <- get_test_model()
  model <- lavaanify(test_cfa_model)
  cat(convert_to_lavaan_syntax(model))
}
