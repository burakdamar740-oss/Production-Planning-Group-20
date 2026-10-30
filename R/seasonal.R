# =============================================================================
# SEASONAL METHODS - SHINY MODULE
#   Method A: Classical seasonal decomposition with moving averages
#             (additive or multiplicative), plus a regression trend line on the
#             deseasonalized data for fitted values and forecasts.
#   Method B: Winters' method (triple exponential smoothing: level, trend and
#             seasonal factors) with additive or multiplicative seasonality.
#
# Defines pure calculation functions and the Shiny module:
#   seasonal_ui(id)
#   seasonal_server(id, selected_series)
#
# `selected_series` must be a reactive that returns a numeric vector.
# Nothing is executed when this file is sourced.
#
# Notation: m = seasonal period, T = number of observations.
#
# Decomposition
#   additive       : Y_t = Trend_t + Seasonal_t + Random_t
#   multiplicative : Y_t = Trend_t * Seasonal_t * Random_t
#   Forecast       : (regression trend line of the deseasonalized data)
#                    combined with the seasonal index of the period.
#
# Winters' method, multiplicative seasonality (S = seasonal factor)
#   Forecast : Yhat_t = (L_(t-1) + T_(t-1)) * S_(t-m)
#   Level    : L_t = alpha * (Y_t / S_(t-m)) + (1 - alpha) * (L_(t-1) + T_(t-1))
#   Trend    : T_t = beta * (L_t - L_(t-1)) + (1 - beta) * T_(t-1)
#   Seasonal : S_t = gamma * (Y_t / L_t) + (1 - gamma) * S_(t-m)
#   h-step   : Yhat_(T+h) = (L_T + h * T_T) * S_(T+h-m)
# Winters' method, additive seasonality
#   Forecast : Yhat_t = L_(t-1) + T_(t-1) + S_(t-m)
#   Level    : L_t = alpha * (Y_t - S_(t-m)) + (1 - alpha) * (L_(t-1) + T_(t-1))
#   Seasonal : S_t = gamma * (Y_t - L_t) + (1 - gamma) * S_(t-m)
#   h-step   : Yhat_(T+h) = L_T + h * T_T + S_(T+h-m)
# =============================================================================


# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------

# Position of each period inside the season (1..m).
seasonal_season_positions <- function(period_idx, season_length) {
  ((period_idx - 1) %% season_length) + 1
}

# Combines a trend/level part with a seasonal part.
seasonal_combine <- function(base_vec, seasonal_vec, type) {
  if (identical(type, "multiplicative")) base_vec * seasonal_vec else base_vec + seasonal_vec
}

# Removes the seasonal part from the data.
seasonal_deseasonalize <- function(demand_vec, seasonal_vec, type) {
  if (identical(type, "multiplicative")) demand_vec / seasonal_vec else demand_vec - seasonal_vec
}

# Number formatting that tolerates NA.
seasonal_format_number <- function(x, digits = 2) {
  if (length(x) == 0 || is.na(x)) {
    "NA"
  } else {
    format(round(x, digits), nsmall = digits)
  }
}


# ----------------------------------------------------------------------------
# Validation helpers
# ----------------------------------------------------------------------------

seasonal_validate_season_length <- function(season_length) {
  if (!is.numeric(season_length) || length(season_length) != 1 ||
      is.na(season_length) || season_length != floor(season_length) ||
      season_length < 2) {
    stop("The seasonal period must be a whole number of at least 2.")
  }
  invisible(TRUE)
}

# Checks the data and requires at least two full seasons.
seasonal_validate_demand <- function(demand_vec, season_length) {
  if (!is.numeric(demand_vec)) {
    stop("Demand data must be numeric.")
  }
  if (any(!is.finite(demand_vec))) {
    stop("Demand data cannot contain missing or non-finite values.")
  }
  seasonal_validate_season_length(season_length)
  if (length(demand_vec) < 2 * season_length) {
    stop(sprintf(
      "At least 2 full seasons are needed: %d observations for a seasonal period of %d, but the data has only %d.",
      2 * season_length, season_length, length(demand_vec)
    ))
  }
  invisible(TRUE)
}

seasonal_validate_type <- function(type) {
  if (!is.character(type) || length(type) != 1 ||
      !(type %in% c("additive", "multiplicative"))) {
    stop("The seasonality type must be 'additive' or 'multiplicative'.")
  }
  invisible(TRUE)
}

# Checks that a smoothing parameter is a single number with 0 < value <= 1.
seasonal_validate_smoothing <- function(param_value, param_name) {
  if (!is.numeric(param_value) || length(param_value) != 1 || is.na(param_value)) {
    stop(sprintf("%s must be a single number.", param_name))
  }
  if (param_value <= 0 || param_value > 1) {
    stop(sprintf("%s must be greater than 0 and at most 1 (entered: %s).",
                 param_name, param_value))
  }
  invisible(TRUE)
}

seasonal_validate_horizon <- function(horizon) {
  if (!is.numeric(horizon) || length(horizon) != 1 || is.na(horizon) ||
      horizon != floor(horizon) || horizon < 1) {
    stop("The forecast horizon must be a whole number of at least 1.")
  }
  invisible(TRUE)
}


# ----------------------------------------------------------------------------
# Method A: classical decomposition with moving averages
# ----------------------------------------------------------------------------

# Classical decomposition via stats::decompose().
# Returns the moving-average trend, the seasonal component, the random
# component and the m seasonal indices (`figure`, position 1..m).
seasonal_decompose <- function(demand_vec, season_length, type = "additive") {
  seasonal_validate_demand(demand_vec, season_length)
  seasonal_validate_type(type)
  if (identical(type, "multiplicative") && any(demand_vec <= 0)) {
    stop("Multiplicative decomposition needs all demand values to be positive. Use additive decomposition instead.")
  }
  
  ts_obj <- stats::ts(demand_vec, frequency = season_length, start = c(1, 1))
  dec    <- stats::decompose(ts_obj, type = type)
  figure <- as.numeric(dec$figure)
  
  if (length(figure) != season_length || any(!is.finite(figure))) {
    stop("The seasonal indices could not be computed. Check the seasonal period.")
  }
  
  list(
    type     = type,
    trend    = as.numeric(dec$trend),
    seasonal = as.numeric(dec$seasonal),
    random   = as.numeric(dec$random),
    figure   = figure
  )
}

# Decomposition-based fitted values and forecasts:
#   1. decompose the series and take the seasonal indices,
#   2. deseasonalize the data and fit a regression line on time,
#   3. fitted value  = trend line combined with the seasonal index,
#      forecast(T+h) = extended trend line combined with the seasonal index.
seasonal_decomposition_forecast <- function(demand_vec, season_length,
                                            type = "additive", horizon = 1) {
  seasonal_validate_horizon(horizon)
  decomposition <- seasonal_decompose(demand_vec, season_length, type)
  
  data_length  <- length(demand_vec)
  period_idx   <- seq_len(data_length)
  seasonal_vec <- decomposition$figure[seasonal_season_positions(period_idx, season_length)]
  deseason_vec <- seasonal_deseasonalize(demand_vec, seasonal_vec, type)
  
  trend_model <- stats::lm(deseason_vec ~ period_idx)
  intercept   <- unname(stats::coef(trend_model)[1])
  slope       <- unname(stats::coef(trend_model)[2])
  trend_line  <- intercept + slope * period_idx
  fitted_vec  <- seasonal_combine(trend_line, seasonal_vec, type)
  
  future_idx       <- data_length + seq_len(horizon)
  future_seasonal  <- decomposition$figure[seasonal_season_positions(future_idx, season_length)]
  future_forecasts <- seasonal_combine(intercept + slope * future_idx, future_seasonal, type)
  
  list(
    decomposition    = decomposition,
    deseasonalized   = deseason_vec,
    intercept        = intercept,
    slope            = slope,
    trend_line       = trend_line,
    fitted           = fitted_vec,
    next_forecast    = future_forecasts[1],
    future_forecasts = future_forecasts
  )
}


# ----------------------------------------------------------------------------
# Method B: Winters' method
# ----------------------------------------------------------------------------

# Automatic initial values for Winters' method:
#   - classical decomposition of the first `n_seasons` full seasons gives the
#     initial seasonal factors,
#   - a regression of the deseasonalized values on time gives
#     L_0 = intercept (value at time 0) and T_0 = slope.
seasonal_winters_init <- function(demand_vec, season_length, type = "multiplicative",
                                  n_seasons = 2) {
  seasonal_validate_demand(demand_vec, season_length)
  max_seasons <- floor(length(demand_vec) / season_length)
  if (!is.numeric(n_seasons) || length(n_seasons) != 1 || is.na(n_seasons) ||
      n_seasons != floor(n_seasons) || n_seasons < 2 || n_seasons > max_seasons) {
    stop(sprintf("The number of seasons used for initial values must be a whole number between 2 and %d.",
                 max_seasons))
  }
  
  init_length  <- n_seasons * season_length
  init_demand  <- demand_vec[seq_len(init_length)]
  decomposition <- seasonal_decompose(init_demand, season_length, type)
  
  period_idx   <- seq_len(init_length)
  seasonal_vec <- decomposition$figure[seasonal_season_positions(period_idx, season_length)]
  deseason_vec <- seasonal_deseasonalize(init_demand, seasonal_vec, type)
  
  init_model <- stats::lm(deseason_vec ~ period_idx)
  list(
    level    = unname(stats::coef(init_model)[1]),
    trend    = unname(stats::coef(init_model)[2]),
    seasonal = decomposition$figure
  )
}

# Parses manually typed seasonal factors (separated by commas, semicolons or
# spaces; use a dot for decimals). Factor k belongs to season position k.
seasonal_parse_factors <- function(factors_text, season_length, type) {
  if (is.null(factors_text) || length(factors_text) != 1 || is.na(factors_text) ||
      !nzchar(trimws(factors_text))) {
    stop(sprintf("Enter %d initial seasonal factors, separated by commas.", season_length))
  }
  parts <- strsplit(trimws(factors_text), "[,;[:space:]]+")[[1]]
  parts <- parts[nzchar(parts)]
  factor_vec <- suppressWarnings(as.numeric(parts))
  if (any(is.na(factor_vec))) {
    stop("Seasonal factors must be numbers separated by commas (use a dot for decimals).")
  }
  if (length(factor_vec) != season_length) {
    stop(sprintf("You entered %d seasonal factor(s), but the seasonal period is %d.",
                 length(factor_vec), season_length))
  }
  if (identical(type, "multiplicative") && any(factor_vec <= 0)) {
    stop("Multiplicative seasonal factors must be positive.")
  }
  factor_vec
}

# Iterative Winters' method.
# `seasonal_init[k]` is the initial seasonal factor of season position k.
# Returns:
#   forecasts        : one-step-ahead forecasts for periods 1..T
#   levels, trends   : L_t and T_t
#   seasonals        : S_t, the seasonal factor updated in period t
#   final_seasonal   : latest seasonal factor of every season position
#   next_forecast    : forecast for period T + 1
#   future_forecasts : forecasts for periods T + 1 .. T + horizon
seasonal_winters_fit <- function(demand_vec, season_length, type,
                                 alpha_value, beta_value, gamma_value,
                                 l0_value, t0_value, seasonal_init, horizon = 1) {
  seasonal_validate_demand(demand_vec, season_length)
  seasonal_validate_type(type)
  seasonal_validate_smoothing(alpha_value, "Alpha")
  seasonal_validate_smoothing(beta_value, "Beta")
  seasonal_validate_smoothing(gamma_value, "Gamma")
  seasonal_validate_horizon(horizon)
  if (!is.numeric(l0_value) || length(l0_value) != 1 || !is.finite(l0_value) ||
      !is.numeric(t0_value) || length(t0_value) != 1 || !is.finite(t0_value)) {
    stop("The initial level and trend must be single finite numbers.")
  }
  if (!is.numeric(seasonal_init) || length(seasonal_init) != season_length ||
      any(!is.finite(seasonal_init))) {
    stop(sprintf("Exactly %d finite initial seasonal factors are required.", season_length))
  }
  if (identical(type, "multiplicative") && (any(seasonal_init <= 0) || any(demand_vec <= 0))) {
    stop("Multiplicative Winters' method needs positive demand values and positive seasonal factors.")
  }
  
  data_length  <- length(demand_vec)
  forecast_vec <- numeric(data_length)
  level_vec    <- numeric(data_length)
  trend_vec    <- numeric(data_length)
  seasonal_vec <- numeric(data_length)
  
  season_state   <- seasonal_init
  previous_level <- l0_value
  previous_trend <- t0_value
  
  for (period_idx in seq_len(data_length)) {
    position         <- seasonal_season_positions(period_idx, season_length)
    old_seasonal     <- season_state[position]
    base_forecast    <- previous_level + previous_trend
    
    if (identical(type, "multiplicative")) {
      forecast_vec[period_idx] <- base_forecast * old_seasonal
      new_level    <- alpha_value * (demand_vec[period_idx] / old_seasonal) +
        (1 - alpha_value) * base_forecast
      new_seasonal <- gamma_value * (demand_vec[period_idx] / new_level) +
        (1 - gamma_value) * old_seasonal
    } else {
      forecast_vec[period_idx] <- base_forecast + old_seasonal
      new_level    <- alpha_value * (demand_vec[period_idx] - old_seasonal) +
        (1 - alpha_value) * base_forecast
      new_seasonal <- gamma_value * (demand_vec[period_idx] - new_level) +
        (1 - gamma_value) * old_seasonal
    }
    new_trend <- beta_value * (new_level - previous_level) +
      (1 - beta_value) * previous_trend
    
    level_vec[period_idx]    <- new_level
    trend_vec[period_idx]    <- new_trend
    seasonal_vec[period_idx] <- new_seasonal
    season_state[position]   <- new_seasonal
    previous_level <- new_level
    previous_trend <- new_trend
  }
  
  if (any(!is.finite(c(forecast_vec, level_vec, trend_vec, seasonal_vec)))) {
    stop("The smoothing produced non-finite values. Try different parameters or initial values.")
  }
  
  future_idx       <- data_length + seq_len(horizon)
  future_base      <- level_vec[data_length] + seq_len(horizon) * trend_vec[data_length]
  future_seasonal  <- season_state[seasonal_season_positions(future_idx, season_length)]
  future_forecasts <- seasonal_combine(future_base, future_seasonal, type)
  
  list(
    forecasts        = forecast_vec,
    levels           = level_vec,
    trends           = trend_vec,
    seasonals        = seasonal_vec,
    final_seasonal   = season_state,
    next_forecast    = future_forecasts[1],
    future_forecasts = future_forecasts
  )
}


# ----------------------------------------------------------------------------
# Error metrics
# ----------------------------------------------------------------------------

# MAD, MSE and MAPE.
# - Only periods where both actual and forecast are available are used.
# - Periods with zero actual demand are excluded from MAPE only
#   (division by zero), and counted in `n_used_mape`.
seasonal_error_metrics <- function(actual_vec, forecast_vec) {
  if (length(actual_vec) != length(forecast_vec)) {
    stop("Actual demand and forecast vectors must have the same length.")
  }
  
  valid_idx <- !is.na(actual_vec) & !is.na(forecast_vec)
  if (!any(valid_idx)) {
    return(list(MAD = NA_real_, MSE = NA_real_, MAPE = NA_real_,
                n_used = 0L, n_used_mape = 0L))
  }
  
  actual_valid   <- actual_vec[valid_idx]
  forecast_valid <- forecast_vec[valid_idx]
  error_vec      <- actual_valid - forecast_valid
  
  nonzero_idx <- actual_valid != 0
  mape_value <- if (any(nonzero_idx)) {
    mean(abs(error_vec[nonzero_idx]) / abs(actual_valid[nonzero_idx])) * 100
  } else {
    NA_real_
  }
  
  list(
    MAD         = mean(abs(error_vec)),
    MSE         = mean(error_vec^2),
    MAPE        = mape_value,
    n_used      = as.integer(sum(valid_idx)),
    n_used_mape = as.integer(sum(nonzero_idx))
  )
}

# Metrics as a small display table.
seasonal_metrics_table <- function(metrics) {
  data.frame(
    Metric = c("MAD", "MSE", "MAPE (%)"),
    Value  = c(
      seasonal_format_number(metrics$MAD, 4),
      seasonal_format_number(metrics$MSE, 4),
      seasonal_format_number(metrics$MAPE, 4)
    ),
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------------------------
# Automatic selection of alpha, beta and gamma (Winters' method)
# ----------------------------------------------------------------------------

# Finds (alpha, beta, gamma) in [0.01, 0.99]^3 that minimize the chosen
# criterion ("MAD", "MSE" or "MAPE") for the given initial values, using optim()
# with the bounded L-BFGS-B method. optim() is started from a 2 x 2 x 2 grid of
# starting points (the error surface can have several local minima) and the best
# result is kept.
seasonal_optimize_winters <- function(demand_vec, season_length, type,
                                      l0_value, t0_value, seasonal_init,
                                      criterion = "MAD") {
  seasonal_validate_demand(demand_vec, season_length)
  
  criterion_type <- toupper(trimws(as.character(criterion)))
  if (length(criterion_type) != 1 ||
      !(criterion_type %in% c("MAD", "MSE", "MAPE"))) {
    stop("Criterion must be 'MAD', 'MSE' or 'MAPE'.")
  }
  
  run_fit <- function(par) {
    seasonal_winters_fit(demand_vec, season_length, type,
                         par[1], par[2], par[3],
                         l0_value, t0_value, seasonal_init, horizon = 1)
  }
  
  # The criterion must be computable (this also checks the initial values).
  reference_fit     <- run_fit(c(0.3, 0.3, 0.3))
  reference_metrics <- seasonal_error_metrics(demand_vec, reference_fit$forecasts)
  if (is.na(reference_metrics[[criterion_type]])) {
    stop(sprintf("%s could not be computed for this data.", criterion_type))
  }
  
  # Value of the chosen criterion; a large finite penalty replaces failures.
  objective_fn <- function(par) {
    fit_result <- tryCatch(run_fit(par), error = function(e) NULL)
    if (is.null(fit_result)) {
      return(1e15)
    }
    metrics <- seasonal_error_metrics(demand_vec, fit_result$forecasts)
    criterion_value <- metrics[[criterion_type]]
    if (is.na(criterion_value) || !is.finite(criterion_value)) 1e15 else criterion_value
  }
  
  start_values <- expand.grid(alpha = c(0.2, 0.7), beta = c(0.2, 0.7), gamma = c(0.2, 0.7))
  best_value   <- Inf
  best_par     <- NULL
  
  for (start_idx in seq_len(nrow(start_values))) {
    optim_result <- tryCatch(
      stats::optim(
        par    = c(start_values$alpha[start_idx], start_values$beta[start_idx],
                   start_values$gamma[start_idx]),
        fn     = objective_fn,
        method = "L-BFGS-B",
        lower  = rep(0.01, 3),
        upper  = rep(0.99, 3)
      ),
      error = function(e) NULL
    )
    if (!is.null(optim_result) && is.finite(optim_result$value) &&
        optim_result$value < best_value) {
      best_value <- optim_result$value
      best_par   <- optim_result$par
    }
  }
  
  if (is.null(best_par) || best_value >= 1e15) {
    stop("The optimization did not converge. Try manual parameters.")
  }
  
  list(
    alpha     = best_par[1],
    beta      = best_par[2],
    gamma     = best_par[3],
    criterion = criterion_type,
    error     = best_value
  )
}


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------

# Four-panel plot of the decomposition components: observed, trend (moving
# average), seasonal and random (base R graphics).
seasonal_plot_components <- function(demand_vec, decomposition) {
  old_par <- graphics::par(mfrow = c(4, 1), mar = c(3, 4.5, 2.5, 1))
  on.exit(graphics::par(old_par), add = TRUE)
  
  period_idx <- seq_along(demand_vec)
  plot(period_idx, demand_vec, type = "l", col = "blue", lwd = 2,
       xlab = "", ylab = "Observed",
       main = paste0("Classical ", decomposition$type, " decomposition"))
  plot(period_idx, decomposition$trend, type = "l", col = "darkorange", lwd = 2,
       xlab = "", ylab = "Trend (MA)")
  plot(period_idx, decomposition$seasonal, type = "l", col = "darkgreen", lwd = 2,
       xlab = "", ylab = "Seasonal")
  plot(period_idx, decomposition$random, type = "h", col = "red",
       xlab = "Period", ylab = "Random")
  graphics::abline(h = if (identical(decomposition$type, "multiplicative")) 1 else 0,
                   lty = 2, col = "gray40")
}

# Actual demand, fitted values and future forecasts (base R graphics).
seasonal_plot_forecast <- function(demand_vec, fitted_vec, future_forecasts,
                                   plot_title = "Actual vs. Fitted and Forecast",
                                   fitted_label = "Fitted values") {
  data_length <- length(demand_vec)
  horizon     <- length(future_forecasts)
  future_idx  <- data_length + seq_len(horizon)
  y_range <- range(c(demand_vec, fitted_vec, future_forecasts), na.rm = TRUE)
  
  plot(seq_len(data_length), demand_vec,
       type = "b", pch = 16, col = "blue",
       xlim = c(1, data_length + horizon), ylim = y_range,
       xlab = "Period", ylab = "Demand", main = plot_title)
  
  lines(seq_len(data_length), fitted_vec,
        type = "b", pch = 17, col = "red", lty = 2)
  
  lines(c(data_length, future_idx), c(fitted_vec[data_length], future_forecasts),
        col = "darkgreen", lty = 3)
  points(future_idx, future_forecasts, pch = 18, cex = 1.6, col = "darkgreen")
  
  legend("topleft",
         legend = c("Actual demand", fitted_label, "Future forecast"),
         col = c("blue", "red", "darkgreen"),
         pch = c(16, 17, 18), lty = c(1, 2, 3), bty = "n")
}


# ----------------------------------------------------------------------------
# Shiny module: UI
# ----------------------------------------------------------------------------

# Small summary card (used by the server).
seasonal_make_card <- function(label, value) {
  shiny::column(
    4,
    shiny::div(
      class = "well well-sm text-center",
      shiny::tags$small(label),
      shiny::h3(value)
    )
  )
}

seasonal_ui <- function(id) {
  ns <- shiny::NS(id)
  
  is_method <- function(method_value) {
    sprintf("input['%s'] == '%s'", ns("method"), method_value)
  }
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::radioButtons(
        ns("method"), "Method",
        choices = c(
          "Seasonal decomposition (moving averages)"      = "decomposition",
          "Winters' method (triple exponential smoothing)" = "winters"
        ),
        selected = "decomposition"
      ),
      shiny::numericInput(
        ns("season_length"), "Seasonal period m (4 = quarterly, 12 = monthly, 7 = weekly)",
        value = 4, min = 2, step = 1
      ),
      shiny::radioButtons(
        ns("season_type"), "Type of seasonality",
        choices = c("Multiplicative" = "multiplicative", "Additive" = "additive"),
        selected = "multiplicative"
      ),
      shiny::numericInput(
        ns("horizon"), "Forecast horizon (periods ahead)",
        value = 4, min = 1, max = 120, step = 1
      ),
      
      # --- Winters' method settings ------------------------------------------------
      shiny::conditionalPanel(
        condition = is_method("winters"),
        shiny::hr(),
        shiny::radioButtons(
          ns("param_mode"), "How should alpha, beta and gamma be chosen?",
          choices = c(
            "Manual"                                 = "manual",
            "Auto (optimize alpha, beta and gamma)"  = "auto"
          ),
          selected = "manual"
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'manual'", ns("param_mode")),
          shiny::sliderInput(ns("alpha_value"), "Alpha (level smoothing)",
                             min = 0.01, max = 1, value = 0.3, step = 0.01),
          shiny::sliderInput(ns("beta_value"), "Beta (trend smoothing)",
                             min = 0.01, max = 1, value = 0.1, step = 0.01),
          shiny::sliderInput(ns("gamma_value"), "Gamma (seasonal smoothing)",
                             min = 0.01, max = 1, value = 0.2, step = 0.01)
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'auto'", ns("param_mode")),
          shiny::selectInput(ns("criterion"), "Minimize",
                             choices = c("MAD", "MSE", "MAPE"), selected = "MAD"),
          shiny::helpText(
            "Searches alpha, beta and gamma between 0.01 and 0.99 with optim() ",
            "(L-BFGS-B) to minimize the selected error metric."
          )
        ),
        
        shiny::hr(),
        shiny::radioButtons(
          ns("init_method"), "Initial level, trend and seasonal factors",
          choices = c(
            "Automatic (decomposition + regression)" = "auto",
            "Manual input"                           = "manual"
          ),
          selected = "auto"
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'auto'", ns("init_method")),
          shiny::numericInput(
            ns("init_seasons"), "Number of full seasons used for initial values",
            value = 2, min = 2, step = 1
          ),
          shiny::helpText(
            "Decomposes the first seasons to get the seasonal factors, then ",
            "regresses the deseasonalized values on time: the intercept is L0 ",
            "and the slope is T0."
          )
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'manual'", ns("init_method")),
          shiny::numericInput(ns("l0_manual"), "Initial level L0", value = NA, step = 1),
          shiny::numericInput(ns("t0_manual"), "Initial trend T0", value = NA, step = 0.1),
          shiny::textInput(
            ns("factors_manual"), "Initial seasonal factors (m numbers, comma-separated)",
            value = "", placeholder = "e.g. 0.8, 1.1, 1.3, 0.8"
          ),
          shiny::helpText(
            "Factor k belongs to season position k (position 1 = first period of ",
            "the series). Use a dot for decimals."
          )
        )
      ),
      
      shiny::helpText(
        "Results update automatically; the data comes from the upload panel at the top. ",
        "At least 2 full seasons of data are required."
      )
    ),
    
    shiny::mainPanel(
      # --- Decomposition outputs ---------------------------------------------------
      shiny::conditionalPanel(
        condition = is_method("decomposition"),
        shiny::uiOutput(ns("dec_summary")),
        shiny::h4("Seasonal Indices"),
        shiny::tableOutput(ns("dec_indices_table")),
        shiny::h4("Error Metrics"),
        shiny::tableOutput(ns("dec_metrics_table")),
        shiny::plotOutput(ns("dec_components_plot"), height = "520px"),
        shiny::plotOutput(ns("dec_forecast_plot"), height = "340px"),
        shiny::h4("Future Forecasts"),
        shiny::tableOutput(ns("dec_future_table")),
        shiny::h4("Detailed Results"),
        shiny::tableOutput(ns("dec_detail_table"))
      ),
      
      # --- Winters' outputs ----------------------------------------------------------
      shiny::conditionalPanel(
        condition = is_method("winters"),
        shiny::uiOutput(ns("win_summary")),
        shiny::h4("Seasonal Factors"),
        shiny::tableOutput(ns("win_factors_table")),
        shiny::h4("Error Metrics"),
        shiny::tableOutput(ns("win_metrics_table")),
        shiny::plotOutput(ns("win_forecast_plot"), height = "340px"),
        shiny::h4("Future Forecasts"),
        shiny::tableOutput(ns("win_future_table")),
        shiny::h4("Detailed Results"),
        shiny::tableOutput(ns("win_detail_table"))
      )
    )
  )
}


# ----------------------------------------------------------------------------
# Shiny module: server
# ----------------------------------------------------------------------------

# selected_series: a reactive returning a numeric vector (supplied by app.R).
seasonal_server <- function(id, selected_series) {
  stopifnot(shiny::is.reactive(selected_series))
  
  shiny::moduleServer(id, function(input, output, session) {
    
    # Turns any problem (validation message, input error, calculation error)
    # into a message on the screen instead of crashing the session.
    safe_run <- function(expr_fn) {
      tryCatch(
        expr_fn(),
        error = function(e) {
          error_msg <- conditionMessage(e)
          if (!nzchar(error_msg)) {
            error_msg <- "Upload a data file and select a numeric column in the panel above."
          }
          list(ok = FALSE, message = error_msg)
        }
      )
    }
    
    warning_box <- function(outcome) {
      shiny::div(class = "alert alert-warning", outcome$message)
    }
    
    # ---- Data and common settings --------------------------------------------------
    prepared <- shiny::reactive({
      raw_series <- selected_series()
      shiny::req(!is.null(raw_series), length(raw_series) > 0)
      shiny::validate(shiny::need(is.numeric(raw_series), "The data must be numeric."))
      demand_vec <- raw_series[is.finite(raw_series)]
      
      season_length <- input$season_length
      shiny::validate(shiny::need(
        is.numeric(season_length) && length(season_length) == 1 && !is.na(season_length) &&
          season_length == floor(season_length) && season_length >= 2,
        "The seasonal period must be a whole number of at least 2."
      ))
      shiny::validate(shiny::need(
        length(demand_vec) >= 2 * season_length,
        sprintf("At least 2 full seasons are needed: %d valid observations for a seasonal period of %d, but the data has only %d.",
                2 * season_length, season_length, length(demand_vec))
      ))
      
      horizon <- input$horizon
      shiny::validate(shiny::need(
        is.numeric(horizon) && length(horizon) == 1 && !is.na(horizon) &&
          horizon == floor(horizon) && horizon >= 1 && horizon <= 120,
        "The forecast horizon must be a whole number between 1 and 120."
      ))
      
      list(
        demand        = demand_vec,
        season_length = season_length,
        horizon       = horizon,
        type          = input$season_type
      )
    })
    
    # =========================================================================
    # Method A: decomposition
    # =========================================================================
    dec_result <- shiny::reactive({
      safe_run(function() {
        prep   <- prepared()
        result <- seasonal_decomposition_forecast(prep$demand, prep$season_length,
                                                  prep$type, prep$horizon)
        list(
          ok      = TRUE,
          prep    = prep,
          result  = result,
          metrics = seasonal_error_metrics(prep$demand, result$fitted)
        )
      })
    })
    
    output$dec_summary <- shiny::renderUI({
      outcome <- dec_result()
      if (!outcome$ok) {
        return(warning_box(outcome))
      }
      prep   <- outcome$prep
      result <- outcome$result
      
      type_text <- if (prep$type == "multiplicative") "Multiplicative" else "Additive"
      shiny::tagList(
        shiny::h4(paste0(type_text, " decomposition (moving averages), m = ", prep$season_length)),
        shiny::fluidRow(
          seasonal_make_card("Seasonal period (m)", prep$season_length),
          seasonal_make_card("Trend slope per period", seasonal_format_number(result$slope, 4)),
          seasonal_make_card("Next-period forecast", seasonal_format_number(result$next_forecast))
        ),
        shiny::tags$p(
          shiny::strong("Trend line (regression on deseasonalized data): "),
          paste0("Trend(t) = ", seasonal_format_number(result$intercept, 4),
                 " + ", seasonal_format_number(result$slope, 4), " * t")
        ),
        shiny::tags$p(
          shiny::strong("Forecast rule: "),
          if (prep$type == "multiplicative") {
            "trend line value x seasonal index of the period"
          } else {
            "trend line value + seasonal index of the period"
          }
        ),
        if (!is.na(outcome$metrics$n_used_mape) &&
            outcome$metrics$n_used_mape < outcome$metrics$n_used) {
          shiny::tags$p(shiny::em(
            "Note: periods with zero actual demand were excluded from MAPE."
          ))
        }
      )
    })
    
    output$dec_indices_table <- shiny::renderTable({
      outcome <- dec_result()
      shiny::req(outcome$ok)
      figure <- outcome$result$decomposition$figure
      data.frame(
        "Season position" = seq_along(figure),
        "Seasonal index"  = format(round(figure, 4), nsmall = 4),
        check.names       = FALSE,
        stringsAsFactors  = FALSE
      )
    })
    
    output$dec_metrics_table <- shiny::renderTable({
      outcome <- dec_result()
      shiny::req(outcome$ok)
      seasonal_metrics_table(outcome$metrics)
    })
    
    output$dec_components_plot <- shiny::renderPlot({
      outcome <- dec_result()
      shiny::req(outcome$ok)
      seasonal_plot_components(outcome$prep$demand, outcome$result$decomposition)
    })
    
    output$dec_forecast_plot <- shiny::renderPlot({
      outcome <- dec_result()
      shiny::req(outcome$ok)
      seasonal_plot_forecast(
        outcome$prep$demand,
        outcome$result$fitted,
        outcome$result$future_forecasts,
        "Decomposition: Actual vs. Fitted (trend line x seasonal) and Forecast",
        "Fitted (trend line + seasonal)"
      )
    })
    
    output$dec_future_table <- shiny::renderTable({
      outcome <- dec_result()
      shiny::req(outcome$ok)
      future_forecasts <- outcome$result$future_forecasts
      future_idx <- length(outcome$prep$demand) + seq_along(future_forecasts)
      data.frame(
        "Period"          = future_idx,
        "Season position" = seasonal_season_positions(future_idx, outcome$prep$season_length),
        "Forecast"        = future_forecasts,
        check.names       = FALSE
      )
    }, digits = 2)
    
    output$dec_detail_table <- shiny::renderTable({
      outcome <- dec_result()
      shiny::req(outcome$ok)
      result <- outcome$result
      demand_vec <- outcome$prep$demand
      data.frame(
        "Period"         = seq_along(demand_vec),
        "Actual"         = demand_vec,
        "Trend (MA)"     = result$decomposition$trend,
        "Seasonal"       = result$decomposition$seasonal,
        "Random"         = result$decomposition$random,
        "Deseasonalized" = result$deseasonalized,
        "Trend Line"     = result$trend_line,
        "Fitted"         = result$fitted,
        "Error"          = demand_vec - result$fitted,
        check.names      = FALSE
      )
    }, na = "-", digits = 2)
    
    # =========================================================================
    # Method B: Winters' method
    # =========================================================================
    win_result <- shiny::reactive({
      safe_run(function() {
        prep          <- prepared()
        demand_vec    <- prep$demand
        season_length <- prep$season_length
        type          <- prep$type
        
        # Initial level, trend and seasonal factors
        if (identical(input$init_method, "manual")) {
          l0_input <- input$l0_manual
          t0_input <- input$t0_manual
          shiny::validate(
            shiny::need(is.numeric(l0_input) && length(l0_input) == 1 && !is.na(l0_input),
                        "Enter a numeric value for the initial level L0."),
            shiny::need(is.numeric(t0_input) && length(t0_input) == 1 && !is.na(t0_input),
                        "Enter a numeric value for the initial trend T0.")
          )
          init_values <- list(
            level    = l0_input,
            trend    = t0_input,
            seasonal = seasonal_parse_factors(input$factors_manual, season_length, type)
          )
          init_label <- "Manual input"
        } else {
          n_seasons   <- input$init_seasons
          max_seasons <- floor(length(demand_vec) / season_length)
          shiny::validate(shiny::need(
            is.numeric(n_seasons) && length(n_seasons) == 1 && !is.na(n_seasons) &&
              n_seasons == floor(n_seasons) && n_seasons >= 2 && n_seasons <= max_seasons,
            sprintf("The number of seasons used for initial values must be a whole number between 2 and %d.",
                    max_seasons)
          ))
          init_values <- seasonal_winters_init(demand_vec, season_length, type, n_seasons)
          init_label  <- sprintf(
            "Automatic: decomposition and regression on the first %d seasons (%d observations)",
            n_seasons, n_seasons * season_length
          )
        }
        
        # Smoothing parameters
        optimization <- NULL
        if (identical(input$param_mode, "auto")) {
          optimization <- seasonal_optimize_winters(
            demand_vec, season_length, type,
            init_values$level, init_values$trend, init_values$seasonal,
            input$criterion
          )
          alpha_value <- optimization$alpha
          beta_value  <- optimization$beta
          gamma_value <- optimization$gamma
        } else {
          alpha_input <- input$alpha_value
          beta_input  <- input$beta_value
          gamma_input <- input$gamma_value
          is_valid_param <- function(x) {
            is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0 && x <= 1
          }
          shiny::validate(
            shiny::need(is_valid_param(alpha_input), "Alpha must be greater than 0 and at most 1."),
            shiny::need(is_valid_param(beta_input),  "Beta must be greater than 0 and at most 1."),
            shiny::need(is_valid_param(gamma_input), "Gamma must be greater than 0 and at most 1.")
          )
          alpha_value <- alpha_input
          beta_value  <- beta_input
          gamma_value <- gamma_input
        }
        
        fit_result <- seasonal_winters_fit(
          demand_vec, season_length, type,
          alpha_value, beta_value, gamma_value,
          init_values$level, init_values$trend, init_values$seasonal,
          horizon = prep$horizon
        )
        
        list(
          ok           = TRUE,
          prep         = prep,
          alpha        = alpha_value,
          beta         = beta_value,
          gamma        = gamma_value,
          param_mode   = input$param_mode,
          optimization = optimization,
          init         = init_values,
          init_label   = init_label,
          fit          = fit_result,
          metrics      = seasonal_error_metrics(demand_vec, fit_result$forecasts)
        )
      })
    })
    
    output$win_summary <- shiny::renderUI({
      outcome <- win_result()
      if (!outcome$ok) {
        return(warning_box(outcome))
      }
      prep    <- outcome$prep
      is_auto <- outcome$param_mode == "auto"
      
      type_text <- if (prep$type == "multiplicative") "multiplicative" else "additive"
      mode_text <- if (is_auto) {
        paste0("Auto mode: alpha, beta and gamma optimized to minimize ",
               outcome$optimization$criterion)
      } else {
        "Manual mode: alpha, beta and gamma chosen by the user"
      }
      
      third_card <- if (is_auto) {
        seasonal_make_card(paste0("Lowest ", outcome$optimization$criterion),
                           seasonal_format_number(outcome$optimization$error, 4))
      } else {
        seasonal_make_card("Periods used for errors", outcome$metrics$n_used)
      }
      
      shiny::tagList(
        shiny::h4(paste0("Winters' method (", type_text, " seasonality), m = ", prep$season_length)),
        shiny::tags$p(mode_text),
        shiny::fluidRow(
          seasonal_make_card(if (is_auto) "Optimal alpha" else "Alpha",
                             seasonal_format_number(outcome$alpha, 4)),
          seasonal_make_card(if (is_auto) "Optimal beta" else "Beta",
                             seasonal_format_number(outcome$beta, 4)),
          seasonal_make_card(if (is_auto) "Optimal gamma" else "Gamma",
                             seasonal_format_number(outcome$gamma, 4))
        ),
        shiny::fluidRow(
          seasonal_make_card("Initial level (L0)", seasonal_format_number(outcome$init$level)),
          seasonal_make_card("Initial trend (T0)", seasonal_format_number(outcome$init$trend, 4)),
          seasonal_make_card("Next-period forecast", seasonal_format_number(outcome$fit$next_forecast))
        ),
        shiny::fluidRow(third_card),
        shiny::tags$p(shiny::strong("Initial values method: "), outcome$init_label),
        if (!is.na(outcome$metrics$n_used_mape) &&
            outcome$metrics$n_used_mape < outcome$metrics$n_used) {
          shiny::tags$p(shiny::em(
            "Note: periods with zero actual demand were excluded from MAPE."
          ))
        }
      )
    })
    
    output$win_factors_table <- shiny::renderTable({
      outcome <- win_result()
      shiny::req(outcome$ok)
      data.frame(
        "Season position"       = seq_along(outcome$init$seasonal),
        "Initial factor"        = format(round(outcome$init$seasonal, 4), nsmall = 4),
        "Final factor (after fitting)" = format(round(outcome$fit$final_seasonal, 4), nsmall = 4),
        check.names             = FALSE,
        stringsAsFactors        = FALSE
      )
    })
    
    output$win_metrics_table <- shiny::renderTable({
      outcome <- win_result()
      shiny::req(outcome$ok)
      seasonal_metrics_table(outcome$metrics)
    })
    
    output$win_forecast_plot <- shiny::renderPlot({
      outcome <- win_result()
      shiny::req(outcome$ok)
      seasonal_plot_forecast(
        outcome$prep$demand,
        outcome$fit$forecasts,
        outcome$fit$future_forecasts,
        paste0("Winters' Method (alpha = ", seasonal_format_number(outcome$alpha, 3),
               ", beta = ", seasonal_format_number(outcome$beta, 3),
               ", gamma = ", seasonal_format_number(outcome$gamma, 3), ")"),
        "Fitted (one-step forecast)"
      )
    })
    
    output$win_future_table <- shiny::renderTable({
      outcome <- win_result()
      shiny::req(outcome$ok)
      future_forecasts <- outcome$fit$future_forecasts
      future_idx <- length(outcome$prep$demand) + seq_along(future_forecasts)
      data.frame(
        "Period"          = future_idx,
        "Season position" = seasonal_season_positions(future_idx, outcome$prep$season_length),
        "Forecast"        = future_forecasts,
        check.names       = FALSE
      )
    }, digits = 2)
    
    output$win_detail_table <- shiny::renderTable({
      outcome <- win_result()
      shiny::req(outcome$ok)
      demand_vec <- outcome$prep$demand
      data.frame(
        "Period"      = seq_along(demand_vec),
        "Actual"      = demand_vec,
        "Level (L_t)" = outcome$fit$levels,
        "Trend (T_t)" = outcome$fit$trends,
        "Seasonal (S_t)" = outcome$fit$seasonals,
        "Forecast"    = outcome$fit$forecasts,
        "Error (e_t)" = demand_vec - outcome$fit$forecasts,
        check.names   = FALSE
      )
    }, na = "-", digits = 3)
  })
}