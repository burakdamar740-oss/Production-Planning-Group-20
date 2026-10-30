# =============================================================================
# HOLT'S METHOD (DOUBLE EXPONENTIAL SMOOTHING FOR TREND SERIES) - SHINY MODULE
#
# Defines pure calculation functions and the Shiny module:
#   holt_ui(id)
#   holt_server(id, selected_series)
#
# `selected_series` must be a reactive that returns a numeric vector.
# Nothing is executed when this file is sourced.
#
# Model (T = number of observations, L_0 / T_0 = initial level / trend):
#   Forecast for period t :  Yhat_t = L_(t-1) + T_(t-1)
#   Level                 :  L_t = alpha * Y_t + (1 - alpha) * (L_(t-1) + T_(t-1))
#   Trend                 :  T_t = beta * (L_t - L_(t-1)) + (1 - beta) * T_(t-1)
#   Next-period forecast  :  Yhat_(T+1) = L_T + T_T
#   h-step-ahead forecast :  Yhat_(T+h) = L_T + h * T_T
# =============================================================================


# ----------------------------------------------------------------------------
# Validation helpers
# ----------------------------------------------------------------------------

# Checks that the demand data is a usable numeric vector.
holt_validate_demand <- function(demand_vec) {
  if (!is.numeric(demand_vec)) {
    stop("Demand data must be numeric.")
  }
  if (length(demand_vec) < 3) {
    stop("Holt's method needs at least 3 observations.")
  }
  if (any(!is.finite(demand_vec))) {
    stop("Demand data cannot contain missing or non-finite values.")
  }
  invisible(TRUE)
}

# Checks that a smoothing parameter is a single number with 0 < value <= 1.
holt_validate_smoothing <- function(param_value, param_name) {
  if (!is.numeric(param_value) || length(param_value) != 1 || is.na(param_value)) {
    stop(sprintf("%s must be a single number.", param_name))
  }
  if (param_value <= 0 || param_value > 1) {
    stop(sprintf("%s must be greater than 0 and at most 1 (entered: %s).",
                 param_name, param_value))
  }
  invisible(TRUE)
}

# Checks that a value is a single finite number.
holt_validate_finite <- function(value, label) {
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value)) {
    stop(sprintf("%s must be a single finite number.", label))
  }
  invisible(TRUE)
}


# ----------------------------------------------------------------------------
# Initial values
# ----------------------------------------------------------------------------

# Initial level and trend from a linear regression of the demand on time:
#   fit <- lm(demand ~ seq_along(demand))
#   L_0 = intercept (coef[1]),  T_0 = slope (coef[2])
holt_regression_init <- function(demand_vec) {
  holt_validate_demand(demand_vec)
  fit_model <- stats::lm(demand_vec ~ seq_along(demand_vec))
  coefficients_vec <- stats::coef(fit_model)
  list(
    l0 = unname(coefficients_vec[1]),
    t0 = unname(coefficients_vec[2])
  )
}


# ----------------------------------------------------------------------------
# Forecasting
# ----------------------------------------------------------------------------

# Iterative Holt smoothing.
# Returns:
#   forecasts        : one-step-ahead forecasts Yhat_t = L_(t-1) + T_(t-1), t = 1..T
#   levels, trends   : L_t and T_t, t = 1..T
#   next_forecast    : Yhat_(T+1) = L_T + T_T
#   future_forecasts : Yhat_(T+h) = L_T + h * T_T, h = 1..horizon
holt_fit <- function(demand_vec, alpha_value, beta_value, l0_value, t0_value,
                     horizon = 1) {
  holt_validate_demand(demand_vec)
  holt_validate_smoothing(alpha_value, "Alpha")
  holt_validate_smoothing(beta_value, "Beta")
  holt_validate_finite(l0_value, "The initial level")
  holt_validate_finite(t0_value, "The initial trend")
  if (!is.numeric(horizon) || length(horizon) != 1 || is.na(horizon) ||
      horizon != floor(horizon) || horizon < 1) {
    stop("The forecast horizon must be a whole number of at least 1.")
  }
  
  data_length  <- length(demand_vec)
  forecast_vec <- numeric(data_length)
  level_vec    <- numeric(data_length)
  trend_vec    <- numeric(data_length)
  
  previous_level <- l0_value
  previous_trend <- t0_value
  for (period_idx in seq_len(data_length)) {
    forecast_vec[period_idx] <- previous_level + previous_trend
    level_vec[period_idx]    <- alpha_value * demand_vec[period_idx] +
      (1 - alpha_value) * (previous_level + previous_trend)
    trend_vec[period_idx]    <- beta_value * (level_vec[period_idx] - previous_level) +
      (1 - beta_value) * previous_trend
    previous_level <- level_vec[period_idx]
    previous_trend <- trend_vec[period_idx]
  }
  
  future_forecasts <- level_vec[data_length] + seq_len(horizon) * trend_vec[data_length]
  
  list(
    forecasts        = forecast_vec,
    levels           = level_vec,
    trends           = trend_vec,
    next_forecast    = level_vec[data_length] + trend_vec[data_length],
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
holt_error_metrics <- function(actual_vec, forecast_vec) {
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


# ----------------------------------------------------------------------------
# Automatic selection of alpha and beta
# ----------------------------------------------------------------------------

# Finds (alpha, beta) in [0.01, 0.99] x [0.01, 0.99] that minimize the chosen
# criterion ("MAD", "MSE" or "MAPE") for given initial values, using optim()
# with the bounded L-BFGS-B method.
# The error surface can have several local minima, so optim() is started from a
# 3 x 3 grid of starting points and the best result is kept.
holt_optimize_params <- function(demand_vec, l0_value, t0_value, criterion = "MAD") {
  holt_validate_demand(demand_vec)
  
  criterion_type <- toupper(trimws(as.character(criterion)))
  if (length(criterion_type) != 1 ||
      !(criterion_type %in% c("MAD", "MSE", "MAPE"))) {
    stop("Criterion must be 'MAD', 'MSE' or 'MAPE'.")
  }
  
  # Value of the chosen criterion for par = c(alpha, beta).
  # A large finite penalty replaces non-finite values so optim() never fails.
  objective_fn <- function(par) {
    fit_result <- holt_fit(demand_vec, par[1], par[2], l0_value, t0_value, horizon = 1)
    metrics    <- holt_error_metrics(demand_vec, fit_result$forecasts)
    criterion_value <- metrics[[criterion_type]]
    if (is.na(criterion_value) || !is.finite(criterion_value)) 1e15 else criterion_value
  }
  
  # The criterion must be computable (for example, MAPE is not computable when
  # every actual value is zero).
  reference_metrics <- holt_error_metrics(
    demand_vec,
    holt_fit(demand_vec, 0.5, 0.5, l0_value, t0_value)$forecasts
  )
  if (is.na(reference_metrics[[criterion_type]])) {
    stop(sprintf("%s could not be computed for this data.", criterion_type))
  }
  
  start_values <- expand.grid(alpha = c(0.1, 0.5, 0.9), beta = c(0.1, 0.5, 0.9))
  best_value   <- Inf
  best_par     <- NULL
  
  for (start_idx in seq_len(nrow(start_values))) {
    optim_result <- tryCatch(
      stats::optim(
        par    = c(start_values$alpha[start_idx], start_values$beta[start_idx]),
        fn     = objective_fn,
        method = "L-BFGS-B",
        lower  = c(0.01, 0.01),
        upper  = c(0.99, 0.99)
      ),
      error = function(e) NULL
    )
    if (!is.null(optim_result) && is.finite(optim_result$value) &&
        optim_result$value < best_value) {
      best_value <- optim_result$value
      best_par   <- optim_result$par
    }
  }
  
  if (is.null(best_par)) {
    stop("The optimization did not converge. Try manual parameters.")
  }
  
  list(
    alpha     = best_par[1],
    beta      = best_par[2],
    criterion = criterion_type,
    error     = best_value
  )
}


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------

# Historical demand, fitted one-step-ahead forecasts and future forecasts
# (base R graphics).
holt_plot <- function(demand_vec, forecast_vec, future_forecasts,
                      plot_title = "Holt's Method") {
  data_length <- length(demand_vec)
  horizon     <- length(future_forecasts)
  future_idx  <- data_length + seq_len(horizon)
  y_range <- range(c(demand_vec, forecast_vec, future_forecasts), na.rm = TRUE)
  
  plot(seq_len(data_length), demand_vec,
       type = "b", pch = 16, col = "blue",
       xlim = c(1, data_length + horizon), ylim = y_range,
       xlab = "Period", ylab = "Demand", main = plot_title)
  
  # Fitted values (one-step-ahead forecasts)
  lines(seq_len(data_length), forecast_vec,
        type = "b", pch = 17, col = "red", lty = 2)
  
  # Future forecasts, connected to the last fitted value
  lines(c(data_length, future_idx), c(forecast_vec[data_length], future_forecasts),
        col = "darkgreen", lty = 3)
  points(future_idx, future_forecasts, pch = 18, cex = 1.6, col = "darkgreen")
  
  legend("topleft",
         legend = c("Actual demand", "Fitted (one-step forecast)", "Future forecast"),
         col = c("blue", "red", "darkgreen"),
         pch = c(16, 17, 18), lty = c(1, 2, 3), bty = "n")
}


# ----------------------------------------------------------------------------
# Shiny module: UI
# ----------------------------------------------------------------------------

holt_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      # --- Smoothing parameters ----------------------------------------------------
      shiny::radioButtons(
        ns("param_mode"), "How should alpha and beta be chosen?",
        choices = c(
          "Manual"                          = "manual",
          "Auto (optimize alpha and beta)"  = "auto"
        ),
        selected = "manual"
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'manual'", ns("param_mode")),
        shiny::sliderInput(
          ns("alpha_value"), "Alpha (level smoothing)",
          min = 0.01, max = 1, value = 0.3, step = 0.01
        ),
        shiny::sliderInput(
          ns("beta_value"), "Beta (trend smoothing)",
          min = 0.01, max = 1, value = 0.2, step = 0.01
        )
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'auto'", ns("param_mode")),
        shiny::selectInput(
          ns("criterion"), "Minimize",
          choices = c("MAD", "MSE", "MAPE"), selected = "MAD"
        ),
        shiny::helpText(
          "Searches alpha and beta between 0.01 and 0.99 with optim() ",
          "(L-BFGS-B) to minimize the selected error metric."
        )
      ),
      
      shiny::hr(),
      
      # --- Initial values ------------------------------------------------------------
      shiny::radioButtons(
        ns("init_method"), "Initial level (L0) and trend (T0)",
        choices = c(
          "Linear regression fit" = "regression",
          "Manual input"          = "manual"
        ),
        selected = "regression"
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'regression'", ns("init_method")),
        shiny::helpText(
          "Fits lm(demand ~ time) on all observations (time = 1, 2, ...). ",
          "The intercept is used as L0 and the slope as T0."
        )
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'manual'", ns("init_method")),
        shiny::numericInput(ns("l0_manual"), "Initial level L0", value = NA, step = 1),
        shiny::numericInput(ns("t0_manual"), "Initial trend T0", value = NA, step = 0.1)
      ),
      
      shiny::hr(),
      
      shiny::numericInput(
        ns("horizon"), "Forecast horizon (periods ahead)",
        value = 3, min = 1, max = 60, step = 1
      ),
      
      shiny::helpText(
        "Forecast: Y-hat(t) = L(t-1) + T(t-1).  ",
        "Level: L(t) = alpha * Y(t) + (1 - alpha) * (L(t-1) + T(t-1)).  ",
        "Trend: T(t) = beta * (L(t) - L(t-1)) + (1 - beta) * T(t-1).  ",
        "Next period: Y-hat(T+1) = L(T) + T(T). ",
        "Results update automatically; the data comes from the upload panel at the top."
      )
    ),
    
    shiny::mainPanel(
      shiny::uiOutput(ns("summary")),
      shiny::h4("Error Metrics"),
      shiny::tableOutput(ns("metrics_table")),
      shiny::plotOutput(ns("forecast_plot"), height = "340px"),
      shiny::h4("Future Forecasts"),
      shiny::tableOutput(ns("future_table")),
      shiny::h4("Detailed Results"),
      shiny::tableOutput(ns("detail_table"))
    )
  )
}


# ----------------------------------------------------------------------------
# Shiny module: server
# ----------------------------------------------------------------------------

# selected_series: a reactive returning a numeric vector (supplied by app.R).
holt_server <- function(id, selected_series) {
  stopifnot(shiny::is.reactive(selected_series))
  
  shiny::moduleServer(id, function(input, output, session) {
    
    # Number formatting helper that tolerates NA.
    format_number <- function(x, digits = 2) {
      if (length(x) == 0 || is.na(x)) {
        "NA"
      } else {
        format(round(x, digits), nsmall = digits)
      }
    }
    
    # ---- Data check ---------------------------------------------------------------
    # Keeps only valid (finite) observations and requires at least 3 of them.
    checked_series <- shiny::reactive({
      raw_series <- selected_series()
      shiny::req(!is.null(raw_series), length(raw_series) > 0)
      shiny::validate(shiny::need(is.numeric(raw_series), "The data must be numeric."))
      
      demand_vec <- raw_series[is.finite(raw_series)]
      shiny::validate(shiny::need(
        length(demand_vec) >= 3,
        "Holt's method needs at least 3 valid observations."
      ))
      demand_vec
    })
    
    # ---- Main computation ------------------------------------------------------------
    # Any problem (validation message, input error, calculation error) becomes a
    # message on the screen instead of crashing the session.
    holt_result <- shiny::reactive({
      tryCatch({
        demand_vec <- checked_series()
        
        # Forecast horizon
        horizon_input <- input$horizon
        shiny::validate(shiny::need(
          is.numeric(horizon_input) && length(horizon_input) == 1 &&
            !is.na(horizon_input) && horizon_input == floor(horizon_input) &&
            horizon_input >= 1 && horizon_input <= 60,
          "The forecast horizon must be a whole number between 1 and 60."
        ))
        
        # Initial level L0 and trend T0
        if (identical(input$init_method, "manual")) {
          l0_input <- input$l0_manual
          t0_input <- input$t0_manual
          shiny::validate(
            shiny::need(is.numeric(l0_input) && length(l0_input) == 1 && !is.na(l0_input),
                        "Enter a numeric value for the initial level L0."),
            shiny::need(is.numeric(t0_input) && length(t0_input) == 1 && !is.na(t0_input),
                        "Enter a numeric value for the initial trend T0.")
          )
          l0_value   <- l0_input
          t0_value   <- t0_input
          init_label <- "Manual input"
        } else {
          init_values <- holt_regression_init(demand_vec)
          l0_value    <- init_values$l0
          t0_value    <- init_values$t0
          init_label  <- sprintf(
            "Linear regression fit on all %d observations (L0 = intercept, T0 = slope)",
            length(demand_vec)
          )
        }
        
        # Alpha and beta
        optimization <- NULL
        if (identical(input$param_mode, "auto")) {
          optimization <- holt_optimize_params(demand_vec, l0_value, t0_value, input$criterion)
          alpha_value  <- optimization$alpha
          beta_value   <- optimization$beta
        } else {
          alpha_input <- input$alpha_value
          beta_input  <- input$beta_value
          shiny::validate(
            shiny::need(is.numeric(alpha_input) && length(alpha_input) == 1 &&
                          !is.na(alpha_input) && alpha_input > 0 && alpha_input <= 1,
                        "Alpha must be greater than 0 and at most 1."),
            shiny::need(is.numeric(beta_input) && length(beta_input) == 1 &&
                          !is.na(beta_input) && beta_input > 0 && beta_input <= 1,
                        "Beta must be greater than 0 and at most 1.")
          )
          alpha_value <- alpha_input
          beta_value  <- beta_input
        }
        
        fit_result <- holt_fit(demand_vec, alpha_value, beta_value,
                               l0_value, t0_value, horizon = horizon_input)
        metrics    <- holt_error_metrics(demand_vec, fit_result$forecasts)
        
        list(
          ok           = TRUE,
          demand       = demand_vec,
          alpha        = alpha_value,
          beta         = beta_value,
          param_mode   = input$param_mode,
          optimization = optimization,
          l0           = l0_value,
          t0           = t0_value,
          init_label   = init_label,
          fit          = fit_result,
          metrics      = metrics
        )
      }, error = function(e) {
        error_msg <- conditionMessage(e)
        if (!nzchar(error_msg)) {
          error_msg <- "Upload a data file and select a numeric column in the panel above."
        }
        list(ok = FALSE, message = error_msg)
      })
    })
    
    # ---- Summary -------------------------------------------------------------------------
    output$summary <- shiny::renderUI({
      outcome <- holt_result()
      if (!outcome$ok) {
        return(shiny::div(class = "alert alert-warning", outcome$message))
      }
      
      is_auto <- outcome$param_mode == "auto"
      mode_text <- if (is_auto) {
        paste0("Auto mode: alpha and beta optimized to minimize ",
               outcome$optimization$criterion)
      } else {
        "Manual mode: alpha and beta chosen by the user"
      }
      
      make_card <- function(label, value) {
        shiny::column(
          4,
          shiny::div(
            class = "well well-sm text-center",
            shiny::tags$small(label),
            shiny::h3(value)
          )
        )
      }
      
      third_card <- if (is_auto) {
        make_card(paste0("Lowest ", outcome$optimization$criterion),
                  format_number(outcome$optimization$error, 4))
      } else {
        make_card("Periods used for errors", outcome$metrics$n_used)
      }
      
      shiny::tagList(
        shiny::h4(mode_text),
        shiny::fluidRow(
          make_card(if (is_auto) "Optimal alpha" else "Alpha", format_number(outcome$alpha, 4)),
          make_card(if (is_auto) "Optimal beta" else "Beta", format_number(outcome$beta, 4)),
          make_card("Next-period forecast", format_number(outcome$fit$next_forecast))
        ),
        shiny::fluidRow(
          make_card("Initial level (L0)", format_number(outcome$l0)),
          make_card("Initial trend (T0)", format_number(outcome$t0, 4)),
          third_card
        ),
        shiny::tags$p(shiny::strong("Initial values method: "), outcome$init_label),
        if (!is.na(outcome$metrics$n_used_mape) &&
            outcome$metrics$n_used_mape < outcome$metrics$n_used) {
          shiny::tags$p(shiny::em(
            "Note: periods with zero actual demand were excluded from MAPE."
          ))
        }
      )
    })
    
    # ---- Metrics table --------------------------------------------------------------------
    output$metrics_table <- shiny::renderTable({
      outcome <- holt_result()
      shiny::req(outcome$ok)
      metrics <- outcome$metrics
      data.frame(
        Metric = c("MAD", "MSE", "MAPE (%)"),
        Value  = c(
          format_number(metrics$MAD, 4),
          format_number(metrics$MSE, 4),
          format_number(metrics$MAPE, 4)
        ),
        stringsAsFactors = FALSE
      )
    })
    
    # ---- Plot ---------------------------------------------------------------------------------
    output$forecast_plot <- shiny::renderPlot({
      outcome <- holt_result()
      shiny::req(outcome$ok)
      holt_plot(
        outcome$demand,
        outcome$fit$forecasts,
        outcome$fit$future_forecasts,
        paste0("Holt's Method (alpha = ", format_number(outcome$alpha, 3),
               ", beta = ", format_number(outcome$beta, 3), ")")
      )
    })
    
    # ---- Future forecasts ------------------------------------------------------------------------
    output$future_table <- shiny::renderTable({
      outcome <- holt_result()
      shiny::req(outcome$ok)
      future_forecasts <- outcome$fit$future_forecasts
      data.frame(
        "Period"   = length(outcome$demand) + seq_along(future_forecasts),
        "Forecast" = future_forecasts,
        check.names = FALSE
      )
    }, digits = 2)
    
    # ---- Detailed results (period 0 shows the initial values) --------------------------------------
    output$detail_table <- shiny::renderTable({
      outcome <- holt_result()
      shiny::req(outcome$ok)
      data.frame(
        "Period"      = 0:length(outcome$demand),
        "Actual"      = c(NA, outcome$demand),
        "Level (L_t)" = c(outcome$l0, outcome$fit$levels),
        "Trend (T_t)" = c(outcome$t0, outcome$fit$trends),
        "Forecast"    = c(NA, outcome$fit$forecasts),
        "Error (e_t)" = c(NA, outcome$demand - outcome$fit$forecasts),
        check.names   = FALSE
      )
    }, na = "-", digits = 2)
  })
}