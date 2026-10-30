# =============================================================================
# SIMPLE EXPONENTIAL SMOOTHING (SES) - SHINY MODULE
#
# Defines pure calculation functions and the Shiny module:
#   ses_ui(id)
#   ses_server(id, selected_series)
#
# `selected_series` must be a reactive that returns a numeric vector.
# Nothing is executed when this file is sourced.
#
# Model (T = number of observations, L_0 = initial level):
#   Forecast for period t :  Yhat_t = L_(t-1)
#   Level update          :  L_t    = alpha * Y_t + (1 - alpha) * L_(t-1)
#   Next-period forecast  :  Yhat_(T+1) = L_T
# =============================================================================


# ----------------------------------------------------------------------------
# Validation helpers
# ----------------------------------------------------------------------------

# Checks that the demand data is a usable numeric vector.
ses_validate_demand <- function(demand_vec) {
  if (!is.numeric(demand_vec)) {
    stop("Demand data must be numeric.")
  }
  if (length(demand_vec) < 3) {
    stop("Simple exponential smoothing needs at least 3 observations.")
  }
  if (any(is.na(demand_vec))) {
    stop("Demand data cannot contain missing values (NA).")
  }
  invisible(TRUE)
}

# Checks that alpha is a single number with 0 < alpha <= 1.
ses_validate_alpha <- function(alpha_value) {
  if (!is.numeric(alpha_value) || length(alpha_value) != 1 || is.na(alpha_value)) {
    stop("Alpha must be a single number.")
  }
  if (alpha_value <= 0 || alpha_value > 1) {
    stop(sprintf("Alpha must be greater than 0 and at most 1 (entered: %s).", alpha_value))
  }
  invisible(TRUE)
}


# ----------------------------------------------------------------------------
# Initial level
# ----------------------------------------------------------------------------

# Initial level L_0 from a linear regression lm(y ~ time) fitted on the first
# `n_periods` observations (time = 1, 2, ...). The intercept is the fitted value
# at time 0, i.e. the level just before the first period.
ses_regression_l0 <- function(demand_vec, n_periods) {
  ses_validate_demand(demand_vec)
  if (!is.numeric(n_periods) || length(n_periods) != 1 || is.na(n_periods) ||
      n_periods != floor(n_periods) || n_periods < 2 || n_periods > length(demand_vec)) {
    stop(sprintf(
      "The number of regression periods must be a whole number between 2 and %d.",
      length(demand_vec)
    ))
  }
  
  period_idx <- seq_len(n_periods)
  fit_data   <- data.frame(y = demand_vec[period_idx], time = period_idx)
  fit_model  <- stats::lm(y ~ time, data = fit_data)
  unname(stats::coef(fit_model)[["(Intercept)"]])
}


# ----------------------------------------------------------------------------
# Forecasting
# ----------------------------------------------------------------------------

# Iterative SES for a given alpha and initial level.
# Returns:
#   forecasts     : one-step-ahead forecasts, Yhat_t = L_(t-1), t = 1..T
#   levels        : smoothed levels L_t, t = 1..T
#   next_forecast : Yhat_(T+1) = L_T
ses_fit <- function(demand_vec, alpha_value, l0_value) {
  ses_validate_demand(demand_vec)
  ses_validate_alpha(alpha_value)
  if (!is.numeric(l0_value) || length(l0_value) != 1 || !is.finite(l0_value)) {
    stop("The initial level must be a single finite number.")
  }
  
  data_length  <- length(demand_vec)
  forecast_vec <- numeric(data_length)
  level_vec    <- numeric(data_length)
  
  previous_level <- l0_value
  for (period_idx in seq_len(data_length)) {
    forecast_vec[period_idx] <- previous_level
    level_vec[period_idx]    <- alpha_value * demand_vec[period_idx] +
      (1 - alpha_value) * previous_level
    previous_level <- level_vec[period_idx]
  }
  
  list(
    forecasts     = forecast_vec,
    levels        = level_vec,
    next_forecast = level_vec[data_length]
  )
}


# ----------------------------------------------------------------------------
# Error metrics
# ----------------------------------------------------------------------------

# MAD, MSE and MAPE.
# - Only periods where both actual and forecast are available are used.
# - Periods with zero actual demand are excluded from MAPE only
#   (division by zero), and counted in `n_used_mape`.
ses_error_metrics <- function(actual_vec, forecast_vec) {
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
# Automatic alpha selection
# ----------------------------------------------------------------------------

# Finds the alpha in (0, 1] that minimizes the chosen criterion ("MAD", "MSE"
# or "MAPE") for a given initial level.
# The error surface is not guaranteed to have a single minimum, so a coarse grid
# (0.01, 0.02, ..., 1) is evaluated first and optimize() then refines the search
# in the neighbourhood of the best grid point. The better of the two results is
# returned.
ses_optimize_alpha <- function(demand_vec, l0_value, criterion = "MAD") {
  ses_validate_demand(demand_vec)
  
  criterion_type <- toupper(trimws(as.character(criterion)))
  if (length(criterion_type) != 1 ||
      !(criterion_type %in% c("MAD", "MSE", "MAPE"))) {
    stop("Criterion must be 'MAD', 'MSE' or 'MAPE'.")
  }
  
  # Value of the chosen criterion for a given alpha (Inf if not computable).
  objective_fn <- function(alpha_value) {
    fit_result <- ses_fit(demand_vec, alpha_value, l0_value)
    metrics    <- ses_error_metrics(demand_vec, fit_result$forecasts)
    criterion_value <- metrics[[criterion_type]]
    if (is.na(criterion_value)) Inf else criterion_value
  }
  
  alpha_grid <- round(seq(0.01, 1, by = 0.01), 2)
  grid_vals  <- vapply(alpha_grid, objective_fn, numeric(1))
  if (all(!is.finite(grid_vals))) {
    stop(sprintf("%s could not be computed for any alpha.", criterion_type))
  }
  
  best_grid_alpha <- alpha_grid[which.min(grid_vals)]
  search_lower    <- max(0.001, best_grid_alpha - 0.02)
  search_upper    <- min(1, best_grid_alpha + 0.02)
  refined         <- stats::optimize(objective_fn,
                                     interval = c(search_lower, search_upper),
                                     tol = 1e-6)
  
  candidate_alphas <- c(best_grid_alpha, refined$minimum)
  candidate_vals   <- c(min(grid_vals), refined$objective)
  best_idx         <- which.min(candidate_vals)
  
  list(
    alpha     = candidate_alphas[best_idx],
    criterion = criterion_type,
    error     = candidate_vals[best_idx]
  )
}


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------

# Actual demand, smoothed level, one-step-ahead forecasts, initial level and the
# next-period forecast (base R graphics).
ses_plot <- function(demand_vec, forecast_vec, level_vec, l0_value, next_forecast,
                     plot_title = "Simple Exponential Smoothing") {
  data_length <- length(demand_vec)
  y_range <- range(c(demand_vec, forecast_vec, level_vec, l0_value, next_forecast),
                   na.rm = TRUE)
  
  plot(seq_len(data_length), demand_vec,
       type = "b", pch = 16, col = "blue",
       xlim = c(0, data_length + 1), ylim = y_range,
       xlab = "Period (0 = initial level)", ylab = "Demand",
       main = plot_title)
  
  # Smoothed level L_t, starting from L_0 at period 0
  lines(0:data_length, c(l0_value, level_vec), col = "darkorange", lwd = 2)
  points(0, l0_value, pch = 15, cex = 1.3, col = "darkorange")
  
  # One-step-ahead forecasts for periods 1..T and the forecast for T + 1
  lines(seq_len(data_length + 1), c(forecast_vec, next_forecast),
        type = "b", pch = 17, col = "red", lty = 2)
  points(data_length + 1, next_forecast, pch = 18, cex = 1.8, col = "darkgreen")
  
  legend("topleft",
         legend = c("Actual demand", "Smoothed level (L_t)", "One-step forecast",
                    "Initial level (L_0)", "Next-period forecast"),
         col = c("blue", "darkorange", "red", "darkorange", "darkgreen"),
         pch = c(16, NA, 17, 15, 18),
         lty = c(1, 1, 2, NA, NA),
         lwd = c(1, 2, 1, 1, 1),
         bty = "n")
}


# ----------------------------------------------------------------------------
# Shiny module: UI
# ----------------------------------------------------------------------------

ses_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      # --- Alpha -------------------------------------------------------------
      shiny::radioButtons(
        ns("alpha_mode"), "How should alpha be chosen?",
        choices = c(
          "Manual"                    = "manual",
          "Auto (optimize alpha)"     = "auto"
        ),
        selected = "manual"
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'manual'", ns("alpha_mode")),
        shiny::numericInput(
          ns("alpha_num"), "Alpha (0 < alpha <= 1)",
          value = 0.3, min = 0.001, max = 1, step = 0.01
        ),
        shiny::sliderInput(
          ns("alpha_slider"), "Or slide to choose alpha",
          min = 0.001, max = 1, value = 0.3, step = 0.001
        )
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'auto'", ns("alpha_mode")),
        shiny::selectInput(
          ns("criterion"), "Minimize",
          choices = c("MAD", "MSE", "MAPE"), selected = "MAD"
        ),
        shiny::helpText(
          "Searches alpha between 0 and 1 with optimize() to minimize the ",
          "selected error metric."
        )
      ),
      
      shiny::hr(),
      
      # --- Initial level -------------------------------------------------------
      shiny::radioButtons(
        ns("init_method"), "Initial level (L0)",
        choices = c(
          "Linear regression fit (intercept)" = "regression",
          "Manual input"                      = "manual"
        ),
        selected = "regression"
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'regression'", ns("init_method")),
        shiny::numericInput(
          ns("reg_periods"), "Number of first periods used in the regression",
          value = 6, min = 2, step = 1
        ),
        shiny::helpText(
          "Fits lm(demand ~ time) on the first periods (time = 1, 2, ...). ",
          "The intercept, i.e. the fitted value at time 0, is used as L0."
        )
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'manual'", ns("init_method")),
        shiny::numericInput(
          ns("l0_manual"), "Initial level L0",
          value = NA, step = 1
        )
      ),
      
      shiny::helpText(
        "Forecast: Y-hat(t) = L(t-1).  Level: L(t) = alpha * Y(t) + (1 - alpha) * L(t-1).  ",
        "Next period: Y-hat(T+1) = L(T). ",
        "Results update automatically; the data comes from the upload panel at the top."
      )
    ),
    
    shiny::mainPanel(
      shiny::uiOutput(ns("summary")),
      shiny::h4("Error Metrics"),
      shiny::tableOutput(ns("metrics_table")),
      shiny::plotOutput(ns("forecast_plot"), height = "340px"),
      shiny::h4("Detailed Results"),
      shiny::tableOutput(ns("detail_table"))
    )
  )
}


# ----------------------------------------------------------------------------
# Shiny module: server
# ----------------------------------------------------------------------------

# selected_series: a reactive returning a numeric vector (supplied by app.R).
ses_server <- function(id, selected_series) {
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
    
    # ---- Keep the alpha number box and the alpha slider in sync ---------------
    shiny::observeEvent(input$alpha_slider, {
      slider_value <- input$alpha_slider
      number_value <- input$alpha_num
      if (is.null(number_value) || is.na(number_value) ||
          abs(number_value - slider_value) > 1e-9) {
        shiny::updateNumericInput(session, "alpha_num", value = slider_value)
      }
    }, ignoreInit = TRUE)
    
    shiny::observeEvent(input$alpha_num, {
      number_value <- input$alpha_num
      if (!is.null(number_value) && !is.na(number_value) &&
          number_value >= 0.001 && number_value <= 1 &&
          abs(number_value - input$alpha_slider) > 1e-9) {
        shiny::updateSliderInput(session, "alpha_slider", value = round(number_value, 3))
      }
    }, ignoreInit = TRUE)
    
    # ---- Data check ---------------------------------------------------------------
    # Shows a friendly message instead of an error when the data is not usable.
    checked_series <- shiny::reactive({
      demand_vec <- selected_series()
      shiny::validate(
        shiny::need(is.numeric(demand_vec), "The data must be numeric."),
        shiny::need(length(demand_vec) >= 3,
                    "Simple exponential smoothing needs at least 3 observations."),
        shiny::need(!anyNA(demand_vec), "The data cannot contain missing values.")
      )
      demand_vec
    })
    
    # ---- Main computation ---------------------------------------------------------
    # Any problem (validation message, input error, calculation error) becomes a
    # message on the screen instead of crashing the session.
    ses_result <- shiny::reactive({
      tryCatch({
        demand_vec <- checked_series()
        n_obs      <- length(demand_vec)
        
        # Initial level L0
        if (identical(input$init_method, "manual")) {
          l0_input <- input$l0_manual
          shiny::validate(shiny::need(
            is.numeric(l0_input) && length(l0_input) == 1 && !is.na(l0_input),
            "Enter a numeric value for the initial level L0."
          ))
          l0_value     <- l0_input
          init_label   <- "Manual input"
        } else {
          reg_periods <- input$reg_periods
          shiny::validate(shiny::need(
            is.numeric(reg_periods) && length(reg_periods) == 1 && !is.na(reg_periods) &&
              reg_periods == floor(reg_periods) && reg_periods >= 2 && reg_periods <= n_obs,
            sprintf("The number of regression periods must be a whole number between 2 and %d.",
                    n_obs)
          ))
          l0_value   <- ses_regression_l0(demand_vec, reg_periods)
          init_label <- sprintf("Linear regression fit (first %d periods)", reg_periods)
        }
        
        # Alpha
        optimization <- NULL
        if (identical(input$alpha_mode, "auto")) {
          optimization <- ses_optimize_alpha(demand_vec, l0_value, input$criterion)
          alpha_value  <- optimization$alpha
        } else {
          alpha_input <- input$alpha_num
          shiny::validate(shiny::need(
            is.numeric(alpha_input) && length(alpha_input) == 1 && !is.na(alpha_input) &&
              alpha_input > 0 && alpha_input <= 1,
            "Alpha must be greater than 0 and at most 1."
          ))
          alpha_value <- alpha_input
        }
        
        fit_result <- ses_fit(demand_vec, alpha_value, l0_value)
        metrics    <- ses_error_metrics(demand_vec, fit_result$forecasts)
        
        list(
          ok           = TRUE,
          demand       = demand_vec,
          alpha        = alpha_value,
          alpha_mode   = input$alpha_mode,
          optimization = optimization,
          l0           = l0_value,
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
    
    # ---- Summary ----------------------------------------------------------------------
    output$summary <- shiny::renderUI({
      outcome <- ses_result()
      if (!outcome$ok) {
        return(shiny::div(class = "alert alert-warning", outcome$message))
      }
      
      mode_text <- if (outcome$alpha_mode == "auto") {
        paste0("Auto mode: alpha optimized to minimize ", outcome$optimization$criterion)
      } else {
        "Manual mode: alpha chosen by the user"
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
      
      shiny::tagList(
        shiny::h4(mode_text),
        shiny::fluidRow(
          make_card(
            if (outcome$alpha_mode == "auto") "Optimized alpha" else "Alpha",
            format_number(outcome$alpha, 4)
          ),
          make_card("Initial level (L0)", format_number(outcome$l0)),
          make_card("Next-period forecast", format_number(outcome$fit$next_forecast))
        ),
        shiny::tags$p(shiny::strong("Initial level method: "), outcome$init_label),
        if (outcome$alpha_mode == "auto") {
          shiny::tags$p(
            shiny::strong(paste0("Lowest ", outcome$optimization$criterion, ": ")),
            format_number(outcome$optimization$error, 4)
          )
        },
        if (!is.na(outcome$metrics$n_used_mape) &&
            outcome$metrics$n_used_mape < outcome$metrics$n_used) {
          shiny::tags$p(shiny::em(
            "Note: periods with zero actual demand were excluded from MAPE."
          ))
        }
      )
    })
    
    # ---- Metrics table ------------------------------------------------------------------
    output$metrics_table <- shiny::renderTable({
      outcome <- ses_result()
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
    
    # ---- Plot -----------------------------------------------------------------------------
    output$forecast_plot <- shiny::renderPlot({
      outcome <- ses_result()
      shiny::req(outcome$ok)
      ses_plot(
        outcome$demand,
        outcome$fit$forecasts,
        outcome$fit$levels,
        outcome$l0,
        outcome$fit$next_forecast,
        paste0("Simple Exponential Smoothing (alpha = ", format_number(outcome$alpha, 3), ")")
      )
    })
    
    # ---- Detailed results --------------------------------------------------------------------
    output$detail_table <- shiny::renderTable({
      outcome <- ses_result()
      shiny::req(outcome$ok)
      data.frame(
        "Period"         = seq_along(outcome$demand),
        "Actual Demand"  = outcome$demand,
        "Forecast"       = outcome$fit$forecasts,
        "Smoothed Level" = outcome$fit$levels,
        "Error"          = outcome$demand - outcome$fit$forecasts,
        check.names      = FALSE
      )
    }, na = "-", digits = 2)
  })
}