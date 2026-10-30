# =============================================================================
# MOVING AVERAGE FORECASTING - SHINY MODULE
#
# Defines pure calculation functions and the Shiny module:
#   moving_average_ui(id)
#   moving_average_server(id, selected_series)
#
# `selected_series` must be a reactive that returns a numeric vector.
# Nothing is executed when this file is sourced.
# =============================================================================


# ----------------------------------------------------------------------------
# Validation helpers
# ----------------------------------------------------------------------------

# Checks that the demand data is a usable numeric vector.
ma_validate_demand <- function(demand_vec) {
  if (!is.numeric(demand_vec)) {
    stop("Demand data must be numeric.")
  }
  if (length(demand_vec) < 2) {
    stop("Demand data must contain at least 2 observations.")
  }
  if (all(is.na(demand_vec))) {
    stop("Demand data cannot be entirely missing (NA).")
  }
  invisible(TRUE)
}

# Checks that a window size is a whole number, >= 1 and < number of observations.
ma_validate_n <- function(n_value, data_length, arg_name = "N") {
  if (!is.numeric(n_value) || length(n_value) != 1 || is.na(n_value)) {
    stop(sprintf("%s must be a single number.", arg_name))
  }
  if (n_value != floor(n_value)) {
    stop(sprintf("%s must be a whole number (entered: %s).", arg_name, n_value))
  }
  if (n_value < 1) {
    stop(sprintf("%s must be at least 1 (entered: %s).", arg_name, n_value))
  }
  if (n_value >= data_length) {
    stop(sprintf(
      "%s must be smaller than the number of observations (%d); entered: %s.",
      arg_name, data_length, n_value
    ))
  }
  invisible(TRUE)
}


# ----------------------------------------------------------------------------
# Forecasting
# ----------------------------------------------------------------------------

# MA(N) forecasts for every period.
# Forecast for period t = mean of actual demand in periods t-N, ..., t-1.
# The demand of period t itself is never used. The first N periods have no
# forecast (NA).
ma_forecast <- function(demand_vec, n_value) {
  ma_validate_demand(demand_vec)
  data_length <- length(demand_vec)
  ma_validate_n(n_value, data_length)
  
  forecast_vec <- rep(NA_real_, data_length)
  for (period_idx in (n_value + 1):data_length) {
    window_vals <- demand_vec[(period_idx - n_value):(period_idx - 1)]
    forecast_vec[period_idx] <- mean(window_vals)
  }
  forecast_vec
}

# One-step-ahead forecast: mean of the last N observations.
ma_next_forecast <- function(demand_vec, n_value) {
  ma_validate_demand(demand_vec)
  ma_validate_n(n_value, length(demand_vec))
  mean(utils::tail(demand_vec, n_value))
}


# ----------------------------------------------------------------------------
# Error metrics
# ----------------------------------------------------------------------------

# MAD, MSE and MAPE.
# - Only periods where both actual and forecast are available are used.
# - Periods with zero actual demand are excluded from MAPE only
#   (division by zero), and counted in `n_used_mape`.
ma_error_metrics <- function(actual_vec, forecast_vec) {
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
# Manual and automatic modes
# ----------------------------------------------------------------------------

# Manual mode: the user chooses N.
ma_manual <- function(demand_vec, n_value) {
  forecast_vec <- ma_forecast(demand_vec, n_value)
  list(
    mode          = "manual",
    n             = n_value,
    forecasts     = forecast_vec,
    next_forecast = ma_next_forecast(demand_vec, n_value),
    metrics       = ma_error_metrics(demand_vec, forecast_vec)
  )
}

# Automatic mode: tries N = 1, ..., max_n and keeps the N with the lowest value
# of the chosen criterion ("MAD", "MSE" or "MAPE"). Ties go to the smaller N.
ma_auto <- function(demand_vec, max_n, criterion = "MAD") {
  ma_validate_demand(demand_vec)
  ma_validate_n(max_n, length(demand_vec), arg_name = "Maximum N")
  
  criterion_type <- toupper(trimws(as.character(criterion)))
  if (length(criterion_type) != 1 ||
      !(criterion_type %in% c("MAD", "MSE", "MAPE"))) {
    stop("Criterion must be 'MAD', 'MSE' or 'MAPE'.")
  }
  
  n_candidates <- seq_len(max_n)
  mad_vals  <- numeric(length(n_candidates))
  mse_vals  <- numeric(length(n_candidates))
  mape_vals <- numeric(length(n_candidates))
  
  for (i in seq_along(n_candidates)) {
    candidate_forecast <- ma_forecast(demand_vec, n_candidates[i])
    candidate_metrics  <- ma_error_metrics(demand_vec, candidate_forecast)
    mad_vals[i]  <- candidate_metrics$MAD
    mse_vals[i]  <- candidate_metrics$MSE
    mape_vals[i] <- candidate_metrics$MAPE
  }
  
  comparison_table <- data.frame(
    N    = n_candidates,
    MAD  = mad_vals,
    MSE  = mse_vals,
    MAPE = mape_vals
  )
  
  criterion_vals <- comparison_table[[criterion_type]]
  if (all(is.na(criterion_vals))) {
    stop(sprintf("%s could not be computed for any N.", criterion_type))
  }
  
  best_row      <- which.min(criterion_vals)
  best_n        <- comparison_table$N[best_row]
  best_forecast <- ma_forecast(demand_vec, best_n)
  
  list(
    mode             = "auto",
    criterion        = criterion_type,
    best_n           = best_n,
    best_error       = criterion_vals[best_row],
    comparison_table = comparison_table,
    forecasts        = best_forecast,
    next_forecast    = ma_next_forecast(demand_vec, best_n),
    metrics          = ma_error_metrics(demand_vec, best_forecast)
  )
}


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------

# Actual demand vs. moving average forecast (base R graphics).
plot_ma <- function(demand_vec, forecast_vec,
                    plot_title = "Actual Demand vs. Moving Average Forecast") {
  y_range <- range(c(demand_vec, forecast_vec), na.rm = TRUE)
  plot(seq_along(demand_vec), demand_vec,
       type = "b", pch = 16, col = "blue", ylim = y_range,
       xlab = "Period", ylab = "Demand", main = plot_title)
  lines(seq_along(forecast_vec), forecast_vec,
        type = "b", pch = 17, col = "red", lty = 2)
  legend("topleft", legend = c("Actual demand", "Forecast"),
         col = c("blue", "red"), pch = c(16, 17), lty = c(1, 2), bty = "n")
}


# ----------------------------------------------------------------------------
# Shiny module: UI
# ----------------------------------------------------------------------------

moving_average_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::radioButtons(
        ns("mode"), "How should N be chosen?",
        choices  = c("Manual" = "manual", "Auto (search best N)" = "auto"),
        selected = "manual"
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'manual'", ns("mode")),
        shiny::numericInput(
          ns("n_value"), "N (number of past periods averaged)",
          value = 3, min = 1, step = 1
        )
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'auto'", ns("mode")),
        shiny::numericInput(
          ns("max_n"), "Maximum N to try",
          value = 5, min = 1, step = 1
        ),
        shiny::selectInput(
          ns("criterion"), "Select the best N by",
          choices = c("MAD", "MSE", "MAPE"), selected = "MAD"
        )
      ),
      shiny::helpText(
        "Results update automatically. The data comes from the upload ",
        "panel at the top of the page."
      )
    ),
    shiny::mainPanel(
      shiny::uiOutput(ns("summary")),
      shiny::h4("Error Metrics"),
      shiny::tableOutput(ns("metrics_table")),
      shiny::plotOutput(ns("forecast_plot"), height = "320px"),
      shiny::h4("Detailed Results"),
      shiny::tableOutput(ns("detail_table")),
      shiny::uiOutput(ns("comparison_ui"))
    )
  )
}


# ----------------------------------------------------------------------------
# Shiny module: server
# ----------------------------------------------------------------------------

# selected_series: a reactive returning a numeric vector (supplied by app.R).
moving_average_server <- function(id, selected_series) {
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
    
    # Main computation. Any error (including validation messages coming from
    # the data panel) is turned into a message instead of crashing the session.
    ma_result <- shiny::reactive({
      tryCatch({
        demand_vec <- selected_series()
        result <- if (identical(input$mode, "auto")) {
          ma_auto(demand_vec, input$max_n, input$criterion)
        } else {
          ma_manual(demand_vec, input$n_value)
        }
        result$demand <- demand_vec
        list(ok = TRUE, result = result)
      }, error = function(e) {
        error_msg <- conditionMessage(e)
        if (!nzchar(error_msg)) {
          error_msg <- "Upload a data file and select a numeric column in the panel above."
        }
        list(ok = FALSE, message = error_msg)
      })
    })
    
    # Summary of results
    output$summary <- shiny::renderUI({
      outcome <- ma_result()
      if (!outcome$ok) {
        return(shiny::div(class = "alert alert-warning", outcome$message))
      }
      
      result <- outcome$result
      if (result$mode == "manual") {
        heading <- "Manual mode"
        n_used  <- result$n
        n_label <- "Selected N: "
      } else {
        heading <- paste0("Auto mode (criterion: ", result$criterion, ")")
        n_used  <- result$best_n
        n_label <- "Best N: "
      }
      
      shiny::tagList(
        shiny::h4(heading),
        shiny::tags$p(shiny::strong(n_label), n_used),
        if (result$mode == "auto") {
          shiny::tags$p(
            shiny::strong(paste0("Lowest ", result$criterion, ": ")),
            format_number(result$best_error, 4)
          )
        },
        shiny::tags$p(
          shiny::strong("Next-period forecast: "),
          format_number(result$next_forecast)
        ),
        if (!is.na(result$metrics$n_used_mape) &&
            result$metrics$n_used_mape < result$metrics$n_used) {
          shiny::tags$p(shiny::em(
            "Note: periods with zero actual demand were excluded from MAPE."
          ))
        }
      )
    })
    
    # Metrics table for the chosen / best N
    output$metrics_table <- shiny::renderTable({
      outcome <- ma_result()
      shiny::req(outcome$ok)
      metrics <- outcome$result$metrics
      data.frame(
        Metric = c("MAD", "MSE", "MAPE (%)", "Periods used"),
        Value  = c(
          format_number(metrics$MAD, 4),
          format_number(metrics$MSE, 4),
          format_number(metrics$MAPE, 4),
          as.character(metrics$n_used)
        ),
        stringsAsFactors = FALSE
      )
    })
    
    # Actual vs. forecast plot
    output$forecast_plot <- shiny::renderPlot({
      outcome <- ma_result()
      shiny::req(outcome$ok)
      result <- outcome$result
      n_used <- if (result$mode == "manual") result$n else result$best_n
      plot_ma(result$demand, result$forecasts,
              paste0("Actual Demand vs. MA(", n_used, ") Forecast"))
    })
    
    # Detailed results: Period, Demand, Forecast, Error
    output$detail_table <- shiny::renderTable({
      outcome <- ma_result()
      shiny::req(outcome$ok)
      result <- outcome$result
      data.frame(
        Period   = seq_along(result$demand),
        Demand   = result$demand,
        Forecast = result$forecasts,
        Error    = result$demand - result$forecasts
      )
    }, na = "-", digits = 2)
    
    # Comparison of all tried N values (auto mode only)
    output$comparison_ui <- shiny::renderUI({
      outcome <- ma_result()
      if (!outcome$ok || outcome$result$mode != "auto") {
        return(NULL)
      }
      shiny::tagList(
        shiny::h4("All Tried Values of N"),
        shiny::tableOutput(session$ns("comparison_table"))
      )
    })
    
    output$comparison_table <- shiny::renderTable({
      outcome <- ma_result()
      shiny::req(outcome$ok, outcome$result$mode == "auto")
      outcome$result$comparison_table
    }, digits = 4)
  })
}