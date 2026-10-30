# =============================================================================
# WEIGHTED MOVING AVERAGE (WMA) - SHINY MODULE
#
# Defines pure calculation functions and the Shiny module:
#   wma_ui(id)
#   wma_server(id, selected_series)
#
# `selected_series` must be a reactive that returns a numeric vector.
# Nothing is executed when this file is sourced.
#
# Weight convention used throughout this file ("lag order"):
#   weights_lag[1] is applied to the most recent observation (period t-1),
#   weights_lag[2] to period t-2, ..., weights_lag[N] to period t-N.
# Weights are always normalized so that they sum to 1.
# =============================================================================


# ----------------------------------------------------------------------------
# Validation helpers
# ----------------------------------------------------------------------------

# Checks that the demand data is a usable numeric vector.
wma_validate_demand <- function(demand_vec) {
  if (!is.numeric(demand_vec)) {
    stop("Demand data must be numeric.")
  }
  if (length(demand_vec) < 2) {
    stop("Demand data must contain at least 2 observations.")
  }
  if (any(is.na(demand_vec))) {
    stop("Demand data cannot contain missing values (NA).")
  }
  invisible(TRUE)
}

# Checks that N is a whole number, >= 1 and smaller than the number of observations.
wma_validate_n <- function(n_value, data_length) {
  if (!is.numeric(n_value) || length(n_value) != 1 || is.na(n_value)) {
    stop("N must be a single number.")
  }
  if (n_value != floor(n_value)) {
    stop(sprintf("N must be a whole number (entered: %s).", n_value))
  }
  if (n_value < 1) {
    stop(sprintf("N must be at least 1 (entered: %s).", n_value))
  }
  if (n_value >= data_length) {
    stop(sprintf(
      "N must be smaller than the number of observations (%d); entered: %s.",
      data_length, n_value
    ))
  }
  invisible(TRUE)
}

# Validates raw weights and rescales them so that they sum to 1.
wma_normalize_weights <- function(raw_weights) {
  if (!is.numeric(raw_weights) || length(raw_weights) == 0) {
    stop("Weights must be numbers.")
  }
  if (any(is.na(raw_weights)) || any(!is.finite(raw_weights))) {
    stop("Weights must be finite numbers.")
  }
  if (any(raw_weights < 0)) {
    stop("Weights cannot be negative.")
  }
  weight_total <- sum(raw_weights)
  if (weight_total <= 0) {
    stop("The weights must not all be zero.")
  }
  raw_weights / weight_total
}


# ----------------------------------------------------------------------------
# Weight construction
# ----------------------------------------------------------------------------

# Choices for the predefined weighting schemes (label = value).
wma_scheme_choices <- function() {
  c(
    "Equal weights (simple moving average)" = "equal",
    "Linear declining (rank-sum)"           = "linear",
    "Reciprocal rank"                       = "reciprocal",
    "Geometric decay"                       = "geometric"
  )
}

# Parses weights typed by the user (separated by commas, semicolons or spaces;
# use a dot for decimals) and returns normalized weights in lag order.
#   most_recent_first = TRUE : the first number belongs to period t-1
#   most_recent_first = FALSE: the first number belongs to period t-N
wma_parse_weights <- function(weights_text, n_value, most_recent_first = TRUE) {
  if (is.null(weights_text) || length(weights_text) != 1 || is.na(weights_text) ||
      !nzchar(trimws(weights_text))) {
    stop("Enter the weights, separated by commas (for example: 0.5, 0.3, 0.2).")
  }
  
  parts <- strsplit(trimws(weights_text), "[,;[:space:]]+")[[1]]
  parts <- parts[nzchar(parts)]
  raw_weights <- suppressWarnings(as.numeric(parts))
  if (any(is.na(raw_weights))) {
    stop("Weights must be numbers separated by commas (use a dot for decimals).")
  }
  if (length(raw_weights) != n_value) {
    stop(sprintf("You entered %d weight(s), but N is %d. Enter exactly N weights.",
                 length(raw_weights), n_value))
  }
  
  normalized <- wma_normalize_weights(raw_weights)
  if (most_recent_first) normalized else rev(normalized)
}

# Builds normalized weights (lag order) for a predefined scheme.
#   equal      : 1, 1, ..., 1
#   linear     : N, N-1, ..., 1                (rank-sum weights)
#   reciprocal : 1, 1/2, ..., 1/N              (reciprocal rank weights)
#   geometric  : 1, r, r^2, ..., r^(N-1)       (r = decay ratio, 0 < r <= 1)
wma_scheme_weights <- function(n_value, scheme, decay_ratio = 0.7) {
  if (!is.numeric(n_value) || length(n_value) != 1 || is.na(n_value) || n_value < 1) {
    stop("N must be a whole number of at least 1.")
  }
  positions <- seq_len(n_value)
  
  raw_weights <- switch(
    scheme,
    equal      = rep(1, n_value),
    linear     = rev(positions),
    reciprocal = 1 / positions,
    geometric  = {
      if (!is.numeric(decay_ratio) || length(decay_ratio) != 1 || is.na(decay_ratio) ||
          decay_ratio <= 0 || decay_ratio > 1) {
        stop("The decay ratio must be greater than 0 and at most 1.")
      }
      decay_ratio^(positions - 1)
    },
    stop("Unknown weighting scheme.")
  )
  wma_normalize_weights(raw_weights)
}


# ----------------------------------------------------------------------------
# Forecasting
# ----------------------------------------------------------------------------

# WMA forecasts for every period.
# Forecast for period t = sum of weight * demand over periods t-N, ..., t-1.
# The demand of period t itself is never used. The first N periods have no
# forecast (NA).
wma_forecast <- function(demand_vec, weights_lag) {
  n_value <- length(weights_lag)
  wma_validate_demand(demand_vec)
  wma_validate_n(n_value, length(demand_vec))
  
  data_length          <- length(demand_vec)
  weights_oldest_first <- rev(weights_lag)
  forecast_vec         <- rep(NA_real_, data_length)
  
  for (period_idx in (n_value + 1):data_length) {
    window_vals <- demand_vec[(period_idx - n_value):(period_idx - 1)]
    forecast_vec[period_idx] <- sum(weights_oldest_first * window_vals)
  }
  forecast_vec
}

# One-step-ahead forecast: weighted average of the last N observations.
wma_next_forecast <- function(demand_vec, weights_lag) {
  n_value <- length(weights_lag)
  wma_validate_demand(demand_vec)
  wma_validate_n(n_value, length(demand_vec))
  sum(rev(weights_lag) * utils::tail(demand_vec, n_value))
}


# ----------------------------------------------------------------------------
# Error metrics
# ----------------------------------------------------------------------------

# MAD, MSE and MAPE.
# - Only periods where both actual and forecast are available are used.
# - Periods with zero actual demand are excluded from MAPE only
#   (division by zero), and counted in `n_used_mape`.
wma_error_metrics <- function(actual_vec, forecast_vec) {
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
# Running a WMA and searching for the best decay ratio
# ----------------------------------------------------------------------------

# Runs a WMA with the given weights (lag order, normalized).
wma_run <- function(demand_vec, weights_lag) {
  forecast_vec <- wma_forecast(demand_vec, weights_lag)
  list(
    weights       = weights_lag,
    forecasts     = forecast_vec,
    next_forecast = wma_next_forecast(demand_vec, weights_lag),
    metrics       = wma_error_metrics(demand_vec, forecast_vec)
  )
}

# Parameter search: for a fixed N, tries geometric weights
# 1, r, r^2, ..., r^(N-1) for every decay ratio r in `decay_grid` and keeps the
# ratio with the lowest value of the chosen criterion ("MAD", "MSE" or "MAPE").
# Ties go to the smallest ratio. r = 1 is equivalent to a simple moving average.
wma_search_decay <- function(demand_vec, n_value, criterion = "MAD",
                             decay_grid = round(seq(0.05, 1, by = 0.05), 2)) {
  wma_validate_demand(demand_vec)
  wma_validate_n(n_value, length(demand_vec))
  
  criterion_type <- toupper(trimws(as.character(criterion)))
  if (length(criterion_type) != 1 ||
      !(criterion_type %in% c("MAD", "MSE", "MAPE"))) {
    stop("Criterion must be 'MAD', 'MSE' or 'MAPE'.")
  }
  
  mad_vals  <- numeric(length(decay_grid))
  mse_vals  <- numeric(length(decay_grid))
  mape_vals <- numeric(length(decay_grid))
  
  for (i in seq_along(decay_grid)) {
    candidate_weights  <- wma_scheme_weights(n_value, "geometric", decay_grid[i])
    candidate_forecast <- wma_forecast(demand_vec, candidate_weights)
    candidate_metrics  <- wma_error_metrics(demand_vec, candidate_forecast)
    mad_vals[i]  <- candidate_metrics$MAD
    mse_vals[i]  <- candidate_metrics$MSE
    mape_vals[i] <- candidate_metrics$MAPE
  }
  
  search_table <- data.frame(
    "Decay ratio" = decay_grid,
    MAD           = mad_vals,
    MSE           = mse_vals,
    MAPE          = mape_vals,
    check.names   = FALSE
  )
  
  criterion_vals <- search_table[[criterion_type]]
  if (all(is.na(criterion_vals))) {
    stop(sprintf("%s could not be computed for any decay ratio.", criterion_type))
  }
  
  best_row   <- which.min(criterion_vals)
  best_decay <- decay_grid[best_row]
  
  list(
    criterion    = criterion_type,
    best_decay   = best_decay,
    best_error   = criterion_vals[best_row],
    weights      = wma_scheme_weights(n_value, "geometric", best_decay),
    search_table = search_table
  )
}


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------

# Actual demand vs. WMA forecast, with the next-period forecast marked
# (base R graphics).
wma_plot <- function(demand_vec, forecast_vec, next_forecast,
                     plot_title = "Actual Demand vs. Weighted Moving Average Forecast") {
  data_length <- length(demand_vec)
  y_range <- range(c(demand_vec, forecast_vec, next_forecast), na.rm = TRUE)
  
  plot(seq_len(data_length), demand_vec,
       type = "b", pch = 16, col = "blue",
       xlim = c(1, data_length + 1), ylim = y_range,
       xlab = "Period", ylab = "Demand", main = plot_title)
  lines(seq_len(data_length), forecast_vec,
        type = "b", pch = 17, col = "red", lty = 2)
  points(data_length + 1, next_forecast, pch = 18, cex = 1.8, col = "darkgreen")
  legend("topleft",
         legend = c("Actual demand", "WMA forecast", "Next-period forecast"),
         col = c("blue", "red", "darkgreen"),
         pch = c(16, 17, 18), lty = c(1, 2, NA), bty = "n")
}


# ----------------------------------------------------------------------------
# Shiny module: UI
# ----------------------------------------------------------------------------

wma_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::numericInput(
        ns("n_value"), "Window size N (number of past periods)",
        value = 3, min = 1, step = 1
      ),
      shiny::radioButtons(
        ns("weight_mode"), "How should the weights be chosen?",
        choices = c(
          "Manual weights"                       = "manual",
          "Predefined scheme"                    = "predefined",
          "Auto (search best decay ratio)"       = "auto"
        ),
        selected = "manual"
      ),
      
      # --- Manual weights ---------------------------------------------------
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'manual'", ns("weight_mode")),
        shiny::textInput(
          ns("weights_text"), "Weights (N numbers, comma-separated)",
          value = "0.5, 0.3, 0.2"
        ),
        shiny::radioButtons(
          ns("weight_order"), "Order of the weights you entered",
          choices = c(
            "Most recent period first (t-1, t-2, ...)" = "recent_first",
            "Oldest period first (t-N, ..., t-1)"      = "oldest_first"
          ),
          selected = "recent_first"
        ),
        shiny::helpText(
          "Weights are normalized automatically so that they sum to 1. ",
          "Use a dot for decimals."
        )
      ),
      
      # --- Predefined scheme ------------------------------------------------
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'predefined'", ns("weight_mode")),
        shiny::selectInput(
          ns("scheme"), "Weighting scheme",
          choices = wma_scheme_choices(), selected = "linear"
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'geometric'", ns("scheme")),
          shiny::numericInput(
            ns("decay_ratio"), "Decay ratio (0 < r <= 1)",
            value = 0.7, min = 0.01, max = 1, step = 0.05
          )
        ),
        shiny::helpText("Every scheme gives the most recent period the largest weight.")
      ),
      
      # --- Auto search ------------------------------------------------------
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'auto'", ns("weight_mode")),
        shiny::selectInput(
          ns("criterion"), "Select the best decay ratio by",
          choices = c("MAD", "MSE", "MAPE"), selected = "MAD"
        ),
        shiny::helpText(
          "For the chosen N, tries geometric weights 1, r, r^2, ... with r from ",
          "0.05 to 1 (step 0.05) and keeps the best one. r = 1 equals a simple ",
          "moving average."
        )
      ),
      
      shiny::helpText(
        "Results update automatically. The data comes from the upload panel ",
        "at the top of the page."
      )
    ),
    
    shiny::mainPanel(
      shiny::uiOutput(ns("summary")),
      shiny::h4("Normalized Weights"),
      shiny::tableOutput(ns("weights_table")),
      shiny::h4("Error Metrics"),
      shiny::tableOutput(ns("metrics_table")),
      shiny::plotOutput(ns("forecast_plot"), height = "320px"),
      shiny::h4("Detailed Results"),
      shiny::tableOutput(ns("detail_table")),
      shiny::uiOutput(ns("search_ui"))
    )
  )
}


# ----------------------------------------------------------------------------
# Shiny module: server
# ----------------------------------------------------------------------------

# selected_series: a reactive returning a numeric vector (supplied by app.R).
wma_server <- function(id, selected_series) {
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
    wma_result <- shiny::reactive({
      tryCatch({
        demand_vec <- selected_series()
        n_value    <- input$n_value
        wma_validate_demand(demand_vec)
        wma_validate_n(n_value, length(demand_vec))
        
        weight_mode  <- input$weight_mode
        search       <- NULL
        scheme_label <- NULL
        
        if (identical(weight_mode, "auto")) {
          search      <- wma_search_decay(demand_vec, n_value, input$criterion)
          weights_lag <- search$weights
        } else if (identical(weight_mode, "predefined")) {
          weights_lag  <- wma_scheme_weights(n_value, input$scheme, input$decay_ratio)
          scheme_choices <- wma_scheme_choices()
          scheme_label   <- names(scheme_choices)[scheme_choices == input$scheme]
        } else {
          weights_lag <- wma_parse_weights(
            input$weights_text, n_value,
            most_recent_first = identical(input$weight_order, "recent_first")
          )
        }
        
        run <- wma_run(demand_vec, weights_lag)
        list(
          ok           = TRUE,
          demand       = demand_vec,
          n            = n_value,
          mode         = weight_mode,
          scheme_label = scheme_label,
          search       = search,
          run          = run
        )
      }, error = function(e) {
        error_msg <- conditionMessage(e)
        if (!nzchar(error_msg)) {
          error_msg <- "Upload a data file and select a numeric column in the panel above."
        }
        list(ok = FALSE, message = error_msg)
      })
    })
    
    # Summary cards
    output$summary <- shiny::renderUI({
      outcome <- wma_result()
      if (!outcome$ok) {
        return(shiny::div(class = "alert alert-warning", outcome$message))
      }
      
      mode_text <- switch(
        outcome$mode,
        manual     = "Manual weights",
        predefined = paste0("Predefined scheme: ", outcome$scheme_label),
        auto       = paste0("Auto search (criterion: ", outcome$search$criterion, ")")
      )
      
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
      
      cards <- list(
        make_card("Window size N", outcome$n),
        make_card("Next-period forecast", format_number(outcome$run$next_forecast))
      )
      if (outcome$mode == "auto") {
        cards[[3]] <- make_card("Best decay ratio", format_number(outcome$search$best_decay))
      } else {
        cards[[3]] <- make_card("Periods used for errors", outcome$run$metrics$n_used)
      }
      
      shiny::tagList(
        shiny::h4(mode_text),
        do.call(shiny::fluidRow, cards),
        if (!is.na(outcome$run$metrics$n_used_mape) &&
            outcome$run$metrics$n_used_mape < outcome$run$metrics$n_used) {
          shiny::tags$p(shiny::em(
            "Note: periods with zero actual demand were excluded from MAPE."
          ))
        }
      )
    })
    
    # Normalized weights (always shown by lag, so there is no ambiguity)
    output$weights_table <- shiny::renderTable({
      outcome <- wma_result()
      shiny::req(outcome$ok)
      weights_lag <- outcome$run$weights
      data.frame(
        Period = paste0("t-", seq_along(weights_lag)),
        Weight = format(round(weights_lag, 4), nsmall = 4),
        stringsAsFactors = FALSE
      )
    })
    
    # Metrics table
    output$metrics_table <- shiny::renderTable({
      outcome <- wma_result()
      shiny::req(outcome$ok)
      metrics <- outcome$run$metrics
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
    
    # Actual vs. forecast plot
    output$forecast_plot <- shiny::renderPlot({
      outcome <- wma_result()
      shiny::req(outcome$ok)
      wma_plot(
        outcome$demand,
        outcome$run$forecasts,
        outcome$run$next_forecast,
        paste0("Actual Demand vs. WMA(", outcome$n, ") Forecast")
      )
    })
    
    # Detailed results: Period, Actual, Forecast, Error
    output$detail_table <- shiny::renderTable({
      outcome <- wma_result()
      shiny::req(outcome$ok)
      data.frame(
        Period   = seq_along(outcome$demand),
        Actual   = outcome$demand,
        Forecast = outcome$run$forecasts,
        Error    = outcome$demand - outcome$run$forecasts
      )
    }, na = "-", digits = 2)
    
    # Table of all decay ratios tried (auto mode only)
    output$search_ui <- shiny::renderUI({
      outcome <- wma_result()
      if (!outcome$ok || outcome$mode != "auto") {
        return(NULL)
      }
      shiny::tagList(
        shiny::h4("All Tried Decay Ratios"),
        shiny::tableOutput(session$ns("search_table"))
      )
    })
    
    output$search_table <- shiny::renderTable({
      outcome <- wma_result()
      shiny::req(outcome$ok, outcome$mode == "auto")
      outcome$search$search_table
    }, digits = 4)
  })
}
