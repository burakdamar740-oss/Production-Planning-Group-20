# =============================================================================
# Production Planning & Forecasting
# Main application: navigation, shared data upload, and module wiring.
#
# Required packages: shiny, readxl, forecast, ggplot2
#   install.packages(c("shiny", "readxl", "forecast", "ggplot2"))
#
# Module files. Keep them in the "R/" subfolder next to app.R: Shiny sources
# every file in "R/" automatically (Shiny 1.5.0 or newer).
#   - Learning Curve module   -> learning_ui(), learning_server()
#   - R/moving_average.R      -> moving_average_ui(), moving_average_server()
#   - R/wma.R                 -> wma_ui(), wma_server()
#   - R/ses.R                 -> ses_ui(), ses_server()
#   - R/holt.R                -> holt_ui(), holt_server()
#   - R/seasonal.R            -> seasonal_ui(), seasonal_server()
#   - R/arima.R               -> arima_ui(), arima_server()
# =============================================================================

library(shiny)
library(readxl)
library(forecast)
library(ggplot2)

options(shiny.maxRequestSize = 30 * 1024^2)   # allow uploads up to 30 MB

# ---- Load module files -----------------------------------------------------------
# Files in "R/" are loaded automatically by Shiny. As a fallback, module files
# that sit next to app.R are sourced here.
for (module_file in c("learning_module.R", "moving_average.R", "wma.R",
                      "ses.R", "holt.R", "seasonal.R", "arima.R")) {
  if (file.exists(module_file)) source(module_file, local = TRUE)
}

required_functions <- c("learning_ui", "learning_server",
                        "moving_average_ui", "moving_average_server",
                        "wma_ui", "wma_server",
                        "ses_ui", "ses_server",
                        "holt_ui", "holt_server",
                        "seasonal_ui", "seasonal_server",
                        "arima_ui", "arima_server")
missing_functions <- required_functions[
  !vapply(required_functions, function(fn) exists(fn, mode = "function"), logical(1))
]
if (length(missing_functions) > 0) {
  stop("Missing module function(s): ", paste(missing_functions, collapse = ", "),
       ". Make sure the module files are in the same folder as app.R ",
       "(or in an 'R/' subfolder) and that the file names match the 'source' list above.")
}


# =============================================================================
# DATA-READING HELPERS
# =============================================================================

# Converts character columns that actually hold numbers (including decimal
# commas such as "12,5") into numeric columns. Other columns are left untouched.
coerce_numeric_columns <- function(df) {
  for (col_name in names(df)) {
    col_vals <- df[[col_name]]
    if (!is.character(col_vals)) next
    
    trimmed    <- trimws(col_vals)
    non_empty  <- !is.na(trimmed) & trimmed != ""
    if (!any(non_empty)) next
    
    direct <- suppressWarnings(as.numeric(trimmed))
    if (all(!is.na(direct[non_empty]))) {
      df[[col_name]] <- direct
    } else if (all(grepl("^-?[0-9]+(,[0-9]+)?$", trimmed[non_empty]))) {
      df[[col_name]] <- suppressWarnings(as.numeric(sub(",", ".", trimmed, fixed = TRUE)))
    }
  }
  df
}

# Reads a CSV file. Detects the separator (comma, semicolon or tab) and falls back
# to Latin-1 encoding if the file is not valid UTF-8.
read_csv_file <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  if (length(first_line) == 0 || !nzchar(gsub("[[:space:]]", "", first_line, useBytes = TRUE))) {
    stop("the file is empty.")
  }
  separators <- c(",", ";", "\t")
  counts <- vapply(separators, function(s) {
    nchar(gsub(paste0("[^", s, "]"), "", first_line, useBytes = TRUE), type = "bytes")
  }, numeric(1))
  sep <- if (max(counts) == 0) "," else separators[which.max(counts)]
  
  read_with_encoding <- function(enc) {
    utils::read.csv(path, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = c("", "NA", "N/A"),
                    fileEncoding = enc, strip.white = TRUE)
  }
  
  df <- tryCatch(
    withCallingHandlers(
      read_with_encoding("UTF-8-BOM"),
      warning = function(w) {
        if (grepl("invalid input|multibyte|embedded nul", conditionMessage(w))) {
          stop("encoding problem")
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) read_with_encoding("latin1")
  )
  coerce_numeric_columns(df)
}

# Reads one sheet of an Excel file (.xlsx or .xls).
read_excel_file <- function(path, sheet = NULL) {
  sheets <- readxl::excel_sheets(path)
  if (length(sheets) == 0) stop("the workbook has no sheets.")
  if (is.null(sheet) || !(sheet %in% sheets)) sheet <- sheets[1]
  df <- suppressMessages(suppressWarnings(
    readxl::read_excel(path, sheet = sheet)
  ))
  coerce_numeric_columns(as.data.frame(df, stringsAsFactors = FALSE))
}


# =============================================================================
# USER INTERFACE
# =============================================================================

# Shared data panel: shown above every tab except "Learning Curve".
data_panel <- conditionalPanel(
  condition = "input.main_nav !== undefined && input.main_nav !== 'learning'",
  div(
    class = "container-fluid",
    wellPanel(
      h4("Data"),
      fluidRow(
        column(
          4,
          fileInput("data_file", "Upload a data file (.csv, .xlsx, .xls)",
                    accept = c(".csv", ".xlsx", ".xls"))
        ),
        column(3, uiOutput("sheet_ui")),
        column(3, uiOutput("column_ui"))
      ),
      uiOutput("data_status")
    )
  )
)

ui <- navbarPage(
  title  = "Production Planning & Forecasting",
  id     = "main_nav",
  header = data_panel,
  
  tabPanel("Learning Curve", learning_ui("learning"), value = "learning"),
  
  navbarMenu(
    "Forecasting",
    tabPanel("Moving Average", moving_average_ui("moving_average")),
    tabPanel("Weighted Moving Average", wma_ui("wma")),
    tabPanel("Simple Exponential Smoothing", ses_ui("ses")),
    tabPanel("Holt's Trend Method", holt_ui("holt")),
    tabPanel("Seasonal & Winters' Methods", seasonal_ui("seasonal")),
    tabPanel("Box-Jenkins (ARIMA)", arima_ui("arima"))
  )
)


# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {
  
  # ---- File type ------------------------------------------------------------
  file_extension <- reactive({
    req(input$data_file)
    tolower(tools::file_ext(input$data_file$name))
  })
  
  # ---- Sheet selector (Excel files only) -------------------------------------
  output$sheet_ui <- renderUI({
    req(input$data_file)
    req(file_extension() %in% c("xlsx", "xls"))
    sheets <- tryCatch(readxl::excel_sheets(input$data_file$datapath),
                       error = function(e) character(0))
    req(length(sheets) > 0)
    selectInput("sheet", "Sheet", choices = sheets, selected = sheets[1])
  })
  
  # ---- Uploaded data ---------------------------------------------------------
  uploaded_data <- reactive({
    req(input$data_file)
    ext <- file_extension()
    validate(need(
      ext %in% c("csv", "xlsx", "xls"),
      "Unsupported file type. Please upload a .csv, .xlsx or .xls file."
    ))
    
    path <- input$data_file$datapath
    result <- tryCatch(
      list(data = if (ext == "csv") read_csv_file(path) else read_excel_file(path, input$sheet)),
      error = function(e) list(error = conditionMessage(e))
    )
    
    validate(need(
      is.null(result$error),
      paste0("The file could not be read: ", result$error)
    ))
    df <- result$data
    validate(need(
      is.data.frame(df) && nrow(df) > 0 && ncol(df) > 0,
      "The uploaded file contains no data rows."
    ))
    df
  })
  
  # ---- Numeric columns only --------------------------------------------------
  numeric_columns <- reactive({
    df <- uploaded_data()
    names(df)[vapply(df, is.numeric, logical(1))]
  })
  
  output$column_ui <- renderUI({
    cols <- tryCatch(numeric_columns(), error = function(e) NULL)
    req(length(cols) > 0)
    selectInput("series_col", "Column to forecast (numeric only)",
                choices = cols, selected = cols[1])
  })
  
  # ---- Selected time series (reactive numeric vector) -------------------------
  selected_series <- reactive({
    df <- uploaded_data()
    validate(need(
      length(numeric_columns()) > 0,
      "The file has no numeric columns. Please upload a file with at least one numeric column."
    ))
    req(input$series_col)
    req(input$series_col %in% names(df))   # ignore a stale selection after a new upload
    validate(need(
      is.numeric(df[[input$series_col]]),
      "The selected column is not numeric. Please choose a numeric column."
    ))
    
    series_vec <- as.numeric(df[[input$series_col]])
    series_vec <- series_vec[is.finite(series_vec)]   # drops NA, NaN and Inf
    validate(need(
      length(series_vec) >= 3,
      "The selected column needs at least 3 valid numeric values."
    ))
    series_vec
  })
  
  # ---- Status line under the upload controls ---------------------------------
  output$data_status <- renderUI({
    req(input$data_file)
    df         <- uploaded_data()
    series_vec <- selected_series()
    
    n_removed <- nrow(df) - length(series_vec)
    removed_note <- if (n_removed > 0) {
      paste0(" (", n_removed, " missing or non-numeric value(s) removed)")
    } else {
      ""
    }
    
    div(
      class = "text-success",
      paste0("Loaded ", nrow(df), " rows and ", ncol(df), " columns. ",
             "Forecasting '", input$series_col, "' with ",
             length(series_vec), " observations", removed_note, ".")
    )
  })
  
  # ---- Notifications on upload ------------------------------------------------
  observeEvent(input$data_file, {
    outcome <- tryCatch({
      uploaded_data()
      NULL
    }, error = function(e) conditionMessage(e))
    
    if (is.null(outcome)) {
      showNotification("File uploaded successfully.", type = "message", duration = 4)
    } else if (nzchar(outcome)) {
      showNotification(outcome, type = "error", duration = 8)
    }
  })
  
  # ---- Module servers -----------------------------------------------------------
  learning_server("learning")
  moving_average_server("moving_average", selected_series)
  wma_server("wma", selected_series)
  ses_server("ses", selected_series)
  holt_server("holt", selected_series)
  seasonal_server("seasonal", selected_series)
  arima_server("arima", selected_series)
}

shinyApp(ui, server)