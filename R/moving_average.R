# =============================================================================
# MOVING AVERAGE (HAREKETLI ORTALAMA) TAHMIN MODULU  -  Saf Base R
# =============================================================================
# Icerik:
#   1) ma_validate_demand()   : talep vektoru kontrolu
#   2) ma_validate_n()        : N degeri kontrolu
#   3) ma_forecast()          : MA(N) tahminleri (ilk N periyot NA)
#   4) ma_next_forecast()     : bir sonraki periyodun tahmini
#   5) ma_error_metrics()     : MAD, MSE, MAPE
#   6) ma_manual()            : Manuel mod (kullanici N secer)
#   7) ma_auto()              : Otomatik mod (1..max_n denenir, en iyi N secilir)
#   8) Ornek calistirma blogu
# =============================================================================


# -----------------------------------------------------------------------------
# 1) Talep vektoru kontrolu
# -----------------------------------------------------------------------------
ma_validate_demand <- function(demand_vec) {
  if (!is.numeric(demand_vec)) {
    stop("Hata: Talep verisi sayisal (numeric) bir vektor olmalidir.")
  }
  if (length(demand_vec) < 2) {
    stop("Hata: Talep verisi en az 2 gozlem icermelidir.")
  }
  if (all(is.na(demand_vec))) {
    stop("Hata: Talep verisinin tamami NA olamaz.")
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# 2) N degeri kontrolu
#    N: tam sayi, >= 1 ve veri sayisindan kucuk olmali
#    (en az bir periyot icin tahmin-gercek karsilastirmasi yapilabilsin diye)
# -----------------------------------------------------------------------------
ma_validate_n <- function(n_value, data_length, arg_name = "N") {
  if (!is.numeric(n_value) || length(n_value) != 1 || is.na(n_value)) {
    stop(sprintf("Hata: %s tek bir sayisal deger olmalidir.", arg_name))
  }
  if (n_value != floor(n_value)) {
    stop(sprintf("Hata: %s bir tam sayi olmalidir (girilen: %s).", arg_name, n_value))
  }
  if (n_value < 1) {
    stop(sprintf("Hata: %s en az 1 olmalidir (girilen: %s).", arg_name, n_value))
  }
  if (n_value >= data_length) {
    stop(sprintf(
      "Hata: %s, veri sayisindan (%d) kucuk olmalidir (girilen: %s).",
      arg_name, data_length, n_value
    ))
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# 3) MA(N) tahmin fonksiyonu
#    Periyot p icin tahmin = mean(demand[(p-N):(p-1)])
#    Gercek talep demand[p] tahmine KATILMAZ. Ilk N periyot NA olur.
# -----------------------------------------------------------------------------
ma_forecast <- function(demand_vec, n_value) {
  ma_validate_demand(demand_vec)
  data_length <- length(demand_vec)
  ma_validate_n(n_value, data_length)
  
  forecast_vec <- rep(NA_real_, data_length)
  
  for (period_idx in (n_value + 1):data_length) {
    window_vals <- demand_vec[(period_idx - n_value):(period_idx - 1)]
    forecast_vec[period_idx] <- mean(window_vals)   # pencerede NA varsa sonuc NA olur
  }
  
  forecast_vec
}


# -----------------------------------------------------------------------------
# 4) Bir sonraki periyodun (veri sonrasi) tahmini = son N gercek talebin ortalamasi
# -----------------------------------------------------------------------------
ma_next_forecast <- function(demand_vec, n_value) {
  ma_validate_demand(demand_vec)
  ma_validate_n(n_value, length(demand_vec))
  mean(utils::tail(demand_vec, n_value))
}


# -----------------------------------------------------------------------------
# 5) Hata metrikleri: MAD, MSE, MAPE
#    - Sadece gercek VE tahmin degeri NA olmayan periyotlar kullanilir.
#    - MAPE'de gercek talebi 0 olan periyotlar (sifira bolme) disarida birakilir.
# -----------------------------------------------------------------------------
ma_error_metrics <- function(actual_vec, forecast_vec) {
  if (length(actual_vec) != length(forecast_vec)) {
    stop("Hata: Gercek talep ve tahmin vektorlerinin uzunlugu esit olmalidir.")
  }
  
  valid_idx <- !is.na(actual_vec) & !is.na(forecast_vec)
  
  if (!any(valid_idx)) {
    return(list(MAD = NA_real_, MSE = NA_real_, MAPE = NA_real_,
                n_used = 0L, n_used_mape = 0L))
  }
  
  actual_valid   <- actual_vec[valid_idx]
  forecast_valid <- forecast_vec[valid_idx]
  error_vec      <- actual_valid - forecast_valid
  
  mad_value <- mean(abs(error_vec))
  mse_value <- mean(error_vec^2)
  
  # MAPE icin sifir talepli periyotlari ele
  nonzero_idx <- actual_valid != 0
  if (any(nonzero_idx)) {
    mape_value <- mean(abs(error_vec[nonzero_idx]) / abs(actual_valid[nonzero_idx])) * 100
  } else {
    mape_value <- NA_real_
  }
  
  if (any(!nonzero_idx)) {
    warning("Uyari: Gercek talebi 0 olan periyotlar MAPE hesabina dahil edilmedi.")
  }
  
  list(
    MAD         = mad_value,
    MSE         = mse_value,
    MAPE        = mape_value,
    n_used      = as.integer(sum(valid_idx)),
    n_used_mape = as.integer(sum(nonzero_idx))
  )
}


# -----------------------------------------------------------------------------
# 6) MANUEL MOD
#    Kullanici N'i belirler; tahminler, sonraki periyot tahmini ve hatalar doner.
# -----------------------------------------------------------------------------
ma_manual <- function(demand_vec, n_value) {
  forecast_vec  <- ma_forecast(demand_vec, n_value)
  next_forecast <- ma_next_forecast(demand_vec, n_value)
  metrics_list  <- ma_error_metrics(demand_vec, forecast_vec)
  
  result_table <- data.frame(
    Period   = seq_along(demand_vec),
    Demand   = demand_vec,
    Forecast = forecast_vec,
    Error    = demand_vec - forecast_vec
  )
  
  list(
    mode          = "manual",
    n             = n_value,
    forecasts     = forecast_vec,
    next_forecast = next_forecast,
    metrics       = metrics_list,
    table         = result_table
  )
}


# -----------------------------------------------------------------------------
# 7) OTOMATIK MOD
#    1..max_n arasindaki tum N degerleri denenir; secilen kriteri (MAD/MSE/MAPE)
#    minimize eden N bulunur.
# -----------------------------------------------------------------------------
ma_auto <- function(demand_vec, max_n, criterion = "MAD") {
  ma_validate_demand(demand_vec)
  data_length <- length(demand_vec)
  ma_validate_n(max_n, data_length, arg_name = "max_n")
  
  criterion_type <- toupper(trimws(as.character(criterion)))
  allowed_criteria <- c("MAD", "MSE", "MAPE")
  if (length(criterion_type) != 1 || !(criterion_type %in% allowed_criteria)) {
    stop("Hata: criterion 'MAD', 'MSE' veya 'MAPE' olmalidir.")
  }
  
  n_candidates <- seq_len(max_n)
  mad_vals  <- numeric(length(n_candidates))
  mse_vals  <- numeric(length(n_candidates))
  mape_vals <- numeric(length(n_candidates))
  
  for (i in seq_along(n_candidates)) {
    fc_tmp <- ma_forecast(demand_vec, n_candidates[i])
    m_tmp  <- suppressWarnings(ma_error_metrics(demand_vec, fc_tmp))
    mad_vals[i]  <- m_tmp$MAD
    mse_vals[i]  <- m_tmp$MSE
    mape_vals[i] <- m_tmp$MAPE
  }
  
  comparison_table <- data.frame(
    N    = n_candidates,
    MAD  = mad_vals,
    MSE  = mse_vals,
    MAPE = mape_vals
  )
  
  criterion_vals <- comparison_table[[criterion_type]]
  if (all(is.na(criterion_vals))) {
    stop(sprintf("Hata: %s tum N degerleri icin hesaplanamadi.", criterion_type))
  }
  
  best_row   <- which.min(criterion_vals)   # esitlikte en kucuk N secilir
  best_n     <- comparison_table$N[best_row]
  best_error <- criterion_vals[best_row]
  
  # En iyi N icin tahminler ve metrikler
  best_forecasts <- ma_forecast(demand_vec, best_n)
  best_metrics   <- suppressWarnings(ma_error_metrics(demand_vec, best_forecasts))
  next_forecast  <- ma_next_forecast(demand_vec, best_n)
  
  list(
    mode             = "auto",
    criterion        = criterion_type,
    best_n           = best_n,
    best_error       = best_error,
    comparison_table = comparison_table,
    forecasts        = best_forecasts,
    next_forecast    = next_forecast,
    metrics          = best_metrics
  )
}


# -----------------------------------------------------------------------------
# 7b) Grafik fonksiyonu: Gercek talep ve MA tahmini
# -----------------------------------------------------------------------------
plot_ma <- function(demand_vec, forecast_vec, plot_title = "Gercek Talep vs MA Tahmini") {
  y_range <- range(c(demand_vec, forecast_vec), na.rm = TRUE)
  plot(seq_along(demand_vec), demand_vec, type = "b", pch = 16, col = "blue",
       ylim = y_range, xlab = "Periyot", ylab = "Talep", main = plot_title)
  lines(seq_along(forecast_vec), forecast_vec, type = "b", pch = 17, col = "red", lty = 2)
  legend("topleft", legend = c("Gercek", "Tahmin"),
         col = c("blue", "red"), pch = c(16, 17), lty = c(1, 2), bty = "n")
}


# =============================================================================
# 8) ORNEK CALISTIRMA BLOGU
# =============================================================================

demand_data <- c(120, 130, 125, 140, 135, 150, 145, 160, 155, 170, 165, 180)

cat("\n=====================================================\n")
cat(" MANUEL MOD  (N = 3)\n")
cat("=====================================================\n")

manual_result <- ma_manual(demand_data, n_value = 3)

cat("Secilen N          :", manual_result$n, "\n")
cat("Sonraki periyot tahmini:", round(manual_result$next_forecast, 2), "\n")
cat("MAD                :", round(manual_result$metrics$MAD, 4), "\n")
cat("MSE                :", round(manual_result$metrics$MSE, 4), "\n")
cat("MAPE (%)           :", round(manual_result$metrics$MAPE, 4), "\n")
cat("Hesaba katilan periyot sayisi:", manual_result$metrics$n_used, "\n\n")
cat("Talep / Tahmin tablosu:\n")
print(manual_result$table, row.names = FALSE)


cat("\n=====================================================\n")
cat(" OTOMATIK MOD  (max_n = 5, kriter = MAD)\n")
cat("=====================================================\n")

auto_result <- ma_auto(demand_data, max_n = 5, criterion = "MAD")

cat("Kriter             :", auto_result$criterion, "\n")
cat("En iyi N           :", auto_result$best_n, "\n")
cat("Minimum hata       :", round(auto_result$best_error, 4), "\n")
cat("Sonraki periyot tahmini (en iyi N):", round(auto_result$next_forecast, 2), "\n\n")
cat("Tum N denemeleri:\n")
print(auto_result$comparison_table, row.names = FALSE)


cat("\n=====================================================\n")
cat(" OTOMATIK MOD  (max_n = 5, kriter = MAPE)\n")
cat("=====================================================\n")

auto_result_mape <- ma_auto(demand_data, max_n = 5, criterion = "MAPE")
cat("En iyi N           :", auto_result_mape$best_n, "\n")
cat("Minimum MAPE (%)   :", round(auto_result_mape$best_error, 4), "\n")
cat("Sonraki periyot tahmini:", round(auto_result_mape$next_forecast, 2), "\n")

# =============================================================================
# 9) SHINY MODUL TANIMLARI
# =============================================================================

# UI Fonksiyonu
moving_average_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        textAreaInput(
          ns("demand_input"), 
          "Talep Verileri (Virgülle ayırın):", 
          value = "120, 130, 125, 140, 135, 150, 145, 160, 155, 170, 165, 180",
          rows = 3
        ),
        radioButtons(
          ns("mode"), 
          "Mod Seçimi:", 
          choices = c("Manuel" = "manual", "Otomatik" = "auto")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'manual'", ns("mode")),
          numericInput(ns("n_val"), "N Değeri:", value = 3, min = 1, step = 1)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'auto'", ns("mode")),
          numericInput(ns("max_n"), "Maksimum N:", value = 5, min = 2, step = 1),
          selectInput(ns("criterion"), "Kriter:", choices = c("MAD", "MSE", "MAPE"))
        ),
        actionButton(ns("calc_btn"), "Hesapla", class = "btn-primary")
      ),
      mainPanel(
        h4("Tahmin Sonuçları"),
        verbatimTextOutput(ns("summary_text")),
        plotOutput(ns("forecast_plot")),
        h4("Tablo"),
        tableOutput(ns("result_table"))
      )
    )
  )
}

# Server Fonksiyonu
moving_average_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    parsed_demand <- reactive({
      req(input$demand_input)
      vals <- as.numeric(trimws(unlist(strsplit(input$demand_input, ","))))
      validate(need(!any(is.na(vals)) && length(vals) >= 2, "Lütfen geçerli sayısal talep verileri girin."))
      vals
    })
    
    calc_results <- eventReactive(input$calc_btn, {
      demands <- parsed_demand()
      if (input$mode == "manual") {
        ma_manual(demands, input$n_val)
      } else {
        ma_auto(demands, input$max_n, input$criterion)
      }
    }, ignoreNULL = FALSE)
    
    output$summary_text <- renderPrint({
      res <- calc_results()
      if (res$mode == "manual") {
        cat("Mod: Manuel | N:", res$n, "\n")
        cat("Sonraki Periyot Tahmini:", round(res$next_forecast, 2), "\n")
        cat("MAD:", round(res$metrics$MAD, 4), "| MSE:", round(res$metrics$MSE, 4), "| MAPE:", paste0(round(res$metrics$MAPE, 2), "%\n"))
      } else {
        cat("Mod: Otomatik | En İyi N:", res$best_n, "(Kriter:", res$criterion, ")\n")
        cat("Sonraki Periyot Tahmini:", round(res$next_forecast, 2), "\n")
        cat("Minimum Hata:", round(res$best_error, 4), "\n")
      }
    })
    
    output$forecast_plot <- renderPlot({
      res <- calc_results()
      demands <- parsed_demand()
      plot_ma(demands, res$forecasts, paste0("Moving Average Tahmini (N = ", if (res$mode == "manual") res$n else res$best_n, ")"))
    })
    
    output$result_table <- renderTable({
      res <- calc_results()
      if (res$mode == "manual") {
        res$table
      } else {
        res$comparison_table
      }
    })
  })
}
# -----------------------------------------------------------------------------
# (Opsiyonel) Grafik: Gercek talep ve MA tahmini
# plot_ma() fonksiyonu yukarida (7b bolumunde) tanimlidir.
# -----------------------------------------------------------------------------
plot_ma(demand_data, auto_result$forecasts,
        paste0("MA(", auto_result$best_n, ") - Kriter: ", auto_result$criterion)) 
