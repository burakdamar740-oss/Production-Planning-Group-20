library(shiny)

# Kendi modul dosyalarini burada yukle (dosya adlarini kendi projene gore duzelt).
# Eger dosyalar projedeki "R/" klasorundeyse Shiny bunlari otomatik yukler,
# bu durumda asagidaki source() satirina gerek yoktur.

# source("learning_module.R")   # Learning Curve modulun dosyasi varsa

ui <- navbarPage("Production Planning",
                 tabPanel("Learning Curve", learning_ui("learning")),
                 tabPanel("Moving Average", moving_average_ui("moving_average"))
)

server <- function(input, output, session) {
  learning_server("learning")
  moving_average_server("moving_average")
}

shinyApp(ui, server)