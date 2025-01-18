install.packages(c("shiny", "clipr")) 
library(shiny)

ui <- fluidPage(
  titlePanel("Обрезка последовательности нуклеотидов"),
  
  textAreaInput("sequence", "Введите последовательность:", "", width = "100%", height = "200px"),
  numericInput("start", "Начальная позиция:", value = 1, min = 1),
  numericInput("end", "Конечная позиция:", value = 1, min = 1),
  actionButton("submit", "Обрезать"),
  
  verbatimTextOutput("trimmed_sequence"),
  actionButton("copy", "Скопировать")
)

server <- function(input, output, session) {
  
  observe({
    seq_length <- nchar(input$sequence)
    updateNumericInput(session, "end", max = seq_length, value = seq_length)
  })
  
  trimmed_seq <- eventReactive(input$submit, {
    start <- input$start
    end <- input$end
    sequence <- input$sequence
    
    if (start > 0 && end > 0 && start <= end && end <= nchar(sequence)) {
      substr(sequence, start, end)
    } else {
      "Некорректные параметры обрезки"
    }
  })
  
  
  output$trimmed_sequence <- renderText({
    trimmed_seq()
  })
  
  observeEvent(input$copy, {
    clipr::write_clip(trimmed_seq())
  })
  
  
}


shinyApp(ui = ui, server = server)

