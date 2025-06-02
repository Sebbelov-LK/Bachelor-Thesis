install.packages(c("shiny", "ggseqlogo"))

library(shiny)
library(ggseqlogo)
library(ggplot2)


# Reuse your function here
compute_logo <- function(foreground_peptides, reference_peptides) {
  clean_peptides <- function(peps) {
    peps <- toupper(trimws(peps))
    peps <- peps[nzchar(peps)]
    stopifnot(length(unique(nchar(peps))) == 1)
    return(peps)
  }
  
  compute_ppm_and_bg <- function(peptides) {
    seq_mat <- do.call(rbind, strsplit(peptides, split = ""))
    alphabet <- sort(unique(as.vector(seq_mat)))
    pfm <- sapply(1:ncol(seq_mat), function(i) table(factor(seq_mat[, i], levels = alphabet)))
    rownames(pfm) <- alphabet
    ppm <- apply(pfm, 2, function(col) col / sum(col))
    all_aa <- unlist(strsplit(peptides, split = ""))
    bg <- table(factor(all_aa, levels = alphabet))
    bg_freq <- as.numeric(bg / sum(bg))
    names(bg_freq) <- alphabet
    list(ppm = ppm, bg = bg_freq)
  }
  
  fg <- compute_ppm_and_bg(clean_peptides(foreground_peptides))
  ref <- compute_ppm_and_bg(clean_peptides(reference_peptides))
  
  fg_log <- log2(sweep(fg$ppm, 1, fg$bg, "/"))
  ref_log <- log2(sweep(ref$ppm, 1, ref$bg, "/"))
  fg_log[!is.finite(fg_log)] <- 0
  ref_log[!is.finite(ref_log)] <- 0
  diff <- fg_log - ref_log
  return(diff)
}

# UI
ui <- fluidPage(
  tags$h1(HTML(
    "<span style='color:red'>B</span>
   <span style='color:orange'>I</span>
   <span style='color:green'>C</span>
   <span style='color:blue'>E</span>
   <span style='color:purple'>L</span>
   <span style='color:brown'>o</span>
   <span style='color:teal'>g</span>
   <span style='color:black'>o</span>"
  )),
  
  tags$p("Paste equal-length peptide sequences, one per line."),
  
  fluidRow(
    column(6, textAreaInput("foreground", "Experimental Peptides", rows = 10,
                            placeholder = "ARNDC\nVRNDC\n...")),
    column(6, textAreaInput("reference", "Reference Peptides", rows = 10,
                            placeholder = "ARNDC\nARNDC\n..."))
  ),
  actionButton("plotbtn", "Generate Logo"),
  tags$hr(),
  plotOutput("logo_plot")
)

# Server
server <- function(input, output, session) {
  observeEvent(input$plotbtn, {
    output$logo_plot <- renderPlot({
      fg <- unlist(strsplit(input$foreground, "\n"))
      ref <- unlist(strsplit(input$reference, "\n"))
      mat <- compute_logo(fg, ref)
      ggseqlogo(mat, method = "custom", seq_type = "aa") +
        ggtitle("Self-Normalized Difference Logo")
    })
  })
}

shinyApp(ui, server)