#!/usr/bin/env Rscript
# converti_pdf.R - Script per convertire GUIDA_BIOCONDUCTOR.md in PDF

cat("=== Conversione Guida Bioconductor in PDF ===\n\n")

# Verifica che rmarkdown sia installato
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  cat("Installazione rmarkdown...\n")
  install.packages("rmarkdown")
}

# Verifica che tinytex sia installato (per LaTeX)
if (!requireNamespace("tinytex", quietly = TRUE)) {
  cat("Installazione tinytex...\n")
  install.packages("tinytex")
  tinytex::install_tinytex()
}

# Converti markdown in PDF
cat("Conversione in corso...\n")

tryCatch({
  rmarkdown::render(
    input = "GUIDA_BIOCONDUCTOR.md",
    output_format = rmarkdown::pdf_document(
      toc = TRUE,
      toc_depth = 3,
      number_sections = TRUE,
      highlight = "tango",
      latex_engine = "xelatex"
    ),
    output_file = "GUIDA_BIOCONDUCTOR.pdf"
  )

  cat("\n✅ Conversione completata con successo!\n")
  cat("📄 File creato: GUIDA_BIOCONDUCTOR.pdf\n")

}, error = function(e) {
  cat("\n❌ Errore durante la conversione:\n")
  cat(paste(e$message, "\n"))
  cat("\nProva uno dei metodi alternativi descritti in CONVERTI_IN_PDF.md\n")
})
