# Come Convertire GUIDA_BIOCONDUCTOR.md in PDF

## Opzione 1: Usando Pandoc (Raccomandato)

```bash
# Installa pandoc (se non l'hai già)
# Su Mac con Homebrew:
brew install pandoc

# Su Ubuntu/Debian:
sudo apt-get install pandoc texlive-latex-base texlive-fonts-recommended

# Converti in PDF con styling avanzato
pandoc GUIDA_BIOCONDUCTOR.md \
  -o GUIDA_BIOCONDUCTOR.pdf \
  --pdf-engine=xelatex \
  -V geometry:margin=1in \
  -V fontsize=11pt \
  -V documentclass=article \
  --toc \
  --toc-depth=2 \
  --highlight-style=tango
```

## Opzione 2: Usando RStudio

1. Apri `GUIDA_BIOCONDUCTOR.md` in RStudio
2. Clicca su "Knit" → "Knit to PDF"
3. Oppure usa:

```r
rmarkdown::render("GUIDA_BIOCONDUCTOR.md", output_format = "pdf_document")
```

## Opzione 3: Usando Visual Studio Code

1. Installa l'estensione "Markdown PDF"
2. Apri `GUIDA_BIOCONDUCTOR.md`
3. Premi `Ctrl+Shift+P` (o `Cmd+Shift+P` su Mac)
4. Digita "Markdown PDF: Export (pdf)" e premi Enter

## Opzione 4: Online (più semplice)

1. Vai su https://www.markdowntopdf.com/
2. Carica `GUIDA_BIOCONDUCTOR.md`
3. Clicca "Convert"
4. Scarica il PDF

## Opzione 5: Usando Chrome/Safari

1. Apri il file markdown in GitHub o in un visualizzatore markdown
2. Usa "Stampa" (Cmd+P / Ctrl+P)
3. Seleziona "Salva come PDF"
