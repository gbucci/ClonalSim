# Guida Completa alla Sottomissione di Pacchetti R su Bioconductor

**Autore:** Gabriele Bucci
**Data:** Gennaio 2026
**Versione:** 1.0

---

## Indice

1. [Introduzione](#introduzione)
2. [Prerequisiti](#prerequisiti)
3. [Preparazione del Pacchetto](#preparazione-del-pacchetto)
4. [Processo di Sottomissione](#processo-di-sottomissione)
5. [Gestione delle Review](#gestione-delle-review)
6. [Push a git.bioconductor.org](#push-a-gitbioconductororg)
7. [Risoluzione Problemi Comuni](#risoluzione-problemi-comuni)
8. [Best Practices](#best-practices)
9. [Risorse Utili](#risorse-utili)

---

## Introduzione

Bioconductor è un progetto open source per l'analisi e la comprensione di dati genomici ad alta produttività. Questa guida ti accompagnerà passo dopo passo nel processo di sottomissione di un pacchetto R a Bioconductor.

### Cos'è Bioconductor?

- Repository specializzato per pacchetti bioinformatici
- Standard di qualità molto elevati
- Review rigorosa del codice
- Sistema di build automatico
- Documentazione obbligatoria e completa

### Timeline Tipica

- **Preparazione:** 2-4 settimane
- **Review iniziale:** 1-2 settimane
- **Iterazioni di review:** 2-8 settimane
- **Accettazione finale:** variabile

---

## Prerequisiti

### 1. Account GitHub

```bash
# Verifica il tuo account GitHub
git config --global user.name
git config --global user.email
```

### 2. Ambiente R Aggiornato

```r
# Verifica la versione di R
R.version.string

# Deve essere almeno R >= 4.4.0 per Bioconductor 3.20
# Installa BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::version()
```

### 3. Tool di Sviluppo

```r
# Installa devtools e BiocCheck
install.packages("devtools")
BiocManager::install("BiocCheck")
BiocManager::install("BiocCredentials")
```

### 4. SSH Key per git.bioconductor.org

```bash
# Genera una chiave SSH (se non ne hai già una)
ssh-keygen -t ed25519 -C "tua.email@esempio.com"

# Copia la chiave pubblica
cat ~/.ssh/id_ed25519.pub

# Aggiungi la chiave su https://git.bioconductor.org/
```

---

## Preparazione del Pacchetto

### Struttura del Pacchetto

Un pacchetto Bioconductor deve avere questa struttura minima:

```
MioPacchetto/
├── DESCRIPTION          # Metadati del pacchetto
├── NAMESPACE           # Funzioni esportate (generato da roxygen2)
├── LICENSE             # File di licenza
├── README.md           # Descrizione del progetto
├── NEWS.md             # Changelog
├── R/                  # Codice sorgente R
│   ├── AllClasses.R    # Definizioni classi S4
│   ├── AllGenerics.R   # Metodi generici S4
│   └── funzioni.R      # Funzioni del pacchetto
├── man/                # Documentazione (generata da roxygen2)
├── vignettes/          # Tutorial e guide
│   └── introduzione.Rmd
├── tests/              # Test unitari
│   └── testthat/
└── inst/               # File addizionali
    ├── CITATION        # Come citare il pacchetto
    └── extdata/        # Dati di esempio
```

### File DESCRIPTION

Esempio completo basato su ClonalSim:

```r
Package: ClonalSim
Type: Package
Title: Simulation of Tumor Clonal Evolution with Realistic Sequencing Noise
Version: 0.99.1
Authors@R: c(
    person("Gabriele", "Bucci",
           email = "bucci.g@gmail.com",
           role = c("aut", "cre"),
           comment = c(ORCID = "0000-0001-9838-7204"))
    )
Description: ClonalSim generates realistic mutational profiles of tumor
    samples with hierarchical clonal structure. It simulates founder, shared,
    and private mutations with biologically realistic noise models including
    intra-tumor heterogeneity (Beta distribution) and technical sequencing
    noise (negative binomial depth variation, binomial read sampling,
    base errors). The package is designed for benchmarking variant callers,
    testing clonal deconvolution algorithms, and teaching tumor heterogeneity
    concepts.
License: MIT + file LICENSE
Encoding: UTF-8
Depends:
    R (>= 4.4.0)
Imports:
    methods,
    stats,
    utils,
    ggplot2,
    tidyr,
    rlang,
    GenomicRanges,
    IRanges,
    S4Vectors,
    VariantAnnotation
Suggests:
    testthat (>= 3.0.0),
    knitr,
    rmarkdown,
    BiocStyle
biocViews: Software, Sequencing, SomaticMutation,
    VariantDetection, Coverage, Visualization, DataImport
VignetteBuilder: knitr
RoxygenNote: 7.3.3
BugReports: https://github.com/tuousername/MioPacchetto/issues
URL: https://github.com/tuousername/MioPacchetto
```

#### Punti Chiave del DESCRIPTION:

1. **Version**: Per nuovi pacchetti, usa `0.99.x`
2. **Authors@R**: Include ORCID se disponibile
3. **Description**: Deve essere di almeno 3 righe
4. **Depends**: Specifica la versione minima di R
5. **biocViews**: Termini controllati da Bioconductor ([lista completa](https://bioconductor.org/packages/release/BiocViews.html))

### Vignette

Le vignette sono **obbligatorie** per Bioconductor. Esempio:

```rmd
---
title: "Introduzione a MioPacchetto"
author: "Il Tuo Nome"
date: "`r Sys.Date()`"
output: BiocStyle::html_document
vignette: >
  %\VignetteIndexEntry{Introduzione a MioPacchetto}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduzione

Questo pacchetto fa...

## Installazione

```{r installation, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MioPacchetto")
```

## Esempio Base

```{r example}
library(MioPacchetto)

# Usa la tua funzione principale
risultato <- funzionePrincipale()
```

## Session Info

```{r sessionInfo}
sessionInfo()
```
```

### Test Unitari

I test sono fortemente raccomandati:

```r
# tests/testthat.R
library(testthat)
library(MioPacchetto)

test_check("MioPacchetto")
```

```r
# tests/testthat/test-funzioni.R
test_that("funzionePrincipale funziona correttamente", {
  risultato <- funzionePrincipale()
  expect_type(risultato, "list")
  expect_true(length(risultato) > 0)
})
```

---

## Processo di Sottomissione

### Passo 1: Verifica Locale con BiocCheck

Prima di sottomettere, esegui BiocCheck:

```r
# Dal terminale R nella directory del pacchetto
library(BiocCheck)

# Check completo
BiocCheck(".")

# Opzionale: build del pacchetto per check più approfonditi
devtools::check()
```

**Risolvi tutti gli errori e avvertimenti prima di procedere.**

### Passo 2: Crea Repository GitHub

```bash
# Inizializza git (se non l'hai già fatto)
cd MioPacchetto
git init
git add .
git commit -m "Initial commit"

# Crea il repository su GitHub e collega
git remote add origin git@github.com:tuousername/MioPacchetto.git
git push -u origin main
```

### Passo 3: Apri Issue su Bioconductor/Contributions

1. Vai su: https://github.com/Bioconductor/Contributions/issues/new/choose
2. Seleziona "New Package Submission"
3. Compila il template:

```markdown
## Package Name
MioPacchetto

## Package Repository
https://github.com/tuousername/MioPacchetto

## Maintainer GitHub Handle
@tuousername

## Package Description
Breve descrizione del pacchetto e del suo scopo principale.

## Type of Package (Software/Annotation/Experiment/Workflow)
Software

## Confirmation
- [x] I understand Bioconductor's package guidelines
- [x] My package passes R CMD check and BiocCheck
- [x] My package includes a vignette
- [x] All dependencies are available on CRAN or Bioconductor
```

### Passo 4: Attendi la Review Automatica

Il sistema Bioconductor eseguirà automaticamente:

- Build del pacchetto
- R CMD check
- BiocCheck
- Test su multiple piattaforme (Linux, macOS, Windows)

Riceverai un commento automatico con i risultati.

---

## Gestione delle Review

### Tipologie di Feedback

#### 1. Errori Critici (da risolvere immediatamente)

```
ERROR: this R is version 4.5.2, package 'MioPacchetto' requires R >= 4.6.0
```

**Soluzione:**
```r
# Aggiorna DESCRIPTION
Depends:
    R (>= 4.4.0)  # Usa versione disponibile su CRAN
```

#### 2. Warning da Risolvere

```
WARNING: fishplot not available in official repositories
```

**Soluzione:**
```r
# Rimuovi dipendenze non ufficiali da DESCRIPTION
# Oppure spostale in Enhances: invece che Suggests:
```

#### 3. Note e Suggerimenti

```
NOTE: Consider using BiocStyle for vignettes
```

**Implementa quando possibile per migliorare la qualità.**

### Workflow di Correzione

```bash
# 1. Crea un branch per le correzioni
git checkout -b fix-bioc-review

# 2. Fai le modifiche necessarie
# ... edita i file ...

# 3. Testa localmente
R CMD check .
BiocCheck::BiocCheck(".")

# 4. Commit e push
git add .
git commit -m "Fix: risolto problema con dipendenze R"
git push origin fix-bioc-review

# 5. Merge su main (tramite PR su GitHub)
# 6. IMPORTANTE: Version bump!
# Aggiorna Version in DESCRIPTION: 0.99.0 -> 0.99.1
git add DESCRIPTION
git commit -m "Bump version to 0.99.1"
git push origin main
```

### Esempio Reale: Caso ClonalSim

**Problema riscontrato:**
```
Error in simulateTumor() : object 'seed' not found
```

**Causa:**
- Il parametro `seed` era usato nel codice ma non definito nella firma della funzione
- La documentazione aveva un valore di default errato

**Soluzione applicata:**

```r
# 1. Aggiornata firma della funzione in R/simulateTumor.R
simulateTumor <- function(
    subclone_freqs = c(0.15, 0.25, 0.30, 0.30),
    n_mut_per_clone = c(20, 25, 30, 15),
    # ... altri parametri ...
    seed = NULL  # ← AGGIUNTO
) {
  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # ... resto del codice ...
}

# 2. Aggiornata documentazione
#' @param seed integer, random seed for reproducibility (default: NULL)

# 3. Corretto man/simulateTumor.Rd
# seed = NULL  (era erroneamente seed = 123)

# 4. Version bump
# Version: 0.99.0 -> 0.99.1
```

**Comunicazione con reviewer:**
```markdown
@lshep Ho corretto il problema:
- Aggiunto parametro `seed = NULL` alla funzione
- Aggiornata documentazione
- Version bump a 0.99.1
- Commit: https://github.com/gbucci/ClonalSim/commit/abc1234
```

---

## Push a git.bioconductor.org

### Setup Iniziale

#### 1. Configura BiocCredentials

```r
# In R
library(BiocCredentials)

# Questo aprirà una pagina web per l'autenticazione
BiocCredentials::setCredentials()

# Verifica le credenziali
BiocCredentials::testCredentials()
```

#### 2. Aggiungi Remote Bioconductor

```bash
# Nella directory del tuo pacchetto
cd MioPacchetto

# Aggiungi il remote Bioconductor
git remote add bioc git@git.bioconductor.org:packages/MioPacchetto

# Verifica i remote
git remote -v
# Output:
# origin  git@github.com:tuousername/MioPacchetto.git (fetch)
# origin  git@github.com:tuousername/MioPacchetto.git (push)
# bioc    git@git.bioconductor.org:packages/MioPacchetto (fetch)
# bioc    git@git.bioconductor.org:packages/MioPacchetto (push)
```

### Push dopo Review

#### 1. Assicurati di avere l'ultimo codice

```bash
# Passa al main branch
git checkout main

# Aggiorna da GitHub
git pull origin main
```

#### 2. Verifica i Branch Disponibili su Bioconductor

```bash
# Controlla quali branch esistono su Bioconductor
git ls-remote bioc

# Output tipico:
# e64ed6229cdbfdac74b7f109b4afab7cd54e181c  HEAD
# e64ed6229cdbfdac74b7f109b4afab7cd54e181c  refs/heads/devel
```

#### 3. Push al Branch Corretto

**Per pacchetti in review o nuovi:**

```bash
# Push al branch devel
git push bioc main:devel
```

**Per pacchetti già accettati:**

```bash
# Push sia a master (release) che devel
git push bioc main:master
git push bioc main:devel
```

#### 4. Verifica il Push

```bash
# Controlla che il push sia andato a buon fine
git ls-remote bioc

# Dovresti vedere il tuo commit hash aggiornato
```

### Esempio Completo: Caso ClonalSim

```bash
# 1. Stato iniziale
$ git status
On branch main
Your branch is up to date with 'origin/main'.

# 2. Pull da GitHub per assicurarsi di avere tutto
$ git pull origin main
Already up to date.

# 3. Aggiungi remote Bioconductor
$ git remote add bioc git@git.bioconductor.org:packages/ClonalSim

# 4. Verifica branch disponibili su Bioconductor
$ git ls-remote bioc
e64ed6229cdbfdac74b7f109b4afab7cd54e181c  HEAD
e64ed6229cdbfdac74b7f109b4afab7cd54e181c  refs/heads/devel

# 5. Push al branch devel (per pacchetti in review)
$ git push bioc main:devel
Enumerating objects: 49, done.
Counting objects: 100% (42/42), done.
Delta compression using up to 8 threads
Compressing objects: 100% (33/33), done.
Writing objects: 100% (33/33), 11.13 KiB | 5.57 MiB/s, done.
Total 33 (delta 20), reused 0 (delta 0)
To git.bioconductor.org:packages/ClonalSim
   e64ed62..cf2b9c8  main -> devel

# 6. Verifica
$ git ls-remote bioc
cf2b9c8... HEAD
cf2b9c8... refs/heads/devel
```

### Notifica su GitHub

Dopo il push, commenta sulla issue di Bioconductor:

```markdown
@lshep Changes have been pushed successfully!

**What was done:**
1. ✅ Version bumped from 0.99.0 to 0.99.1
2. ✅ Pushed to git.bioconductor.org:packages/ClonalSim (devel branch)
3. ✅ All changes synced to GitHub

**Changes included:**
- Updated R dependency to >= 4.4.0
- Fixed seed parameter in simulateTumor() function
- Fixed documentation default values

**Commits:**
- GitHub: https://github.com/gbucci/ClonalSim/commit/cf2b9c8
- Bioconductor devel branch updated from e64ed62 to cf2b9c8

The package should now trigger a new build. Thank you!
```

---

## Risoluzione Problemi Comuni

### Problema 1: Errore di Permessi su Push

```
remote: FATAL: W refs/heads/master packages/MioPacchetto user DENIED by fallthru
remote: error: hook declined to update refs/heads/master
```

**Causa:** Stai cercando di pushare a un branch protetto o non hai i permessi.

**Soluzione:**
```bash
# 1. Verifica quali branch sono disponibili
git ls-remote bioc

# 2. Push al branch corretto (solitamente 'devel' per pacchetti in review)
git push bioc main:devel

# 3. Se ancora non funziona, chiedi supporto sulla GitHub issue
```

### Problema 2: Version Bump Dimenticato

```
Il build non si triggera dopo il push
```

**Causa:** Bioconductor richiede un version bump per triggerare un nuovo build.

**Soluzione:**
```r
# 1. Aggiorna Version in DESCRIPTION
# Da: Version: 0.99.0
# A:  Version: 0.99.1

# 2. Commit
git add DESCRIPTION
git commit -m "Bump version to 0.99.1 for Bioconductor build"

# 3. Push a GitHub E Bioconductor
git push origin main
git push bioc main:devel
```

### Problema 3: Conflitti di Dipendenze

```
Package 'xyz' is not available in CRAN or Bioconductor
```

**Soluzione:**
```r
# Opzione 1: Rimuovi la dipendenza se non essenziale
# Opzione 2: Spostala in Enhances: invece di Suggests:
# Opzione 3: Sostituisci con un pacchetto disponibile ufficialmente

# Esempio: fishplot non è disponibile → rimuovi da Suggests
```

### Problema 4: Errori nella Documentazione

```
WARNING: Rd files with mismatched \alias
```

**Soluzione:**
```r
# Usa roxygen2 per rigenerare la documentazione
library(roxygen2)
roxygenize()

# Commit i file man/ aggiornati
git add man/
git commit -m "Update documentation"
```

### Problema 5: Test Falliti

```
ERROR: Test failures
```

**Soluzione:**
```r
# 1. Esegui i test localmente
devtools::test()

# 2. Debug del test specifico
testthat::test_file("tests/testthat/test-funzione.R")

# 3. Correggi il codice o il test
# 4. Re-run per verificare
```

### Problema 6: Vignette non Compila

```
ERROR: re-building vignettes failed
```

**Soluzione:**
```r
# 1. Build vignette localmente
devtools::build_vignettes()

# 2. Controlla errori specifici
# Spesso sono dovuti a:
# - Pacchetti mancanti in Suggests
# - Code chunk che generano errori
# - Path relativi non corretti

# 3. Testa la vignette
rmarkdown::render("vignettes/introduzione.Rmd")
```

---

## Best Practices

### 1. Versioning

```
0.99.0 → Prima submission
0.99.1 → Prima revisione
0.99.2 → Seconda revisione
...
1.0.0  → Prima release ufficiale su Bioconductor (gestita da loro)
```

### 2. Commit Messages Chiari

**Buono:**
```
Fix: resolve seed parameter issue in simulateTumor()

The function was referencing an undefined 'seed' variable.
Added seed = NULL parameter to function signature and
implemented proper seed handling.

Fixes #123
```

**Evita:**
```
fix stuff
```

### 3. Documentazione Completa

```r
#' Simula l'evoluzione tumorale clonale
#'
#' @description Questa funzione genera profili mutazionali realistici
#'   di campioni tumorali con struttura clonale gerarchica.
#'
#' @param subclone_freqs numeric vector, frequenze di ciascun subclone
#' @param seed integer, seed casuale per riproducibilità (default: NULL)
#'
#' @return Un oggetto ClonalSimData contenente le mutazioni simulate
#'
#' @examples
#' # Esempio base
#' sim <- simulateTumor()
#'
#' # Con seed per riproducibilità
#' sim <- simulateTumor(seed = 123)
#'
#' @export
simulateTumor <- function(subclone_freqs, seed = NULL) {
  # Implementazione
}
```

### 4. Test Completi

```r
# Testa casi normali
test_that("simulateTumor funziona con parametri di default", {
  sim <- simulateTumor()
  expect_s4_class(sim, "ClonalSimData")
})

# Testa casi limite
test_that("simulateTumor gestisce input non validi", {
  expect_error(simulateTumor(subclone_freqs = c(0.6, 0.6)))
})

# Testa riproducibilità
test_that("simulateTumor è riproducibile con seed", {
  sim1 <- simulateTumor(seed = 42)
  sim2 <- simulateTumor(seed = 42)
  expect_equal(getMutations(sim1), getMutations(sim2))
})
```

### 5. Vignette Informative

Includi:
- **Introduzione** chiara al problema
- **Installazione** passo-passo
- **Esempi pratici** funzionanti
- **Use cases** reali
- **Interpretazione** dei risultati
- **Session info** alla fine

### 6. NEWS.md Aggiornato

```markdown
# ClonalSim 0.99.1

## Bug Fixes

* Fixed 'object seed not found' error in simulateTumor() (#7)
* Corrected documentation default values

## Changes

* Updated R dependency to >= 4.4.0 for compatibility with Bioconductor 3.20

# ClonalSim 0.99.0

* Initial Bioconductor submission
```

### 7. Rispondere Prontamente ai Reviewer

- Rispondi entro **1-2 giorni** quando possibile
- Sii professionale e cortese
- Se non capisci un commento, chiedi chiarimenti
- Mostra cosa hai fatto per risolvere il problema

Esempio:
```markdown
@reviewer Thank you for the feedback!

I've addressed all the points:

1. **R dependency**: Updated from 4.6.0 to 4.4.0 ✅
2. **Seed parameter**: Added proper implementation ✅
3. **Documentation**: Fixed default values ✅

Version bumped to 0.99.1 and pushed to both GitHub and git.bioconductor.org.

Let me know if you need any other changes!
```

---

## Risorse Utili

### Documentazione Ufficiale

- **Bioconductor Package Guidelines:** https://contributions.bioconductor.org/
- **Package Submission:** https://github.com/Bioconductor/Contributions
- **BiocCheck:** https://bioconductor.org/packages/BiocCheck/
- **BiocViews Terms:** https://bioconductor.org/packages/release/BiocViews.html

### Tool e Pacchetti

```r
# Essenziali
install.packages("devtools")
install.packages("roxygen2")
BiocManager::install("BiocCheck")
BiocManager::install("BiocCredentials")
BiocManager::install("BiocStyle")

# Per testing
install.packages("testthat")
install.packages("covr")  # Coverage dei test
```

### Comandi Utili

```bash
# Build locale
R CMD build .
R CMD check MioPacchetto_0.99.1.tar.gz

# BiocCheck
R -e "BiocCheck::BiocCheck('.')"

# Test coverage
R -e "covr::package_coverage()"

# Build vignettes
R -e "devtools::build_vignettes()"
```

### Template Checklist Pre-Submission

```markdown
## Pre-Submission Checklist

- [ ] Version è 0.99.0 (o 0.99.x per re-submission)
- [ ] DESCRIPTION completo e corretto
- [ ] biocViews appropriati
- [ ] R dependency corretta (>= 4.4.0)
- [ ] Tutte le dipendenze disponibili su CRAN/Bioc
- [ ] LICENSE file presente
- [ ] README.md informativo
- [ ] NEWS.md con changelog
- [ ] Almeno una vignette completa
- [ ] Documentazione completa (roxygen2)
- [ ] Test unitari presenti
- [ ] R CMD check senza errori/warning
- [ ] BiocCheck passa senza errori
- [ ] Repository GitHub pubblico
- [ ] SSH key configurata per git.bioconductor.org
- [ ] Issue creata su Bioconductor/Contributions
```

### Community e Supporto

- **Bioconductor Support Site:** https://support.bioconductor.org/
- **Bioconductor Slack:** https://community-bioc.slack.com/
- **GitHub Discussions:** https://github.com/Bioconductor/Contributions/discussions

---

## Esempio Completo: Timeline di ClonalSim

### Week 1: Preparazione
- ✅ Struttura pacchetto completa
- ✅ Vignette scritta
- ✅ Test implementati
- ✅ R CMD check passa
- ✅ BiocCheck passa

### Week 2: Submission
- ✅ Repository GitHub creato
- ✅ Issue aperta su Bioconductor/Contributions
- ⏳ Build automatico iniziato

### Week 3-4: Prima Review
- ❌ Errore: R >= 4.6.0 non disponibile
- ❌ Errore: fishplot non in repository ufficiali
- ❌ Errore: set.seed() nel codice del pacchetto
- 🔧 Correzioni applicate
- ✅ Version bump a 0.99.1

### Week 5: Seconda Review
- ❌ Errore: object 'seed' not found
- 🔧 Parametro seed aggiunto alla funzione
- 🔧 Documentazione corretta
- ✅ Version bump a 0.99.2
- ✅ Push a git.bioconductor.org riuscito

### Week 6: Review Finale
- ✅ Build su tutte le piattaforme: SUCCESS
- ✅ BiocCheck: PASS
- 🎉 Pacchetto accettato

### Post-Accettazione
- Pacchetto disponibile in Bioconductor devel
- Alla prossima release diventerà parte di Bioconductor release
- Manutenzione continua richiesta

---

## Conclusione

La sottomissione a Bioconductor richiede:

1. **Preparazione accurata** - Tempo investito in documentazione e test
2. **Pazienza** - Il processo di review può richiedere settimane
3. **Collaborazione** - Lavoro con reviewer per migliorare il pacchetto
4. **Manutenzione** - Impegno continuo dopo l'accettazione

**Benefici:**

- ✅ Visibilità nella comunità bioinformatica
- ✅ Standard di qualità riconosciuti
- ✅ Build automatici su multiple piattaforme
- ✅ Integrazione nell'ecosistema Bioconductor
- ✅ Citazioni e impatto scientifico

**Ricorda:** Ogni errore o richiesta di modifica è un'opportunità per migliorare la qualità del tuo pacchetto!

---

## Appendice: Script Utili

### Script per Check Pre-Submission

```r
#!/usr/bin/env Rscript
# pre_submission_check.R

cat("=== Pre-Submission Check ===\n\n")

# 1. R CMD check
cat("1. Running R CMD check...\n")
devtools::check()

# 2. BiocCheck
cat("\n2. Running BiocCheck...\n")
BiocCheck::BiocCheck(".")

# 3. Test coverage
cat("\n3. Calculating test coverage...\n")
cov <- covr::package_coverage()
print(cov)
covr::report(cov)

# 4. Build vignettes
cat("\n4. Building vignettes...\n")
devtools::build_vignettes()

cat("\n=== Check Complete ===\n")
```

### Script per Version Bump

```bash
#!/bin/bash
# bump_version.sh

# Leggi versione corrente
CURRENT=$(grep "^Version:" DESCRIPTION | sed 's/Version: //')
echo "Current version: $CURRENT"

# Incrementa patch version
IFS='.' read -ra PARTS <<< "$CURRENT"
MAJOR=${PARTS[0]}
MINOR=${PARTS[1]}
PATCH=${PARTS[2]}
NEW_PATCH=$((PATCH + 1))
NEW_VERSION="$MAJOR.$MINOR.$NEW_PATCH"

echo "New version: $NEW_VERSION"

# Aggiorna DESCRIPTION
sed -i '' "s/Version: $CURRENT/Version: $NEW_VERSION/" DESCRIPTION

# Commit
git add DESCRIPTION
git commit -m "Bump version to $NEW_VERSION"

echo "Version bumped! Don't forget to push."
```

### Script per Sincronizzazione

```bash
#!/bin/bash
# sync_bioc.sh

# Assicurati di essere su main
git checkout main
git pull origin main

# Push a GitHub
echo "Pushing to GitHub..."
git push origin main

# Push a Bioconductor
echo "Pushing to Bioconductor..."
git push bioc main:devel

echo "Sync complete!"
```

---

**Buona fortuna con la tua sottomissione a Bioconductor! 🚀**
