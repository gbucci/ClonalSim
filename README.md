# Simulatore di Profili Mutazionali Tumorali

Script R per la simulazione di profili mutazionali di campioni tumorali con struttura clonale gerarchica.

## Descrizione

Questo strumento permette di simulare un campione di DNA tumorale composto da una miscela di cloni e sottocloni con frequenze alleliche (VAF - Variant Allele Frequency) specifiche. Lo script genera dati realistici che includono:

- **Mutazioni founder**: presenti in tutti i sottocloni (eventi iniziali)
- **Mutazioni condivise**: presenti in sottogruppi di cloni secondo una gerarchia evolutiva
- **Mutazioni private**: specifiche di singoli sottocloni
- **Rumore tecnico**: simulazione di errori di sequenziamento

## Caratteristiche

### Struttura Clonale Simulata

Lo script simula un tumore eterogeneo con:
- 4 sottocloni con frequenze configurabili (default: 0.15, 0.25, 0.30, 0.30)
- Gerarchia evolutiva: Clone1 → Clone2,3 → Clone4
- Mutazioni distribuite secondo modello evolutivo realistico

### Output Generati

1. **profilo_mutazionale_simulato.csv**: Dataset completo con tutte le mutazioni
2. **VAF_density_cumulativo.png**: Density plot del campione mischiato (come sequenziamento reale)
3. **VAF_scatter_plot.png**: Scatter plot di tutte le VAF
4. **VAF_distribuzione_per_tipo.png**: Istogrammi per tipo di mutazione
5. **VAF_violin_plot.png**: Violin plot per visualizzare distribuzioni
6. **matrice_clonale.png**: Heatmap presenza/assenza mutazioni nei cloni

## Requisiti

### Software
- R >= 4.0
- Pacchetti R:
  - `ggplot2`
  - `tidyr`

### Installazione dipendenze

```bash
# Ubuntu/Debian
sudo apt-get install r-base r-cran-ggplot2 r-cran-tidyr

# Oppure da R
install.packages(c("ggplot2", "tidyr"))
```

## Utilizzo

### Esecuzione Base

```bash
# Rendi eseguibile lo script
chmod +x simulate_tumor_clones.R

# Esegui
./simulate_tumor_clones.R
```

Oppure da R:

```r
source("simulate_tumor_clones.R")
```

### Personalizzazione Parametri

Modifica i parametri all'inizio dello script:

```r
# Frequenze dei 4 sottocloni (devono sommare <= 1)
freq_sottocloni <- c(0.15, 0.25, 0.30, 0.30)

# Numero di mutazioni per ciascun sottoclone (private)
n_mut_per_clone <- c(20, 25, 30, 15)

# Numero di mutazioni founder
n_mut_founder <- 10

# Struttura mutazioni condivise (gerarchia evolutiva)
n_mut_condivise <- list(
  "2 3 4" = 15,    # Condivise da Clone2, Clone3, Clone4
  "3 4" = 12,      # Condivise da Clone3 e Clone4
  "1 2" = 8        # Condivise da Clone1 e Clone2
)

# Rumore tecnico (deviazione standard)
rumore_sd <- 0.02
```

## Struttura dei Dati di Output

### File CSV

Il file `profilo_mutazionale_simulato.csv` contiene le seguenti colonne:

| Colonna | Descrizione |
|---------|-------------|
| `Mutazione` | Identificatore univoco della mutazione |
| `Cromosoma` | Cromosoma (chr1-chr22) |
| `Posizione` | Posizione genomica |
| `Ref` | Allele di riferimento |
| `Alt` | Allele alternativo |
| `VAF` | Variant Allele Frequency (0-1) |
| `Depth` | Copertura di sequenziamento |
| `Alt_reads` | Numero di read con allele alternativo |
| `Clone` | Clone/i che portano la mutazione |
| `Tipo` | Tipo di mutazione (troncale/condivisa/privata) |
| `Clone_IDs` | ID numerici dei cloni coinvolti |

### Esempio Output

```
Mutazione          Cromosoma  Posizione  Ref  Alt  VAF     Depth  Alt_reads  Clone       Tipo
Founder_1          chr19      26274770   C    T    0.9888  101    100        Founder     founder
Shared_C2_3_4_mut1 chr9       76062670   C    C    0.8745  85     74         Clone2+3+4  condivisa
Clone1_mut1        chr15      45123456   A    G    0.1498  95     14         Clone1      privata
```

## Interpretazione dei Risultati

### Density Plot Cumulativo

Il grafico **VAF_density_cumulativo.png** rappresenta la distribuzione delle frequenze alleliche come se avessi sequenziato il DNA del tumore reale (miscela di tutti i cloni). 

**Picchi attesi:**
- **~1.0**: Mutazioni founder (presenti in tutte le cellule tumorali)
- **~0.85**: Mutazioni condivise da Clone2+3+4
- **~0.60**: Mutazioni condivise da Clone3+4
- **~0.40**: Mutazioni condivise da Clone1+2
- **0.15-0.30**: Mutazioni private dei singoli cloni

### Matrice Clonale

L'heatmap **matrice_clonale.png** mostra quali mutazioni sono presenti in quali cloni, ordinata per VAF decrescente. Questo visualizza chiaramente la gerarchia evolutiva del tumore.

## Modello Biologico

### Assunzioni del Modello

1. **Purezza del campione**: La somma delle frequenze clonali rappresenta la frazione tumorale
2. **Ploidia**: Assunzione diploide (VAF max teorica = 0.5 per mutazioni eterozigoti)
3. **Evoluzione clonale**: Struttura gerarchica (albero filogenetico)
4. **Errori tecnici**: Rumore gaussiano con SD configurabile

### Limitazioni

- Non simula Copy Number Alterations (CNA)
- Non considera contaminazione con cellule normali esplicita
- Assunzione di mutazioni heterozygous
- Non include subclonal loss of heterozygosity (LOH)

## Casi d'Uso

### Ricerca

- Testing di algoritmi di deconvoluzione clonale
- Benchmarking di variant caller
- Sviluppo di metodi di analisi filogenetica tumorale
- Studi di eterogeneità intratumorale

### Didattica

- Insegnamento di concetti di evoluzione clonale
- Visualizzazione di eterogeneità tumorale
- Esercitazioni su analisi di NGS data

### Validazione

- Controllo positivo per pipeline di analisi
- Test di sensibilità di metodi di detection
- Validazione di strumenti bioinformatici

## Esempi Avanzati

### Simulare un Tumore con Alta Eterogeneità

```r
# 6 sottocloni con frequenze diverse
freq_sottocloni <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.25)

# Aumenta il numero di mutazioni
n_mut_per_clone <- c(30, 40, 50, 60, 70, 50)

# Struttura gerarchica più complessa
n_mut_condivise <- list(
  "1 2 3 4 5 6" = 20,  # Founder (in alternativa a n_mut_founder)
  "2 3 4 5 6" = 15,
  "3 4 5 6" = 12,
  "4 5 6" = 10,
  "5 6" = 8,
  "1 2" = 5
)
```

### Simulare Tumore con Bassa Purezza

```r
# Frequenze basse per simulare contaminazione normale
freq_sottocloni <- c(0.05, 0.10, 0.12, 0.13)  # Somma = 0.40 (40% cellule tumorali)

# Aumenta il rumore per simulare errori maggiori
rumore_sd <- 0.03
```

### Esportare per Analisi Downstream

```r
# Dopo aver generato i dati
dati <- read.csv("profilo_mutazionale_simulato.csv")

# Converti in formato VCF-like
vcf_like <- dati[, c("Cromosoma", "Posizione", "Ref", "Alt", "VAF", "Depth", "Alt_reads")]

# Salva
write.table(vcf_like, "mutations.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

## Troubleshooting

### Errore: "package 'ggplot2' is not available"

```bash
# Installa i pacchetti mancanti
R -e "install.packages(c('ggplot2', 'tidyr'), repos='https://cloud.r-project.org')"
```

### Warning: "removed rows containing missing values"

Normale se alcune mutazioni hanno VAF molto basse o alte che vengono troncate a 0.01-0.99.

### Le frequenze non sommano a 1

Questo è intenzionale! Rappresenta la possibilità di contaminazione con cellule normali. Se vuoi purezza 100%, assicurati che `sum(freq_sottocloni) = 1.0`.

## Estensioni Future

- [ ] Supporto per Copy Number Variations (CNV)
- [ ] Simulazione di Loss of Heterozygosity (LOH)
- [ ] Modelli di evoluzione non-gerarchici (reticolare)
- [ ] Export diretto in formato VCF
- [ ] Interfaccia Shiny per uso interattivo
- [ ] Simulazione di dati RNA-seq con allelic expression

## Contribuire

Contributi, bug reports e feature requests sono benvenuti! 

### Come Contribuire

1. Fork del repository
2. Crea un branch per la tua feature (`git checkout -b feature/AmazingFeature`)
3. Commit delle modifiche (`git commit -m 'Add some AmazingFeature'`)
4. Push al branch (`git push origin feature/AmazingFeature`)
5. Apri una Pull Request

## Licenza

Questo progetto è rilasciato sotto licenza MIT. Vedi il file `LICENSE` per i dettagli.

## Autori

Creato per applicazioni bioinformatiche in oncologia computazionale.

## Citazione

Se utilizzi questo strumento nella tua ricerca, per favore cita:

```
Tumor Clonal Simulator - Tool per la simulazione di profili mutazionali tumorali
con struttura clonale gerarchica
```

## Contatti

Per domande, suggerimenti o collaborazioni, apri una issue su GitHub.

## Riferimenti

### Letture Consigliate

1. McGranahan N, Swanton C. **Clonal Heterogeneity and Tumor Evolution: Past, Present, and the Future.** Cell. 2017
2. Dentro SC, Wedge DC, Van Loo P. **Principles of Reconstructing the Subclonal Architecture of Cancers.** Cold Spring Harb Perspect Med. 2017
3. Schwarz RF, et al. **Phylogenetic Quantification of Intra-tumour Heterogeneity.** PLoS Comput Biol. 2014

### Strumenti Correlati

- **PyClone**: Analisi statistica di sottocloni tumorali
- **SciClone**: Identificazione di cluster di mutazioni
- **ABSOLUTE**: Stima di purezza e ploidia tumorale
- **THetA**: Inferenza di copy number e purezza

---

**Note**: Questo è uno strumento di simulazione. I dati generati sono sintetici e per scopi educativi/di ricerca.
