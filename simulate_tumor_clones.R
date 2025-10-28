#!/usr/bin/env Rscript

# Simulazione profilo mutazionale di un campione tumorale con sottocloni
# Script per bioinformatica: simula frequenze alleliche da 4 sottocloni
# Include mutazioni condivise secondo gerarchia evolutiva

library(ggplot2)

set.seed(123)  # Per riproducibilità

# ===== PARAMETRI DI SIMULAZIONE =====

# Frequenze dei 4 sottocloni (devono sommare <= 1)
freq_sottocloni <- c(0.15, 0.25, 0.30, 0.30)  # Sottoclone 1, 2, 3, 4
names(freq_sottocloni) <- paste0("Clone", 1:4)

# Numero di mutazioni per ciascun sottoclone (private mutations)
n_mut_per_clone <- c(20, 25, 30, 15)  # Mutazioni uniche di ogni clone

# Numero di mutazioni founder (presenti in tutti i cloni)
n_mut_founder <- 10

# Numero di mutazioni condivise tra sottogruppi di cloni
# Struttura gerarchica: Clone1 -> Clone2,3 -> Clone4
n_mut_condivise <- list(
  "2 3 4" = 15,     # Condivise da Clone2, Clone3, Clone4
  "3 4" = 12,        # Condivise da Clone3 e Clone4
  "1 2" = 8          # Condivise da Clone1 e Clone2
)

# Rumore tecnico (errore di sequenziamento)
rumore_sd <- 0.02  # Deviazione standard del rumore gaussiano

# ===== FUNZIONI =====

# Funzione per generare mutazioni con frequenza allelica
genera_mutazioni <- function(n_mut, freq_base, clone_ids, tipo = "privata", rumore_sd = 0.02) {
  
  if (n_mut == 0) return(NULL)
  
  # Genera nomi mutazioni
  if (tipo == "founder") {
    mut_names <- paste0("Founder_", 1:n_mut)
    clone_label <- "Founder"
  } else if (tipo == "condivisa") {
    clone_str <- paste(clone_ids, collapse = "_")
    mut_names <- paste0("Shared_C", clone_str, "_mut", 1:n_mut)
    clone_label <- paste0("Clone", paste(clone_ids, collapse = "+"))
  } else {
    mut_names <- paste0("Clone", clone_ids, "_mut", 1:n_mut)
    clone_label <- paste0("Clone", clone_ids)
  }
  
  # Frequenza allelica base (con piccola variazione biologica)
  vaf <- rnorm(n_mut, mean = freq_base, sd = rumore_sd)
  
  # Limita VAF tra 0 e 1
  vaf <- pmax(0.01, pmin(0.99, vaf))
  
  # Crea dataframe
  data.frame(
    Mutazione = mut_names,
    VAF = vaf,
    Clone = clone_label,
    Tipo = tipo,
    Clone_IDs = paste(clone_ids, collapse = ","),
    stringsAsFactors = FALSE
  )
}

# ===== GENERAZIONE DATI =====

# Lista per raccogliere tutte le mutazioni
lista_mutazioni <- list()
idx <- 1

# 1. Mutazioni founder (presenti in tutti i cloni)
freq_founder <- sum(freq_sottocloni)  # Somma di tutte le frequenze
lista_mutazioni[[idx]] <- genera_mutazioni(
  n_mut_founder, 
  freq_founder, 
  clone_ids = 1:4, 
  tipo = "founder",
  rumore_sd = rumore_sd
)
idx <- idx + 1

# 2. Mutazioni condivise tra sottogruppi di cloni
for (i in 1:length(n_mut_condivise)) {
  # Estrai i cloni che condividono le mutazioni
  cloni_condivisi <- as.numeric(strsplit(names(n_mut_condivise)[i], " ")[[1]])
  
  # Calcola la frequenza come somma delle frequenze dei cloni coinvolti
  freq_condivisa <- sum(freq_sottocloni[cloni_condivisi])
  
  lista_mutazioni[[idx]] <- genera_mutazioni(
    n_mut_condivise[[i]],
    freq_condivisa,
    clone_ids = cloni_condivisi,
    tipo = "condivisa",
    rumore_sd = rumore_sd
  )
  idx <- idx + 1
}

# 3. Mutazioni private per ogni sottoclone
for (i in 1:length(freq_sottocloni)) {
  lista_mutazioni[[idx]] <- genera_mutazioni(
    n_mut_per_clone[i],
    freq_sottocloni[i],
    clone_ids = i,
    tipo = "privata",
    rumore_sd = rumore_sd
  )
  idx <- idx + 1
}

# Combina tutti i dati (rimuovi eventuali NULL)
lista_mutazioni <- lista_mutazioni[!sapply(lista_mutazioni, is.null)]
dati_completi <- do.call(rbind, lista_mutazioni)

# Aggiungi colonne aggiuntive per realismo
dati_completi$Cromosoma <- sample(paste0("chr", 1:22), nrow(dati_completi), replace = TRUE)
dati_completi$Posizione <- sample(1e6:2e8, nrow(dati_completi), replace = TRUE)
dati_completi$Ref <- sample(c("A", "T", "C", "G"), nrow(dati_completi), replace = TRUE)
dati_completi$Alt <- sample(c("A", "T", "C", "G"), nrow(dati_completi), replace = TRUE)

# Simula copertura di sequenziamento (depth)
dati_completi$Depth <- rpois(nrow(dati_completi), lambda = 100)
dati_completi$Alt_reads <- round(dati_completi$VAF * dati_completi$Depth)

# Riordina colonne
dati_completi <- dati_completi[, c("Mutazione", "Cromosoma", "Posizione", 
                                    "Ref", "Alt", "VAF", "Depth", "Alt_reads",
                                    "Clone", "Tipo")]

# ===== OUTPUT =====

# Salva in file CSV
write.csv(dati_completi, "profilo_mutazionale_simulato.csv", row.names = FALSE)

# Stampa sommario
cat("\n===== SOMMARIO SIMULAZIONE =====\n")
cat("Frequenze sottocloni:\n")
print(freq_sottocloni)
cat("\nNumero totale mutazioni:", nrow(dati_completi), "\n")
cat("  - Mutazioni founder:", n_mut_founder, "\n")
cat("  - Mutazioni condivise:\n")
for (i in 1:length(n_mut_condivise)) {
  cat("    Cloni", names(n_mut_condivise)[i], ":", n_mut_condivise[[i]], "\n")
}
cat("  - Mutazioni private:\n")
for (i in 1:length(n_mut_per_clone)) {
  cat("    Clone", i, ":", n_mut_per_clone[i], "\n")
}

# Statistiche per clone
cat("\n===== VAF MEDIE PER TIPO DI MUTAZIONE =====\n")
print(aggregate(VAF ~ Tipo + Clone, data = dati_completi, FUN = mean))

# Conta mutazioni per tipo
cat("\n===== CONTEGGIO MUTAZIONI PER TIPO =====\n")
print(table(dati_completi$Tipo))

# ===== VISUALIZZAZIONE =====

# Definisci colori per i tipi di mutazione
colori_tipo <- c("founder" = "#E41A1C", 
                 "condivisa" = "#377EB8", 
                 "privata" = "#4DAF4A")

# Grafico 1: Distribuzione VAF per tipo di mutazione
p1 <- ggplot(dati_completi, aes(x = VAF, fill = Tipo)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~Tipo, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = colori_tipo) +
  theme_minimal() +
  labs(title = "Distribuzione VAF per tipo di mutazione",
       x = "Variant Allele Frequency (VAF)",
       y = "Numero di mutazioni") +
  theme(legend.position = "none")

ggsave("VAF_distribuzione_per_tipo.png", p1, width = 10, height = 8)

# Grafico 2: Scatter plot VAF colorato per tipo
p2 <- ggplot(dati_completi, aes(x = 1:nrow(dati_completi), y = VAF, color = Tipo)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = colori_tipo) +
  theme_minimal() +
  labs(title = "Profilo mutazionale: VAF di tutte le mutazioni",
       x = "Indice mutazione",
       y = "Variant Allele Frequency (VAF)",
       color = "Tipo") +
  geom_hline(yintercept = freq_sottocloni, linetype = "dashed", alpha = 0.3, color = "gray50") +
  annotate("text", x = nrow(dati_completi) * 0.95, 
           y = freq_sottocloni, 
           label = paste0("Clone", 1:4), 
           hjust = 1, vjust = -0.5, size = 3, color = "gray30")

ggsave("VAF_scatter_plot.png", p2, width = 12, height = 6)

# Grafico 3: DENSITY PLOT CUMULATIVO - campione mischiato
# Questo simula quello che vedresti sequenziando il DNA del tumore reale
p3 <- ggplot(dati_completi, aes(x = VAF)) +
  geom_density(fill = "#984EA3", alpha = 0.6, color = "#984EA3", linewidth = 1) +
  geom_rug(alpha = 0.3, color = "#984EA3") +
  theme_minimal() +
  labs(title = "Density Plot Cumulativo - Campione Tumorale Mischiato",
       subtitle = paste0("Simulazione di ", nrow(dati_completi), 
                        " mutazioni da 4 sottocloni con frequenze: ",
                        paste(freq_sottocloni, collapse = ", ")),
       x = "Variant Allele Frequency (VAF)",
       y = "Densità") +
  # Aggiungi linee verticali per le frequenze attese dei cloni
  geom_vline(xintercept = freq_sottocloni, 
             linetype = "dashed", 
             color = "red", 
             alpha = 0.5,
             linewidth = 0.8) +
  geom_vline(xintercept = sum(freq_sottocloni), 
             linetype = "dashed", 
             color = "darkred", 
             alpha = 0.7,
             linewidth = 1) +
  annotate("text", 
           x = c(freq_sottocloni, sum(freq_sottocloni)), 
           y = max(density(dati_completi$VAF)$y) * 0.9,
           label = c(paste0("C", 1:4), "Founder"),
           angle = 90, 
           vjust = -0.5, 
           size = 3.5,
           color = "red")

ggsave("VAF_density_cumulativo.png", p3, width = 12, height = 7)

# Grafico 4: Violin plot per tipo di mutazione
p4 <- ggplot(dati_completi, aes(x = Tipo, y = VAF, fill = Tipo)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 1) +
  scale_fill_manual(values = colori_tipo) +
  theme_minimal() +
  labs(title = "Distribuzione VAF per tipo di mutazione (Violin Plot)",
       x = "Tipo di mutazione",
       y = "Variant Allele Frequency (VAF)") +
  theme(legend.position = "none")

ggsave("VAF_violin_plot.png", p4, width = 10, height = 7)

# Grafico 5: Heatmap-like per visualizzare la gerarchia clonale
# Crea una matrice presenza/assenza per ogni clone
# Prima assicuriamoci che Clone_IDs esista
if (!"Clone_IDs" %in% colnames(dati_completi)) {
  # Se non esiste, creiamo Clone_IDs basandoci sul campo Clone
  dati_completi$Clone_IDs <- sapply(dati_completi$Clone, function(x) {
    if (x == "Founder") return("1,2,3,4")
    # Estrai numeri dal nome del clone
    nums <- gsub("[^0-9,+]", "", x)
    nums <- gsub("\\+", ",", nums)
    return(nums)
  })
}

matrice_clonale <- data.frame(
  Mutazione = dati_completi$Mutazione,
  VAF = dati_completi$VAF,
  Tipo = dati_completi$Tipo,
  Clone1 = grepl("1", dati_completi$Clone_IDs),
  Clone2 = grepl("2", dati_completi$Clone_IDs),
  Clone3 = grepl("3", dati_completi$Clone_IDs),
  Clone4 = grepl("4", dati_completi$Clone_IDs)
)

# Ordina per VAF decrescente
matrice_clonale <- matrice_clonale[order(-matrice_clonale$VAF), ]
matrice_clonale$Indice <- 1:nrow(matrice_clonale)

# Reshape per ggplot
library(tidyr)
matrice_long <- pivot_longer(matrice_clonale, 
                              cols = starts_with("Clone"),
                              names_to = "Clone",
                              values_to = "Presente")

p5 <- ggplot(matrice_long, aes(x = Clone, y = Indice, fill = Presente)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#756bb1"),
                    labels = c("Assente", "Presente")) +
  theme_minimal() +
  labs(title = "Matrice di presenza mutazioni nei sottocloni",
       subtitle = "Ordinato per VAF decrescente",
       x = "Sottoclone",
       y = "Mutazione (ordinata per VAF)",
       fill = "Stato") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())

ggsave("matrice_clonale.png", p5, width = 8, height = 10)

cat("\n===== FILE GENERATI =====\n")
cat("1. profilo_mutazionale_simulato.csv - Tabella completa delle mutazioni\n")
cat("2. VAF_distribuzione_per_tipo.png - Istogrammi VAF per tipo di mutazione\n")
cat("3. VAF_scatter_plot.png - Scatter plot di tutte le VAF\n")
cat("4. VAF_density_cumulativo.png - DENSITY PLOT del campione mischiato\n")
cat("5. VAF_violin_plot.png - Violin plot per tipo di mutazione\n")
cat("6. matrice_clonale.png - Heatmap presenza mutazioni nei cloni\n\n")

# Mostra prime righe
cat("Prime righe del dataset:\n")
print(head(dati_completi, 15))
