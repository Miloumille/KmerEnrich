#!/usr/bin/env Rscript
# KMerEnrich v2.0 — analyse de 2 virus SANS background (null statistique par défaut).
#
# Prérequis : le package KMerEnrich installé
#   devtools::install("chemin/vers/2026_KMerEnrich_2.0", build_vignettes = FALSE)
#
# Données attendues : un dossier `virus/` avec UN sous-dossier par virus
#   virus/
#   ├── virusA/  (une ou plusieurs souches .fasta)
#   └── virusB/  (une ou plusieurs souches .fasta)
#
# Lancement :  Rscript analyze_two_viruses.R /chemin/vers/virus /chemin/vers/sorties

suppressMessages(library(KMerEnrich))
suppressMessages(library(Biostrings))

# ---- Paramètres (à calibrer sur tes données) --------------------------------
args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[[1]] else "virus"      # dossier des virus
OUT_DIR  <- if (length(args) >= 2) args[[2]] else "kmerenrich_out"
K              <- 9L     # longueur du primer (SMART-9N)
MAX_DEGENERATE <- 3L     # positions dégénérées max par primer
READ_LEN       <- 200L   # empreinte de la RT (bases couvertes par primer)
BIN_SIZE       <- 100L   # finesse de la couverture
FLOOR          <- 0.5    # couverture minimale garantie PAR virus (équité panel)
MIN_ENRICH     <- 2      # enrichissement minimal admis
TM_MIN <- 0; TM_MAX <- 60; GC_MIN <- 0.2; GC_MAX <- 0.8        # filtre thermo
MAX_CANDIDATES <- 150L   # primers max envoyés à l'ILP (au-delà de ~250, l'optim explose)
KMAX           <- 20L    # taille de pool max évaluée par la courbe de saturation
MIN_SITES      <- 2L     # écarte les primers ne se fixant qu'à un seul locus (peu utiles)
# -----------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Cibles (2 virus) — pas de background : null statistique uniforme
targets <- load_targets(DATA_DIR)
stopifnot(length(targets) == 2)
model   <- default_background()          # 25 % par base, indépendance (order-0)

# 2. Candidats : k-mers conservés de CHAQUE virus, fusionnés en primers IUPAC
kmers   <- unique(unlist(lapply(targets, enumerate_conserved_kmers, k = K), use.names = FALSE))
primers <- collapse_to_iupac(kmers, max_degenerate = MAX_DEGENERATE, anchor = 2L)

# 3. Enrichissement (vs null) + filtre thermo/enrichissement
viral  <- do.call(c, unname(targets))    # toutes les souches des 2 virus
scored <- score_enrichment(primers, viral, model)
admit  <- filter_candidates(scored, min_enrichment = MIN_ENRICH,
                            tm_min = TM_MIN, tm_max = TM_MAX, gc_min = GC_MIN, gc_max = GC_MAX)
stopifnot(nrow(admit) > 0)

# Plafond : l'ILP de couverture devient intraitable au-delà de ~250 primers.
if (nrow(admit) > MAX_CANDIDATES) {
  message(sprintf("[info] %d primers admis -> on garde les %d meilleurs par enrichissement (au-dela de ~250, l'optimisation ILP peut prendre des minutes/heures).",
                  nrow(admit), MAX_CANDIDATES))
  admit <- head(admit[order(-admit$enrichment), ], MAX_CANDIDATES)
}

# Écarte les primers qui ne se fixent pas >= MIN_SITES fois sur CHAQUE génome de
# référence : un primer utile au panel doit viser plusieurs loci sur tous les virus
# (sinon il ne sera pas retenu et pollue les résultats).
if (MIN_SITES > 1L) {
  refs  <- lapply(targets, function(x) x[1])
  nsite <- vapply(admit$pattern, function(p)
    min(vapply(refs, function(r)
      Biostrings::vcountPattern(p, r, fixed = c(pattern = FALSE, subject = TRUE)), integer(1))),
    integer(1))
  message(sprintf("[info] %d primers ecartes (< %d sites sur au moins un genome).",
                  sum(nsite < MIN_SITES), MIN_SITES))
  admit <- admit[nsite >= MIN_SITES, ]
  stopifnot(nrow(admit) > 0)
}

# 4. Couverture, saturation (N*), pool optimal avec plancher par virus
cov  <- coverage_matrix(admit$pattern, targets, read_len = READ_LEN, bin_size = BIN_SIZE)
sb   <- saturation_balanced(cov, k_max = min(nrow(admit), KMAX))
best <- sb$pools[[sb$n_star]]

# 5. Score global + séquences à commander (complément inverse)
enr    <- mean(admit$enrichment[admit$pattern %in% best$patterns])
score  <- desirability(enr, best$coverage, best$regularity)
orders <- as.character(reverseComplement(DNAStringSet(best$patterns)))
pool   <- data.frame(primer = best$patterns, a_commander = orders,
                     enrichissement = round(admit$enrichment[match(best$patterns, admit$pattern)], 1))

# 6. Sortie : un unique rapport HTML interactif (export CSV intégré dans la page)
report <- kmer_report(targets, admit, cov = cov, read_len = READ_LEN, floor = FLOOR,
                      k_max = KMAX, min_sites = 1L, file = file.path(OUT_DIR, "rapport.html"))

# 7. Résumé console
cat("\n=== KMerEnrich — 2 virus, sans background ===\n")
cat("Virus            :", paste(names(targets), collapse = ", "), "\n")
cat("Candidats / admis:", length(primers), "/", nrow(admit), "\n")
cat("Pool recommandé  : N* =", sb$n_star, "primers\n")
cat("Couverture globale:", round(100 * best$coverage), "%  | score:", round(score, 1), "/100\n")
cat("Couverture par virus:\n"); print(round(best$per_target, 3))
cat("\nSéquences à commander:\n"); print(pool)
cat("\nRapport interactif : ", file.path(normalizePath(OUT_DIR), "rapport.html"), "\n")
cat("Autres sorties dans :", normalizePath(OUT_DIR), "\n")
