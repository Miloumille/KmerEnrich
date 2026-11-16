#!/usr/bin/env Rscript
# KMerEnrich v2.0 — analyse d'un panel de virus (nombre quelconque), background optionnel.
#
# Prérequis : le package KMerEnrich installé
#   devtools::install("chemin/vers/2026_KMerEnrich_2.0", build_vignettes = FALSE)
#
# Données attendues : un dossier `virus/` avec UN sous-dossier par virus
#   virus/
#   ├── virusA/  (une ou plusieurs souches .fasta)
#   ├── virusB/  ...
#   └── virusC/  ...
#
# Lancement :
#   Rscript analyze_viruses.R <dossier_virus> <dossier_sorties> [<dossier_background>]
#   (sans 3e argument -> null statistique par défaut, aucun background requis)

suppressMessages(library(KMerEnrich))
suppressMessages(library(Biostrings))

# ---- Paramètres (à calibrer sur tes données) --------------------------------
args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[[1]] else "virus"
OUT_DIR  <- if (length(args) >= 2) args[[2]] else "kmerenrich_out"
BG_DIR   <- if (length(args) >= 3) args[[3]] else ""     # optionnel : dossier background
K              <- 9L
MAX_DEGENERATE <- 3L
READ_LEN       <- 200L
BIN_SIZE       <- 100L
FLOOR          <- 0.5    # couverture minimale garantie PAR virus (équité panel)
MIN_ENRICH     <- 2
TM_MIN <- 0; TM_MAX <- 60; GC_MIN <- 0.2; GC_MAX <- 0.8
MAX_CANDIDATES <- 150L   # primers max envoyés à l'ILP (au-delà de ~250, l'optim explose)
KMAX           <- 6L     # taille de pool max évaluée (courbe de saturation + curseur N : 1..6)
MIN_SITES      <- { v <- Sys.getenv("KE_MIN_SITES"); if (nzchar(v)) as.integer(v) else 2L }
                         # nb minimal de sites de fixation exigé sur CHAQUE génome de
                         # référence (surchargeable : KE_MIN_SITES=5). Écarte les primers
                         # à trop peu de loci (peu utiles).
# -----------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Cibles (1..N virus)
targets <- load_targets(DATA_DIR)
stopifnot(length(targets) >= 1)
n_str <- sum(vapply(targets, length, integer(1)))
if (n_str > 500)
  message(sprintf("[info] %d souches au total : le calcul de couverture prend ~1 min par ~1000 souches.", n_str))

# 2. Background : dossier fourni -> fit en streaming ; sinon null statistique
bg_files <- if (nzchar(BG_DIR)) list.files(BG_DIR, recursive = TRUE, full.names = TRUE) else character(0)
model <- if (length(bg_files)) markov_background_files(bg_files) else default_background()

# 3. Candidats : k-mers conservés de CHAQUE virus, fusionnés en primers IUPAC
kmers   <- unique(unlist(lapply(targets, enumerate_conserved_kmers, k = K), use.names = FALSE))
primers <- collapse_to_iupac(kmers, max_degenerate = MAX_DEGENERATE, anchor = 2L)

# 4. Enrichissement + filtre thermo/enrichissement
viral  <- do.call(c, unname(targets))
scored <- score_enrichment(primers, viral, model)
admit  <- filter_candidates(scored, min_enrichment = MIN_ENRICH,
                            tm_min = TM_MIN, tm_max = TM_MAX, gc_min = GC_MIN, gc_max = GC_MAX)
stopifnot(nrow(admit) > 0)

# Contrainte panel (appliquée à TOUT le vivier, avant toute présélection) : chaque primer
# retenu doit se fixer sur >= MIN_SITES loci de CHAQUE génome de référence. Écarte ainsi
# (i) les primers absents d'un des virus et (ii) les primers à site unique.
if (MIN_SITES > 1L) {
  refs  <- lapply(targets, function(x) x[1])
  nsite <- vapply(admit$pattern, function(p)
    min(vapply(refs, function(r)
      as.integer(sum(Biostrings::vcountPattern(p, r, fixed = c(pattern = FALSE, subject = TRUE)))),
      integer(1))), integer(1))
  message(sprintf("[info] %d/%d primers ecartes (< %d sites sur au moins un genome de reference).",
                  sum(nsite < MIN_SITES), nrow(admit), MIN_SITES))
  admit <- admit[nsite >= MIN_SITES, ]
  stopifnot(nrow(admit) > 0)
}

# Sélection des candidats pour le solveur max-min.
# L'enrichissement est dominé par le virus le + représenté (Z : 498 souches vs 45) ;
# garder le top-N par enrichissement AFFAME les petits virus -> leur couverture plafonne.
# On prend donc un large vivier enrichi, on calcule sa couverture, puis on choisit
# gloutonnement les primers qui PAVENT le mieux TOUS les génomes (priorité au virus le
# moins couvert). Le solveur max-min peut alors réellement équilibrer les deux.
PRESELECT <- 400L
pre     <- head(admit[order(-admit$enrichment), ], min(PRESELECT, nrow(admit)))
cov_all <- coverage_matrix(pre$pattern, targets, read_len = READ_LEN, bin_size = BIN_SIZE)

balanced_greedy <- function(m, bt, enrich, k) {
  vs <- sort(unique(bt)); covered <- rep(FALSE, nrow(m))
  avail <- seq_len(ncol(m)); chosen <- integer(0)
  per <- function() vapply(vs, function(t) { idx <- which(bt == t); mean(covered[idx]) }, numeric(1))
  repeat {
    if (length(chosen) >= k || length(avail) == 0) break
    v   <- vs[which.min(per())]                       # virus actuellement le moins couvert
    idx <- which(bt == v & !covered)
    gain <- if (length(idx)) colSums(m[idx, avail, drop = FALSE]) else rep(0, length(avail))
    if (max(gain) == 0) { gain <- colSums(m[!covered, avail, drop = FALSE]); if (max(gain) == 0) break }
    top <- avail[gain == max(gain)]; j <- top[which.max(enrich[top])]
    chosen <- c(chosen, j); covered <- covered | (m[, j] > 0); avail <- setdiff(avail, j)
  }
  if (length(chosen) < k) chosen <- c(chosen, head(setdiff(order(-enrich), chosen), k - length(chosen)))
  sort(unique(chosen))
}
enr_pre <- pre$enrichment[match(cov_all$patterns, pre$pattern)]
sel     <- balanced_greedy(cov_all$matrix, cov_all$bin_target, enr_pre, MAX_CANDIDATES)
admit   <- pre[match(cov_all$patterns[sel], pre$pattern), ]
cov     <- list(matrix = cov_all$matrix[, sel, drop = FALSE], bin_target = cov_all$bin_target,
                patterns = cov_all$patterns[sel], targets = cov_all$targets)
ceil_per <- vapply(sort(unique(cov$bin_target)), function(t) {
  idx <- which(cov$bin_target == t); mean(rowSums(cov$matrix[idx, , drop = FALSE]) > 0) }, numeric(1))
message(sprintf("[info] vivier %d -> %d candidats couvrants (plafond par virus atteignable: %s)",
                nrow(pre), length(sel),
                paste(sprintf("%s=%.0f%%", cov$targets[sort(unique(cov$bin_target))], 100*ceil_per), collapse=" ")))

# 5. Pool ÉQUILIBRÉ optimal pour chaque N (max-min, toujours faisable)
sb   <- saturation_balanced(cov, k_max = min(ncol(cov$matrix), KMAX))
best <- sb$pools[[sb$n_star]]

# 6. Score global + séquences à commander (complément inverse)
enr    <- mean(admit$enrichment[admit$pattern %in% best$patterns])
score  <- desirability(enr, best$coverage, best$regularity)
orders <- as.character(reverseComplement(DNAStringSet(best$patterns)))
pool   <- data.frame(primer = best$patterns, a_commander = orders,
                     enrichissement = round(admit$enrichment[match(best$patterns, admit$pattern)], 1))

# 7. Sortie : un unique rapport HTML interactif (export CSV intégré dans la page)
report <- kmer_report(targets, admit, cov = cov, read_len = READ_LEN, floor = FLOOR,
                      k_max = KMAX, min_sites = 1L, file = file.path(OUT_DIR, "rapport.html"))

# 8. Résumé console
cat("\n=== KMerEnrich — panel de", length(targets), "virus",
    if (length(bg_files)) "(avec background)" else "(sans background, null par défaut)", "===\n")
cat("Virus            :", paste(names(targets), collapse = ", "), "\n")
cat("Candidats / admis:", length(primers), "/", nrow(admit), "\n")
cat("Pool recommandé  : N* =", sb$n_star, "primers\n")
cat("Couverture globale:", round(100 * best$coverage), "%  | score:", round(score, 1), "/100\n")
cat("Couverture par virus:\n"); print(round(best$per_target, 3))
cat("\nSéquences à commander:\n"); print(pool)
cat("\nRapport (fichier unique) :", normalizePath(report), "\n")
