#!/usr/bin/env Rscript
# symbol_alias.R — resolve a reference set's HGNC symbols into this compartment's vocabulary.
# =============================================================================
# The defect this exists for. GSE160097 was quantified against a CellRanger hg19 reference,
# so the count matrix is frozen to that build's HGNC vintage: cGAS is MB21D1, STING is
# TMEM173, MARCHF5 is MARCH5, MRE11 is MRE11A. Reference gene sets ship current symbols.
# Every match in this compartment is an exact string match, so a renamed gene is dropped from
# a set silently and the loss reads as biological absence. In the Treg synovial-vs-blood
# ranked list TMEM173 sits at rank 265 and MB21D1 at 458 of 13,999 — the two genes that name
# the cGAS-STING axis are its two strongest members here, and every STING gene set misses
# them.
#
# The direction of resolution. The reference symbol is NEWER and the matrix symbol is OLDER,
# so a reference symbol is resolved DOWN into the vocabulary the data carries. The opposite
# traversal (stale matrix symbol -> current, for a tool that keys on current symbols) means
# something different under the one-to-one safety condition and belongs in its own function.
#
# The hazard that shapes the guard. Many retired symbols were reassigned as the official
# symbol of a DIFFERENT gene. PGF carries the alias PIGF, and PIGF now names a GPI-anchor
# biosynthesis gene; THPO carries TPO, and TPO now names thyroid peroxidase; ACOD1 carries
# CAD, which now names carbamoyl-phosphate synthetase. Accepting any of those attaches one
# gene's expression to another gene's set membership. So a candidate that is the official
# symbol of any other Entrez id is rejected, counted, and published beside the accepted pairs.
#
# Alias resolution is a correctness fix. Acceptance requires surviving the ownership guard,
# and the reporting contract is a ledger with a bucket per cause — matched, matched-via-alias,
# expression-filtered, below-detection, absent-from-reference.
#
# Provides:
#   build_alias_map(reference_symbols, matrix_vocabulary, db, flagged_pairs)
#   resolve_sets(sets, matrix_vocabulary, alias_map)
#   symbol_ledger(sets, alias_map, ranked_vocabulary, matrix_vocabulary, reference_vocabulary)
#   accepted_pairs(alias_map)   the reference_symbol -> matrix_symbol lookup, flagged excluded
#
# Source from the compartment root:
#   source("02_analysis/helpers/symbol_alias.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# The resolutions a candidate pair can carry. Every value except `accepted` withholds the
# pair, and `flagged_for_review` is the one that is withheld by human decision rather than
# by the automated guard.
ALIAS_RESOLUTIONS <- c(
  "accepted",
  "flagged_for_review",
  "rejected_symbol_belongs_to_another_gene",
  "rejected_multiple_aliases_in_vocabulary",
  "rejected_reference_symbol_ambiguous_in_org_db")

.alias_map_empty <- function() {
  tibble(reference_symbol = character(), matrix_symbol = character(),
         entrez_id = character(), n_aliases_in_vocabulary = integer(),
         resolution = character())
}

#' Candidate pairs mapping a reference symbol onto the vocabulary the data carries.
#'
#' A pair is a candidate when the reference symbol is absent from `matrix_vocabulary`,
#' resolves to exactly one Entrez id, and exactly one alias of that same Entrez id is
#' present in the vocabulary. Candidates are then filtered by the ownership guard.
#'
#' @param reference_symbols character  symbols as the reference set ships them.
#' @param matrix_vocabulary character  the target vocabulary to resolve into.
#' @param db  an org.*.eg.db object; org.Hs.eg.db for this compartment.
#' @param flagged_pairs character  "REF->MATRIX" strings withheld for human review.
#' @return tibble(reference_symbol, matrix_symbol, entrez_id, n_aliases_in_vocabulary,
#'   resolution), one row per candidate pair, carrying attributes `summary` (the
#'   no-candidate counts) and `n_rejected_ambiguous` (kept for the caller that reports it).
build_alias_map <- function(reference_symbols, matrix_vocabulary,
                            db = NULL, flagged_pairs = character()) {
  reference_symbols <- unique(reference_symbols[!is.na(reference_symbols) &
                                                  nzchar(reference_symbols)])
  matrix_vocabulary <- unique(matrix_vocabulary)
  missing_symbols <- sort(setdiff(reference_symbols, matrix_vocabulary))
  empty <- .alias_map_empty()
  summ <- list(n_reference_symbols = length(reference_symbols),
               n_absent_from_vocabulary = length(missing_symbols),
               n_not_in_org_db = NA_integer_,
               n_no_alias_in_vocabulary = NA_integer_,
               n_candidates = 0L, n_accepted = 0L, n_flagged = 0L, n_rejected = 0L)
  finish <- function(map) {
    summ$n_candidates <- nrow(map)
    summ$n_accepted   <- sum(map$resolution == "accepted")
    summ$n_flagged    <- sum(map$resolution == "flagged_for_review")
    summ$n_rejected   <- sum(!map$resolution %in% c("accepted", "flagged_for_review"))
    structure(map, summary = summ,
              n_rejected_ambiguous = sum(
                map$resolution == "rejected_symbol_belongs_to_another_gene"))
  }

  if (is.null(db)) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
      message("  org.Hs.eg.db unavailable, alias recovery skipped and 0 symbols recovered.")
      return(finish(empty))
    }
    db <- org.Hs.eg.db::org.Hs.eg.db
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(finish(empty))
  if (!length(missing_symbols)) return(finish(empty))

  # AnnotationDbi::select() throws outright when NONE of its keys is valid for the
  # keytype, so every call is prefiltered against the live keyspace and returns early
  # when nothing survives. With ~1,200 candidates this never fires; on a 7-gene gene set
  # it does, and TMEM173 arrives as an ALIAS where a SYMBOL is expected, which is what takes the
  # ownership guard down on GOBP_POSITIVE_REGULATION_OF_CGAS_STING_SIGNALING_PATHWAY.
  symbol_keys <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  sel <- function(keys, keytype, columns) {
    keys <- unique(keys[!is.na(keys) & nzchar(keys)])
    if (identical(keytype, "SYMBOL")) keys <- intersect(keys, symbol_keys)
    if (!length(keys)) return(NULL)
    suppressMessages(AnnotationDbi::select(db, keys = keys, keytype = keytype,
                                           columns = columns))
  }

  eg <- sel(missing_symbols, "SYMBOL", "ENTREZID")
  summ$n_not_in_org_db <- length(setdiff(missing_symbols, symbol_keys))
  if (is.null(eg)) return(finish(empty))
  eg <- eg[!is.na(eg$ENTREZID), , drop = FALSE]
  # A reference symbol carrying more than one Entrez id names more than one gene here,
  # so it is withheld — but only reported once a candidate for it actually exists.
  ambiguous <- unique(eg$SYMBOL[duplicated(eg$SYMBOL)])
  eg <- eg[!(eg$SYMBOL %in% ambiguous), , drop = FALSE]
  if (!nrow(eg)) return(finish(empty))

  al <- sel(unique(eg$ENTREZID), "ENTREZID", "ALIAS")
  if (is.null(al)) return(finish(empty))
  al <- al[al$ALIAS %in% matrix_vocabulary, , drop = FALSE]
  if (!nrow(al)) return(finish(empty))
  # Two aliases of one gene present in the vocabulary leaves no unique target, so the
  # pair is withheld with both names recorded, leaving no silent NA.
  cand <- al %>% group_by(.data$ENTREZID) %>%
    summarise(matrix_symbol = paste(sort(unique(.data$ALIAS)), collapse = "/"),
              n_aliases_in_vocabulary = n_distinct(.data$ALIAS), .groups = "drop")
  hits <- eg %>% inner_join(cand, by = "ENTREZID") %>%
    transmute(reference_symbol = .data$SYMBOL, matrix_symbol = .data$matrix_symbol,
              entrez_id = .data$ENTREZID,
              n_aliases_in_vocabulary = as.integer(.data$n_aliases_in_vocabulary))
  summ$n_no_alias_in_vocabulary <- length(setdiff(eg$SYMBOL, hits$reference_symbol))
  if (!nrow(hits)) return(finish(empty))

  multi <- hits %>% filter(.data$n_aliases_in_vocabulary > 1L) %>%
    mutate(resolution = "rejected_multiple_aliases_in_vocabulary")
  hits <- hits %>% filter(.data$n_aliases_in_vocabulary == 1L)
  if (!nrow(hits)) return(finish(multi))

  # The ownership guard: is the candidate the official symbol of some other gene?
  own <- sel(unique(hits$matrix_symbol), "SYMBOL", "ENTREZID")
  taken <- rep(FALSE, nrow(hits))
  if (!is.null(own)) {
    own <- own[!is.na(own$ENTREZID), , drop = FALSE]
    owner <- setNames(own$ENTREZID, own$SYMBOL)
    taken <- !is.na(owner[hits$matrix_symbol]) &
      owner[hits$matrix_symbol] != hits$entrez_id
    taken[is.na(taken)] <- FALSE
  }
  hits$resolution <- ifelse(taken, "rejected_symbol_belongs_to_another_gene", "accepted")
  # Pairs a human has to decide on are withheld at this point, so the
  # exclusion travels with the map and is visible in every consumer's ledger.
  flagged <- paste0(hits$reference_symbol, "->", hits$matrix_symbol) %in% flagged_pairs
  hits$resolution[flagged & hits$resolution == "accepted"] <- "flagged_for_review"

  amb <- if (length(ambiguous))
    tibble(reference_symbol = intersect(ambiguous, missing_symbols),
           matrix_symbol = NA_character_, entrez_id = NA_character_,
           n_aliases_in_vocabulary = NA_integer_,
           resolution = "rejected_reference_symbol_ambiguous_in_org_db")
  else .alias_map_empty()

  finish(bind_rows(hits, multi, amb) %>% arrange(.data$reference_symbol))
}

#' The reference_symbol -> matrix_symbol lookup an alias map licenses.
#'
#' Only `accepted` rows are returned; `flagged_for_review` is withheld by construction so
#' no consumer can apply a pair a human has not signed off on.
accepted_pairs <- function(alias_map) {
  if (is.null(alias_map) || !nrow(alias_map)) return(setNames(character(), character()))
  a <- alias_map[alias_map$resolution == "accepted", , drop = FALSE]
  setNames(a$matrix_symbol, a$reference_symbol)
}

#' Resolve every set in an fgsea-shaped `pathways` list into one vocabulary.
#'
#' @param sets named list of character vectors.
#' @param matrix_vocabulary character  the vocabulary the sets are being matched against.
#' @param alias_map tibble from build_alias_map(), or a named reference->matrix vector.
#' @return list(sets = resolved and de-duplicated, applied = per set x pair,
#'   collapsed = sets that gained fewer genes than pairs applied,
#'   many_to_one = two reference symbols resolving onto one matrix symbol)
resolve_sets <- function(sets, matrix_vocabulary, alias_map) {
  pairs <- if (is.data.frame(alias_map)) accepted_pairs(alias_map) else alias_map
  matrix_vocabulary <- unique(matrix_vocabulary)
  applied <- list(); collapsed <- list(); many_to_one <- list()

  out <- lapply(names(sets), function(nm) {
    g <- unique(sets[[nm]])
    hit <- names(pairs)[names(pairs) %in% g]
    hit <- hit[!hit %in% matrix_vocabulary]           # a set may carry both vintages
    tgt <- unname(pairs[hit])
    hit <- hit[tgt %in% matrix_vocabulary]
    tgt <- tgt[tgt %in% matrix_vocabulary]
    if (length(hit)) {
      applied[[nm]] <<- tibble(gene_set = nm, reference_symbol = hit, matrix_symbol = tgt)
      # Two reference symbols landing on one matrix symbol merges two set members into
      # one measurement. It is reported and kept separate.
      dup_tgt <- unique(tgt[duplicated(tgt)])
      if (length(dup_tgt))
        many_to_one[[nm]] <<- tibble(
          gene_set = nm, matrix_symbol = dup_tgt,
          reference_symbols = vapply(dup_tgt,
                                     function(s) paste(sort(hit[tgt == s]), collapse = "/"),
                                     character(1)))
    }
    res <- unique(c(g, tgt))
    # A set that does not grow by the number of pairs applied already carried the matrix
    # symbol under its other vintage, or two pairs collapsed onto one target.
    if (length(hit) && length(res) - length(g) != length(hit))
      collapsed[[nm]] <<- tibble(gene_set = nm, n_genes_before = length(g),
                                 n_genes_after = length(res),
                                 n_pairs_applied = length(hit),
                                 n_collapsed = length(hit) - (length(res) - length(g)))
    res
  })
  names(out) <- names(sets)
  list(sets = out, applied = bind_rows(applied), collapsed = bind_rows(collapsed),
       many_to_one = bind_rows(many_to_one))
}

#' The per-set symbol ledger: one bucket per cause, and they are kept apart.
#'
#' Conflating vocabulary loss with never-detected and with expression-filtered is the
#' whole defect, so each unmatched gene lands in exactly one bucket and the buckets close
#' against the set's size. The vocabulary layers are nested, and each is a different kind
#' of statement:
#'   reference_vocabulary  the CellRanger feature union — the only true "absent" is here
#'   matrix_vocabulary     post-QC gene_symbols.csv — below it is a detection fact
#'   ranked_vocabulary     post-filterByExpr — below it is a power statement about the contrast
#'
#' @return tibble, one row per set, closing on
#'   n_unique_set_genes == n_matched + n_matched_via_alias + n_alias_flagged_for_review +
#'     n_alias_rejected_ambiguous + n_expression_filtered + n_below_detection +
#'     n_absent_from_reference
#'
#' The buckets count REFERENCE genes, so they close against the set as the reference ships
#' it. `set_size_resolved` counts MATRIX symbols and is therefore smaller than
#' n_matched + n_matched_via_alias wherever pairs collapsed; `n_alias_collapsed` carries
#' that difference so the two never have to be reconciled by hand.
symbol_ledger <- function(sets, alias_map, ranked_vocabulary, matrix_vocabulary,
                          reference_vocabulary = NULL) {
  pairs <- accepted_pairs(alias_map)
  flag  <- alias_map[alias_map$resolution == "flagged_for_review", , drop = FALSE]
  rej   <- alias_map[alias_map$resolution ==
                       "rejected_symbol_belongs_to_another_gene", , drop = FALSE]
  ranked_vocabulary <- unique(ranked_vocabulary)
  matrix_vocabulary <- unique(matrix_vocabulary)
  have_ref <- !is.null(reference_vocabulary) && length(reference_vocabulary) > 0
  ref_vocab <- if (have_ref) unique(reference_vocabulary) else character()

  bind_rows(lapply(names(sets), function(nm) {
    g <- unique(sets[[nm]])
    matched <- intersect(g, ranked_vocabulary)
    rest <- setdiff(g, matched)
    # Precedence: a gene recoverable through an accepted pair is a vocabulary result and
    # is counted as such before any statement about detection or power is made.
    via <- rest[rest %in% names(pairs) & unname(pairs[rest]) %in% ranked_vocabulary]
    rest <- setdiff(rest, via)
    flagged <- rest[rest %in% flag$reference_symbol]
    rest <- setdiff(rest, flagged)
    ambig <- rest[rest %in% rej$reference_symbol]
    rest <- setdiff(rest, ambig)
    # The remaining buckets are read off the symbol the DATA would carry, so a retired
    # reference name whose matrix twin was dropped by the expression filter is reported as
    # expression-filtered, which is a different cause from never detected. Bucketing on the reference name
    # there would put a power statement in a detection bucket.
    eff <- ifelse(rest %in% names(pairs), unname(pairs[rest]), rest)
    expr_filtered <- rest[eff %in% matrix_vocabulary]
    rest <- setdiff(rest, expr_filtered)
    eff <- ifelse(rest %in% names(pairs), unname(pairs[rest]), rest)
    below <- if (have_ref) rest[eff %in% ref_vocab] else character()
    absent_ref <- if (have_ref) setdiff(rest, below) else rest
    tibble(
      gene_set = nm, n_unique_set_genes = length(g),
      n_matched = length(matched),
      n_matched_via_alias = length(via),
      n_alias_flagged_for_review = length(flagged),
      n_alias_rejected_ambiguous = length(ambig),
      n_expression_filtered = length(expr_filtered),
      n_below_detection = length(below),
      n_absent_from_reference = length(absent_ref),
      # The true resolved size, de-duplicated. Two things collapse here and both are
      # reported as its own count: a set carrying both vintages of one gene, and
      # several reference paralogs whose current names all resolve onto one matrix row
      # (NOTCH2NLA/B/C -> NOTCH2NL is the real case). Either way the set grows by fewer
      # genes than pairs applied, and `n_alias_collapsed` is that shortfall.
      set_size_resolved = length(unique(c(matched, unname(pairs[via])))),
      n_alias_collapsed = length(via) -
        (length(unique(c(matched, unname(pairs[via])))) - length(matched)),
      reference_vocabulary_available = have_ref,
      alias_pairs_applied = paste(sort(paste0(via, "->", unname(pairs[via]))),
                                  collapse = "/"),
      alias_pairs_flagged = paste(sort(paste0(
        flagged, "->", flag$matrix_symbol[match(flagged, flag$reference_symbol)])),
        collapse = "/"),
      alias_pairs_rejected = paste(sort(paste0(
        ambig, "->", rej$matrix_symbol[match(ambig, rej$reference_symbol)])),
        collapse = "/"))
  }))
}

#' Hard closure check on a ledger, asserted in-script on every run.
assert_ledger_closes <- function(ledger, label = "symbol ledger") {
  s <- with(ledger, n_matched + n_matched_via_alias + n_alias_flagged_for_review +
              n_alias_rejected_ambiguous + n_expression_filtered + n_below_detection +
              n_absent_from_reference)
  bad <- which(s != ledger$n_unique_set_genes)
  if (length(bad))
    stop(sprintf(paste0("[symbol_alias] %s does not close for %d set(s), e.g. %s ",
                        "(%d genes, buckets sum to %d). Every unmatched gene must land ",
                        "in exactly one cause bucket."),
                 label, length(bad), ledger$gene_set[bad[1]],
                 ledger$n_unique_set_genes[bad[1]], s[bad[1]]))
  invisible(TRUE)
}
