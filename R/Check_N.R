#' Check proportion of N bases in each sequence.
#'
#' Calculates the proportion of "N" bases (ambiguous bases) in each sequence and flags if above the given threshold.
#'
#' @param rep_seqs DNAStringSet or character vector of representative sequences.
#' @param cutoff Numeric, percent threshold (default 1).
#' @return Data frame with columns: qseqid, N_percent, N_flag.
#' @examples
#' seqs <- Biostrings::DNAStringSet(c(seq1 = "ATGCNNNN", seq2 = "NNNNATGC"))
#' check_N(seqs)
#' check_N(seqs, cutoff = 10)
#' @export
check_N <- function(rep_seqs, cutoff = 1) {
  seqs <- as.character(rep_seqs)
  n_perc <- sapply(seqs, function(s) sum(strsplit(s, NULL)[[1]] == "N") / nchar(s) * 100)
  data.frame(qseqid = names(rep_seqs), N_percent = n_perc, N_flag = n_perc > cutoff)
}