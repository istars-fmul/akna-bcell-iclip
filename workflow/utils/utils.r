# Shared helpers for the downstream R analyses.

# Sum a run-length encoded BigWig track over each supplied GRanges interval.
# Callers select the appropriate strand before calling this function. Return
# one numeric value per interval in the same order, or numeric() for no sites.
# Input signal is already library-normalized; this does not perform rescaling.
sum_track_over_ranges <- function(track, ranges) {
  if (length(ranges) == 0) return(numeric())
  vapply(as.list(track[ranges]), sum, numeric(1), na.rm = TRUE)
}

# Convert list-valued metadata (e.g. multiple transcript-region boundaries) to
# semicolon-separated strings for CSV/Excel export. Scalar columns are retained.
flatten_list_columns <- function(data) {
  data[] <- lapply(data, function(column) {
    if (is.list(column)) {
      vapply(column, paste, collapse = ";", FUN.VALUE = character(1))
    } else {
      column
    }
  })
  data
}
