# Shared helpers for the downstream R analyses.

sum_track_over_ranges <- function(track, ranges) {
  if (length(ranges) == 0) return(numeric())
  vapply(as.list(track[ranges]), sum, numeric(1), na.rm = TRUE)
}

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
