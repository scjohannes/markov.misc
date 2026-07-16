sample_categorical_rows <- function(probabilities) {
  probabilities <- as.matrix(probabilities)
  storage.mode(probabilities) <- "double"
  cpp_sample_categorical_rows(probabilities)
}
