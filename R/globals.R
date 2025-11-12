# This tells R CMD check that these variables are used intentionally
# inside dplyr/ggplot2 pipes and are not global variables.
utils::globalVariables(c(
  "id",
  "y",
  "tx",
  "time",
  "yprev",
  "last_not_home",
  "dead",
  "drs"
))
