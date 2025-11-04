#' Global variables
#'
#' @description Declare global variables to avoid R CMD check NOTEs
#' for ggplot2 NSE (non-standard evaluation) and tidyverse-style programming
#'
#' @noRd
#' @importFrom rlang .data
utils::globalVariables(c(
  # ggplot2 aesthetics
  "VAF",
  "Type",
  "Depth",
  "Clone",
  "Index",
  "Present",
  "vaf",
  "label"
))
