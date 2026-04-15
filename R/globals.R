# Suppress R CMD check NOTEs for variables that are genuinely used before
# assignment in conditional branches (checked with exists()) across
# pathway analysis helpers.
#
# All dplyr/ggplot2 NSE column-name NOTEs have been resolved by replacing
# bare names with .data$col (dplyr) or ggplot2::aes(.data$col) (ggplot2).
#
# The remaining entries below are true forward-references: the variable
# may or may not exist at runtime (guarded by exists()), so R CMD check
# cannot verify the binding.
utils::globalVariables(c(
  "pathway_report"   # read from disk in pathway helpers, guarded by exists()
))
