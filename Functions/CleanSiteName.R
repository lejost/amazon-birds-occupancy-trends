# --------------------------------------
# FUNCTION: CleanSiteName
# required packages: none
# description: this function cleans and standardizes sites names.
# required objects: none
########################################
CleanSiteName <- function(x) {
  # First two letters should be uppercase
  first_two <- toupper(substr(x, 1, 2))
  # Third letter should be lowercase
  third_letter <- tolower(substr(x, 3, 3))
  # Middle part (from 4th character to second-to-last)
  middle_part <- substr(x, 4, nchar(x) - 1)
  # Last letter: lowercase if 'x', uppercase otherwise
  last_letter <- substr(x, nchar(x), nchar(x))
  if (tolower(last_letter) == "x") {
    last_letter <- tolower(last_letter)
  } else {
    last_letter <- toupper(last_letter)
  }
  # Combine parts to form the cleaned site name
  paste0(first_two, third_letter, middle_part, last_letter)
}
