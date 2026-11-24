# Generate the exact cols4all friendly13 palette for this static website.
# Run once from the website repository root:
#   Rscript scripts/export_friendly13_palette.R
#
# The generated JSON is committed with the site. The browser does not run R.

if (!requireNamespace("cols4all", quietly = TRUE)) {
  install.packages(
    "cols4all",
    repos = c("https://cols4all.r-universe.dev", "https://cloud.r-project.org")
  )
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}

palette <- unname(as.character(
  cols4all::c4a(
    palette = "friendly13",
    n = 13,
    type = "cat",
    format = "hex",
    verbose = FALSE
  )
))

if (
  length(palette) != 13L ||
  anyNA(palette) ||
  any(!grepl("^#[0-9A-Fa-f]{6}$", palette)) ||
  length(unique(toupper(palette))) != 13L
) {
  stop("Expected exactly 13 unique HEX colors from cols4all::c4a('friendly13', 13).")
}

payload <- list(
  palette = "cols4all::friendly13",
  palette_type = "cat",
  color_count = 13L,
  colors = as.list(toupper(palette)),
  generated = TRUE,
  package = "cols4all",
  package_version = as.character(utils::packageVersion("cols4all")),
  source_call = 'cols4all::c4a(palette = "friendly13", n = 13, type = "cat", format = "hex")'
)

jsonlite::write_json(
  payload,
  path = "data/friendly13_palette.json",
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("Generated cols4all friendly13 palette\n")
cat("cols4all version: ", payload$package_version, "\n", sep = "")
cat(paste(sprintf("%02d  %s", seq_along(palette), toupper(palette)), collapse = "\n"), "\n")
