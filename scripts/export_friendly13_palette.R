# Generate the exact ordered cols4all::friendly13 palette for GitHub Pages.
# Package version is pinned to 0.11 so the published colors are reproducible.

required_version <- "0.11"
repos <- c("https://cols4all.r-universe.dev", "https://cloud.r-project.org")

if (
  !requireNamespace("cols4all", quietly = TRUE) ||
  as.character(utils::packageVersion("cols4all")) != required_version
) {
  install.packages("cols4all", repos = repos)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = repos)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  install.packages("digest", repos = repos)
}

installed_version <- as.character(utils::packageVersion("cols4all"))
if (!identical(installed_version, required_version)) {
  stop("Expected cols4all ", required_version, "; found ", installed_version, ".")
}

colors <- unname(as.character(
  cols4all::c4a(
    palette = "friendly13",
    n = 13,
    type = "cat",
    format = "hex",
    verbose = FALSE
  )
))

if (
  length(colors) != 13L ||
  anyNA(colors) ||
  any(!grepl("^#[0-9A-Fa-f]{6}$", colors)) ||
  length(unique(toupper(colors))) != 13L
) {
  stop("cols4all did not return exactly 13 unique HEX colors for friendly13.")
}

colors <- toupper(colors)
palette_hash <- digest::digest(
  paste(colors, collapse = "|"),
  algo = "sha256",
  serialize = FALSE
)

payload <- list(
  palette = "cols4all::friendly13",
  palette_type = "cat",
  color_count = 13L,
  colors = as.list(colors),
  generated = TRUE,
  package = "cols4all",
  package_version = installed_version,
  source_call = 'cols4all::c4a(palette = "friendly13", n = 13, type = "cat", format = "hex")',
  ordered_hex_sha256 = palette_hash
)

jsonlite::write_json(
  payload,
  path = "data/friendly13_palette.json",
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("Generated exact cols4all friendly13 palette, version ", installed_version, "\n", sep = "")
cat(paste(sprintf("%02d  %s", seq_along(colors), colors), collapse = "\n"), "\n")
cat("Ordered HEX SHA256: ", palette_hash, "\n", sep = "")
