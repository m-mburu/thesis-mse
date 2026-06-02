remove_utf8_bom <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  size <- file.info(path)$size
  bytes <- readBin(path, what = "raw", n = size)
  bom <- as.raw(c(0xEF, 0xBB, 0xBF))

  has_bom <- length(bytes) >= 3 && identical(bytes[1:3], bom)

  if (has_bom) {
    con <- file(path, "wb")
    on.exit(close(con), add = TRUE)
    writeBin(bytes[-(1:3)], con)
    message("Removed UTF-8 BOM from: ", path)
  } else {
    message("No UTF-8 BOM found in: ", path)
  }
}

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) > 0) args[1] else here::here("report/methodology.qmd")

remove_utf8_bom(target)


#quarto::quarto_render(target)