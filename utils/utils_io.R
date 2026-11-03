#util_io.R

library(fs)
library(readr)

write_csv_safe <- function(df, path) {
  dir_create(path_dir(path))
  readr::write_csv(df, path)
}

save_rds_safe <- function(obj, path) {
  dir_create(path_dir(path))
  saveRDS(obj, path)
}

write_text_safe <- function(lines, path) {
  dir_create(path_dir(path))
  writeLines(lines, con = path, sep = "\n")
}