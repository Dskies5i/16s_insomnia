#utils_plot.R
library(fs)
library(ggplot2)

save_plot <- function(p, path, w = 8, h = 5, dpi = 300) {
  dir_create(path_dir(path))
  ggsave(filename = path, plot = p, width = w, height = h, dpi = dpi)
}