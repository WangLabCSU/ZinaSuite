# Create clean hexsticker logo for ZinaSuite
# Minimalist design with clear visual identity

if (!requireNamespace("hexSticker", quietly = TRUE)) {
  install.packages("hexSticker")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(hexSticker)
library(ggplot2)

output_dir <- "man/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create minimalist gene/DNA icon using simple geometric shapes
create_gene_icon <- function() {
  # Central horizontal bar (representing a gene/chromosome)
  center_y <- 0.5
  bar_height <- 0.08
  
  # Main horizontal bar
  segments <- data.frame(
    x = c(0.15, 0.85),
    xend = c(0.15, 0.85),
    y = c(center_y, center_y),
    yend = c(center_y, center_y),
    width = c(6, 6)
  )
  
  # Vertical tick marks representing exons/gene features
  ticks <- data.frame(
    x = seq(0.25, 0.75, length.out = 5),
    y = center_y - 0.15,
    yend = center_y + 0.15,
    width = 4
  )
  
  # Circular nodes at ends (representing start/end codons or data points)
  nodes <- data.frame(
    x = c(0.15, 0.85),
    y = c(center_y, center_y),
    size = c(8, 8)
  )
  
  # Data visualization dots (representing analysis results)
  data_points <- data.frame(
    x = c(0.35, 0.50, 0.65),
    y = c(0.75, 0.82, 0.78),
    size = c(3, 4, 3)
  )
  
  p <- ggplot() +
    # Main gene bar
    geom_segment(data = segments,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "#FFFFFF", linewidth = 6, lineend = "round") +
    # Exon/feature ticks
    geom_segment(data = ticks,
                 aes(x = x, y = y, xend = x, yend = yend),
                 color = "#00D4AA", linewidth = 4, lineend = "round") +
    # End nodes
    geom_point(data = nodes,
               aes(x = x, y = y),
               color = "#FFFFFF", size = 8) +
    geom_point(data = nodes,
               aes(x = x, y = y),
               color = "#FF6B6B", size = 5) +
    # Data visualization points
    geom_point(data = data_points,
               aes(x = x, y = y),
               color = "#FFD93D", size = data_points$size) +
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  return(p)
}

gene_plot <- create_gene_icon()

# Clean minimalist sticker design with GitHub URL
sticker(
  gene_plot,
  package = "ZinaSuite",
  p_size = 18,
  p_color = "#FFFFFF",
  p_family = "sans",
  p_fontface = "bold",
  p_y = 1.45,
  s_x = 1.0,
  s_y = 0.9,
  s_width = 0.9,
  s_height = 0.75,
  h_fill = "#1A365D",
  h_color = "#00D4AA",
  h_size = 1.5,
  spotlight = FALSE,
  url = "github.com/WangLabCSU/ZinaSuite",
  u_size = 3.5,
  u_color = "#00D4AA",
  u_x = 1.0,
  u_y = 0.06,
  filename = file.path(output_dir, "logo.png"),
  dpi = 300,
  white_around_sticker = FALSE
)

message("Logo created successfully at: ", file.path(output_dir, "logo.png"))
