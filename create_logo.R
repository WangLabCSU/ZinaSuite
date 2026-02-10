# Create ZinaSuite Hex Sticker Logo
# Uses hexSticker package

# Install required packages if needed
if (!requireNamespace("hexSticker", quietly = TRUE)) {
  install.packages("hexSticker")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("showtext", quietly = TRUE)) {
  install.packages("showtext")
}

library(hexSticker)
library(ggplot2)
library(showtext)

# Create output directory
output_dir <- "man/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create a DNA helix visualization for the sticker
set.seed(42)
n_points <- 100
t <- seq(0, 4*pi, length.out = n_points)
x <- cos(t)
y <- sin(t) + seq(-2, 2, length.out = n_points)
z <- t

# Create data frame for DNA strands
dna_data <- data.frame(
  x = c(x, x + 0.3),
  y = c(y, y),
  strand = rep(c("Strand1", "Strand2"), each = n_points),
  base_pair = rep(1:n_points, 2)
)

# Create the plot
p <- ggplot(dna_data, aes(x = x, y = y, color = strand)) +
  geom_line(size = 1.2) +
  geom_point(data = dna_data[seq(1, nrow(dna_data), 5), ], 
             aes(x = x, y = y), size = 2) +
  scale_color_manual(values = c("Strand1" = "#2E86AB", "Strand2" = "#A23B72")) +
  theme_void() +
  theme_transparent() +
  theme(legend.position = "none")

# Create the hex sticker
sticker(
  p,
  package = "ZinaSuite",
  p_size = 22,
  p_color = "#FFFFFF",
  p_family = "sans",
  p_fontface = "bold",
  s_x = 1.0,
  s_y = 0.75,
  s_width = 1.3,
  s_height = 1.0,
  h_fill = "#1B4965",
  h_color = "#5FA8D3",
  h_size = 1.5,
  spotlight = TRUE,
  l_x = 0.8,
  l_y = 0.7,
  l_width = 2,
  l_height = 2,
  l_alpha = 0.3,
  url = "github.com/WangLabCSU/ZinaSuite",
  u_color = "#BEE9E8",
  u_size = 3.5,
  filename = file.path(output_dir, "logo.png"),
  dpi = 300
)

message("Logo created successfully at: ", file.path(output_dir, "logo.png"))
