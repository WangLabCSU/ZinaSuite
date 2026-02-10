# Create ZinaSuite Hex Sticker Logo
# Uses hexSticker package with improved design

# Install required packages if needed
if (!requireNamespace("hexSticker", quietly = TRUE)) {
  install.packages("hexSticker")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(hexSticker)
library(ggplot2)

# Create output directory
output_dir <- "man/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create a simple, clean DNA helix visualization
set.seed(42)
n_points <- 50
t <- seq(0, 3*pi, length.out = n_points)

# First strand
x1 <- cos(t) * 0.8
y1 <- seq(-1.2, 1.2, length.out = n_points)

# Second strand (offset)
x2 <- cos(t + pi) * 0.8
y2 <- y1

# Create data frame
dna_data <- data.frame(
  x = c(x1, x2),
  y = c(y1, y2),
  strand = rep(c("A", "B"), each = n_points)
)

# Create base pairs (horizontal lines)
base_pairs <- data.frame(
  x = rep(NA, n_points),
  y = rep(NA, n_points),
  xend = rep(NA, n_points),
  yend = rep(NA, n_points)
)

for(i in seq(1, n_points, 3)) {
  base_pairs$x[i] <- x1[i]
  base_pairs$y[i] <- y1[i]
  base_pairs$xend[i] <- x2[i]
  base_pairs$yend[i] <- y2[i]
}
base_pairs <- base_pairs[!is.na(base_pairs$x), ]

# Create the plot with clean design
p <- ggplot() +
  # Draw strands
  geom_line(data = dna_data, aes(x = x, y = y, color = strand), 
            size = 1.5, lineend = "round") +
  # Draw base pairs
  geom_segment(data = base_pairs, 
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "#4ECDC4", size = 0.8, alpha = 0.7) +
  # Color scheme
  scale_color_manual(values = c("A" = "#FF6B6B", "B" = "#4ECDC4")) +
  # Clean theme
  theme_void() +
  theme_transparent() +
  theme(legend.position = "none") +
  coord_fixed(ratio = 0.6)

# Create the hex sticker with improved parameters
sticker(
  p,
  package = "ZinaSuite",
  p_size = 24,
  p_color = "#FFFFFF",
  p_family = "sans",
  p_fontface = "bold",
  s_x = 1.0,
  s_y = 0.78,
  s_width = 1.4,
  s_height = 1.1,
  h_fill = "#1A535C",
  h_color = "#4ECDC4",
  h_size = 1.8,
  spotlight = FALSE,
  url = "github.com/WangLabCSU/ZinaSuite",
  u_color = "#BEE9E8",
  u_size = 3.2,
  filename = file.path(output_dir, "logo.png"),
  dpi = 300
)

message("✓ Logo created successfully at: ", file.path(output_dir, "logo.png"))

# Verify the file was created
if (file.exists(file.path(output_dir, "logo.png"))) {
  file_info <- file.info(file.path(output_dir, "logo.png"))
  message("✓ File size: ", round(file_info$size / 1024, 2), " KB")
  message("✓ Dimensions: Check the file to verify it looks correct")
} else {
  message("✗ Error: Logo file was not created!")
}
