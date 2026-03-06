library(ggplot2)
library(hexSticker)
library(showtext)

font_add_google("Lato", "lato")
showtext_auto()

# --- DAG subplot ---
# Node positions (scaled to fit nicely in subplot area)
nodes <- data.frame(
  label = c("D", "Y", "X"),
  x     = c(-1.1, 1.1, 0),
  y     = c(0.6, 0.6, -0.7),
  size  = c(14, 14, 16)
)

# Edge: D->Y (normal white)
# Edge: D->X (red border, white center — the "bad" connection)

p <- ggplot() +
  # D->Y edge (normal)
  geom_segment(aes(x = -1.1, y = 0.6, xend = 1.1, yend = 0.6),
               color = "white", linewidth = 0.8, alpha = 0.8) +
  # X->Y edge (normal)
  geom_segment(aes(x = 0, y = -0.7, xend = 1.1, yend = 0.6),
               color = "white", linewidth = 0.8, alpha = 0.8) +
  # D->X edge: red border (thick red underneath)
  geom_segment(aes(x = -1.1, y = 0.6, xend = 0, yend = -0.7),
               color = "#FF2D20", linewidth = 3.5) +
  # D->X edge: white center (thinner white on top)
  geom_segment(aes(x = -1.1, y = 0.6, xend = 0, yend = -0.7),
               color = "white", linewidth = 1.2) +
  # Nodes (circles)
  geom_point(data = nodes, aes(x = x, y = y),
             color = "#E05A47", fill = "#E05A47", shape = 21, size = 12, stroke = 1.2) +
  # Node labels
  geom_text(data = nodes, aes(x = x, y = y, label = label),
            color = "white", size = 8, fontface = "bold", family = "lato") +
  coord_fixed(xlim = c(-2, 2), ylim = c(-1.5, 1.3)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# --- Hex sticker ---
sticker(
  subplot   = p,
  package   = "badcontrols",
  # Package name: bigger, inside the hex
  p_size    = 18,
  p_color   = "#E8E8E8",
  p_family  = "lato",
  p_fontface = "bold",
  p_y       = 1.42,
  # Subplot position and size
  s_x       = 1,
  s_y       = 0.72,
  s_width   = 1.6,
  s_height  = 1.2,
  # Hex colors
  h_fill    = "#2C3E50",
  h_color   = "#E05A47",
  h_size    = 2.5,
  # White border
  white_around_sticker = FALSE,
  # Output
  filename  = "inst/figures/hexsticker.png",
  dpi       = 300
)
