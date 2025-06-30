library(hexSticker)
library(ggplot2)

# Create a sample graph layout representing DIF "response network"
set.seed(123)
g <- sample_pa(30, directed = FALSE)
V(g)$group <- sample(c("reference", "focal"), 30, replace = TRUE)
V(g)$color <- ifelse(V(g)$group == "reference", "#0072B2", "#F0E442")

# Generate the graph plot
data(GMAT)
fit <- difNLR::difNLR(GMAT[, 1:20], GMAT$group, focal.name = 1, model = "3PL")
(g <- plot(fit, item = 2)[[1]] +
    scale_size(breaks = c(0, 15, 30, 45, 60, 150), range = c(0.5, 3.5)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        rect = element_blank(),
        axis.line = element_blank()
        ))


for (i in seq_along(g$layers)) {
  if (inherits(g$layers[[i]]$geom, "GeomLine")) {
    g$layers[[i]]$aes_params$linewidth <- 0.6  # Change this value to your desired width
  }
}


getAnywhere(plot.difNLR)


library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Comfortaa", "Comfortaa")
theme_get()

# Build the sticker
sticker(
  subplot = g,
  package = "difNLR",
  p_x = 1.33,
  p_y = 0.66,
  p_size = 35,
  p_color = "black",
  h_fill = "white",
  h_color = "black",
  p_family = "Comfortaa",
  s_x = 0.95,
  s_y = 1,
  s_width = 2,
  s_height = 1.45,
  dpi = 600,
  filename = "difNLR_hex_sticker.png"
)
