plot.prediction2 <- function(input, ...){
  ggplot() +
    geom_polygon(data = input, aes(x=x, y), colour = "white", fill = "white", size = 0.25) +
    geom_raster(input,mapping = aes(x = x, y = y, fill = maxent)) +
    scale_fill_viridis(option = "B", guide = guide_colourbar(title = "Suitability"),na.value = "white") +
    coord_fixed() + 
    xlab("Longitude") +
    ylab("Latitude") +
    theme(
      axis.line = element_line(),
      axis.text.y = element_text(colour = "black", size = 10),
      axis.text.x = element_text(colour = "black", size = 10),
      legend.justification = c(0, 1), legend.position = c(.90, .35),
      panel.background = element_blank()
    )
}

plot.enmtools.maxent2 <- function(x, ...){
  suit.points <- data.frame(rasterToPoints(x$suitability))
  colnames(suit.points) <- c("Longitude", "Latitude", "Suitability")
  suit.plot <- ggplot(data = suit.points, aes_string(y = "Latitude", x = "Longitude")) +
    geom_raster(aes_string(fill = "Suitability")) +
    scale_fill_viridis(option = "B", guide = guide_colourbar(title = "Suitability")) +
    coord_fixed() + theme_classic()
}