#######################################################################################################
# MAP FIGURE SYSTEM
#######################################################################################################

# Define patch type colors
patch_colours_map <- c(
  "Hedge" = "#FFA500",        # Hedge
  "Grass" = "#00FF00",        # Grass (lawn)
  "Tree" = "#8000FF",         # Trees
  "Water" = "#0000FF",        # Water
  "Roof (user mapped)" = "#FF0000",         # Roof (building)
  "Roof (imported)" = "#FF0000",         # Roof (building)
  "Artificial surface" = "#FF00FF",        # Artificial surface
  "Bush" = "#FFFF00",         # Bush
  "Cultivated bed" = "#00FFFF",    # Cultivated bed
  "Sparrow observation" = "#000000"
)

plot_map_figure <- function(figure_data,
                           column_widths = c(1, 1, 0.7)) {

  # Load required libraries
  library(ggplot2)
  library(png)
  library(cowplot)
  library(magick)

  # Order data by discrimination score
  figure_data <- figure_data[order(figure_data$discrimination_score), ]

  # Create legend data with linetype and linewidth as columns
  legend_df <- data.frame(
    x = 1:length(patch_colours_map),
    y = 1:length(patch_colours_map),
    category = factor(names(patch_colours_map), levels = names(patch_colours_map)),
    stringsAsFactors = FALSE
  )

  # Add linetype and linewidth columns
  legend_df$linetype <- "solid"
  legend_df$linetype[legend_df$category == "Roof (imported)"] <- "dotted"
  legend_df$linetype[legend_df$category == "Sparrow observation"] <- "dashed"

  legend_df$linewidth_val <- 0.5
  legend_df$linewidth_val[legend_df$category == "Sparrow observation"] <- 1.2

  legend_plot <- ggplot(legend_df, aes(x = x, y = y)) +
    geom_tile(aes(fill = category, color = category),
              alpha = 0.3,
              linetype = legend_df$linetype,
              linewidth = legend_df$linewidth_val) +
    scale_fill_manual(values = patch_colours_map, name = "Land Cover Type") +
    scale_color_manual(values = patch_colours_map, name = "Land Cover Type") +
    scale_linetype_manual(
      values = c("solid" = "solid", "dotted" = "dotted", "dashed" = "22"),
      guide = "none"
    ) +
    guides(fill = guide_legend(override.aes = list(
      alpha = 0.3,
      linetype = c("solid", "solid", "solid", "solid", "solid",
                   "dotted", "solid", "solid", "solid", "22"),
      linewidth = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.2)
    ))) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 13),
      legend.key.size = unit(1.5, "cm"),
      legend.key.spacing.y = unit(0.5, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )

  # Extract and save legend
  legend_grob <- suppressWarnings(get_legend(legend_plot))
  ggsave("./maps/legend_temp_original.png",
         plot = wrap_elements(legend_grob),
         width = 5, height = 10, units = "in", dpi = 300, bg = "white")

  # Crop the legend and add asymmetric border (more on left)
  legend_img_magick <- image_read("./maps/legend_temp_original.png")
  legend_img_magick <- image_trim(legend_img_magick)
  legend_img_magick <- image_border(legend_img_magick, "white", "20x20")
  legend_info <- image_info(legend_img_magick)
  legend_img_magick <- image_extent(legend_img_magick,
                                    geometry = paste0((legend_info$width + 80), "x", legend_info$height, "+80+0"),
                                    gravity = "east")
  image_write(legend_img_magick, "./maps/legend_temp.png")

  # Read the cropped legend
  legend_img <- readPNG("./maps/legend_temp.png")
  legend_height <- nrow(legend_img)
  legend_width <- ncol(legend_img)
  legend_aspect <- legend_height / legend_width

  # Calculate positioning
  total_width <- sum(column_widths)
  col1_width <- column_widths[1] / total_width
  col2_width <- column_widths[2] / total_width
  legend_width_prop <- column_widths[3] / total_width

  col1_end <- col1_width
  col2_end <- col1_width + col2_width

  # Function to calculate centered aspect-ratio-preserving bounds
  calc_bounds <- function(xmin, xmax, ymin, ymax, aspect_ratio, scale = 0.9) {
    available_width <- (xmax - xmin) * scale
    available_height <- (ymax - ymin) * scale
    available_aspect <- available_height / available_width

    x_center <- (xmin + xmax) / 2
    y_center <- (ymin + ymax) / 2

    if (available_aspect > aspect_ratio) {
      # Too tall, adjust height
      new_height <- available_width * aspect_ratio
      list(xmin = x_center - available_width/2, xmax = x_center + available_width/2,
           ymin = y_center - new_height/2, ymax = y_center + new_height/2)
    } else {
      # Too wide, adjust width
      new_width <- available_height / aspect_ratio
      list(xmin = x_center - new_width/2, xmax = x_center + new_width/2,
           ymin = y_center - available_height/2, ymax = y_center + available_height/2)
    }
  }

  # Read all map images and calculate bounds
  row_height <- 1/3
  maps <- list()
  text_positions <- list()

  # Fine-tuned positions for bottom row
  positions <- list(
    list(row = 3, col = 1, xmin = 0, xmax = col1_end, scale = 0.9),
    list(row = 3, col = 2, xmin = col1_end, xmax = col2_end, scale = 0.9),
    list(row = 2, col = 1, xmin = 0, xmax = col1_end, scale = 0.9),
    list(row = 2, col = 2, xmin = col1_end, xmax = col2_end, scale = 0.9),
    list(row = 1, col = 1, xmin = 0, xmax = col1_end * 0.88, scale = 0.88),  # Slightly bigger
    list(row = 1, col = 2, xmin = col1_end * 0.84, xmax = col2_end, scale = 1.02)  # Smaller than before
  )

  for (i in 1:6) {
    img <- readPNG(paste0("./maps/", figure_data$garden[i], ".png"))
    img_aspect <- nrow(img) / ncol(img)

    pos <- positions[[i]]
    ymin <- (pos$row - 1) * row_height
    ymax <- pos$row * row_height

    bounds <- calc_bounds(pos$xmin, pos$xmax, ymin, ymax, img_aspect, scale = pos$scale)
    maps[[i]] <- list(img = img, bounds = bounds)

    text_positions[[i]] <- list(
      x = (bounds$xmin + bounds$xmax) / 2,
      y = bounds$ymax + 0.005,
      label = sprintf("Brier discrimination: %.2f", figure_data$discrimination_score[i])
    )
  }

  # Calculate legend bounds with aspect ratio
  legend_bounds <- calc_bounds(col2_end, 1, 0, 1, legend_aspect, scale = 0.98)

  # Create single plot with all images positioned
  final_plot <- ggplot() +
    annotation_raster(maps[[1]]$img,
                     xmin = maps[[1]]$bounds$xmin, xmax = maps[[1]]$bounds$xmax,
                     ymin = maps[[1]]$bounds$ymin, ymax = maps[[1]]$bounds$ymax) +
    annotation_raster(maps[[2]]$img,
                     xmin = maps[[2]]$bounds$xmin, xmax = maps[[2]]$bounds$xmax,
                     ymin = maps[[2]]$bounds$ymin, ymax = maps[[2]]$bounds$ymax) +
    annotation_raster(maps[[3]]$img,
                     xmin = maps[[3]]$bounds$xmin, xmax = maps[[3]]$bounds$xmax,
                     ymin = maps[[3]]$bounds$ymin, ymax = maps[[3]]$bounds$ymax) +
    annotation_raster(maps[[4]]$img,
                     xmin = maps[[4]]$bounds$xmin, xmax = maps[[4]]$bounds$xmax,
                     ymin = maps[[4]]$bounds$ymin, ymax = maps[[4]]$bounds$ymax) +
    annotation_raster(maps[[5]]$img,
                     xmin = maps[[5]]$bounds$xmin, xmax = maps[[5]]$bounds$xmax,
                     ymin = maps[[5]]$bounds$ymin, ymax = maps[[5]]$bounds$ymax) +
    annotation_raster(maps[[6]]$img,
                     xmin = maps[[6]]$bounds$xmin, xmax = maps[[6]]$bounds$xmax,
                     ymin = maps[[6]]$bounds$ymin, ymax = maps[[6]]$bounds$ymax) +
    annotation_raster(legend_img,
                     xmin = legend_bounds$xmin, xmax = legend_bounds$xmax,
                     ymin = legend_bounds$ymin, ymax = legend_bounds$ymax) +
    annotate("text", x = text_positions[[1]]$x, y = text_positions[[1]]$y,
             label = text_positions[[1]]$label, size = 3, fontface = "bold", vjust = 0) +
    annotate("text", x = text_positions[[2]]$x, y = text_positions[[2]]$y,
             label = text_positions[[2]]$label, size = 3, fontface = "bold", vjust = 0) +
    annotate("text", x = text_positions[[3]]$x, y = text_positions[[3]]$y,
             label = text_positions[[3]]$label, size = 3, fontface = "bold", vjust = 0) +
    annotate("text", x = text_positions[[4]]$x, y = text_positions[[4]]$y,
             label = text_positions[[4]]$label, size = 3, fontface = "bold", vjust = 0) +
    annotate("text", x = text_positions[[5]]$x, y = text_positions[[5]]$y,
             label = text_positions[[5]]$label, size = 3, fontface = "bold", vjust = 0) +
    annotate("text", x = text_positions[[6]]$x, y = text_positions[[6]]$y,
             label = text_positions[[6]]$label, size = 3, fontface = "bold", vjust = 0) +
    xlim(0, 1) +
    ylim(0, 1) +
    coord_fixed(clip = "off") +
    theme_void()

  # Save the final plot
  ggsave("./maps/map_figure_temp.png", plot = final_plot, bg = "white",
         width = 6.27, height = 6.27, units = "in", dpi = 300)

  # Trim the final output and add back 20px border
  final_img <- image_read("./maps/map_figure_temp.png")
  final_img <- image_trim(final_img)
  final_img <- image_border(final_img, "white", "20x20")
  image_write(final_img, "./maps/map_figure.png")

  return(final_plot)
}