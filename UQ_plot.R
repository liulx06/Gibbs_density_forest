rm(list=ls())

library(ggplot2)

## Split the sample
# forest temperture 2 or 10.
setwd("/Users/linxiliu/Documents/1Dbeta_2026Feb24/UQ/1Dbeta1e4_ratio6_inferP2/UQ/")
setwd("/Users/linxiliu/Documents/1Dtrig_2026Feb24/UQ/1Dtrig1e4_ratio6_inferP2/UQ/")
setwd("/Users/linxiliu/Documents/1Dcircle_2026Feb24/UQ/1Dcircle1e4_ratio6_inferP2/UQ/")

total_num = 1024*4
x <- seq.int(from = 1, to= total_num, by=1)/ total_num + 1/(2*total_num) 

# mixture-of-beta distribution
# mixture weights and component densities
w1 <- 0.7
w2 <- 0.3

d1 <- dbeta(x, shape1 = 3,   shape2 = 12)
d2 <- dbeta(x, shape1 = 300, shape2 = 10)
dmix <- w1 * d1 + w2 * d2

density_df <- data.frame(
  x = x,
  y0 = dmix
)

# trigonometric distribution
dtrig <- 1+ sin(2*pi * x - pi/2)
density_df <- data.frame(
  x = x,
  y0 = dtrig
)

# semi-circle distribution
dcircle <- ifelse( (x > 0.5 & x <1),  (8.0/pi) * sqrt( 1 - (4*x -3) *(4*x-3)),0)
density_df <- data.frame(
  x = x,
  y0 = dcircle
)

### Read in data
K=200

for(i in 1:K) {
  filename <- paste0("./UQforest_", i-1, "/forest_density.txt")
  varname <- paste0("y", i)
  density_data <- read.table(filename)
  density_df[[varname]] = density_data[,2]
}

density_df$ymean <- apply(density_df[,3:(K+2)], 1, mean)


library(tidyr)

df_long <- density_df |>
  pivot_longer(
    cols = starts_with("y"),
    names_to = "curve",
    values_to = "value"
  )

df_long$col_index <- ifelse(df_long$curve == "y0" | df_long$curve == "ymean", df_long$curve, "samples")

df_long$col_index [df_long$curve == "y0"] <- "true"
df_long$col_index [df_long$curve == "ymean"] <- "posterior_mean"

# For the beta distribution example
density_plot <- ggplot(df_long, aes(x, value, color = col_index)) +
  scale_color_manual(
    values = c(
      "true" = "violetred3",
      "samples" = "gray70",
      "posterior_mean" = "royalblue"
    ),
    labels = c(
      "true" = "True density",
      "samples" = "Posterior samples",
      "posterior_mean" = "Posterior mean"
    ))+
  geom_line(
    data = subset(df_long, col_index == "samples")
  ) +
  geom_line(
    data = subset(df_long, col_index == "posterior_mean")
  ) +
  geom_line(
    data = subset(df_long, col_index == "true")
  )+
  coord_cartesian(ylim = c(0, 12.5)) + #, xlim = c(0.75, 1)
  #coord_cartesian(xlim = c(0.75, 1)) +
  labs(
    x = "x",
    y = "Density",
    color =" ",
    title = "UQ of GDF for the mixture-of-beta distribution"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(size = 18),      # title size
    legend.title = element_text(size = 14),    # legend title size
    legend.text = element_text(size = 14),     # legend labels size
    axis.title = element_text(size = 14),      # x/y axis title
    axis.text = element_text(size = 12)        # axis tick labels
  )

density_plot
  
ggsave(
  filename = "UQ_plot.pdf",
  plot = density_plot,
  width = 10,
  height = 6
)
  
  
# For the trigonometric distribution example
density_plot <- ggplot(df_long, aes(x, value, color = col_index)) +
  scale_color_manual(
    values = c(
      "true" = "violetred3",
      "samples" = "gray70",
      "posterior_mean" = "royalblue"
    ),
    labels = c(
      "true" = "True density",
      "samples" = "Posterior samples",
      "posterior_mean" = "Posterior mean"
    ))+
  geom_line(
    data = subset(df_long, col_index == "samples")
  ) +
  geom_line(
    data = subset(df_long, col_index == "posterior_mean")
  ) +
  geom_line(
    data = subset(df_long, col_index == "true")
  )+
  coord_cartesian(ylim = c(0, 2.5)) + #, xlim = c(0.75, 1)
  #coord_cartesian(xlim = c(0.75, 1)) +
  labs(
    x = "x",
    y = "Density",
    color =" ",
    title = "UQ of GDF for the trigonometric distribution"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(size = 18),      # title size
    legend.title = element_text(size = 14),    # legend title size
    legend.text = element_text(size = 14),     # legend labels size
    axis.title = element_text(size = 14),      # x/y axis title
    axis.text = element_text(size = 12)        # axis tick labels
  )

density_plot

ggsave(
  filename = "UQ_plot.pdf",
  plot = density_plot,
  width = 10,
  height = 6
)

# For the circular distribution example
density_plot <- ggplot(df_long, aes(x, value, color = col_index)) +
  scale_color_manual(
    values = c(
      "true" = "violetred3",
      "samples" = "gray70",
      "posterior_mean" = "royalblue"
    ),
    labels = c(
      "true" = "True density",
      "samples" = "Posterior samples",
      "posterior_mean" = "Posterior mean"
    ))+
  geom_line(
    data = subset(df_long, col_index == "samples")
  ) +
  geom_line(
    data = subset(df_long, col_index == "posterior_mean")
  ) +
  geom_line(
    data = subset(df_long, col_index == "true")
  )+
  coord_cartesian(ylim = c(0, 4)) + 
  labs(
    x = "x",
    y = "Density",
    color =" ",
    title = "UQ of GDF for the semi-circle distribution"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(size = 18),      # title size
    legend.title = element_text(size = 14),    # legend title size
    legend.text = element_text(size = 14),     # legend labels size
    axis.title = element_text(size = 14),      # x/y axis title
    axis.text = element_text(size = 12)        # axis tick labels
  )

density_plot

ggsave(
  filename = "UQ_plot.pdf",
  plot = density_plot,
  width = 10,
  height = 6
)
