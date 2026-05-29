## ============================================================
## FINAL PhyloBayes chronogram with correctly mapped node-age HPD bars
## Uses *_sample.labels as labelled tree
## 0 Ma on RIGHT + geological timescale
## Fixes:
##   1. correct node mapping
##   2. bottom tip not hidden
##   3. tip labels visible
##   4. geological timescale kept clean
## ============================================================

## ---- Packages ----
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggtree", quietly = TRUE)) BiocManager::install("ggtree")
if (!requireNamespace("deeptime", quietly = TRUE)) install.packages("deeptime")

library(ape)
library(ggplot2)
library(dplyr)
library(ggtree)
library(deeptime)
library(grid)

## ---- Directory ----
base_dir <- "C:\Users\33651\Documents\Astage M1\methode_final\tree_Nalda_5N\arbre final\PhyloBayes chrono_both2"
setwd(base_dir)

## ---- Find files ----
pick_one <- function(pattern) {
  x <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  if (length(x) == 0) stop("No file matching: ", pattern)
  if (length(x) > 1) message("More than one file matched ", pattern, ". Using: ", basename(x[1]))
  x[1]
}

dates_file  <- pick_one("_sample\\.dates$")
labels_file <- pick_one("_sample\\.labels$")

message("Dates file:  ", basename(dates_file))
message("Labels tree: ", basename(labels_file))

## ---- Read .labels file as tree ----
tr <- read.tree(labels_file)

if (inherits(tr, "multiPhylo")) {
  message("Multiple trees detected in labels file. Using the first one.")
  tr <- tr[[1]]
}

message("Number of tips: ", Ntip(tr))
message("Number of internal nodes: ", tr$Nnode)

## ---- Read .dates summary table ----
dates <- read.table(
  dates_file,
  header = TRUE,
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

colnames(dates)[1:9] <- c(
  "pb_node",
  "meandate",
  "date_stderr",
  "inf95",
  "sup95",
  "instant_rate",
  "instant_rate_stderr",
  "average_rate",
  "average_rate_stderr"
)

dates_internal <- dates %>%
  mutate(
    pb_node = suppressWarnings(as.integer(pb_node)),
    meandate = suppressWarnings(as.numeric(meandate)),
    date_stderr = suppressWarnings(as.numeric(date_stderr)),
    inf95 = suppressWarnings(as.numeric(inf95)),
    sup95 = suppressWarnings(as.numeric(sup95))
  ) %>%
  filter(!is.na(pb_node), !is.na(meandate), meandate > 0)

## ---- Get ggtree coordinates ----
p0 <- ggtree(tr, size = 0.45)
tree_df <- p0$data

## Convert x to age before present
root_age <- max(tree_df$x[tree_df$isTip], na.rm = TRUE)

tree_df <- tree_df %>%
  mutate(age = root_age - x)

tip_df <- tree_df %>%
  filter(isTip) %>%
  mutate(label_x = -4)   ## extra position just to the right of 0 Ma

## ---- Extract PhyloBayes internal node labels ----
node_xy <- tree_df %>%
  filter(!isTip) %>%
  mutate(
    pb_node = suppressWarnings(as.integer(as.character(label)))
  ) %>%
  filter(!is.na(pb_node)) %>%
  select(
    pb_node,
    ape_node = node,
    y,
    tree_age = age
  )

message("Number of labelled internal nodes found: ", nrow(node_xy))

if (nrow(node_xy) == 0) {
  stop("No internal PhyloBayes node labels were found after reading the .labels tree.")
}

## ---- Join .dates to correct labelled tree nodes ----
ci_df <- dates_internal %>%
  left_join(node_xy, by = "pb_node") %>%
  filter(!is.na(y), !is.na(tree_age))

missing_nodes <- setdiff(dates_internal$pb_node, ci_df$pb_node)

if (length(missing_nodes) > 0) {
  warning(
    "Some .dates nodes were not found in the labelled tree: ",
    paste(missing_nodes, collapse = ", ")
  )
}

if (nrow(ci_df) == 0) {
  stop("No .dates nodes could be mapped to the labelled tree.")
}

## ---- Check mapping quality ----
ci_df <- ci_df %>%
  mutate(age_minus_meandate = tree_age - meandate)

message("Median |tree age - .dates meandate| = ",
        round(median(abs(ci_df$age_minus_meandate), na.rm = TRUE), 4), " Ma")

message("Maximum |tree age - .dates meandate| = ",
        round(max(abs(ci_df$age_minus_meandate), na.rm = TRUE), 4), " Ma")

write.csv(
  ci_df,
  file = file.path(base_dir, "chrono_both2_node_mapping_FROM_LABEL_TREE.csv"),
  row.names = FALSE
)

## ---- Center HPD bars on correct plotted nodes ----
ci_df <- ci_df %>%
  mutate(
    lower_width = meandate - inf95,
    upper_width = sup95 - meandate,
    bar_center = tree_age,
    bar_xmin = pmax(0, bar_center - lower_width),
    bar_xmax = bar_center + upper_width
  )

## ---- Axis range ----
max_age_data <- max(c(root_age, ci_df$bar_xmax, ci_df$sup95), na.rm = TRUE)
max_age <- ceiling(max_age_data / 50) * 50

## Extra right-hand space for tip labels.
## Because the x-axis is reversed, values below 0 are to the right of present.
right_label_space <- 35

x_breaks <- seq(0, max_age, by = 50)

## ---- Extra vertical space to prevent WSSV from being hidden ----
y_min <- min(tree_df$y, na.rm = TRUE) - 1.4
y_max <- max(tree_df$y, na.rm = TRUE) + 1.0

## ---- Alternating vertical background bands ----
band_width <- 50

band_df <- data.frame(
  xmin = seq(0, max_age - band_width, by = band_width),
  xmax = seq(band_width, max_age, by = band_width)
)

band_df <- band_df[seq_len(nrow(band_df)) %% 2 == 0, ]

## ============================================================
## Final plot
## ============================================================

p <- ggplot() +
  
  ## Alternating vertical time bands
  geom_rect(
    data = band_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey95",
    colour = NA
  ) +
  
  ## Tree
  geom_tree(
    data = tree_df,
    aes(x = age, y = y),
    colour = "grey10",
    linewidth = 0.45
  ) +
  
  ## Correctly mapped node-age HPD bars
  geom_segment(
    data = ci_df,
    aes(x = bar_xmin, xend = bar_xmax, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 1.6,
    colour = "grey45",
    alpha = 0.95,
    lineend = "round"
  ) +
  
  ## Node centers
  geom_point(
    data = ci_df,
    aes(x = bar_center, y = y),
    inherit.aes = FALSE,
    size = 0.8,
    colour = "black"
  ) +
  
  ## Tip labels, placed in the extra right-hand space
  geom_text(
    data = tip_df,
    aes(x = label_x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.1,
    colour = "grey10"
  ) +
  
  ## Reverse axis: older left, 0 Ma right
  ## Limit extends slightly past 0 to make room for labels.
  scale_x_reverse(
    limits = c(max_age, -right_label_space),
    breaks = x_breaks,
    labels = x_breaks,
    expand = expansion(mult = c(0.01, 0.00)),
    name = "Age before present (Ma)",
    sec.axis = dup_axis(
      name = NULL,
      breaks = x_breaks,
      labels = x_breaks
    )
  ) +
  
  ## Geological timescale
  coord_geo(
    xlim = c(max_age, -right_label_space),
    ylim = c(y_min, y_max),
    dat = list("periods", "eras"),
    pos = list("bottom", "bottom"),
    abbrv = list(TRUE, FALSE),
    rot = 0,
    neg = FALSE,
    clip = "on",
    height = unit(c(0.65, 0.55), "cm"),
    center_end_labels = TRUE
  ) +
  
  labs(
    title = "PhyloBayes chronogram with node-age uncertainty",
    y = NULL
  ) +
  
  theme_tree2() +
  
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0,
      margin = margin(b = 6)
    ),
    axis.title.x = element_text(
      size = 12,
      face = "bold",
      margin = margin(t = 22)
    ),
    axis.text.x = element_text(
      size = 9,
      colour = "grey15"
    ),
    axis.text.x.top = element_text(
      size = 9,
      colour = "grey15"
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid.major.x = element_line(
      colour = "grey86",
      linewidth = 0.35
    ),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(12, 130, 105, 12)
  )

print(p)

## ---- Save final figure ----
ggsave(
  filename = file.path(base_dir, "chrono_both2_fancy_geological_timescale_FINAL_TIP_LABELS_FIXED.pdf"),
  plot = p,
  width = 12,
  height = 10
)

