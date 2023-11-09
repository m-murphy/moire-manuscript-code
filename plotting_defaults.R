library(ggplot2)
library(ggpubfigs)

plot_crop <- function(x, margins = c(0, 0, 0, 0)) {
  cmd <- sprintf(
    "--margins '%s %s %s %s' %s %s",
    margins[1], margins[2], margins[3], margins[4], x, x
  )
  # system2("/usr/local/texlive/2022/bin/x86_64-linux/pdfcrop", cmd, stdout = FALSE)
  system2("pdfcrop", cmd, stdout = FALSE)
}

convert_img <- function(x, from = ".pdf", to = ".eps") {
  out <- gsub(stringr::str_c(from, "$"), to, x)
  cmd <- sprintf("%s %s", x, out)
  system2("convert", cmd, stdout = FALSE)
}

save_and_crop <- function(fname, g, dir = NULL, width = NULL, height = NULL) {
  parent_dir <- figures_dir
  if (!is.null(dir)) {
    parent_dir <- file.path(figures_dir, dir)
  }

  dir.create(parent_dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(parent_dir, sprintf("%s.%s", fname, figure_fmt))
  ggsave(out,
    plot = g, dpi = "retina", device = figure_fmt,
    width = width, height = height, units = "in"
  )
  plot_crop(out)
}

plot_chain_swap_dist <- function(moire_results) {
  swap_dist <- moire_results$chains[[1]]$swap_acceptances / (moire_results$args$samples_per_chain)
  temps <- moire_results$chains[[1]]$temp_gradient
  swap_idx <- (temps[1:length(temps) - 1] + temps[2:length(temps)]) / 2
  dat <- data.frame(swap_rate = swap_dist, temp = swap_idx)
  g <- ggplot(dat, aes(x = temp, y = swap_rate)) +
    geom_point() +
    geom_vline(data = data.frame(x = temps), aes(xintercept = x), linetype = "dashed", alpha = 0.25) +
    coord_cartesian(ylim = c(0, 1))
  g
}


theme_set(theme_classic())

black <- "#000000"
grey <- "#808080"

palette <- friendly_pal("bright_seven")
palette2 <- friendly_pal("nickel_five")

# naive estimator color
naive_color <- palette[6]

theme_and_axis_nolegend <- theme(
  legend.position = "None",
  text = element_text(size = 25, face = "bold"),
  axis.text = element_text(size = 18, face = "bold", color = black),
  axis.line = element_line(color = black, linewidth = 0.6),
  panel.border = element_rect(fill = NA)
)

theme_no_axis_no_legend <- theme(
  legend.position = "None",
  text = element_text(size = 25, face = "bold"),
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  strip.background = element_blank()
)

theme_and_axis_legend <- theme(
  text = element_text(size = 25, face = "bold"),
  axis.text = element_text(size = 18, face = "bold", color = black),
  axis.line = element_line(color = black, linewidth = 0.6),
  panel.border = element_rect(fill = NA),
  strip.background = element_blank()
)

rotate_axis_text <- theme(
  axis.text.x = element_text(angle = 45, vjust = .75),
)

remove_facet_label <- theme(
  strip.background = element_blank(),
  strip.text = element_blank()
)

remove_x_facet_label <- theme(
  strip.background.x = element_blank(),
  strip.text.x = element_blank()
)

remove_y_facet_label <- theme(
  strip.background.y = element_blank(),
  strip.text.y = element_blank()
)

remove_y_axis <- theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank()
)

remove_x_axis <- theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank()
)

remove_border <- theme(
  panel.border = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(vjust = 3)
)

increase_legend_size <- guides(
  color = guide_legend(override.aes = list(size = 5)),
  shape = guide_legend(override.aes = list(size = 5))
)

figure_fmt <- "pdf"
figures_dir <- "figures_output/output"

plot_labeler <- c(
  "10" = "10 Alleles",
  "5" = "5 Alleles",
  "2" = "100 SNPs",
  "R: 0" = "r: 0",
  "R: 2" = "r: 2",
  "R: 4" = "r: 4",
  "R: 6" = "r: 6",
  "R: 8" = "r: 8",
  "R: 10" = "r: 10"
)

panel_order <- c("100*SNPs", "24*SNP", "101*SNP", "Low~Div.", "Med.~Div.", "High~Div.", "MaD^{4}*HatTeR", "Ampliseq", "Amplseq")

my_labeller <- function(values) {
  l1 <- list(
    "2" = "Synthetic",
    "5" = "Synthetic",
    "10" = "Synthetic",
    "20" = "Synthetic",
    "madhatter" = "MaD^{4}*HatTeR",
    "amplseq" = "AMPLseq",
    "ampliseq" = "AmpliSeq",
    "broad" = "24*SNP",
    "sanger" = "101*SNP"
  )

  l2 <- list(
    "2" = "Loci:~100",
    "5" = "Loci:~30",
    "10" = "Loci:~30",
    "20" = "Loci:~30",
    "madhatter" = "Loci:~165",
    "amplseq" = "Loci:~128",
    "ampliseq" = "Loci:~233",
    "broad" = "Loci:~24",
    "sanger" = "Loci:~101"
  )

  l3 <- list(
    "2" = "Alleles:~2",
    "5" = "Alleles:~5",
    "10" = "Alleles:~10",
    "20" = "Alleles:~20",
    "madhatter" = "Alleles:Mixed",
    "amplseq" = "Alleles:Mixed",
    "ampliseq" = "Alleles:Mixed",
    "broad" = "Alleles:~2",
    "sanger" = "Alleles:~2"
  )

  l1_entries <- lapply(unname(l1[as.character(values$panel)]), function(x) (parse(text = x)))
  l2_entries <- lapply(unname(l2[as.character(values$panel)]), function(x) (parse(text = x)))
  l3_entries <- lapply(unname(l3[as.character(values$panel)]), function(x) (parse(text = x)))

  list(
    l1_entries,
    l2_entries,
    l3_entries
  )
}
