library(ggplot2)
library(moire)
library(dplyr)

source("plotting_defaults.R")


eda <- function(afs, n) {
    sum(1 - (1 - afs)^n)
}

moire::calculate_he(c(.25, .75))

eda(c(.5, .5), 2)

afs_files <- list.files("afs", full.names = TRUE)
region_order <- c(
    "Central Africa",
    "Central West Africa",
    "West Africa",
    "East Africa",
    "South Africa",
    "South East Africa",
    "South Asia",
    "South East Asia - East",
    "South East Asia - West",
    "Papua New Guinea",
    "South America - Central",
    "South America - North"
)

panel_labels <- c("24*SNP", "101*SNP", "MaD^{4}*HatTeR", "AMPLseq", "AmpliSeq")

allele_freqs <- purrr::map(afs_files, function(x) {
    afs <- readRDS(x)
    panel <- strsplit(basename(x), "_")[[1]][1]
    purrr::imap(afs, function(af, region) {
        purrr::imap(unname(af), function(x, i) {
            dplyr::tibble(
                region = region,
                panel = panel,
                locus = paste0("L", i),
                allele_frequency = x
            )
        }) |> dplyr::bind_rows()
    }) |> dplyr::bind_rows()
}) |>
    dplyr::bind_rows() |>
    dplyr::mutate(region = factor(region, levels = region_order)) |>
    dplyr::mutate(panel = dplyr::case_when(
        panel == "madhatter-full" ~ "MaD^{4}*HatTeR",
        panel == "amplseq" ~ "AMPLseq",
        panel == "ampliseq" ~ "AmpliSeq",
        panel == "broad" ~ "24*SNP",
        panel == "sanger" ~ "101*SNP"
    )) |>
    dplyr::mutate(
        panel = factor(panel, levels = panel_labels)
    )



allele_freq_summary <- allele_freqs |>
    dplyr::group_by(region, panel, locus) |>
    dplyr::summarise(
        he = moire::calculate_he(allele_frequency),
        eda2 = eda(allele_frequency, 2),
        eda3 = eda(allele_frequency, 3),
        eda4 = eda(allele_frequency, 4),
        eda5 = eda(allele_frequency, 5),
        eda6 = eda(allele_frequency, 6),
        eda7 = eda(allele_frequency, 7),
        eda8 = eda(allele_frequency, 8),
        eda9 = eda(allele_frequency, 9),
        eda10 = eda(allele_frequency, 10)
    )

eda_summary <- allele_freq_summary |>
    tidyr::pivot_longer(
        cols = tidyselect::starts_with("eda"),
        names_to = "strains",
        values_to = "eda_value"
    ) |>
    dplyr::mutate(strains = as.numeric(stringr::str_replace(strains, "eda", ""))) |>
    dplyr::group_by(region, panel, strains) |>
    dplyr::summarize(mean_eda = mean(eda_value), max_eda = max(eda_value), q75_eda = quantile(eda_value, 0.75), q25_eda = quantile(eda_value, 0.25)) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(
        cols = c(mean_eda, q75_eda),
        names_to = "eda_type",
        values_to = "eda_value"
    ) |>
    dplyr::mutate(eda_type = dplyr::case_when(
        eda_type == "mean_eda" ~ "Mean EDA",
        eda_type == "q75_eda" ~ "75th Percentile EDA"
    )) |>
    dplyr::mutate(eda_type = factor(eda_type, levels = c("Mean EDA", "75th Percentile EDA"))) |>
    dplyr::mutate(
        continent = dplyr::case_when(
            region %in% c("Central Africa", "Central West Africa", "West Africa", "East Africa", "South Africa", "South East Africa") ~ "Africa",
            region %in% c("South Asia", "South East Asia - East", "South East Asia - West", "Papua New Guinea") ~ "Asia/Oceania",
            region %in% c("South America - Central", "South America - North") ~ "Americas"
        )
    )


africa_eda_summary <- eda_summary |>
    dplyr::filter(region %in% c("Central Africa", "Central West Africa", "West Africa", "East Africa", "South Africa", "South East Africa"))

asia_eda_summary <- eda_summary |>
    dplyr::filter(region %in% c("South Asia", "South East Asia - East", "South East Asia - West", "Papua New Guinea"))

americas_eda_summary <- eda_summary |>
    dplyr::filter(region %in% c("South America - Central", "South America - North"))


g_afr <- ggplot(
    africa_eda_summary,
    aes(
        x = strains, y = eda_value,
        group = interaction(panel, eda_type),
        shape = eda_type, fill = panel,
        color = panel, linetype = eda_type
    )
) +
    geom_line() +
    geom_point() +
    facet_wrap(~region, nrow = 1, labeller = label_wrap_gen(25)) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = .5),
    ) +
    labs(
        x = "Number of Strains",
        y = "Mean EDA",
        color = "Panel",
        linetype = "EDA"
    ) +
    theme_and_axis_legend +
    scale_color_discrete(labels = scales::parse_format()) +
    scale_shape_manual(
        values = c(21, 1),
        labels = c("Mean EDA", "75th Percentile EDA")
    ) +
    lims(y = c(1, 5)) +
    guides(
        shape = "none",
        fill = "none"
    )

g_asia <- ggplot(
    asia_eda_summary,
    aes(
        x = strains, y = eda_value,
        group = interaction(panel, eda_type),
        shape = eda_type, fill = panel,
        color = panel, linetype = eda_type
    )
) +
    geom_line() +
    geom_point() +
    facet_wrap(~region, nrow = 1, labeller = label_wrap_gen(25)) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = .5),
    ) +
    labs(
        x = "Number of Strains",
        y = "Mean EDA",
        color = "Panel",
        linetype = "EDA"
    ) +
    theme_and_axis_legend +
    scale_color_discrete(labels = scales::parse_format()) +
    scale_shape_manual(
        values = c(21, 1),
        labels = c("Mean EDA", "75th Percentile EDA")
    ) +
    lims(y = c(1, 5)) +
    guides(
        shape = "none",
        fill = "none"
    ) +
    theme(
        axis.title.x = element_blank()
    )

g_americas <- ggplot(
    americas_eda_summary,
    aes(
        x = strains, y = eda_value,
        group = interaction(panel, eda_type),
        shape = eda_type, fill = panel,
        color = panel, linetype = eda_type
    )
) +
    geom_line() +
    geom_point() +
    facet_wrap(~region, nrow = 1, labeller = label_wrap_gen(25)) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = .5),
    ) +
    labs(
        x = "Number of Strains",
        y = "Mean EDA",
        color = "Panel",
        linetype = "EDA"
    ) +
    theme_and_axis_legend +
    scale_color_discrete(labels = scales::parse_format()) +
    scale_shape_manual(
        values = c(21, 1),
        labels = c("Mean EDA", "75th Percentile EDA")
    ) +
    lims(y = c(1, 5)) +
    guides(
        shape = "none",
        fill = "none"
    ) +
    theme(
        axis.title.x = element_blank()
    )

layout <- "
AA####
BBBB##
CCCCCC
"

g <- g_americas / g_asia / g_afr + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "A")


ggsave("figures/regional_eda.pdf", g, width = 25, height = 20)
