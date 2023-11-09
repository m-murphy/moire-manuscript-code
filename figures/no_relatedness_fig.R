library(ggplot2)
library(patchwork)

source("plotting_defaults.R")
source("load_data.R")


plot_allele_frequencies <- function(df) {
    ggplot(
        df,
        aes(x = true_allele_frequency, y = value, group = Estimator, color = Estimator)
    ) +
        geom_point(alpha = .5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        theme(axis.title = element_blank()) +
        facet_wrap(~panel, nrow = 1) +
        labs(color = "Allele Frequency") +
        coord_fixed() +
        theme_and_axis_legend +
        rotate_axis_text
}

plot_he <- function(df) {
    ggplot(
        df,
        aes(x = true_he, y = value, group = Estimator, color = Estimator)
    ) +
        geom_point(alpha = .5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        theme(axis.title = element_blank()) +
        facet_wrap(~panel, nrow = 1) +
        labs(color = "Heterozygosity") +
        coord_fixed() +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}

plot_moi <- function(df) {
    ggplot(
        df,
        aes(x = true_coi, y = value, group = Estimator, color = Estimator)
    ) +
        geom_boxplot(aes(group = interaction(true_coi, Estimator)), alpha = .5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        theme(axis.title = element_blank()) +
        facet_wrap(~panel, nrow = 1) +
        labs(color = "MOI") +
        coord_fixed() +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}

plot_relatedness <- function(df) {
    ggplot(
        df,
        aes(y = value, group = Estimator, color = Estimator)
    ) +
        geom_boxplot() +
        theme(axis.title = element_blank()) +
        facet_wrap(~panel, nrow = 1) +
        labs(color = "Relatedness") +
        coord_fixed() +
        remove_facet_label +
        theme_and_axis_legend
}

plot_emoi <- function(df) {
    ggplot(
        df,
        aes(x = true_effective_coi, y = value, group = Estimator, color = Estimator)
    ) +
        geom_point(alpha = .5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        theme(axis.title = element_blank()) +
        facet_wrap(~panel, nrow = 1) +
        labs(color = "eMOI") +
        coord_fixed() +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}


emoi_to_plot <- ecoi_results_plot |>
    dplyr::filter(
        expected_fp == 1 | is.na(expected_fp),
        expected_fn == 10 | is.na(expected_fn),
        !(panel %in% c("Mad4hatter-pool3", "Mad4hatter-pool4")),
        relatedness == "0",
    ) |>
    dplyr::mutate(
        panel = case_when(
            panel == "Mad4hatter-full" ~ "Mad4hatter",
            TRUE ~ panel
        ),
        panel = factor(
            panel,
            levels = c("100 SNPs", "Low Div.", "Med. Div.", "High Div.", "Mad4hatter", "Broad", "Sanger")
        )
    )


af_to_plot <- allele_results_plot |>
    dplyr::filter(
        expected_fp == 1 | is.na(expected_fp),
        expected_fn == 10 | is.na(expected_fn),
        !(panel %in% c("Mad4hatter-pool3", "Mad4hatter-pool4")),
        relatedness == "0",
    ) |>
    dplyr::mutate(
        panel = case_when(
            panel == "Mad4hatter-full" ~ "Mad4hatter",
            TRUE ~ panel
        ),
        panel = factor(
            panel,
            levels = c("100 SNPs", "Low Div.", "Med. Div.", "High Div.", "Mad4hatter", "Broad", "Sanger")
        )
    )


he_to_plot <- he_results_plot |>
    dplyr::filter(
        expected_fp == 1 | is.na(expected_fp),
        expected_fn == 10 | is.na(expected_fn),
        # mean_coi == 3 | is.na(mean_coi),
        !(panel %in% c("Mad4hatter-pool3", "Mad4hatter-pool4")),
        relatedness == "0",
        # relatedness == t_relatedness,
        # region == "East-Africa" | is.na(region)
    ) |>
    dplyr::mutate(
        panel = case_when(
            panel == "Mad4hatter-full" ~ "Mad4hatter",
            TRUE ~ panel
        ),
        panel = factor(
            panel,
            levels = c("100 SNPs", "Low Div.", "Med. Div.", "High Div.", "Mad4hatter", "Broad", "Sanger")
        )
    )


moi_to_plot <- coi_results_plot |>
    dplyr::filter(
        expected_fp == 1 | is.na(expected_fp),
        expected_fn == 10 | is.na(expected_fn),
        # mean_coi == 3 | is.na(mean_coi),
        relatedness == "0",
        # relatedness == t_relatedness,
        # region == "East-Africa" | is.na(region),
        !(panel %in% c("Mad4hatter-pool3", "Mad4hatter-pool4")),
        Estimator %in% c("Posterior Median", "Naive Offset")
    ) |>
    dplyr::mutate(
        panel = case_when(
            panel == "Mad4hatter-full" ~ "Mad4hatter",
            TRUE ~ panel
        ),
        panel = factor(
            panel,
            levels = c("100 SNPs", "Low Div.", "Med. Div.", "High Div.", "Mad4hatter", "Broad", "Sanger")
        )
    )



relatedness_to_plot <- relatedness_results_plot |>
    dplyr::filter(
        expected_fp == 1 | is.na(expected_fp),
        expected_fn == 10 | is.na(expected_fn),
        # mean_coi == 3 | is.na(mean_coi),
        relatedness == "0",
        true_effective_coi > 1,
        !(panel %in% c("Mad4hatter-pool3", "Mad4hatter-pool4")),
        # relatedness == t_relatedness,
        # region == "East-Africa" | is.na(region),
        post_coi_med > 1
    )
    # dplyr::mutate(
    #     panel = case_when(
    #         panel == "Mad4hatter-full" ~ "Mad4hatter",
    #         TRUE ~ panel
    #     ),
    #     panel = factor(
    #         panel,
    #         levels = c("100 SNPs", "Low Div.", "Med. Div.", "High Div.", "Mad4hatter", "Broad", "Sanger")
    #     )
    # )

design <- "
A
B
C
D
E
"


g <- plot_allele_frequencies(af_to_plot) +
    plot_he(he_to_plot) +
    plot_moi(moi_to_plot) +
    plot_relatedness(relatedness_to_plot) +
    plot_emoi(emoi_to_plot) +
    plot_layout(design = design)
g
ggsave("figures/no_relatedness_fig.png", g, width = 24, height = 20)
