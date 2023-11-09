library(ggplot2)
library(patchwork)
library(latex2exp)

setwd("/home/mmurphy/Workspace/moire-manuscript/R/0moire_v20230328")
source("plotting_defaults.R")
source("load_data.R")

mean_he <- he_results_df |>
    dplyr::filter(fws_sim) |>
    dplyr::group_by(simulation) |>
    dplyr::summarise(
        true_mean_he = mean(true_he),
        naive_mean_he = mean(naive_he)
    )


naive_fws_summary <- sample_results_df |>
    dplyr::filter(fws_sim, relatedness_str %in% c("Low", "High")) |>
    dplyr::group_by(simulation, region, panel, relatedness, relatedness_str, expected_fp, expected_fn, mean_coi) |>
    dplyr::summarize(
        mean_fws = mean(estimated_fws),
        mean_true_fws = mean(true_fws),
        mean_naive_coi = mean(offset_naive_coi),
        true_mean_coi = mean(true_coi),
        n = n()
    ) |>
    dplyr::left_join(mean_he) |>
    dplyr::arrange(relatedness, expected_fp, expected_fn) |>
    dplyr::mutate(latex_panel = dplyr::case_when(
        panel == "MaD^{4}*HatTeR" ~ "$\\text{MaD}^{4}\\text{HatTeR}$",
        panel == "24*SNP" ~ "$24\\text{SNP}$",
        panel == "101*SNP" ~ "$101\\text{SNP}$",
        panel == "Ampliseq" ~ "Ampliseq",
        panel == "Amplseq" ~ "Amplseq"
    )) |>
    dplyr::mutate(relatedness_str = dplyr::case_when(
        relatedness_str == "Low" ~ "Low relatedness",
        relatedness_str == "High" ~ "High relatedness"
    )) |>
    dplyr::mutate(
        relatedness_str = factor(relatedness_str, levels = c("Low relatedness", "High relatedness"))
    )


to_plot <- pop_results_df |>
    dplyr::filter(fws_sim, relatedness_str %in% c("Low", "High")) |>
    dplyr::left_join(mean_he) |>
    dplyr::arrange(relatedness, expected_fp, expected_fn) |>
    dplyr::mutate(latex_panel = dplyr::case_when(
        panel == "MaD^{4}*HatTeR" ~ "$\\text{MaD}^{4}\\text{HatTeR}$",
        panel == "24*SNP" ~ "$24\\text{SNP}$",
        panel == "101*SNP" ~ "$101\\text{SNP}$",
        panel == "Ampliseq" ~ "Ampliseq",
        panel == "Amplseq" ~ "Amplseq"
    )) |>
    dplyr::mutate(relatedness_str = dplyr::case_when(
        relatedness_str == "Low" ~ "Low relatedness",
        relatedness_str == "High" ~ "High relatedness"
    )) |>
    dplyr::mutate(
        relatedness_str = factor(relatedness_str, levels = c("Low relatedness", "High relatedness"))
    )


fws_plt <- ggplot(naive_fws_summary, aes(
    x = true_mean_he,
    y = mean_fws,
    color = relatedness_str,
    linetype = as.character(true_mean_coi),
    group = interaction(relatedness_str, mean_coi)
)) +
    geom_point(
        size = 1
    ) +
    geom_smooth(method = "lm") +
    facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
    labs(x = "True mean heterozygosity", y = TeX("\\textbf{Mean $F_{WS}$}"), color = "Relatedness", linetype = "Mean MOI") +
    guides(color = "none") +
    theme_minimal() +
    remove_x_axis

nmoi_plt <- ggplot(naive_fws_summary, aes(
    x = true_mean_he,
    y = mean_naive_coi,
    color = relatedness_str,
    linetype = as.character(true_mean_coi),
    group = interaction(relatedness_str, mean_coi)
)) +
    geom_point(
        size = 1
    ) +
    geom_smooth(method = "lm") +
    geom_hline(
        aes(
            yintercept = true_mean_coi
        ),
        color = "black",
        linetype = "dashed"
    ) +
    facet_wrap(~panel, nrow = 1) +
    labs(x = "True mean heterozygosity", y = "Mean naive MOI", color = "Relatedness", linetype = "Mean MOI") +
    guides(color = "none") +
    theme_minimal() +
    remove_x_axis +
    remove_facet_label

ecoi_plt <- ggplot(
    to_plot,
    aes(
        y = mean_mean_ecoi,
        x = true_mean_he,
        ymin = lower_mean_ecoi,
        ymax = upper_mean_ecoi,
        color = relatedness_str,
        linetype = as.character(true_mean_coi),
        group = interaction(relatedness_str, mean_coi)
    )
) +
    geom_point() +
    geom_errorbar(linetype = "solid") +
    geom_smooth(method = "lm") +
    geom_hline(aes(
        yintercept = true_mean_ecoi
    ), color = "black", linetype = "dashed") +
    facet_wrap(~panel, nrow = 1) +
    labs(x = "True mean heterozygosity", y = "Mean eMOI", color = "Relatedness", linetype = "Mean MOI") +
    rotate_axis_text +
    theme_minimal() +
    remove_facet_label

relatedness_plt <- ggplot(
    to_plot,
    aes(
        y = mean_mean_rel_ignore_monoclonal,
        x = true_mean_he,
        ymin = lower_mean_rel_ignore_monoclonal,
        ymax = upper_mean_rel_ignore_monoclonal,
        color = relatedness_str,
        linetype = as.character(true_mean_coi),
        group = interaction(relatedness_str, mean_coi)
    )
) +
    geom_point() +
    geom_errorbar(linetype = "solid") +
    geom_smooth(method = "lm") +
    geom_hline(aes(
        yintercept = true_mean_relatedness
    ), color = "black", linetype = "dashed") +
    facet_wrap(~panel, nrow = 1) +
    labs(x = "True mean heterozygosity", y = "Mean within-host relatedness", color = "Relatedness") +
    lims(y = c(0, 1)) +
    remove_facet_label +
    theme_minimal() +
    remove_x_axis

moi_plt <- ggplot(
    to_plot,
    aes(
        y = mean_mean_coi,
        x = true_mean_he,
        ymin = lower_mean_coi,
        ymax = upper_mean_coi,
        color = relatedness_str,
        linetype = as.character(round(true_mean_coi, 2)),
        group = interaction(relatedness_str, mean_coi)
    )
) +
    geom_point() +
    geom_errorbar(linetype = "solid") +
    geom_smooth(method = "lm") +
    geom_hline(aes(
        yintercept = true_mean_coi
    ), color = "black", linetype = "dashed") +
    facet_wrap(~panel, nrow = 1) +
    labs(x = "True mean heterozygosity", y = "Mean MOI", color = "Relatedness", linetype = "Mean MOI") +
    remove_facet_label +
    theme_minimal() +
    remove_x_axis


combo_plt <- fws_plt / nmoi_plt / ecoi_plt + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme_and_axis_legend & theme(panel.spacing.x = unit(15, "pt"))
ggsave("figures/diversity_metrics.pdf", combo_plt, height = 12, width = 16)
