library(ggplot2)
library(patchwork)

source("load_data.R")
source("plotting_defaults.R")

panel_order <- c("100 SNPs", "Low Div.", "Med. Div.", "High Div.", "Mad4hatter", "Broad", "Sanger")
priors <- pop_results_df$prior |> unique()

for (p in priors) {
    to_plot <- pop_results_df |>
        dplyr::filter(
            expected_fp == 1 | is.na(expected_fp),
            expected_fn == 10 | is.na(expected_fn),
            prior == p
        ) |>
        dplyr::mutate(
            mean_coi = case_when(
                is.na(mean_coi) ~ "3",
                TRUE ~ mean_coi
            ),
            mean_coi_val = case_when(
                mean_coi == "1" ~ "1.58",
                mean_coi == "3" ~ "3.16",
                mean_coi == "5" ~ "5.03",
            )
        ) |>
        dplyr::select(
            panel, relatedness_str, true_mean_relatedness, 
            inv_weighted_mean_relatedness, mean_mean_relatedness, 
            true_mean_rel_ignore_monoclonal, mean_mean_rel_ignore_monoclonal, 
            upper_mean_rel_ignore_monoclonal, lower_mean_rel_ignore_monoclonal,
            mean_mean_rel_polyclonal_only, upper_mean_rel_polyclonal_only, 
            lower_mean_rel_polyclonal_only,
            prior, mean_coi_val, expected_fp, expected_fn
        )

    g <- to_plot |>
        ggplot(aes(x = true_mean_rel_ignore_monoclonal, y = mean_mean_rel_ignore_monoclonal, color =)) +
        geom_point() +
        geom_errorbar(aes(ymin = lower_mean_rel_ignore_monoclonal, ymax = upper_mean_rel_ignore_monoclonal)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        facet_grid(mean_coi_val ~ panel) +
        coord_fixed() +
        scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 1)) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1)
        ) +
        labs(
            x = "True mean relatedness",
            y = "Post. mean relatedness"
        ) +
        theme_and_axis_legend

    path <- paste0("figures/mean_relatedness_figs/mean_relatedness_", p, ".png")
    ggsave(path, g, width = 15, height = 15)
}

