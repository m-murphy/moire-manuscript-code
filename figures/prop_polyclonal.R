library(ggplot2)
library(patchwork)

source("load_data.R")
source("plotting_defaults.R")

to_plot <- pop_results_df |>
    dplyr::filter(
        expected_fp == 1 | is.na(expected_fp),
        expected_fn == 10 | is.na(expected_fn),
        prior == "Flat"
    ) |>
    dplyr::mutate(
        panel = case_when(
            panel == "Mad4hatter-full" ~ "Mad4hatter",
            TRUE ~ panel
        ),
        panel = factor(
            panel,
            levels = panel_order
        ),
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
        true_prop_polyclonal, mean_prop_polyclonal, lower_prop_polyclonal, upper_prop_polyclonal, panel, prior, mean_coi_val
    )

g <- to_plot |>
    ggplot(aes(x = true_prop_polyclonal, y = mean_prop_polyclonal, color = panel)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_prop_polyclonal, ymax = upper_prop_polyclonal)) +
    facet_grid(
        ~ panel
    ) +
    theme_bw() +
    xlab("True proportion polyclonal") +
    ylab("Post. prop. polyclonal") +
    coord_fixed() +
    theme_and_axis_legend +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

ggsave("figures/prop_polyclonal.pdf", g, width = 12, height = 12)
