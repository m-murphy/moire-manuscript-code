library(ggplot2)
library(patchwork)

setwd("/home/mmurphy/Workspace/moire-manuscript/R/0moire_v20230328")
source("plotting_defaults.R")
source("load_data.R")

# performance of population parameter estimation

## mean moi estimation
mean_moi_plot <- function(df) {
    ggplot(df, aes(x = true_mean_coi, y = mean_mean_coi, ymin = lower_mean_coi, ymax = upper_mean_coi, color = mean_coi_val)) +
        geom_point() +
        geom_errorbar() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        labs(x = "True Mean MOI", y = "Estimate", color = "Mean MOI") +
        lims(
            x = c(0, max(c(df$true_mean_coi, df$upper_mean_coi), na.rm = TRUE)),
            y = c(0, max(c(df$mean_mean_coi, df$upper_mean_coi), na.rm = TRUE))
        ) +
        theme_and_axis_legend
}

## mean eMOI estimation
mean_emoi_plot <- function(df) {
    ggplot(df, aes(x = true_mean_ecoi, y = mean_mean_ecoi, ymin = lower_mean_ecoi, ymax = upper_mean_ecoi, color = mean_coi_val)) +
        geom_point() +
        geom_errorbar() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        labs(x = "True Mean eMOI", y = "Estimate", color = "Mean MOI") +
        remove_facet_label +
        theme_and_axis_legend
}

## mean within-host relatedness estimation
mean_relatedness_plot <- function(df) {
    ggplot(df, aes(x = true_mean_rel_ignore_monoclonal, y = mean_mean_rel_ignore_monoclonal, ymin = lower_mean_rel_ignore_monoclonal, ymax = upper_mean_rel_ignore_monoclonal, color = mean_coi_val)) +
        geom_point() +
        geom_errorbar() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        labs(x = "True Mean Relatedness", y = "Estimate", color = "Mean MOI") +
        lims(x = c(0, 1), y = c(0, 1)) +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}

prop_polyclonal_plot <- function(df) {
    ggplot(df, aes(x = true_prop_polyclonal, y = mean_prop_polyclonal, ymin = lower_prop_polyclonal, ymax = upper_prop_polyclonal, color = mean_coi_val)) +
        geom_point() +
        geom_errorbar() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        labs(x = "True Proportion Polyclonal", y = "Estimate", color = "Mean MOI") +
        lims(x = c(0, 1), y = c(0, 1)) +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}

df <- pop_results_df |>
    dplyr::filter(
        prior == "None",
        !fws_sim,
        expected_fp == 1 | is.na(expected_fp),
        expected_fn == 10 | is.na(expected_fn)
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
    )

moi_plot <- mean_moi_plot(df)
emoi_plot <- mean_emoi_plot(df)
rel_plot <- mean_relatedness_plot(df)
# prop_polyclonal_plot <- prop_polyclonal_plot(df)

g <- moi_plot / rel_plot / emoi_plot + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave("figures/population_fig.pdf", g, width = 20, height = 12)
