library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

source("plotting_defaults.R")

parse_sample_results <- function(mcmc_res, metadata) {
    naive_coi <- moire::calculate_naive_coi_offset(mcmc_res$args$data$data, 2)
    coi <- moire::summarize_coi(mcmc_res)
    coi$naive_coi <- naive_coi
    ecoi <- moire::summarize_effective_coi(mcmc_res)
    relatedness <- moire::summarize_relatedness(mcmc_res)
    eps_pos <- moire::summarize_epsilon_pos(mcmc_res)
    eps_neg <- moire::summarize_epsilon_neg(mcmc_res)

    sample_metadata <- data.frame(sample_id = mcmc_res$args$data$sample_ids) |>
        dplyr::left_join(metadata, by = "sample_id")

    return(
        coi |>
            dplyr::left_join(ecoi, by = "sample_id") |>
            dplyr::left_join(relatedness, by = "sample_id") |>
            dplyr::left_join(eps_pos, by = "sample_id") |>
            dplyr::left_join(eps_neg, by = "sample_id") |>
            dplyr::left_join(sample_metadata, by = "sample_id")
    )
}

summarize_pop_params_by_facility <- function(mcmc_results, metadata) {
    facilities <- metadata |>
        dplyr::filter(sample_id %in% mcmc_results$args$data$sample_ids) |>
        dplyr::pull(HealthFacility) |>
        unique()

    pop_summary <- data.frame()

    for (f in facilities) {
        pop_summary <- rbind(
            pop_summary,
            data.frame(
                HealthFacility = f,
                summarize_pop_params(mcmc_results, subset = metadata |> dplyr::filter(HealthFacility == f) |> dplyr::pull(sample_id))
            )
        )
    }

    return(pop_summary)
}

summarize_pop_params <- function(mcmc_results, subset = NULL) {
    if (is.null(subset)) {
        sample_ids <- mcmc_results$args$data$sample_ids
        sample_idxs <- 1:length(sample_ids)
    } else {
        sample_ids <- subset
        sample_idxs <- which(mcmc_results$args$data$sample_ids %in% subset)
    }

    coi <- matrix(
        unlist(mcmc_results$chains[[1]]$coi[sample_idxs]),
        nrow = length(sample_idxs),
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )
    relatedness <- matrix(
        unlist(mcmc_results$chains[[1]]$relatedness[sample_idxs]),
        nrow = length(sample_idxs),
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )

    eps_pos <- matrix(
        unlist(mcmc_results$chains[[1]]$eps_pos[sample_idxs]),
        nrow = length(sample_idxs),
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )

    eps_neg <- matrix(
        unlist(mcmc_results$chains[[1]]$eps_neg[sample_idxs]),
        nrow = length(sample_idxs),
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )

    naive_cois <- moire::calculate_naive_coi_offset(mcmc_results$args$data$data, 2)[sample_idxs]
    naive_afs <- moire::calculate_naive_allele_frequencies(lapply(mcmc_results$args$data$data, function(x) x[sample_idxs]))
    naive_hes <- sapply(naive_afs, moire::calculate_he)

    post_allele_freqs <- lapply(1:length(mcmc_results$args$data$loci), function(x) c())

    for (l in 1:length(mcmc_results$chains[[1]]$allele_freqs)) {
        post_allele_freqs[[l]] <- c(post_allele_freqs[[l]], mcmc_results$chains[[1]]$allele_freqs[[l]])
    }

    post_allele_freqs_mats <- lapply(
        post_allele_freqs,
        function(x) matrix(unlist(x), ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin, byrow = FALSE)
    )

    he_dists <- lapply(post_allele_freqs_mats, function(x) {
        apply(x, 2, moire::calculate_he)
    })

    he_dists_combined <- matrix(unlist(he_dists), ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin, byrow = TRUE)
    mean_he_dist <- colMeans(he_dists_combined)


    rel_ignore_monoclonal <- relatedness
    rel_ignore_monoclonal[coi == 1] <- NA

    ecoi <- (coi - 1) * (1 - relatedness) + 1

    mean_coi_dist <- colMeans(coi)
    mean_ecoi_dist <- colMeans(ecoi)
    mean_relatedness_dist <- colMeans(relatedness)
    mean_rel_ignore_monoclonal_dist <- colMeans(rel_ignore_monoclonal, na.rm = TRUE)
    prop_polyclonal_dist <- colMeans(coi > 1)
    mean_eps_pos_dist <- colMeans(eps_pos)
    mean_eps_neg_dist <- colMeans(eps_neg)

    pop_summary <- data.frame(
        naive_mean_coi = mean(naive_cois),
        naive_prop_polyclonal = mean(naive_cois > 1),
        lower_mean_coi = quantile(mean_coi_dist, 0.025),
        upper_mean_coi = quantile(mean_coi_dist, 0.975),
        mean_mean_coi = mean(mean_coi_dist),
        med_mean_coi = median(mean_coi_dist),
        lower_mean_ecoi = quantile(mean_ecoi_dist, 0.025),
        upper_mean_ecoi = quantile(mean_ecoi_dist, 0.975),
        mean_mean_ecoi = mean(mean_ecoi_dist),
        med_mean_ecoi = median(mean_ecoi_dist),

        # lower_mean_relatedness = quantile(mean_relatedness_dist, 0.025),
        # upper_mean_relatedness = quantile(mean_relatedness_dist, 0.975),
        # mean_mean_relatedness = mean(mean_relatedness_dist),
        # med_mean_relatedness = median(mean_relatedness_dist),

        lower_mean_rel_ignore_monoclonal = quantile(mean_rel_ignore_monoclonal_dist, 0.025, na.rm = TRUE),
        upper_mean_rel_ignore_monoclonal = quantile(mean_rel_ignore_monoclonal_dist, 0.975, na.rm = TRUE),
        mean_mean_rel_ignore_monoclonal = mean(mean_rel_ignore_monoclonal_dist, na.rm = TRUE),
        med_mean_rel_ignore_monoclonal = median(mean_rel_ignore_monoclonal_dist, na.rm = TRUE),
        lower_mean_eps_pos = quantile(mean_eps_pos_dist, 0.025),
        upper_mean_eps_pos = quantile(mean_eps_pos_dist, 0.975),
        mean_mean_eps_pos = mean(mean_eps_pos_dist),
        med_mean_eps_pos = median(mean_eps_pos_dist),
        lower_mean_eps_neg = quantile(mean_eps_neg_dist, 0.025),
        upper_mean_eps_neg = quantile(mean_eps_neg_dist, 0.975),
        mean_mean_eps_neg = mean(mean_eps_neg_dist),
        med_mean_eps_neg = median(mean_eps_neg_dist),
        lower_prop_polyclonal = quantile(prop_polyclonal_dist, 0.025),
        upper_prop_polyclonal = quantile(prop_polyclonal_dist, 0.975),
        mean_prop_polyclonal = mean(prop_polyclonal_dist),
        med_prop_polyclonal = median(prop_polyclonal_dist),
        lower_mean_he = quantile(mean_he_dist, 0.025),
        upper_mean_he = quantile(mean_he_dist, 0.975),
        mean_mean_he = mean(mean_he_dist),
        med_mean_he = median(mean_he_dist),
        naive_mean_he = mean(naive_hes),
        n = length(sample_idxs)
    )

    return(pop_summary)
}


summarize_mean_relatedness <- function(mcmc_results) {
    relatedness <- unlist(mcmc_results$chains[[1]]$relatedness) * as.integer(unlist(mcmc_results$chains[[1]]$coi) > 1)
    relatedness_mcmc_chain <- matrix(
        relatedness,
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )
    mean_relatedness_dist <- colMeans(relatedness_mcmc_chain)
    return(mean_relatedness_dist)
}

hd_order <- c(
    "Rundu",
    "Nyangana",
    "Andara",
    "Zambezi"
)

hf_order <- c(
    "Sauyemwa C", "Ndama C", "Rundu DH", "Nkarapamwe C", "Kaisosi C",
    "Kayengona C", "Sambyu HC", "Mashare C", "Mabushe", "Ndonga Linenea",
    "Karukuta C", "Nyangana DH", "Katere C", "Shinyungwe C", "Mbambi East",
    "Kangongo C", "Mayara C", "Biro C", "Shadikongoro C", "Andara DH",
    "Divundu C", "Mutjiku C", "Old Bagani C", "Choi C", "Sesheke C",
    "Sibbinda HC", "Kasheshe C"
)
excluded_hfs <- c("Chinchimani_C.rds", "Kanono_C.rds")

basedir <- "namibia_20230726"
mcmc_dir <- file.path(basedir, "mcmc_output")
hf_files <- list.files(path = mcmc_dir, pattern = "*.rds")
hf_files <- hf_files[!hf_files %in% excluded_hfs]
metadata <- readxl::read_xlsx(file.path(basedir, "namibia_data.xlsx"), skip = 1) |>
    dplyr::rename(sample_id = ID) |>
    dplyr::select(1:4)

all_sample_res <- do.call(rbind, lapply(hf_files, function(x) {
    mcmc_results <- readRDS(file.path(mcmc_dir, x))
    data.frame(
        parse_sample_results(mcmc_results, metadata)
    )
})) |>
    dplyr::mutate(
        HealthFacility = factor(HealthFacility, levels = hf_order),
        HealthDistrict = factor(HealthDistrict, levels = hd_order)
    )



he_res <- do.call(rbind, lapply(hf_files, function(x) {
    mcmc_results <- readRDS(file.path(mcmc_dir, x))
    HealthFacility <- strsplit(x, "\\.")[[1]][1] |> stringr::str_replace_all("_", " ")
    HealthDistrict <- unique(metadata$HealthDistrict[metadata$HealthFacility == HealthFacility])[1]
    data.frame(
        HealthFacility = HealthFacility,
        HealthDistrict = HealthDistrict,
        moire::summarize_he(mcmc_results)
    )
})) |>
    dplyr::mutate(
        HealthFacility = factor(HealthFacility, levels = hf_order),
        HealthDistrict = factor(HealthDistrict, levels = hd_order)
    )

pop_res <- do.call(rbind, lapply(hf_files, function(x) {
    mcmc_results <- readRDS(file.path(mcmc_dir, x))
    HealthFacility <- strsplit(x, "\\.")[[1]][1] |> stringr::str_replace_all("_", " ")
    HealthDistrict <- unique(metadata$HealthDistrict[metadata$HealthFacility == HealthFacility])[1]
    data.frame(
        HealthFacility = HealthFacility,
        HealthDistrict = HealthDistrict,
        summarize_pop_params(mcmc_results)
    )
})) |>
    dplyr::mutate(
        HealthFacility = factor(HealthFacility, levels = hf_order),
        HealthDistrict = factor(HealthDistrict, levels = hd_order)
    )

color_scale <- scale_color_brewer(palette = "Set1")
fill_scale <- scale_fill_brewer(palette = "Set1")
estimated_annotation_color <- "black"
estimated_annotation_fill <- "black"
estimated_annotation_shape <- 21


naive_annotation_color <- "black"
naive_annotation_fill <- "black"
naive_annotation_shape <- 24
region_alpha <- 0.15

rectangle_annotations <- all_sample_res |>
    dplyr::select(HealthFacility, HealthDistrict) |>
    dplyr::distinct()

mean_moi <- ggplot() +
    geom_tile(
        data = rectangle_annotations,
        aes(
            x = HealthFacility,
            y = 15 / 2 + .5,
            height = 15,
            fill = HealthDistrict
        ),
        alpha = region_alpha,
    ) +
    geom_jitter(
        data = all_sample_res,
        aes(
            y = post_coi_med,
            x = HealthFacility,
            color = HealthDistrict
        ),
        width = 0.1,
        height = 0,
        alpha = 0.25
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = mean_mean_coi,
            x = HealthFacility
        ),
        size = 4,
        color = estimated_annotation_color,
        shape = estimated_annotation_shape,
        fill = estimated_annotation_fill
    ) +
    geom_errorbar(
        data = pop_res,
        aes(
            ymin = lower_mean_coi,
            ymax = upper_mean_coi,
            x = HealthFacility
        ),
        width = .5,
        color = estimated_annotation_color
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = naive_mean_coi,
            x = HealthFacility,
        ),
        size = 4,
        color = naive_annotation_color,
        shape = naive_annotation_shape,
        fill = naive_annotation_fill
    ) +
    labs(
        x = "Health Facility",
        y = "MOI",
        fill = "Health District"
    ) +
    theme_and_axis_legend +
    coord_cartesian(ylim = c(0.5, 15.5), expand = FALSE) +
    guides(color = "none") +
    fill_scale +
    color_scale +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
mean_moi

mean_emoi <- ggplot() +
    geom_tile(
        data = rectangle_annotations,
        aes(
            x = HealthFacility,
            y = 13 / 2 + .5,
            height = 13,
            fill = HealthDistrict
        ),
        alpha = region_alpha,
    ) +
    # geom_boxplot(
    #     data = all_sample_res,
    #     aes(
    #         y = post_effective_coi_mean,
    #         x = HealthFacility,
    #         fill = HealthDistrict,
    #     ),
    #     alpha = .1
    # ) +
    geom_jitter(
        data = all_sample_res,
        aes(
            y = post_effective_coi_mean,
            x = HealthFacility,
            color = HealthDistrict
        ),
        width = 0.1,
        height = 0,
        alpha = 0.25
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = mean_mean_ecoi,
            x = HealthFacility
        ),
        size = 4,
        color = estimated_annotation_color,
        shape = estimated_annotation_shape,
        fill = estimated_annotation_fill
    ) +
    geom_errorbar(
        data = pop_res,
        aes(
            ymin = lower_mean_ecoi,
            ymax = upper_mean_ecoi,
            x = HealthFacility
        ),
        width = .5,
        color = estimated_annotation_color
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = naive_mean_coi,
            x = HealthFacility,
        ),
        size = 4,
        color = naive_annotation_color,
        shape = naive_annotation_shape,
        fill = naive_annotation_fill
    ) +
    labs(
        x = "Health Facility",
        y = "eMOI",
        fill = "Health District"
    ) +
    theme_and_axis_legend +
    coord_cartesian(ylim = c(0.5, 13), expand = FALSE) +
    guides(color = "none", fill = "none") +
    fill_scale +
    color_scale +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
mean_emoi


all_sample_res |>
    dplyr::filter(is.na(post_relatedness_mean)) |>
    dplyr::group_by(HealthFacility) |>
    dplyr::summarise(mean(post_effective_coi_mean))


mean_relatedness <- ggplot() +
    geom_tile(
        data = rectangle_annotations,
        aes(
            x = HealthFacility,
            y = .5,
            height = 1,
            fill = HealthDistrict
        ),
        alpha = region_alpha,
    ) +
    # geom_violin(
    #     data = all_sample_res,
    #     aes(
    #         y = post_relatedness_mean,
    #         x = HealthFacility,
    #         fill = HealthDistrict,
    #     ),
    #     alpha = .1
    # ) +
    # geom_boxplot(
    #     data = all_sample_res,
    #     aes(
    #         y = post_relatedness_mean,
    #         x = HealthFacility,
    #         fill = HealthDistrict,
    #     ),
    #     alpha = .1
    # ) +
    geom_jitter(
        data = all_sample_res,
        aes(
            y = post_relatedness_mean,
            x = HealthFacility,
            alpha = prob_polyclonal,
            color = HealthDistrict
        ),
        width = 0.1,
        height = 0
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = mean_mean_rel_ignore_monoclonal,
            x = HealthFacility
        ),
        size = 4,
        color = estimated_annotation_color,
        shape = estimated_annotation_shape,
        fill = estimated_annotation_fill
    ) +
    geom_errorbar(
        data = pop_res,
        aes(
            ymin = lower_mean_rel_ignore_monoclonal,
            ymax = upper_mean_rel_ignore_monoclonal,
            x = HealthFacility
        ),
        width = .5,
        color = estimated_annotation_color
    ) +
    labs(
        x = "Health Facility",
        y = "Relatedness",
        fill = "Health District"
    ) +
    theme_and_axis_legend +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    guides(color = "none", alpha = "none") +
    fill_scale +
    color_scale +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
mean_relatedness


mean_he <- ggplot() +
    geom_tile(
        data = rectangle_annotations,
        aes(
            x = HealthFacility,
            y = .5,
            height = 1,
            fill = HealthDistrict
        ),
        alpha = region_alpha,
    ) +
    geom_line(
        data = he_res,
        aes(
            y = post_stat_mean,
            x = HealthFacility,
            group = locus
        ),
        alpha = .15
    ) +
    geom_point(
        data = he_res,
        aes(
            y = post_stat_mean,
            x = HealthFacility,
            color = HealthDistrict
        ),
        alpha = .25,
        size = 3
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = mean_mean_he,
            x = HealthFacility
        ),
        size = 5,
        color = estimated_annotation_color,
        shape = estimated_annotation_shape,
        fill = estimated_annotation_fill
    ) +
    geom_errorbar(
        data = pop_res,
        aes(
            ymin = lower_mean_he,
            ymax = upper_mean_he,
            x = HealthFacility
        ),
        width = .5,
        color = estimated_annotation_color
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = naive_mean_he,
            x = HealthFacility
        ),
        size = 5,
        color = naive_annotation_color,
        shape = naive_annotation_shape,
        fill = naive_annotation_fill
    ) +
    guides(color = "none") +
    labs(
        x = "Health Facility",
        y = "Heterozygosity",
        fill = "Health District"
    ) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    fill_scale +
    color_scale +
    theme_and_axis_legend +
    # rotate the x-axis labels 45 degrees
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
mean_he

nam_plot <- mean_moi / mean_relatedness / mean_emoi / mean_he + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") & guides(fill = guide_legend(override.aes = list(alpha = .5)), color = "none")
nam_plot

ggsave("figures/namibia_figures/namibia_by_hf_fig.pdf", nam_plot, width = 24, height = 15)




# moi_by_relatedness <- ggplot(
#     pop_res,
#     aes(
#         y = mean_mean_rel_ignore_monoclonal,
#         x = mean_mean_coi,
#         color = HealthDistrict,
#         fill = HealthDistrict
#     )
# ) +
#     geom_point() +
#     geom_errorbar(
#         aes(
#             xmin = lower_mean_coi,
#             xmax = upper_mean_coi
#         )
#     ) +
#     geom_errorbar(
#         aes(
#             ymin = lower_mean_rel_ignore_monoclonal,
#             ymax = upper_mean_rel_ignore_monoclonal
#         )
#     ) +
#     lims(y = c(0, 1))


ggplot(
    all_sample_res,
    aes(
        x = post_relatedness_mean,
        y = post_coi_mean,
        alpha = prob_polyclonal,
    )
) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~HealthFacility) +
    lims(x = c(0, 1))

# all_sample_res |>
#     dplyr::group_by(HealthFacility) |>
#     dplyr::summarise(
#         mean_relatedness_iqr = IQR(post_relatedness_mean, na.rm = TRUE),
#         mean_coi_iqr = IQR(post_coi_mean),
#         mean_ecoi_iqr = IQR(post_effective_coi_mean),
#         n = n()
#     ) |> View()

pop_res |>
    dplyr::select(HealthFacility, HealthDistrict, mean_mean_coi, mean_mean_ecoi, mean_mean_rel_ignore_monoclonal, n) |>
    View()

pop_res |>
    ggplot(aes(x = n, y = mean_mean_rel_ignore_monoclonal)) +
    geom_point() +
    geom_smooth(method = "lm")

pop_res |>
    ggplot(aes(x = n, y = mean_prop_polyclonal)) +
    geom_point() +
    geom_smooth(method = "lm")

pop_res |>
    ggplot(aes(x = mean_mean_rel_ignore_monoclonal, y = mean_prop_polyclonal)) +
    geom_point() +
    geom_smooth(method = "lm")

pop_res |>
    ggplot(aes(x = mean_mean_rel_ignore_monoclonal, y = mean_mean_coi)) +
    geom_point() +
    geom_smooth(method = "lm")

cor.test(pop_res$n, pop_res$mean_mean_coi)
cor.test(pop_res$n, pop_res$mean_mean_rel_ignore_monoclonal)
