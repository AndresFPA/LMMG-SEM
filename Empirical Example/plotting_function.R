library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)

#' Generate Alluvial and Geographical Transition Plots for LMM Fit
#'
#' @param lmm_fit An LMM model object containing a 3D posterior array `$post` [groups, clusters, times].
#' @param country_codes Vector of unique country ISO-2 codes (length G).
#' @param custom_palette Optional named or unnamed character vector of hex colors (length K).
#' @param crop_bbox Named numeric vector for bounding box (xmin, xmax, ymin, ymax).
#' @param sem_model_syntax Optional character string with lavaan syntax for the structural model.
#' @param cluster_betas Optional list of matrices containing cluster-specific regression coefficients.
#'
#' @return A named list containing `alluvial_plot`, `map_plot`, and `path_plots`.
plot_lmm_results <- function(lmm_fit,
                             country_codes,
                             custom_palette = NULL,
                             crop_bbox = c(xmin = -11, xmax = 36, ymin = 31, ymax = 71),
                             sem_model_syntax = NULL,
                             cluster_betas = NULL) {
  
  # ---------------------------------------------------------------------------
  # 1. Dimension & Metadata Extraction
  # ---------------------------------------------------------------------------
  posteriors <- round(lmm_fit$post, 4)
  ngroups    <- dim(posteriors)[1]
  k_clusters <- dim(posteriors)[2]
  ntimes     <- dim(posteriors)[3]
  
  if (length(country_codes) != ngroups) {
    stop(sprintf("Length of country_codes (%d) must match number of groups (%d).", 
                 length(country_codes), ngroups))
  }
  
  time_cols   <- paste0("Time_", 1:ntimes)
  time_labels <- paste0("Time ", 1:ntimes)
  cluster_lvls <- as.character(1:k_clusters)
  cluster_lbls <- paste0("Cluster ", 1:k_clusters)
  
  # Color Palette Resolution
  if (is.null(custom_palette)) {
    if (k_clusters <= 8) {
      cluster_colors <- RColorBrewer::brewer.pal(max(3, k_clusters), "Dark2")[1:k_clusters]
    } else {
      cluster_colors <- scales::hue_pal()(k_clusters)
    }
  } else {
    if (length(custom_palette) < k_clusters) {
      stop(sprintf("custom_palette must contain at least %d colors.", k_clusters))
    }
    cluster_colors <- custom_palette[1:k_clusters]
  }
  names(cluster_colors) <- cluster_lvls
  
  # ---------------------------------------------------------------------------
  # 2. Extract Discrete Cluster Assignments
  # ---------------------------------------------------------------------------
  assignments_mat <- matrix(NA, nrow = ngroups, ncol = ntimes)
  for (t in 1:ntimes) {
    assignments_mat[, t] <- apply(posteriors[, , t, drop = FALSE], 1, which.max)
  }
  
  alluvial_data <- as.data.frame(assignments_mat)
  colnames(alluvial_data) <- time_cols
  alluvial_data$Country   <- country_codes
  
  for (col in time_cols) {
    alluvial_data[[col]] <- factor(alluvial_data[[col]], levels = cluster_lvls)
  }
  
  # ---------------------------------------------------------------------------
  # 3. Build Alluvial Plot
  # ---------------------------------------------------------------------------
  data_lodes <- ggalluvial::to_lodes_form(
    alluvial_data,
    axes = time_cols,
    id = "Country_id",
    key = "Wave",
    value = "Cluster",
    discern = FALSE
  )
  
  data_lodes <- data_lodes %>%
    group_by(Country) %>%
    mutate(Baseline_Cluster = factor(Cluster[Wave == "Time_1"], levels = cluster_lvls)) %>%
    ungroup()
  
  p_alluvial <- ggplot(data_lodes, aes(x = Wave, stratum = Cluster, alluvium = Country, y = 1)) +
    geom_alluvium(
      aes(fill = Baseline_Cluster),
      width = 0.22,
      alpha = 0.75,
      knot.pos = 0.35,
      color = "white",
      linewidth = 0.35
    ) +
    geom_text(
      stat = "alluvium",
      aes(label = Country),
      size = 3,
      fontface = "plain",
      color = "white",
      nudge_x = 0.16
    ) +
    geom_stratum(
      width = 0.22,
      fill = "white",
      color = "grey35",
      linewidth = 0.6
    ) +
    geom_text(
      stat = "stratum",
      aes(label = paste0("Cluster ", after_stat(stratum))),
      size = 3,
      fontface = "bold",
      color = "grey20"
    ) +
    scale_x_discrete(
      limits = time_cols,
      labels = time_labels,
      expand = c(0.12, 0.08)
    ) +
    scale_fill_manual(
      values = cluster_colors,
      labels = cluster_lbls,
      name = "Baseline Cluster (T1)",
      drop = FALSE
    ) +
    labs(
      title = sprintf("Latent Cluster Transitions (N = %d)", ngroups),
      x = "Wave",
      y = "Number of Countries"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 13)),
      axis.text.y = element_text(color = "grey30")
    )
  
  # ---------------------------------------------------------------------------
  # 4. Build Spatial/Geographical Map Plot (Updated & Robust)
  # ---------------------------------------------------------------------------
  df_map_long <- alluvial_data %>%
    pivot_longer(
      cols = starts_with("Time_"),
      names_to = "Wave",
      values_to = "Cluster"
    ) %>%
    mutate(
      Wave = factor(Wave, levels = time_cols, labels = time_labels),
      Cluster = factor(Cluster, levels = cluster_lvls)
    )
  
  # Load and sanitize geometries
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
    sf::st_make_valid() %>%
    filter(!sf::st_is_empty(geometry)) %>%
    mutate(
      iso2 = case_when(
        admin == "Portugal"       ~ "PT",
        admin == "France"         ~ "FR",
        admin == "Norway"         ~ "NO",
        admin == "United Kingdom" ~ "GB",
        admin == "Israel"         ~ "IL",
        TRUE                      ~ coalesce(iso_a2_eh, iso_a2)
      )
    ) %>%
    select(iso2, name = admin, geometry)
  
  # Define bbox explicitly with WGS84 CRS
  bbox_crop <- sf::st_bbox(
    c(
      xmin = as.numeric(crop_bbox["xmin"]),
      ymin = as.numeric(crop_bbox["ymin"]),
      xmax = as.numeric(crop_bbox["xmax"]),
      ymax = as.numeric(crop_bbox["ymax"])
    ),
    crs = sf::st_crs(4326)
  )
  
  # Safely crop and heal boundary lines
  europe_base <- suppressWarnings(
    sf::st_crop(world, bbox_crop)
  ) %>%
    sf::st_make_valid() %>%
    filter(!sf::st_is_empty(geometry))
  
  map_all_waves <- tidyr::crossing(
    Wave = factor(time_labels, levels = time_labels),
    iso2 = europe_base$iso2
  ) %>%
    left_join(europe_base, by = "iso2") %>%
    left_join(df_map_long, by = c("iso2" = "Country", "Wave" = "Wave")) %>%
    sf::st_as_sf()
  
  p_map <- ggplot(map_all_waves) +
    geom_sf(aes(fill = Cluster), color = "white", linewidth = 0.25) +
    scale_fill_manual(
      values = cluster_colors,
      labels = cluster_lbls,
      na.value = "#ECECEC",
      name = "Latent Cluster",
      drop = FALSE
    ) +
    coord_sf(
      xlim = c(crop_bbox["xmin"] + 0.5, crop_bbox["xmax"] - 0.5),
      ylim = c(crop_bbox["ymin"] + 0.5, crop_bbox["ymax"] - 3),
      expand = FALSE
    ) +
    facet_wrap(~Wave, nrow = 1) +
    labs(
      title = "Geographic Distribution of Structural Clusters Across Waves",
      caption = "Gray indicates countries not included in the sample."
    ) +
    theme_void(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 11, margin = margin(b = 6)),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13, margin = margin(b = 10)),
      plot.caption = element_text(size = 8, color = "grey40", hjust = 0.5, margin = margin(t = 8)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      panel.spacing = unit(0.5, "lines")
    )
  
  # --------------------------------------------------------------------------
  # Path Diagram Plots (using DiagrammeR)
  # --------------------------------------------------------------------------
  cluster_betas <- lmm_fit$param$beta_ks
  path_plot_html <- NULL
  path_plots <- list()
  
  if (!is.null(cluster_betas) && !is.null(sem_model_syntax)) {
    structural_pars <- lavaan::lavaanify(sem_model_syntax)
    regression_pars <- structural_pars[structural_pars$op == "~", c("lhs", "rhs")]
    model_vars <- unique(c(regression_pars$lhs, regression_pars$rhs))
    sem_models <- vector("list", length(cluster_betas))
    names(sem_models) <- paste("Cluster", seq_along(cluster_betas))

    for (i in 1:length(cluster_betas)) {
      beta_mat <- cluster_betas[[i]]

      dummy_df <- as.data.frame(
        setNames(replicate(length(model_vars), rnorm(20), simplify = FALSE), model_vars)
      )

      fit <- lavaan::sem(
        sem_model_syntax,
        data = dummy_df,
        do.fit = FALSE
      )

      mod <- semPlot::semPlotModel(fit)

      for (path_index in seq_len(nrow(regression_pars))) {
        lhs <- regression_pars$lhs[path_index]
        rhs <- regression_pars$rhs[path_index]

        if (lhs %in% rownames(beta_mat) && rhs %in% colnames(beta_mat)) {
          mod@Pars$est[mod@Pars$lhs == rhs & mod@Pars$rhs == lhs & mod@Pars$edge == "~>"] <-
            beta_mat[lhs, rhs]
        }
      }

      sem_models[[i]] <- mod
    }

    n_columns <- min(2, length(sem_models))
    n_rows <- ceiling(length(sem_models) / n_columns)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(n_rows, n_columns))

    for (i in seq_along(sem_models)) {
      semPlot::semPaths(
        sem_models[[i]],
        whatLabels = "est",
        layout = "tree",
        shapeMan = "circle",
        shapeLat = "circle",
        sizeMan = 20,
        sizeLat = 20,
        nCharNodes = 0,
        intercepts = FALSE,
        residuals = FALSE,
        edge.label.cex = 1.2,
        title = FALSE
      )
      graphics::title(main = names(sem_models)[i])
    }
    path_plots[[1]] <- grDevices::recordPlot()
  }

  if (!is.null(cluster_betas)) {
    full_graph_def <- "
    digraph SEM_All_Clusters {
      graph [rankdir = LR, layout = dot, nodesep = 0.5, ranksep = 0.5]
      node [shape = circle, style = solid, fontname = Helvetica, fontsize = 10, width = 0.6, height = 0.6]
      edge [fontname = Helvetica, fontsize = 8, color = black, arrowsize = 0.6]
    "
    
    for (i in 1:length(cluster_betas)) {
      beta_mat <- cluster_betas[[i]]
      
      full_graph_def <- paste0(full_graph_def, "  subgraph cluster_", i, " { \n")
      full_graph_def <- paste0(full_graph_def, "    label = 'Cluster ", i, "'; \n")
      
      # Only iterate over the paths you care about:
      # e.g., rows/cols that are not 0 in your structural model
      for (rhs in colnames(beta_mat)) {
        for (lhs in rownames(beta_mat)) {
          val <- beta_mat[lhs, rhs]
          # Only add edge if there is a substantive path (non-zero)
          if (!is.na(val) && abs(val) > 0.001) {
            full_graph_def <- paste0(full_graph_def, "    ", rhs, "_cl", i, " -> ", lhs, "_cl", i, " [label='", round(val, 2), "'];\n")
          }
        }
      }
      full_graph_def <- paste0(full_graph_def, "  }\n")
    }
    full_graph_def <- paste0(full_graph_def, "}")
    
    path_plot_html <- DiagrammeR::grViz(full_graph_def)
  }
  
  return(list(
    alluvial_plot   = p_alluvial,
    map_plot        = p_map,
    path_plots      = path_plots,
    path_plots_html = path_plot_html
  ))
}