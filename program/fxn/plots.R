#' create_box_plot_cx
#'   Create box plot of the signal (module score) of selected single gene in selected sample, grouped by cluster.
#'   Clusters are on the Y-axis, signal on the X-axis.
#'
#' @param genes_data      - data.frame containing passed genes signal data
#' @param cluster_values  - Selected object cluster data.frame
#' @param cluster_names   - Selected object user defined cluster names vector
#' @param color_scheme    - Selected object cluster color scheme vector
#' @param title           - plot title
#' @param file_name       - file name to use when downloading the plot as an image or as json
#'
#' @return canvasXpress plot object or NULL
create_box_plot_cx <- function(genes_data,
                               cluster_values,
                               cluster_names,
                               color_scheme,
                               title,
                               file_name) {
    plot <- NULL

    if (NROW(cluster_values) > 0) {
        clusters <- cluster_values %>%
            mutate(cluster = as.character(cluster))

        colors            <- color_scheme$colors
        samplesClustered  <- !color_scheme$object_colors

        if (!is.null(cluster_names)) {
            clusters <- clusters %>%
                mutate(cluster = case_when(cluster %in% names(cluster_names) ~ cluster_names[as.character(cluster)],
                                           TRUE                              ~ cluster))
            # sort cluster names based on color scheme
            colors <- assign_metadata_to_colors(colors   = colors,
                                                metadata = cluster_names)
        }

        # adjust the overlay widths based on name length
        overlay_width <- 10 + (10 * max(nchar(unique(clusters$cluster))))

        # order the boxes based on the specified color order
        clusters <- sort_plot_data_by_colors(data   = clusters,
                                             colors = colors)

        if (!is.null(genes_data)) {
            data <- bind_cols(
                genes_data[rownames(clusters), , drop = F],
                Cluster = clusters$cluster) %>%
                mutate(Cluster = as.character(Cluster))
            cx.data <- data %>%
                select(-Cluster) %>%
                rename(!!title := all_of(names(genes_data))) %>%
                as.data.frame()
            rownames(cx.data) <- rownames(data)
            cx.data           <- cx.data %>%
                t() %>%
                as.data.frame()

            smp.annot <- data[, "Cluster", drop = FALSE]

            plot <- canvasXpress(
                data                       = cx.data,
                smpAnnot                   = smp.annot,
                backgroundType             = "solid", # overrides common config
                stringSampleFactors        = list("Cluster"),
                groupingFactors            = list("Cluster"),
                graphType                  = "Boxplot", # overrides common config
                printMagnification         = 1, # overrides common config
                showLegend                 = FALSE,
                colorBy                    = "Cluster",
                colorKey                   = list(Cluster = as.list(colors)),
                smpTitle                   = "Cluster",
                xAxisTitle                 = NULL,
                xAxisShow                  = FALSE,
                xAxis2Show                 = TRUE,
                title                      = glue("{title} Normalized Expression"),
                titleScaleFontFactor       = 0.45,
                titleFontStyle             = "bold",
                smpOverlays                = list("Cluster"),
                showSampleNames            = FALSE,
                smpOverlayProperties       = list(Cluster = list(thickness = overlay_width,
                                                                 rotate    = 90,
                                                                 showBox   = FALSE)),
                showNameOverlays           = FALSE,
                smpTitleScaleFontFactor    = 1.3,
                overlayTextColor           = "white",
                overlayTextFontStyle       = "bold",
                overlayTextScaleFontFactor = 1.3,
                overlaysThickness          = 45,
                boxplotTransparency        = 1,
                boxplotOutliersRatio       = 4,
                xAxisTicks                 = 10,
                xAxisTextScaleFontFactor   = 0.7,
                boxplotBorderColor         = "black",
                decorations                = list(line = list(list(color = "gray", x = 0))),
                broadcast                  = FALSE,
                saveFilename               = file_name) %>% add_cx_common_config()
        }
    }

    plot
}


#' create_cluster_plot_cx
#'   Create spatial plot of clusters for selected sample
#'
#' @param plot_img_id           - id of tissue background image for selected sample
#' @param plot_img_dim          - dimensions of the background image plot (result of dim())
#' @param coordinates           - image spot coordinates (data.frame of imagerow, imagecol)
#' @param cluster_values        - Selected object cluster data.frame
#' @param cluster_names         - Selected object user defined cluster names vector
#' @param color_scheme          - Selected object cluster color scheme vector
#' @param file_name             - file name to use when downloading the plot as an image or as json
#' @param title                 - title text to use
#' @param subtitle              - subtitle text to use (default NULL)
#' @param plot_size             - Format the plot for large, small, xsmall sizes
#' @param selected_points_input - name of the input variable to store CX selected points into it (default = NULL)
#'
#' @return canvasXpress plot object or NULL
create_cluster_plot_cx <- function(plot_img_id,
                                   plot_img_dim,
                                   coordinates,
                                   cluster_values,
                                   cluster_names,
                                   color_scheme,
                                   file_name,
                                   title,
                                   subtitle              = NULL,
                                   plot_size             = "small",
                                   selected_points_input = NULL) {
    buffer_space <- 10 # plot padding around points
    plot         <- NULL

    if (NROW(cluster_values) > 0) {
        var_Annot <- cluster_values %>%
            mutate(cluster = as.character(cluster))

        colors <- color_scheme$colors

        if (!is.null(cluster_names)) {
            var_Annot <- var_Annot %>%
                mutate(cluster = case_when(cluster %in% names(cluster_names) ~ cluster_names[as.character(cluster)],
                                           TRUE                              ~ cluster))
            # sort pathology names based on color scheme
            colors <- assign_metadata_to_colors(colors   = colors,
                                                metadata = cluster_names)
        }

        # scaling and trimming
        spot_scale_factor <- get_spot_scaling(plot_img_dim)
        legend_text_max   <- 40

        # sort the legend by sorting the data
        var_Annot   <- sort_plot_data_by_colors(data = var_Annot, colors = colors)
        coordinates <- coordinates[rownames(var_Annot), ]

        mouse_event <- "'mousemove' : function(o, e, t) {
                                 if (o != null && o != false &&
                                     o.z != null && o.z.cluster != null &&
                                     o.z.cluster.length > 0) {
                                     t.showInfoSpan(e, '<b>Cluster</b>: ' + o.z.cluster);
                                 }
                  }"

        events <- htmlwidgets::JS(glue("{[[mouse_event]]}",
                                       .open = "[[", .close = "]]"))

        if (!is.null(selected_points_input)) {
            # this is to make sure we can store the selections on the chart as
            # a value we can access as a regular shiny input later
            point_selection_event <- glue("'select': function(o, e, t) {
                                           if (typeof o === 'boolean' && o === false) {
                                              Shiny.onInputChange('[[selected_points_input]]', null);
                                           } else {
                                              if (CanvasXpress.selector.selections > 0) {
                                                  Shiny.onInputChange('[[selected_points_input]]', Object.keys(CanvasXpress.selector.vars));
                                              } else {
                                                  Shiny.onInputChange('[[selected_points_input]]', null);
                                              }
                                           };
                                      }", .open = "[[", .close = "]]")
            events <- htmlwidgets::JS(glue("{[[mouse_event]],
                                             [[point_selection_event]]}",
                                           .open = "[[", .close = "]]"))
        }

        plot <- canvasXpress(
            data                      = coordinates,
            varAnnot                  = var_Annot,
            backgroundImage           = glue('javascript://{plot_img_id}'),
            colorBy                   = "cluster",
            colorKey                  = list("cluster" = as.list(colors)),
            stringVariableFactors     = list("cluster"),
            title                     = title,
            subtitle                  = subtitle,
            yAxis                     = list("imagerow"),
            xAxis                     = list("imagecol"),
            setMinX                   = max(0,   floor(min(coordinates$imagecol) - buffer_space)),
            setMaxX                   = min(plot_img_dim[2], ceiling(max(coordinates$imagecol) + buffer_space)),
            setMinY                   = max(0,   floor(min(coordinates$imagerow) - buffer_space)),
            setMaxY                   = min(plot_img_dim[1], ceiling(max(coordinates$imagerow) + buffer_space)),
            dataPointSizeScaleFactor  = spot_scale_factor,
            scatterOutlineThreshold   = 100000,
            outlineWidth              = 0.001,
            transparency              = 1,
            transparencyHidden        = 0.5,
            visiumHideWhenFilter      = TRUE,
            visiumFixedAspectRatio    = FALSE,
            filterMode                = "color",
            missingDataColor          = MISSING_DATA_COLOR,
            legendOrder               = list("cluster" = names(colors)),
            maxTextSize               = legend_text_max,
            movable                   = FALSE,
            broadcastGroup            = plot_img_id,
            events                    = events,
            saveFilename              = file_name) %>% add_cx_common_config()

        if (tolower(plot_size) == "large") {
            plot <- plot %>%
                add_cx_large_chart_config()
        } else if (tolower(plot_size) == "small") {
            plot <- plot %>%
                add_cx_small_chart_config()
        } else if (tolower(plot_size) == "xsmall") {
            plot <- plot %>%
                add_cx_xsmall_chart_config()
        }
    }

    plot
}


#' create_dot_plot_cx
#'   Create heatmap dot plot of normalized expression of each gene in each cluster, for selected gene signature.
#'   Genes are on the X-axis, clusters on the Y-axis. Each dot is colored by the scaled mean of expm1(Expression) of the
#'   gene in that cluster, and sized by the percent expressed. An error is generated if there is zero variance in the cluster values
#'
#' @param genes_data      - data.frame containing passed genes signal data
#' @param cluster_values  - Selected object cluster data.frame
#' @param cluster_names   - Selected object user defined cluster names vector
#' @param color_scheme    - Selected object cluster color scheme vector
#' @param features        - list of individual genes for selected gene signature
#' @param title           - plot title
#' @param file_name       - file name to use when downloading the plot as an image or as json
#'
#' @return canvasXpress plot object, error message or NULL
create_dot_plot_cx <- function(genes_data,
                               cluster_values,
                               cluster_names,
                               color_scheme,
                               features,
                               title,
                               file_name) {
    plot <- NULL

    if (NROW(cluster_values) > 0) {
        clusters <- cluster_values %>%
            mutate(cluster = as.character(cluster))

        colors            <- color_scheme$colors
        samplesClustered  <- !color_scheme$object_colors

        if (!is.null(cluster_names)) {
            clusters <- clusters %>%
                mutate(cluster = case_when(cluster %in% names(cluster_names) ~ cluster_names[as.character(cluster)],
                                           TRUE                              ~ cluster))
            # sort pathology names based on color scheme
            colors <- assign_metadata_to_colors(colors   = colors,
                                                metadata = cluster_names)
        }

        # adjust the overlay widths based on name length
        overlay_width <- 10 + (10 * max(nchar(unique(clusters$cluster))))

        # order the dots based on the specified color order
        clusters <- sort_plot_data_by_colors(data                 = clusters,
                                             colors               = colors,
                                             convert_to_character = FALSE)
        # NOTE: convert_to_character is false here to keep the factor levels

        if (!is.null(genes_data)) {
            data <- bind_cols(
                genes_data[rownames(clusters), , drop = F],
                Cluster = clusters$cluster) %>%
                gather("Gene", "Expression", any_of(features)) %>%
                group_by(Cluster, Gene) %>%
                summarise(mean_exp = mean(expm1(Expression), na.rm = TRUE),
                          pct_exp  = sum(Expression > 0, na.rm = TRUE)/length(Expression) * 100) %>%
                group_by(Gene) %>%
                mutate(mean_exp = scale(mean_exp),
                       mean_exp = case_when(mean_exp < 0 ~ 0,
                                            mean_exp > 1 ~ 1,
                                            TRUE         ~ mean_exp))
            if (all(is.nan(data$mean_exp) | is.na(data$mean_exp))) {
                plot <- tags$p(style = "text-align:center; color: #A94442; font-style: italic;",
                               get_message_text("extended_dot_plot_zero_variance"))
            } else {
                cx.mean <- data %>%
                    select(-pct_exp) %>%
                    spread(Cluster, mean_exp) %>%
                    column_to_rownames("Gene") %>%
                    as.data.frame()

                cx.pct <- data %>%
                    select(-mean_exp) %>%
                    spread(Cluster, pct_exp) %>%
                    column_to_rownames("Gene") %>%
                    as.data.frame()

                smp.annot           <- data.frame("Cluster" = colnames(cx.mean), stringsAsFactors = F)
                rownames(smp.annot) <- smp.annot$Cluster

                events <- htmlwidgets::JS("{'mousemove' : function(o, e, t) {
                                                      if ((o != null) && (o.y != null) && !(o.objectType)) {
                                                          t.showInfoSpan(e,
                                                              '<b>' + o.y.vars[0] + '</b><br/>' +
                                                              '<b>Cluster ' + o.y.smps[0] + '</b><br/>' +
                                                              'Avg Expression: ' + o.y.data[0] + '<br/>' +
                                                              'Pct Expressed: &nbsp;' + o.y.data2[0])
                                                      };
                                                  }}")

                plot <- canvasXpress(
                    data                       = list(y = cx.mean, data2 = cx.pct),
                    smpAnnot                   = smp.annot,
                    graphType                  = "Heatmap", # overrides common config
                    objectBorderColor          = "black",
                    backgroundType             = "solid", # overrides common config
                    printMagnification         = 1, # overrides common config
                    sizeBy                     = "Percent\nExpressed",
                    sizes                      = seq(3, 36, by = 3),
                    sizeByData                 = "data2",
                    sizeByContinuous           = TRUE,
                    colorSpectrum              = list("lightgray", "darkblue"),
                    colorKey                   = list(Cluster = as.list(colors)),
                    yAxisTitle                 = "Cluster",
                    xAxisTitle                 = NULL,
                    title                      = glue("{title} Normalized Expression"),
                    titleScaleFontFactor       = 0.45,
                    titleFontStyle             = "bold",
                    stringSampleFactors        = list("Cluster"),
                    smpOverlays                = list("Cluster"),
                    smpOverlayProperties       = list(Cluster = list(thickness = overlay_width,
                                                                     rotate    = 90,
                                                                     showBox   = FALSE)),
                    showNameOverlays           = FALSE,
                    smpTitleScaleFontFactor    = 0.7,
                    overlayTextColor           = "white",

                    overlayTextFontStyle       = "bold",
                    overlayTextScaleFontFactor = 1.3,
                    overlaysThickness          = 45,
                    smpTitle                   = "Cluster",
                    samplesClustered           = samplesClustered,
                    showSmpDendrogram          = FALSE,
                    variablesClustered         = ifelse(NROW(cx.mean) > 1, TRUE, FALSE),
                    showVarDendrogram          = FALSE,
                    varTextRotate              = 45,
                    varTextFontStyle           = "bold",
                    showSampleNames            = FALSE,
                    showHeatmapIndicator       = TRUE,
                    heatmapIndicatorPosition   = "top",
                    heatmapIndicatorHeight     = 15,
                    heatmapIndicatorWidth      = 600,
                    legendTitleAlign           = "center",
                    legendTitleMargin          = 25,
                    marginBottom               = 30,
                    broadcast                  = FALSE,
                    events                     = events,
                    saveFilename               = file_name) %>% add_cx_common_config()
            }
        }
    }

    plot
}


#' create_metadata_plot_cx
#'   Create spatial plot of passed metadata for selected sample.
#'
#' @param plot_img_id           - id of tissue background image for selected sample
#' @param plot_img_dim          - dimensions of the background image plot (result of dim())
#' @param coordinates           - image spot coordinates (data.frame of imagerow, imagecol)
#' @param meta_data             - dataframe from seurat object
#' @param file_name             - file name to use when downloading the plot as an image or as json
#' @param title                 - title text to use (default NULL)
#' @param subtitle              - subtitle text to use (default NULL)
#' @param plot_size             - Format the plot for "large" or "small" sizes
#' @param selected_points_input - name of the input variable to store CX selected points into it (default = NULL)
#' @param color_key             - color key for categorical metadata (default = NULL)
#' @param color_spectrum_breaks - color spectrum breaks for numeric metadata (default = NULL)
#'
#' @return canvasXpress plot object or NULL
create_metadata_plot_cx <- function(plot_img_id,
                                    plot_img_dim,
                                    coordinates,
                                    meta_data,
                                    file_name,
                                    title                 = NULL,
                                    subtitle              = NULL,
                                    plot_size             = "small",
                                    selected_points_input = NULL,
                                    color_key             = NULL,
                                    color_spectrum_breaks = NULL) {
    plot <- NULL

    if (NROW(meta_data) > 0) {
        buffer_space     <- 10 # plot padding around points
        metadata_clean   <- str_replace_all(names(meta_data), "\\.", "_") # clean the names for JS use
        names(meta_data) <- metadata_clean

        # scaling and trimming
        spot_scale_factor <- get_spot_scaling(plot_img_dim)

        mouse_event <- glue("'mousemove' : function(o, e, t) {
                                 if (o != null && o != false &&
                                     o.z != null && o.z.[metadata_clean] != null &&
                                     o.z.[metadata_clean].length > 0) {
                                     t.showInfoSpan(e, '<b>[metadata_clean]</b>: ' + o.z.[metadata_clean]);
                                 }
                  }", .open = "[", .close = "]")

        events <- htmlwidgets::JS(glue("{[[mouse_event]]}",
                                       .open = "[[", .close = "]]"))

        if (!is.null(selected_points_input)) {
            # this is to make sure we can store the selections on the chart as
            # a value we can access as a regular shiny input later
            point_selection_event <- glue("'select': function(o, e, t) {
                                           if (typeof o === 'boolean' && o === false) {
                                              Shiny.onInputChange('[[selected_points_input]]', null);
                                           } else {
                                              if (CanvasXpress.selector.selections > 0) {
                                                  Shiny.onInputChange('[[selected_points_input]]', Object.keys(CanvasXpress.selector.vars));
                                              } else {
                                                  Shiny.onInputChange('[[selected_points_input]]', null);
                                              }
                                           };
                                      }", .open = "[[", .close = "]]")
            events <- htmlwidgets::JS(glue("{[[mouse_event]],
                                             [[point_selection_event]]}",
                                           .open = "[[", .close = "]]"))
        }

        if (!is.numeric(meta_data[[metadata_clean]])) {
            meta_data[[metadata_clean]] <- meta_data[[metadata_clean]] %>%
                as.character() %>%
                replace_na("NA") %>% # use "NA" string to work with the color key
                as.factor()
        }

        color_key_cx <- list()
        legend_order <- list()
        if (length(color_key) > 0) {
            color_key_cx[[metadata_clean]] <- as.list(color_key)
            legend_order[[metadata_clean]] <- names(color_key)
        }

        var_Annot   <- meta_data
        coordinates <- coordinates[rownames(var_Annot), ]

        plot <- canvasXpress(
            data                     = coordinates,
            varAnnot                 = var_Annot,
            backgroundImage          = glue('javascript://{plot_img_id}'),
            colorBy                  = metadata_clean,
            colorKey                 = color_key_cx,
            colorSpectrum            = COLOR_SCHEME_SPECTRUM_SCALING,
            colorSpectrumBreaks      = color_spectrum_breaks,
            legendOrder              = legend_order,
            yAxis                    = list("imagerow"),
            xAxis                    = list("imagecol"),
            setMinX                  = max(0,   floor(min(coordinates$imagecol) - buffer_space)),
            setMaxX                  = min(plot_img_dim[2], ceiling(max(coordinates$imagecol) + buffer_space)),
            setMinY                  = max(0,   floor(min(coordinates$imagerow) - buffer_space)),
            setMaxY                  = min(plot_img_dim[1], ceiling(max(coordinates$imagerow) + buffer_space)),
            dataPointSizeScaleFactor = spot_scale_factor,
            scatterOutlineThreshold  = 100000,
            outlineWidth             = 0.001,
            transparency             = 1,
            transparencyHidden       = 0.5,
            visiumHideWhenFilter     = TRUE,
            visiumFixedAspectRatio   = FALSE,
            missingDataColor         = MISSING_DATA_COLOR,
            title                    = title,
            subtitle                 = subtitle,
            broadcastGroup           = plot_img_id,
            events                   = events,
            saveFilename             = file_name) %>% add_cx_common_config()

        if (tolower(plot_size) == "large") {
            plot <- plot %>%
                add_cx_large_chart_config()
        } else if (tolower(plot_size) == "xsmall") {
            plot <- plot %>%
                add_cx_xsmall_chart_config()
        }

    }

    plot
}


#' create_pathology_plot_cx
#'   Create spatial plot of pathology groups for selected sample
#'
#' @param plot_img_id           - id of tissue background image for selected sample
#' @param plot_img_dim          - dimensions of the background image plot (result of dim())
#' @param coordinates           - image spot coordinates (data.frame of imagerow, imagecol)
#' @param pathology_values      - Selected object pathology data.frame
#' @param pathology_names       - Selected object user defined pathology names vector
#' @param color_scheme          - Selected object pathology color scheme vector
#' @param file_name             - file name to use when downloading the plot as an image or as json
#' @param title                 - title text to use
#' @param subtitle              - subtitle text to use (default NULL)
#' @param plot_size             - Format the plot for large, small or xsmall sizes
#' @param selected_points_input - name of the input variable to store CX selected points into it (default = NULL)
#'
#' @return canvasXpress plot object or NULL
create_pathology_plot_cx <- function(plot_img_id,
                                     plot_img_dim,
                                     coordinates,
                                     pathology_values,
                                     pathology_names,
                                     color_scheme,
                                     file_name,
                                     title,
                                     subtitle              = NULL,
                                     plot_size,
                                     selected_points_input = NULL) {
    buffer_space <- 10 # plot padding around points
    plot         <- NULL

    if (NROW(pathology_values) > 0) {
        var_Annot <- pathology_values %>%
            mutate(pathology = as.character(pathology))

        colors <- color_scheme$colors

        if (!is.null(pathology_names)) {
            var_Annot <- var_Annot %>%
                mutate(pathology = case_when(pathology %in% names(pathology_names) ~ pathology_names[as.character(pathology)],
                                             TRUE                                  ~ pathology))
            # sort pathology names based on color scheme
            colors <- assign_metadata_to_colors(colors   = colors,
                                                metadata = pathology_names)
        }

        # scaling and trimming
        spot_scale_factor <- get_spot_scaling(plot_img_dim)
        legend_text_max   <- 40

        # sort the legend by sorting the data
        var_Annot    <- sort_plot_data_by_colors(data = var_Annot, colors = colors)
        coordinates  <- coordinates[row.names(var_Annot), ]

        mouse_event  <- "'mousemove' : function(o, e, t) {
                                        if (o != null && o != false &&
                                            o.z != null && o.z.pathology != null &&
                                            o.z.pathology.length > 0) {
                                            t.showInfoSpan(e, '<b>Pathology</b>: ' + o.z.pathology);
                                    }}"

        events <- htmlwidgets::JS(glue("{[[mouse_event]]}",
                                       .open = "[[", .close = "]]"))

        if (!is.null(selected_points_input)) {
            # this is to make sure we can store the selections on the chart as
            # a value we can access as a regular shiny input later
            point_selection_event <- glue("'select': function(o, e, t) {
                                           if (typeof o === 'boolean' && o === false) {
                                              Shiny.onInputChange('[[selected_points_input]]', null);
                                           } else {
                                              if (CanvasXpress.selector.selections > 0) {
                                                  Shiny.onInputChange('[[selected_points_input]]', Object.keys(CanvasXpress.selector.vars));
                                              } else {
                                                  Shiny.onInputChange('[[selected_points_input]]', null);
                                              }
                                           };
                                      }", .open = "[[", .close = "]]")
            events <- htmlwidgets::JS(glue("{[[mouse_event]],
                                             [[point_selection_event]]}",
                                           .open = "[[", .close = "]]"))
        }

        plot <- canvasXpress(
            data                      = coordinates,
            varAnnot                  = var_Annot,
            backgroundImage           = glue('javascript://{plot_img_id}'),
            colorBy                   = "pathology",
            colorKey                  = list("pathology" = as.list(colors)),
            legendOrder               = list("pathology" = names(colors)),
            stringVariableFactors     = list("pathology"),
            title                     = title,
            subtitle                  = subtitle,
            yAxis                     = list("imagerow"),
            xAxis                     = list("imagecol"),
            setMinX                   = max(0,   floor(min(coordinates$imagecol) - buffer_space)),
            setMaxX                   = min(plot_img_dim[2], ceiling(max(coordinates$imagecol) + buffer_space)),
            setMinY                   = max(0,   floor(min(coordinates$imagerow) - buffer_space)),
            setMaxY                   = min(plot_img_dim[1], ceiling(max(coordinates$imagerow) + buffer_space)),
            dataPointSizeScaleFactor  = spot_scale_factor,
            scatterOutlineThreshold   = 100000,
            outlineWidth              = 0.001,
            transparency              = 1,
            transparencyHidden        = 0.5,
            visiumHideWhenFilter      = TRUE,
            visiumFixedAspectRatio    = FALSE,
            colors                    = as.list(colors),
            filterMode                = "color",
            missingDataColor          = MISSING_DATA_COLOR,
            maxTextSize               = legend_text_max,
            movable                   = FALSE,
            broadcastGroup            = plot_img_id,
            events                    = events,
            saveFilename              = file_name) %>% add_cx_common_config()

        if (tolower(plot_size) == "large") {
            plot <- plot %>%
                add_cx_large_chart_config()
        } else if (tolower(plot_size) == "small") {
            plot <- plot %>%
                add_cx_small_chart_config()
        } else if (tolower(plot_size) == "xsmall") {
            plot <- plot %>%
                add_cx_xsmall_chart_config()
        }
    }

    plot
}


#' create_spatial_plot_cx
#'   Create spatial plot of signal (module score) for selected gene signature or single gene in selected sample
#'
#' @param plot_img_id           - id of tissue background image for selected sample
#' @param plot_img_dim          - dimensions of the background image plot (result of dim())
#' @param features              - list of individual genes for selected gene signature or single gene
#' @param coordinates           - image spot coordinates (data.frame of imagerow, imagecol)
#' @param signal                - spot values
#' @param title                 - plot title
#' @param subtitle              - plot subtitle
#' @param top_genes             - data.frame of top represented genes for each spot in the sample
#' @param scaling               - selected option for color scale range: "native", "trim", "center", or list of min and max values
#' @param file_name             - file name to use when downloading the plot as an image or as json
#' @param coloring_signal       - dataframe that used for coloring for some of matrix color scaling options (default = NULL)
#' @param plot_size             - Format the plot for large, small or xsmall sizes
#' @param selected_points_input - name of the input variable to store CX selected points into it (default = NULL)
#'
#' @return canvasXpress plot object or NULL
create_spatial_plot_cx <- function(plot_img_id,
                                   plot_img_dim,
                                   features,
                                   coordinates, # data.frame with imagerow, imagecol
                                   signal,      # must be in the SAME coordinate order as top_genes, due to legacy Seurat fetching
                                   title,
                                   subtitle,
                                   top_genes,
                                   scaling,
                                   file_name,
                                   coloring_signal = NULL,
                                   plot_size,
                                   selected_points_input = NULL) {
    buffer_space <- 10   # plot padding around points
    plot         <- NULL

    if (!is.null(signal)) {
        # get per-spot top genes
        spot_cal <- NULL

        if (NROW(top_genes) > 0) {
            if (length(features) > 1) {
                #filter by signature genes
                top_genes <- top_genes %>%
                    filter(Symbol %in% features)
            }

            spot_cal <- top_genes %>%
                group_by(Coord) %>%
                top_n(g_display_top_n_spot_genes, Value) %>%
                group_by(Coord) %>%
                arrange(desc(Value)) %>%
                summarise(top = glue_collapse(Symbol, sep = ", "))
        }

        if (is.null(spot_cal)) {
            spot_cal <- data.frame(Coord = rownames(coordinates), top = "")
        }

        # ensure every spot has a value - this is needed because we are using raw
        # expression data for genes instead of scores now

        if (NROW(signal) != NROW(coordinates)) {
            signal <- coordinates %>%
                mutate(cellid = as.character(cellid)) %>%
                left_join(signal %>% rownames_to_column("cellid"), by = "cellid") %>%
                column_to_rownames("cellid") %>%
                select(all_of(colnames(signal)))
        }

        if (length(rownames(signal)) == NROW(signal)) {
            signal <- data.frame(Coord = as.character(rownames(signal)),
                                 signal = setNames(signal, NULL))
            var_Annot <- left_join(signal, spot_cal, by = "Coord") %>%
                as.data.frame() %>%
                arrange(ordered(Coord, rownames(coordinates)))
            rownames(var_Annot) <- var_Annot$Coord
            var_Annot$Coord <- NULL
        } else {
            # this happens with seurat objects, we have to assume they are in the same position, we can't get names on both
            var_Annot <- base::merge(signal, spot_cal, by.x = 0, by.y = "Coord", all.x = TRUE) %>%
                column_to_rownames("Row.names") %>%
                rename(signal = 1)
        }

        if (is.null(coloring_signal)) {
            coloring_signal <- var_Annot$signal
        }

        dr_continue  <- quantile(coloring_signal, probs = c(0, 0.01, 0.99, 1), na.rm = TRUE, names = FALSE)
        dr_continue  <- abs(dr_continue[3] - dr_continue[2]) >= g_display_dynamic_range_min

        # scaling and trimming
        if (length(plot_img_dim) < 2) {
            plot_img_dim <- c(600, 600)
        }
        spot_scale_factor <- get_spot_scaling(plot_img_dim)
        legend_text_max   <- 40

        # disallow small DR on single-gene plots only
        if ((length(features) > 1) || (dr_continue)) {
            breaks <- get_color_scale_breaks(signal      = coloring_signal,
                                             scaling     = scaling,
                                             color_count = length(COLOR_SCHEME_SPECTRUM_SCALING))
            mouse_event <- glue("'mousemove' : function(o, e, t) {
                                 if (o != null && o != false &&
                                     o.z != null && o.z.signal != null &&
                                     o.z.signal.length > 0) {
                                     if (o.z.top[0] == null) {
                                       top_genes = '';
                                     } else {
                                       top_genes = o.z.top[0];
                                     }

                                     in_matrix = '[[is.null(top_genes)]]';
                                     signal    = Math.round((o.z.signal[0] + Number.EPSILON) * 10000)/10000;
                                     spot_tip  = '<b>Signal</b>: ' + signal;

                                     if (in_matrix == 'FALSE') {
                                         spot_tip =  spot_tip + '<br/>' + '<b> Top Genes</b>: ' + top_genes;
                                     }

                                     t.showInfoSpan(e, spot_tip);
                                 }
                           }", .open = "[[", .close = "]]")
            events <- htmlwidgets::JS(glue("{[[mouse_event]]}",
                                           .open = "[[", .close = "]]"))

            if (!is.null(selected_points_input)) {
                # this is to make sure we can store the selections on the chart as
                # a value we can access as a regular shiny input later
                point_selection_event <- glue("'select': function(o, e, t) {
                                           if (typeof o === 'boolean' && o === false) {
                                              Shiny.onInputChange('[[selected_points_input]]', null);
                                           } else {
                                              if (CanvasXpress.selector.selections > 0) {
                                                  Shiny.onInputChange('[[selected_points_input]]', Object.keys(CanvasXpress.selector.vars));
                                              } else {
                                                  Shiny.onInputChange('[[selected_points_input]]', null);
                                              }
                                           };
                                      }", .open = "[[", .close = "]]")
                events <- htmlwidgets::JS(glue("{[[mouse_event]],
                                             [[point_selection_event]]}",
                                               .open = "[[", .close = "]]"))
            }

            plot <- canvasXpress(
                data                      = coordinates,
                varAnnot                  = var_Annot,
                backgroundImage           = glue('javascript://{plot_img_id}'),
                colorBy                   = "signal",
                title                     = title,
                subtitle                  = subtitle,
                yAxis                     = list("imagerow"),
                xAxis                     = list("imagecol"),
                setMinX                   = max(0,   floor(min(coordinates$imagecol) - buffer_space)),
                setMaxX                   = min(plot_img_dim[2], ceiling(max(coordinates$imagecol) + buffer_space)),
                setMinY                   = max(0,   floor(min(coordinates$imagerow) - buffer_space)),
                setMaxY                   = min(plot_img_dim[1], ceiling(max(coordinates$imagerow) + buffer_space)),
                dataPointSizeScaleFactor  = spot_scale_factor,
                scatterOutlineThreshold   = 100000,
                outlineWidth              = 0.001,
                transparency              = 1,
                transparencyHidden        = 0.5,
                visiumHideWhenFilter      = TRUE,
                visiumFixedAspectRatio    = FALSE,
                colorSpectrum             = COLOR_SCHEME_SPECTRUM_SCALING,
                colorSpectrumBreaks       = breaks,
                filterMode                = "color",
                missingDataColor          = MISSING_DATA_COLOR,
                maxTextSize               = legend_text_max,   # needed to harmonize the title text sizes in Advanced
                movable                   = FALSE,
                broadcastGroup            = plot_img_id,
                heatmapIndicatorWidth     = 300,
                events                    = events,
                saveFilename              = file_name) %>% add_cx_common_config()

            if (tolower(plot_size) == "large") {
                plot <- plot %>%
                    add_cx_large_chart_config()
            } else if (tolower(plot_size) == "small") {
                plot <- plot %>%
                    add_cx_small_chart_config()
            } else if (tolower(plot_size) == "xsmall") {
                plot <- plot %>%
                    add_cx_xsmall_chart_config()
            }
        }
    }

    plot
}


#' create_tissue_plot_cx
#'   Create spatial plot of tissue background image for selected sample
#'
#' @param plot_img_id     - id of tissue background image for selected sample
#' @param plot_img_dim    - dimensions of the background image plot (result of dim())
#' @param coordinates     - image spot coordinates (data.frame of imagerow, imagecol)
#' @param file_name       - file name to use when downloading the plot as an image or as json
#' @param title           - plot title (default - "Tissue")
#' @param plot_size       - Format the plot for small or xsmall sizes
#'
#' @return canvasXpress plot object or NULL
create_tissue_plot_cx <- function(plot_img_id,
                                  plot_img_dim,
                                  coordinates,
                                  file_name,
                                  title           = "Tissue",
                                  plot_size) {
    buffer_space <- 10 # plot padding around points
    plot         <- NULL

    # scaling and trimming
    legend_text_max <- 40

    plot <- canvasXpress(
        data                      = data.frame(imagerow = -1, imagecol = -1),
        backgroundImage           = glue('javascript://{plot_img_id}'),
        title                     = title,
        setMinX                   = max(0,   floor(min(coordinates$imagecol) - buffer_space)),
        setMaxX                   = min(plot_img_dim[2], ceiling(max(coordinates$imagecol) + buffer_space)),
        setMinY                   = max(0,   floor(min(coordinates$imagerow) - buffer_space)),
        setMaxY                   = min(plot_img_dim[1], ceiling(max(coordinates$imagerow) + buffer_space)),
        broadcast                 = FALSE,
        movable                   = FALSE,
        visiumFixedAspectRatio    = FALSE,
        maxTextSize               = legend_text_max,   # needed to harmonize the title text sizes in Advanced
        saveFilename              = file_name) %>% add_cx_common_config()

    if (tolower(plot_size) == "small") {
        plot <- plot %>%
            add_cx_small_chart_config()
    } else if (tolower(plot_size) == "xsmall") {
        plot <- plot %>%
            add_cx_xsmall_chart_config()
    }

    plot
}
