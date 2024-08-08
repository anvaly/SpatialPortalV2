#' stMatrixPlotServer
#'   Module 7sServer functionality for handling spatial plots matrix cell
#'
#' @param id              - spatial plot id
#' @param sample_name     - selected sample name
#' @param feature_name    - signature name, single gene value or NULL in case of rendering object metadata field
#' @param signature_genes - vector of individual genes for selected signature (NULL if the plot is for single gene or metadata field)
#' @param slice_image     - Tissue background image for selected sample
#' @param coordinates     - image spot coordinates (data.frame of imagerow, imagecol)
#' @param signal          - spot values
#' @param metadata_type   - metadata field that is used to display its plot (i.e Pathology.Group' and 'Clustering', ...)
#' @param meta_names      - Selected metadata user defined names vector (currently NULL for any metadata field except 'Pathology.Group' and 'Clustering')
#' @param meta_colors     - Selected metadata color scheme vector (currently NULL for any metadata field except 'Pathology.Group' and 'Clustering')
#' @param coloring        - Selected color scaling option (NULL if the plot is for metadata field)
#' @param coloring_signal - Calculated coloring signal values based on both spatial plot coloring and spatial plot scaling values
#' @param logger          - reactive logger S4 object
#'
stMatrixPlotServer <- function(id,
                               sample_name,
                               feature_name,
                               signature_genes,
                               slice_image,
                               coordinates,
                               signal,
                               metadata_type,
                               meta_names,
                               meta_colors,
                               coloring,
                               coloring_signal,
                               logger) {
    moduleServer(
        id,
        function(input, output, session) {
            internal <- reactiveValues()
                internal$id              <- id
                internal$sample          <- sample_name
                internal$feature         <- ifelse(!is.null(feature_name), feature_name, metadata_type)
                internal$signature_genes <- signature_genes
                internal$slice           <- slice_image
                internal$coord           <- coordinates
                internal$signal          <- signal
                internal$meta_type       <- metadata_type
                internal$meta_names      <- meta_names
                internal$meta_colors     <- meta_colors
                internal$coloring        <- coloring
                internal$coloring_signal <- coloring_signal

            output$st_plot <- renderUI({
                plot_id      <- glue("{internal$id}_plot")
                plot_img_id  <- glue("{plot_id}_{SUFFIX_SPATIAL_PLOT_ID}")

                bkg_img_file <- tempfile(pattern = plot_img_id, fileext = CX_BKG_IMG_FILE_EXT)
                write_png_bkg_image(raw_img = internal$slice, target_file = bkg_img_file)

                if (is.null(internal$meta_type)) {
                    if (is.null(internal$signature_genes)) {
                        genes <- internal$feature
                    } else {
                        genes <- internal$signature_genes
                    }

                    plot <- create_spatial_plot_cx(
                        plot_img_id     = session$ns(plot_img_id),
                        plot_img_dim    = dim(internal$slice),
                        coordinates     = internal$coord,
                        signal          = internal$signal,
                        features        = genes,
                        title           = "",
                        subtitle        = "",
                        top_genes       = NULL,
                        scaling         = internal$coloring,
                        file_name       = get_download_filename(text = c(internal$sample, internal$feature)),
                        coloring_signal = internal$coloring_signal,
                        plot_size       = "xsmall")
                } else {
                    if (internal$meta_type == METADATA_ITEM_TISSUE) {
                        plot <- create_tissue_plot_cx(
                            plot_img_id     = session$ns(plot_img_id),
                            plot_img_dim    = dim(internal$slice),
                            coordinates     = internal$coord,
                            file_name       = get_download_filename(text = c(internal$sample, internal$feature)),
                            title           = "",
                            plot_size       = "xsmall")
                    } else if (internal$meta_type == METADATA_ITEM_CLUSTERING) {
                        plot <- create_cluster_plot_cx(
                            plot_img_id     = session$ns(plot_img_id),
                            plot_img_dim    = dim(internal$slice),
                            coordinates     = internal$coord,
                            cluster_values  = internal$signal,
                            cluster_names   = internal$meta_names,
                            color_scheme    = internal$meta_colors,
                            file_name       = get_download_filename(text = c(internal$sample, internal$feature)),
                            title           = "",
                            plot_size       = "xsmall")
                    } else if (internal$meta_type == METADATA_ITEM_PATHOLOGY) {
                        plot <- create_pathology_plot_cx(
                            plot_img_id      = session$ns(plot_img_id),
                            plot_img_dim     = dim(internal$slice),
                            coordinates      = internal$coord,
                            pathology_values = internal$signal,
                            pathology_names  = internal$meta_names,
                            color_scheme     = internal$meta_colors,
                            file_name        = get_download_filename(text = c(internal$sample, internal$feature)),
                            title            = "",
                            plot_size        = "xsmall")
                    } else {
                        color_key <- NULL
                        breaks    <- NULL
                        vals      <- internal$signal[[internal$feature]]

                        # setup the color spectrum breaks or color key to be the same as what is used in Pathology tab
                        if (is.numeric(vals)) {
                            breaks <- get_color_scale_breaks(signal      = vals,
                                                             scaling     = "native",
                                                             color_count = length(COLOR_SCHEME_SPECTRUM_SCALING))
                        } else {
                            color_scheme <- build_plot_color_scheme(color_string = "",
                                                                    data_groups  = unique(vals))
                            color_key    <- color_scheme$colors
                        }

                        plot <- create_metadata_plot_cx(
                            plot_img_id           = session$ns(plot_img_id),
                            plot_img_dim          = dim(internal$slice),
                            coordinates           = internal$coord,
                            meta_data             = internal$signal,
                            file_name             = get_download_filename(text = c(internal$sample, internal$feature)),
                            plot_size             = "xsmall",
                            color_key             = color_key,
                            color_spectrum_breaks = breaks)
                    }
                }

                if (is.null(plot)) {
                    tagList(tags$div(style = "text-align:center;",
                                     tags$h6(class = "text-danger",
                                             "Plot Unavailable")))
                } else {
                    output[[plot_img_id]] <- renderImage({
                        list(src         = bkg_img_file,
                             contentType = "image/png",
                             style       = "display:none;")
                    },
                    deleteFile = FALSE)
                    output[[plot_id]] <- renderCanvasXpress(plot)

                    aspect_dim <- c(plot$x$config$setMaxY - plot$x$config$setMinY,
                                    plot$x$config$setMaxX - plot$x$config$setMinX)

                    # NOTE: no adjustment is made for layout of extremely unbalanced images in the matrix

                    # fixed width layout for matrix
                    dim_width <- 165
                    list(imageOutput(session$ns(plot_img_id), height = "0px"),
                         canvasXpressOutput(session$ns(plot_id),
                                            height = glue("{dim_width * (aspect_dim[1]/aspect_dim[2])}px"),
                                            width  = glue("{dim_width}px")))
                }
            })
        })
}
