#' get_ui_addl_metadata_table
#'     Helper function to format the headings and additional metadata data tables consistently
#'
#' @param samples_data     - dataframe containing all the sample data
#' @param metadata_cols    - metadata column data frame or NULL if not found
#' @param smp_id           - character sample identifier
#' @param type             - type of table to show (ie. Secondary, Technical, Other, ...)
#'
#' @return shiny.tag or NULL
get_ui_addl_metadata_table <- function(samples_data, metadata_cols, smp_id, type) {
    result         <- NULL
    metadata_table <- data.frame()

    if (identical(tolower(get_env_value("data_mode")), "file")) {
        metadata_table <- get_file_metadata_table(samples_data  = samples_data,
                                                  metadata_cols = metadata_cols,
                                                  smp_id        = smp_id,
                                                  type          = type)
    } else {
        metadata_table <- get_api_metadata_table(samples_data  = samples_data,
                                                 metadata_cols = metadata_cols,
                                                 smp_id        = smp_id,
                                                 type          = type)
    }

    # format the table
    if (NROW(metadata_table) > 0) {
        formatted_tbl <- datatable(metadata_table,
                                   rownames  = FALSE,
                                   colnames  = "",
                                   selection = "none",
                                   options   = list(dom            = "t",
                                                    paging         = FALSE,
                                                    scrollY        = "75vh",
                                                    scrollCollapse = TRUE,
                                                    bSort          = FALSE,
                                                    autoWidth      = TRUE,
                                                    columnDefs     = list(list(width = "30%", targets = 0)))) %>%
            formatStyle("key", fontWeight = "bold")

        result <- tagList(renderUI(tags$div(
            style = "text-align:center; font-size:large; font-weight:bold",
            tags$br(),
            type)),
            renderDataTable(formatted_tbl))
    }

    result
}


#' get_ui_matrix_plot_table
#'   Build matrix html table along with ST plots
#'
#' @param selected_samples    - sample objects selected
#' @param selected_signatures - selected list of signatures
#' @param selected_genes      - selected list of genes
#' @param selected_metadata   - selected list of additional metadata to be displayed
#' @param coloring            - selected option for color scale range
#' @param coloring_scaling    - selected color scaling option
#' @param samples_signals     - signals data for the selected samples and signatures/genes
#' @param selectedMatrixCol   - reactive variable to track which column or row is selected
#' @param selectedMatrixRow   - reactive variable to track which column or row is selected
#' @param samples             - the cached sample table
#' @param logger              - reactive logger S4 object
#'
#' @return html table
get_ui_matrix_plot_table <- function(selected_samples,
                                     selected_signatures,
                                     selected_genes,
                                     selected_metadata,
                                     coloring,
                                     coloring_scaling,
                                     samples_signals,
                                     selectedMatrixCol,
                                     selectedMatrixRow,
                                     samples,
                                     logger) {
    samples_rows              <- tagList()
    table_header_row          <- tagList()
    header_columns            <- tagList(tags$td())
    gene_sig_coloring_signal  <- NULL

    sample_names_by_id <- setNames(sapply(selected_samples, function(x) { x@name }),
                                   nm = sapply(selected_samples, function(x) { x@st_id }))

    if ((coloring_scaling == "bySignatureGene") || (coloring_scaling == "allTogether")) {
        coloring_signals <- bind_rows(samples_signals)
    } else {
        coloring_signals <- samples_signals
    }

    if (NROW(selected_signatures) > 0) {
        for (signature_index in 1:NROW(selected_signatures)) {
            sample_uuid    <- glue("sample-{UUIDgenerate()}")
            signature_name <- selected_signatures$Name[signature_index]

            header_columns <- tagList(
                header_columns,
                tags$th(matrixGridtoSingleUI(
                    id = sample_uuid,
                    display_name = signature_name)))

            matrixGridtoSingleServer(id                 = sample_uuid,
                                     sample_names_by_id = sample_names_by_id,
                                     signatures         = signature_name,
                                     genes              = NULL,
                                     coloring           = coloring,
                                     selectedMatrix     = selectedMatrixCol)
        }
    }

    if (length(selected_genes) > 0) {
        for (gene_index in 1:length(selected_genes)) {
            sample_uuid <- glue("sample-{UUIDgenerate()}")
            gene        <- str_trim(selected_genes[gene_index])

            header_columns <- tagList(
                header_columns,
                tags$th(matrixGridtoSingleUI(
                    id = sample_uuid,
                    display_name = gene)))

            matrixGridtoSingleServer(id                 = sample_uuid,
                                     sample_names_by_id = sample_names_by_id,
                                     signatures         = NULL,
                                     genes              = gene,
                                     coloring           = coloring,
                                     selectedMatrix     = selectedMatrixCol)
        }
    }

    if (length(selected_metadata) > 0) {
        for (meta in selected_metadata) {
            header_columns <- tagList(
                header_columns,
                tags$th(meta))
        }
    }

    table_header_row <- tags$tr(header_columns)

    for (s_obj in selected_samples) {
        if (NROW(selected_signatures) > 0 || length(selected_genes) > 0) {
            sample_uuid  <- glue("sample-{UUIDgenerate()}")

            sample_row <- tagList(
                tags$th(matrixGridtoSingleUI(
                    id           = sample_uuid,
                    display_name = s_obj@name)))

            matrixGridtoSingleServer(id                 = sample_uuid,
                                     sample_names_by_id = setNames(s_obj@name, nm = s_obj@st_id),
                                     signatures         = selected_signatures$Name,
                                     genes              = selected_genes,
                                     coloring           = coloring,
                                     selectedMatrix     = selectedMatrixRow)
        } else {
            sample_row <- tagList(
                tags$th(s_obj@name))
        }

        if (NROW(selected_signatures) > 0) {
            for (signature_index in 1:NROW(selected_signatures)) {
                plot_id          <- glue("matrix-{UUIDgenerate()}")
                signature_name   <- selected_signatures$Name[signature_index]
                obj_signal       <- samples_signals[[s_obj@st_id]]

                if (s_obj@is_api_based) {
                    item_name <- glue("{signature_name}")
                    coord     <- s_obj@api_data$coordinates
                } else {
                    item_name <- glue("sig.{signature_name}")
                    coord     <- get_seurat_coordinates(value(s_obj@file_f_seurat))
                }

                if (!is.null(obj_signal) && (item_name %in% colnames(obj_signal))) {
                    obj_signal <- obj_signal[, item_name, drop = F]
                } else {
                    obj_signal <- NULL
                }

                if (coloring_scaling == "bySignatureGene") {
                    gene_sig_coloring_signal <- coloring_signals %>%
                        select(any_of(item_name))
                } else if (coloring_scaling == "bySample") {
                    gene_sig_coloring_signal <- coloring_signals[[s_obj@st_id]]
                } else if (coloring_scaling == "allTogether") {
                    gene_sig_coloring_signal <- coloring_signals
                }

                stMatrixPlotServer(
                    id                  = plot_id,
                    sample_name         = s_obj@name,
                    feature_name        = signature_name,
                    signature_genes     = unlist(selected_signatures$Genes[signature_index]),
                    slice_image         = s_obj@slice_image,
                    coordinates         = coord,
                    signal              = obj_signal,
                    metadata_type       = NULL,
                    meta_names          = NULL,
                    meta_colors         = NULL,
                    coloring            = coloring,
                    coloring_signal     = gene_sig_coloring_signal,
                    logger              = logger)

                sample_row <- tagList(
                    sample_row,
                    tags$td(stMatrixPlotUI(plot_id)))
            }
        }

        if (length(selected_genes) > 0) {
            for (gene_index in 1:length(selected_genes)) {
                plot_id    <- UUIDgenerate()
                gene       <- str_trim(selected_genes[gene_index])
                obj_signal <- samples_signals[[s_obj@st_id]]

                if (s_obj@is_api_based) {
                    item_name <- glue("{gene}")
                    coord     <- s_obj@api_data$coordinates
                } else {
                    item_name <- gene
                    coord     <- get_seurat_coordinates(value(s_obj@file_f_seurat))
                }

                if (!is.null(obj_signal) && (item_name %in% colnames(obj_signal))) {
                    obj_signal <- obj_signal[, item_name, drop = F]
                } else {
                    obj_signal <- NULL
                }

                if (coloring_scaling == "bySignatureGene") {
                    gene_sig_coloring_signal <- coloring_signals %>%
                        select(any_of(item_name))
                } else if (coloring_scaling == "bySample") {
                    gene_sig_coloring_signal <- coloring_signals[[s_obj@st_id]]
                } else if (coloring_scaling == "allTogether") {
                    gene_sig_coloring_signal <- coloring_signals
                }

                stMatrixPlotServer(
                    id                  = plot_id,
                    sample_name         = s_obj@name,
                    feature_name        = gene,
                    signature_genes     = NULL,
                    slice_image         = s_obj@slice_image,
                    coordinates         = coord,
                    signal              = obj_signal,
                    metadata_type       = NULL,
                    meta_names          = NULL,
                    meta_colors         = NULL,
                    coloring            = coloring,
                    coloring_signal     = gene_sig_coloring_signal,
                    logger              = logger)

                sample_row <- tagList(
                    sample_row,
                    tags$td(stMatrixPlotUI(plot_id)))
            }
        }

        if (length(selected_metadata) > 0) {
            for (meta in selected_metadata) {
                plot_id  <- UUIDgenerate()
                m_names  <- NULL
                m_colors <- NULL

                if (s_obj@is_api_based) {
                    coord    <- s_obj@api_data$coordinates

                    if (meta == METADATA_ITEM_CLUSTERING) {
                        data    <- get_api_cell_tag_column(s_obj, "cluster")
                        m_names <- get_api_user_defined_names(
                            sample    = s_obj,
                            slot_name = get_env_value("clustering_name_field"))
                        m_colors <- get_api_color_scheme(
                            sample      = s_obj,
                            data_groups = unique(data$cluster),
                            color_slot  = get_env_value("clustering_color_field"))
                    } else if (meta == METADATA_ITEM_PATHOLOGY) {
                        data    <- get_api_cell_tag_column(s_obj, "pathology")
                        m_names <- get_api_user_defined_names(
                            sample     = s_obj,
                            slot_name  = get_env_value("pathology_name_field"),
                            lower_case = FALSE)
                        m_colors <- get_api_color_scheme(
                            sample      = s_obj,
                            data_groups = unique(data$pathology),
                            color_slot  = get_env_value("pathology_color_field"))
                    } else {
                        data <- get_api_cell_tag_column(s_obj, meta)
                    }
                } else {
                    sample_seurat <- value(s_obj@file_f_seurat)
                    coord         <- get_seurat_coordinates(sample_seurat)
                    smp_sel       <- samples %>%
                        filter(ST_ID == s_obj@st_id)

                    if (meta == METADATA_ITEM_CLUSTERING) {
                        data     <- get_seurat_clusters(seurat = sample_seurat, cluster_field = smp_sel$CLUSTER_FIELD)
                        m_names  <- build_plot_custom_names(naming_string = smp_sel[[get_env_value("clustering_name_field")]],
                                                            lower_case    = TRUE) # set to lower case on purpose for clusters
                        m_colors <- build_plot_color_scheme(color_string  = smp_sel[[get_env_value("clustering_color_field")]],
                                                            data_groups   = unique(data$cluster))
                    } else if (meta == METADATA_ITEM_PATHOLOGY) {
                        data     <- get_seurat_pathology(seurat = sample_seurat, pathology_field = smp_sel$PATHOLOGY_FIELD)
                        m_names  <- build_plot_custom_names(naming_string = smp_sel[[get_env_value("pathology_name_field")]],
                                                            lower_case    = FALSE)
                        m_colors <- build_plot_color_scheme(color_string  = smp_sel[[get_env_value("pathology_color_field")]],
                                                            data_groups   = unique(data$pathology))
                    } else {
                        data <- get_seurat_data(seurat = sample_seurat,
                                                items  = meta,
                                                slot   = SEURAT_DATA_SLOT_NAME)
                    }
                }

                stMatrixPlotServer(
                    id                  = plot_id,
                    sample_name         = s_obj@name,
                    feature_name        = NULL,
                    signature_genes     = NULL,
                    slice_image         = s_obj@slice_image,
                    coordinates         = coord,
                    signal              = data,
                    metadata_type       = meta,
                    meta_names          = m_names,
                    meta_colors         = m_colors,
                    coloring            = NULL,
                    coloring_signal     = NULL,
                    logger              = logger)

                sample_row <- tagList(
                    sample_row,
                    tags$td(stMatrixPlotUI(plot_id)))
            }
        }

        samples_rows <- tagList(
            samples_rows,
            tags$tr(sample_row))
    }

    tagList(
        tags$table(class = "paddingBetweenCols",
                   table_header_row,
                   samples_rows),
        shinyjs::runjs(
            glue("dismissMatrixBusyModal({length(selected_samples) * (NROW(selected_signatures) + length(selected_genes) + length(selected_metadata))})")))
}


#' get_ui_matrix_table_button
#'     a button to trigger the modal popup with additional metadata fields
#'
#' @param button_id  - an id to identify this button
#' @param sample_ids - a unique identifier
#'
#' @return character string of a action button that can be added to a datatable
get_ui_matrix_table_button <- function(button_id, sample_ids) {
    matrix_btns <- character(length(sample_ids))

    for (i in seq_len(length(sample_ids))) {
        matrix_btns[i] <- glue("<div style = 'display:inline-flex; float:right;'>",
                               as.character(shiny::actionButton(
                                   inputId     = glue("extended_{button_id}_{sample_ids[i]}"),
                                   label       = NULL,
                                   onclick     = "Shiny.setInputValue('extended_button',  this.id, {priority: 'event'})",
                                   onmousedown = I("stopClickPropagation(event)"), # requires custom js function in ui_body
                                   icon        = icon("plus"),
                                   title       = "Extended Attributes")),
                               "</div>",
                               .open = "{{", .close = "}}")
    }

    matrix_btns
}


