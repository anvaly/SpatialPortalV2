#' stSingleBoxServer
#'   Server functionality for spatial transcriptomics visualization box module
#'
#' @param id                    - box id
#' @param genes_choices         - list of genes for each available organism or NULL
#' @param signature_choices     - character vector of all gene signature names
#' @param signatures            - data.frame of genes for each gene signature
#' @param samples               - data.frame of metadata for all available samples
#' @param st_objects            - list of st_data_class objects for all available samples reactive variable
#' @param custom_signatures     - a list of custom signatures that have been by the user added so far
#' @param selected_dataset      - name of initial dataset selection
#' @param selected_signature    - name of initial gene signature selection
#' @param selected_sample_st_id - ST_ID of initial sample selection
#' @param selected_singleGene   - name of initial single gene selection
#' @param selected_custom_sig   - name of custom signatures
#' @param selected_coloring     - selected coloring from matrix tab if any
#' @param plot_on_create        - should the plot be created immediately (default FALSE)
#' @param logger                - reactive logger S4 object
#'
#' @return reactive expression containing a logical - TRUE if the UI box is open, FALSE if it has been closed
stSingleBoxServer <- function(id,
                              genes_choices,
                              signature_choices,
                              signatures,
                              samples,
                              st_objects,
                              custom_signatures,
                              selected_dataset      = character(0),
                              selected_signature    = character(0),
                              selected_sample_st_id = character(0),
                              selected_singleGene   = character(0),
                              selected_custom_sig   = character(0),
                              selected_coloring     = character(0),
                              plot_on_create        = FALSE,
                              logger) {
    moduleServer(id, function(input, output, session) {
        internal <- reactiveValues()
            internal$id             <- id
            internal$valid          <- TRUE
            internal$datasetSel     <- selected_dataset
            internal$signatureSel   <- selected_signature
            internal$sampleSel      <- selected_sample_st_id
            internal$sampleAppName  <- NULL
            internal$singleGeneSel  <- selected_singleGene
            internal$customSigSel   <- selected_custom_sig
            internal$plotScaling    <- "trim"
            internal$bkg_img_file   <- tempfile(fileext = CX_BKG_IMG_FILE_EXT)
            internal$plot_on_create <- plot_on_create

        update_dataset_selectized_input(session, samples, internal$datasetSel)
        update_samples_selectized_input(session, samples, internal$datasetSel, internal$sampleSel)
        update_signatures_selectized_input(session, signatures, signature_choices, internal$signatureSel)
        update_customSigSel_selectized_input(session, custom_signatures, internal$customSigSel)

        if (length(selected_coloring) > 0) {
            updateRadioButtons(session  = session,
                               inputId  = "plotScalingRadio",
                               selected = selected_coloring)
        }

        output$boxTitle <- renderText({
            "Create a Spatial Plot"
        })

        output$signatureGeneLbl <- renderText({
            "Signature/Gene/Custom: "
        })

        output$organismTooltip <- renderUI({
            ui_tooltip(session$ns("organismTooltip"),
                       "Dataset",
                       "Choose from the available Organisms and Tissues")
        })

        output$plotScalingRadioTT <- renderUI({
            tags$span(class = "periscope-input-label-with-tt",
                      "Color Scale",
                      shiny::img(id     = session$ns("plotScalingRadioTT"),
                                 src    = isolate(periscope:::.g_opts$tt_image),
                                 height = isolate(periscope:::.g_opts$tt_height),
                                 width  = isolate(periscope:::.g_opts$tt_width)),
                      bsTooltip(id        = session$ns("plotScalingRadioTT"),
                                placement = "left",
                                options   = list("data-html" = "true"),
                                glue('Plots display pre-calculated normalized expression values (SCT). ',
                                     'The options modify the range covered by the colormap:<br/><br/>',
                                     '<b>Native</b>: Full data range &nbsp;&nbsp;&nbsp;<br/>',
                                     '<b>Trimmed</b>: 1st to 99th percentile<br/>',
                                     '<b>Fixed</b>: Fixed range of [0-1]<br/>',
                                     '<b>Custom</b>: User specified min/max')))
        })

        output$sampleTooltip <- renderUI({
            ui_tooltip(session$ns("sampleTooltip"),
                       "Sample",
                       "Choose from the available samples")
        })

        output$customSignatureTooltip <- renderUI({
            ui_tooltip(session$ns("customSignatureTooltip"),
                       "Custom",
                       "Add or Choose a Custom Gene Signature<br/><i>Create a new signature by typing a list of symbols separated by spaces, commas or semicolons and pressing Enter</i><br/>Note: Custom values are available in this and new panels for <b>this app session only</b>.")
        })

        output$signatureTooltip <- renderUI({
            ui_tooltip(session$ns("signatureTooltip"),
                       "Gene Signature",
                       "Choose from custom gene signatures.<br/><i>Search for signatures containing a gene of interest by typing that gene\\'s name.</i>")
        })

        output$geneTooltip <- renderUI({
            ui_tooltip(session$ns("geneTooltip"),
                       "Single Gene",
                       "Choose from genes available in the sample.<br/><i>Search for genes of interest by typing that gene\\'s name.</i>")
        })

        observeEvent(input$closeBoxBtn, {
            internal$valid <- FALSE
        }, ignoreInit = TRUE)

        observeEvent(input$viewToggleBtn_std, {
            if (!is.null(input$viewToggleBtn_std) && !is.null(input$viewTypeTabs)) {
                switch_view_type(session, standard = input$viewToggleBtn_std)
            }
        }, ignoreInit = TRUE)

        observeEvent(input$viewToggleBtn_adv, {
            if (!is.null(input$viewToggleBtn_adv) && !is.null(input$viewTypeTabs)) {
                switch_view_type(session, standard = !(input$viewToggleBtn_adv))
            }
        }, ignoreInit = TRUE)

        observeEvent(input$datasetSel, {
            sample     <- character(0)
            signature  <- character(0)
            singleGene <- character(0)
            customSig  <- character(0)

            if (all(!is.null(internal$datasetSel),
                    length(internal$datasetSel) != 0,
                    internal$datasetSel != "",
                    input$datasetSel == internal$datasetSel)) {
                # this is a cloned box being created, use the preset values
                sample    <- internal$sampleSel
                signature <- internal$signatureSel
                customSig <- internal$customSigSel
            }

            clear_chart_box_inputs_feedback()

            # get selected organism samples
            update_samples_selectized_input(session, samples, input$datasetSel, sample)

            # clear signature input if it has data
            if (!is_valid_selection(input$datasetSel) && is_valid_selection(input$signatureSel)) {
                update_signatures_selectized_input(session, signatures, signature_choices, signature)
            }

            # clear custom signature input if it has data
            if (all(!is_valid_selection(internal$datasetSel),
                    !is_valid_selection(input$datasetSel),
                    is_valid_selection(input$customSigSel))) {
                reset("customSigSel")
            }
        })

        observeEvent(c(input$signatureSel, input$singleGeneSel, input$customSigSel), {
            clear_chart_box_inputs_feedback()
        })

        observeEvent(input$sampleSel, {
            clear_chart_box_inputs_feedback()

            sample_reference <- samples %>%
                filter(ST_ID == input$sampleSel) %>%
                pull(Reference)

            s_gene               <- internal$singleGeneSel
            sample_genes_choices <- get_sample_genes(dataset_sel        = input$datasetSel,
                                                     sample_sel         = input$sampleSel,
                                                     st_objects         = st_objects,
                                                     initial_genes_list = genes_choices,
                                                     sample_reference   = sample_reference)

            if (((length(internal$datasetSel) > 0) &&
                 (internal$datasetSel != "") &&
                 (input$datasetSel != internal$datasetSel)) ||
                (length(s_gene) == 0)                       ||
                !(s_gene %in% sample_genes_choices)) {
                # clear invalid choice
                s_gene <- character(0)
            }

            updateSelectizeInput(session,
                                 "singleGeneSel",
                                 server   = TRUE,
                                 selected = s_gene,
                                 choices  = sample_genes_choices,
                                 options  = list(
                                     placeholder = PLACE_HOLDER_SELECTIZE_INPUTS,
                                     render      = I("{ option: function(item, escape) {
                                           return '<div><b>' + item.value + '</b></div>'; }}")))

        })

        observeEvent(input$plotScalingRadio, {
            ns <- NS(id)
            if (is_valid_selection(input$plotScalingRadio) &&
                (internal$plotScaling != input$plotScalingRadio)) {
                if (input$plotScalingRadio == "custom") {
                    shinyjs::show(ns('customScalingLayout'), asis = T)
                }
                else {
                    shinyjs::hide(ns('customScalingLayout'), asis = T)
                    internal$plotScaling <- input$plotScalingRadio
                }
            }
        })

        observeEvent(input$customScalingBtn, {
            if (is_valid_selection(input$plotScalingRadio) &&
                (input$plotScalingRadio == "custom")) {
                feedbackDanger("customScaling1",
                               (!(is.numeric(type.convert(input$customScaling1, as.is = TRUE))) ||
                                    as.double(input$customScaling1) >= as.double(input$customScaling2)),
                               text = NULL)

                feedbackDanger("customScaling2",
                               (!(is.numeric(type.convert(input$customScaling2, as.is = TRUE))) ||
                                    as.double(input$customScaling2) <= as.double(input$customScaling1)),
                               text = NULL)

                if (is.numeric(type.convert(input$customScaling1, as.is = TRUE)) &&
                    is.numeric(type.convert(input$customScaling2, as.is = TRUE)) &&
                    (as.double(input$customScaling2) > as.double(input$customScaling1))) {
                    internal$plotScaling <- list(input$customScaling1, input$customScaling2)
                }
            }
        })

        observeEvent(c(input$getPlotBtn, internal$plotScaling, internal$plot_on_create), {
            closeAlert(session, "plotsAlert")

            if (!is.null(internal$plot_on_create) && internal$plot_on_create) {
                # submit on initialization for matrix -> single creation of boxes
                in_dataset    <- internal$datasetSel
                in_sample     <- internal$sampleSel
                in_signature  <- internal$signatureSel
                in_singleGene <- internal$singleGeneSel
                in_customSig  <- internal$customSigSel

                isolate({
                    # remove the flag that indicates this should be processed on create
                    internal$plot_on_create <- FALSE
                })
            } else {
                # standard processing of plot button/scaling
                in_dataset    <- input$datasetSel
                in_sample     <- input$sampleSel
                in_signature  <- input$signatureSel
                in_singleGene <- input$singleGeneSel
                in_customSig  <- input$customSigSel
            }

            feedbackDanger("datasetSel",
                           !is_valid_selection(in_dataset),
                           get_message_text("dataset_required"))
            req(in_dataset)

            feedbackDanger("sampleSel",
                           !is_valid_selection(in_sample),
                           get_message_text("sample_required"))
            req(in_sample)
            internal$datasetSel    <- in_dataset
            internal$sampleSel     <- in_sample
            internal$sampleAppName <- samples %>%
                filter(ST_ID == in_sample) %>%
                pull(Sample_App_Name)

            only_signatures_selected  <- all(is_valid_selection(in_signature),
                                             !is_valid_selection(in_singleGene),
                                             !is_valid_selection(in_customSig))

            only_single_gene_selected <- all(!is_valid_selection(in_signature),
                                             is_valid_selection(in_singleGene),
                                             !is_valid_selection(in_customSig))

            only_custom_sig_selected  <- all(!is_valid_selection(in_signature),
                                             !is_valid_selection(in_singleGene),
                                             is_valid_selection(in_customSig))

            is_only_one_input         <- any(only_signatures_selected,
                                             only_single_gene_selected,
                                             only_custom_sig_selected)

            if (!is_only_one_input) {
                shinyjs::show(id = "selectionFeedbackDiv", anim = TRUE, animType = "fade")
            }

            feedbackDanger("signatureSel",
                           !is_only_one_input, "")
            feedbackDanger("singleGeneSel",
                           !is_only_one_input, "")
            feedbackDanger("customSigSel",
                           !is_only_one_input, "")

            req(is_only_one_input)

            gene_title <- ""

            if (only_signatures_selected) {
                internal$signatureSel  <- in_signature
                internal$singleGeneSel <- NULL
                internal$customSigSel  <- NULL
                gene_title             <- internal$signatureSel
                selected_signature     <- signatures[(signatures$Name == internal$signatureSel), , drop = F]
                signature_gene_lbl     <- "Signature: "
            } else if (only_single_gene_selected) {
                internal$singleGeneSel <- in_singleGene
                internal$signatureSel <- NULL
                internal$customSigSel <- NULL
                gene_title            <- internal$singleGeneSel
                signature_gene_lbl    <- "Gene: "
            } else {
                custom_signature <- str_split_1(in_customSig, " |,|;")
                feedbackDanger("customSigSel",
                               length(custom_signature) == 0,
                               get_message_text("enter_one_gene"))
                req(length(custom_signature) > 0)

                internal$singleGeneSel <- NULL
                internal$signatureSel  <- NULL
                internal$customSigSel  <- custom_signature
                gene_title             <- glue_collapse(custom_signature, ", ")
                signature_gene_lbl     <- "Custom: "
            }

            if (only_custom_sig_selected) {
                custom_signature  <- toupper(custom_signature)
            }

            output$signatureGeneLbl <- renderText({
                signature_gene_lbl
            })

            output$boxTitle <- renderText({
                glue("{internal$sampleAppName}: {gene_title}")
            })

            output$infoSample <- renderText({
                internal$sampleAppName
            })

            output$infoSignature <- renderText({
                gene_title
            })

            selected_sample         <- get_full_sample_data_objects(internal$sampleSel, st_objects)[[1]]
            spatial_plot_id         <- UUIDgenerate()
            spatial_plot_img_id     <- glue("{spatial_plot_id}_{SUFFIX_SPATIAL_PLOT_ID}")
            spatial_plot_adv_img_id <- glue("{spatial_plot_id}_{SUFFIX_ADVANCED_SPATIAL_PLOT_ID}")
            features                <- list()

            # gene, signature and custom signature cannot have values at the same time
            if (only_signatures_selected) {
                features <- list(selected_signature$Genes[[1]])
            } else if (only_single_gene_selected) {
                features <- list(str_trim(internal$singleGeneSel))
            } else {
                features <- list(custom_signature)
            }

            scores            <- NULL
            coords            <- NULL
            img_dim           <- c(600, 600)
            overall_top_genes <- ""
            cluster_data      <- NULL
            cluster_names     <- NULL
            c_color_scheme    <- NULL
            pathology_data    <- NULL
            pathology_names   <- NULL
            p_color_scheme    <- NULL
            missing_item      <- NULL # if any required data item is missing we don't try to create plots
            valid_genes       <- character(0)

            if (selected_sample@is_api_based) {
                # API-Based App
                tryCatch({
                    if (is.null(selected_sample@api_data$features)) {
                        missing_item <- "features list"
                    } else if (is.null(selected_sample@api_data$cells)) {
                        missing_item <- "cells"
                    } else if (is.null(selected_sample@api_data$coordinates)) {
                        missing_item <- "coordinates"
                    } else if (is.null(selected_sample@slice_image)) {
                        missing_item <- "slice image"
                    } else if (is.null(selected_sample@top_genes)) {
                        missing_item <- "top genes"
                    } else {
                        valid_genes <- intersect(unlist(features), unique(selected_sample@api_data$features$gene_symbol))
                    }

                    if (is.null(missing_item) && (length(valid_genes) > 0)) {
                        img_dim <- dim(selected_sample@slice_image)
                        write_png_bkg_image(raw_img = selected_sample@slice_image, target_file = internal$bkg_img_file)

                        overall_top_genes <- get_overall_top_genes(valid_genes, selected_sample@top_genes)

                        coords <- selected_sample@api_data$coordinates

                        scores <- get_api_ucell_score(st_id    = selected_sample@st_id,
                                                      features = list(valid_genes),
                                                      st_objects)
                        if (is.null(scores)) {
                            # if the returned data is null even though we only requested valid genes for this sample,
                            # then the API retrieval failed
                            missing_item <- "scores"
                        } else {
                            genes_data <- get_api_feature_data(selected_sample, valid_genes)

                            if (is.null(genes_data)) {
                                # if the returned data is null even though we only requested valid genes for this sample,
                                # then the API retrieval failed
                                missing_item <- "feature data"
                            }
                        }

                        if (is.null(missing_item)) {
                            cluster_data  <- get_api_cell_tag_column(selected_sample, "cluster")
                            cluster_names <- get_api_user_defined_names(
                                sample    = selected_sample,
                                slot_name = get_env_value("clustering_name_field"))
                            c_color_scheme <- get_api_color_scheme(
                                sample      = selected_sample,
                                data_groups = unique(cluster_data$cluster),
                                color_slot  = get_env_value("clustering_color_field"))

                            pathology_data  <- get_api_cell_tag_column(selected_sample, "pathology")
                            pathology_names <- get_api_user_defined_names(
                                sample     = selected_sample,
                                slot_name  = get_env_value("pathology_name_field"),
                                lower_case = FALSE)
                            p_color_scheme <- get_api_color_scheme(
                                sample      = selected_sample,
                                data_groups = unique(pathology_data$pathology),
                                color_slot  = get_env_value("pathology_color_field"))
                        }
                    }
                },
                warning = function(w) {
                    warning(w)
                },
                error = function(e) {
                    warning(e)
                })
            } else {
                # Seurat-Based App
                tryCatch({
                    seurat        <- value(selected_sample@file_f_seurat)
                    sig_gene_name <- NULL

                    # Setup seurat score names
                    if (only_signatures_selected) {
                        sig_gene_name <- glue("sig.{selected_signature$Name}")
                    } else if (only_single_gene_selected) {
                        # make.names is used here because Seurat changes non-compliant names when
                        # adding on the module scores, however holds the gene names correctly on the
                        # object so in order to retrieve genes we have to not munge the gene names
                        sig_gene_name <- glue("sig.{make.names(str_trim(internal$singleGeneSel))}")
                    } else {
                        sig_gene_name <- "sig.custom"
                    }

                    single_gene_valid <- TRUE
                    sample_reference  <- samples %>%
                        filter(ST_ID == internal$sampleSel) %>%
                        pull(Reference)

                    if (!is.na(sample_reference) &&
                        (sample_reference %in% names(genes_choices)) &&
                        ((length(in_singleGene) == 1)) &&
                        (!(in_singleGene %in% genes_choices[[sample_reference]]))) {
                        single_gene_valid <- FALSE
                    }

                    if (only_single_gene_selected && !single_gene_valid) {
                        logwarn(glue("Gene: {internal$singleGeneSel} does not exist in sample: {internal$sampleAppName} genes list"))
                    } else {
                        if (only_single_gene_selected) {
                            overall_top_genes <- get_overall_top_genes(features, selected_sample@top_genes)
                            scores            <- get_seurat_data(seurat = seurat,
                                                                 items  = unlist(features),
                                                                 assay  = SEURAT_ASSAY)
                            valid_genes <- unlist(features)
                        } else {
                            seurat.scored <- add_ucell_module_score(seurat   = seurat,
                                                                    features = features,
                                                                    name     = sig_gene_name)
                            invalid_genes <- setdiff(unlist(features), rownames(seurat.scored))
                            valid_genes   <- setdiff(unlist(features), invalid_genes)
                            if (length(valid_genes) > 0) {
                                overall_top_genes <- get_overall_top_genes(valid_genes, selected_sample@top_genes)
                                scores     <- get_seurat_data(seurat = seurat.scored,
                                                              items  = glue("{sig_gene_name}"),
                                                              slot   = SEURAT_DATA_SLOT_NAME)
                            }
                        }

                        img_dim <- dim(selected_sample@slice_image)
                        write_png_bkg_image(raw_img = selected_sample@slice_image, target_file = internal$bkg_img_file)

                        coords     <- get_seurat_coordinates(seurat)
                        genes_data <- get_seurat_data(seurat = seurat,
                                                      items  = valid_genes,
                                                      slot   = SEURAT_DATA_SLOT_NAME)

                        smp_sel <- samples %>%
                            filter(ST_ID == in_sample)

                        cluster_data   <- get_seurat_clusters(seurat = seurat, cluster_field = smp_sel$CLUSTER_FIELD)
                        cluster_names  <- build_plot_custom_names(naming_string = smp_sel[[get_env_value("clustering_name_field")]],
                                                                  lower_case    = TRUE) # set to lower case on purpose for clusters
                        c_color_scheme <- build_plot_color_scheme(color_string  = smp_sel[[get_env_value("clustering_color_field")]],
                                                                  data_groups   = unique(cluster_data$cluster))

                        pathology_data  <- get_seurat_pathology(seurat = seurat, pathology_field = smp_sel$PATHOLOGY_FIELD)
                        pathology_names <- build_plot_custom_names(naming_string = smp_sel[[get_env_value("pathology_name_field")]],
                                                                   lower_case    = FALSE)
                        p_color_scheme  <- build_plot_color_scheme(color_string  = smp_sel[[get_env_value("pathology_color_field")]],
                                                                  data_groups    = unique(pathology_data$pathology))
                    }
                },
                warning = function(w) {
                    warning(w)
                },
                error = function(e) {
                    warning(e)
                })
            }

            if (!is.null(missing_item)) {
                logwarn(glue("Unable to retrieve {missing_item} for '{in_dataset}' sample '{internal$sampleAppName}"))
                createAlert(session,
                            "bodyAlert",
                            "plotsAlert",
                            style   = "warning",
                            content = get_message_text("data_item_missing",
                                                       message_parameters = list(item    = missing_item,
                                                                                 dataset = in_dataset,
                                                                                 sample  = internal$sampleAppName)),
                            append  = FALSE)
            }

            output$infoGenes <- renderUI({
                result <- HTML(overall_top_genes)
                if (is.null(overall_top_genes)) {
                    result <- tags$p(style = "text-align: center; color: #A94442; font-style: italic;",
                                     get_message_text("top_genes_unexpected_error"))
                }
                result
            })

            output$spatialPlot <- renderUI({
                spatial_plot <- NULL

                if (!is.null(missing_item)) {
                    sp_msg <- get_message_text("unexpected_data_retrieval_error")
                } else if (is.null(scores)) {
                    sp_msg <- get_message_text("genes_not_represented")
                } else {
                    spatial_plot <- create_spatial_plot_cx(
                        plot_img_id   = session$ns(spatial_plot_img_id),
                        plot_img_dim  = img_dim,
                        coordinates   = coords,
                        signal        = scores,
                        features      = unlist(features),
                        title         = internal$sampleAppName,
                        subtitle      = gene_title,
                        top_genes     = selected_sample@top_genes,
                        scaling       = internal$plotScaling,
                        file_name     = get_download_filename(text = c(internal$sampleAppName, gene_title)),
                        plot_size     = "large")

                    if (is.null(spatial_plot)) {
                        # the plot is empty due to dynamic range limitations
                        sp_msg <- get_message_text("small_dynamic_range")
                    }
                }

                if (is.null(spatial_plot)) {
                    shinyjs::hide('plotScalingRadio')
                    shinyjs::hide('customScalingLayout')

                    tags$div(style = "text-align:center;
                                  white-space:nowrap;
                                  margin-left:150px;
                                  margin-top:150px;",
                             tags$h4(class = "text-danger",
                                     get_message_text("unable_create_spatial_overview")),
                             tags$small(class = "text-muted",
                                        tags$em(sp_msg,
                                                tags$br(),
                                                get_message_text("contact_app_author"))))
                } else {
                    output[[spatial_plot_img_id]] <- renderImage({
                        list(src         = internal$bkg_img_file,
                             contentType = 'image/png',
                             style       = "display:none;")
                    }, deleteFile = FALSE)

                    output[[spatial_plot_id]] <- renderCanvasXpress(spatial_plot)
                    shinyjs::show('plotScalingRadio')
                    if (is_valid_selection(input$plotScalingRadio) &&
                        (input$plotScalingRadio == "custom")) {
                        shinyjs::show('customScalingLayout')
                    }

                    aspect_dim <- c(spatial_plot$x$config$setMaxY - spatial_plot$x$config$setMinY,
                                    spatial_plot$x$config$setMaxX - spatial_plot$x$config$setMinX)

                    # handle extremely unbalanced images
                    dim_factor <- 1
                    if ((aspect_dim[2]/aspect_dim[1]) > 1.75) {
                        dim_factor <- 0.66 #two-thirds of the fixed aspect
                    }

                    # fixed height layout for overview in plot box
                    dim_height       <- 550
                    legend_allowance <- 100
                    list(
                        imageOutput(session$ns(spatial_plot_img_id), height = "0px"),
                        canvasXpressOutput(session$ns(spatial_plot_id),
                                           height = glue("{dim_height * dim_factor}px"),
                                           width  = glue("{(dim_height + legend_allowance) * dim_factor * (aspect_dim[2]/aspect_dim[1])}px"))
                    )
                }
            })

            output$spatialPlotAdvanced <- renderUI({
                spatial_plot   <- NULL
                cluster_plot   <- NULL
                tissue_plot    <- NULL
                pathology_plot <- NULL
                dot_box_plot   <- NULL
                sp_msg         <- NULL

                if (!is.null(missing_item)) {
                    sp_msg <- get_message_text("unexpected_data_retrieval_error")
                } else if (!is.null(scores)) {
                    spatial_plot <- create_spatial_plot_cx(
                        plot_img_id   = session$ns(spatial_plot_adv_img_id),
                        plot_img_dim  = img_dim,
                        coordinates   = coords,
                        signal        = scores,
                        features      = unlist(features),
                        title         = "Spatial",
                        subtitle      = NULL,
                        top_genes     = selected_sample@top_genes,
                        scaling       = internal$plotScaling,
                        file_name     = get_download_filename(text = c(internal$sampleAppName, gene_title)),
                        plot_size     = "small")

                    if (is.null(spatial_plot)) {
                        # the plot is empty due to dynamic range limitations
                        sp_msg <- get_message_text("small_dynamic_range")
                    }
                }

                if (!is.null(spatial_plot)) {
                    cluster_plot <- create_cluster_plot_cx(
                        plot_img_id    = session$ns(spatial_plot_adv_img_id),
                        plot_img_dim   = img_dim,
                        coordinates    = coords,
                        cluster_values = cluster_data,
                        cluster_names  = cluster_names,
                        color_scheme   = c_color_scheme,
                        file_name      = get_download_filename(text = c(internal$sampleAppName, METADATA_ITEM_CLUSTERING)),
                        title          = METADATA_ITEM_CLUSTERING,
                        plot_size      = "small")

                    pathology_plot  <- create_pathology_plot_cx(
                        plot_img_id      = session$ns(spatial_plot_adv_img_id),
                        plot_img_dim     = img_dim,
                        coordinates      = coords,
                        pathology_values = pathology_data,
                        pathology_names  = pathology_names,
                        color_scheme     = p_color_scheme,
                        file_name        = get_download_filename(text = c(internal$sampleAppName, METADATA_ITEM_PATHOLOGY)),
                        title            = METADATA_ITEM_PATHOLOGY,
                        plot_size        = "small")

                    tissue_plot  <- create_tissue_plot_cx(
                        plot_img_id  = session$ns(spatial_plot_adv_img_id),
                        plot_img_dim = img_dim,
                        coordinates  = coords,
                        file_name    = get_download_filename(text = c(internal$sampleAppName, "Tissue")),
                        plot_size    = "small")

                    if (only_custom_sig_selected) {
                        gene_title <- glue_collapse(valid_genes, ", ")
                    }

                    if (length(valid_genes) > 1) {
                        dot_box_plot <- create_dot_plot_cx(
                            genes_data     = genes_data,
                            cluster_values = cluster_data,
                            cluster_names  = cluster_names,
                            color_scheme   = c_color_scheme,
                            features       = valid_genes, # list of genes (used individually)
                            title          = gene_title,
                            file_name      = get_download_filename(text = c(internal$sampleAppName, gene_title)))
                    } else if (length(valid_genes) == 1) {
                        dot_box_plot <- create_box_plot_cx(
                            genes_data     = genes_data,
                            cluster_values = cluster_data,
                            cluster_names  = cluster_names,
                            color_scheme   = c_color_scheme,
                            title          = gene_title,
                            file_name      = get_download_filename(text = c(internal$sampleAppName, gene_title)))
                    }
                }

                if (any(is.null(spatial_plot),
                        is.null(cluster_plot),
                        is.null(dot_box_plot),
                        is.null(pathology_plot),
                        is.null(tissue_plot))) {

                    if (is.null(sp_msg)) {
                        sp_msg <- get_message_text("genes_not_represented_info_unavailable")
                    }

                    tags$div(style = "text-align:center;
                                  margin-top:100px;",
                             tags$h4(class = "text-danger",
                                     style = "text-align:center;",
                                     get_message_text("unable_create_spatial_extended")),
                             tags$small(class = "text-muted",
                                        tags$em(sp_msg,
                                                tags$br(),
                                                get_message_text("contact_app_author"))))
                } else {
                    aspect_dim <- c(spatial_plot$x$config$setMaxY - spatial_plot$x$config$setMinY,
                                    spatial_plot$x$config$setMaxX - spatial_plot$x$config$setMinX)

                    # handle extremely unbalanced images
                    dim_factor <- 1
                    if ((aspect_dim[2]/aspect_dim[1]) > 1.75) {
                        dim_factor <- 0.66 #two-thirds of the fixed aspect
                    }

                    # fixed height layout for advanced view
                    dim_height       <- 225
                    legend_allowance <- 75
                    height           <- glue("{dim_height * dim_factor}px")
                    width            <- glue("{(dim_height + legend_allowance) * dim_factor * (aspect_dim[2]/aspect_dim[1])}px")

                    # update the background image for the spatial plot to be this tab's content
                    adv_spatial_plot <- spatial_plot
                    adv_spatial_plot$x$config$backgroundImage <- glue("javascript://{session$ns(spatial_plot_adv_img_id)}")

                    output[[spatial_plot_adv_img_id]] <- renderImage({
                        list(src         = internal$bkg_img_file,
                             contentType = 'image/png',
                             style       = "display:none;")
                    }, deleteFile = FALSE)

                    adv_spatial_id <- glue("{spatial_plot_id}_adv")
                    output[[adv_spatial_id]] <- renderCanvasXpress(adv_spatial_plot)

                    adv_cluster_id <- glue("{spatial_plot_id}_clust")
                    output[[adv_cluster_id]] <- renderCanvasXpress(cluster_plot)

                    adv_tissue_id <- glue("{spatial_plot_id}_tissue")
                    output[[adv_tissue_id]] <- renderCanvasXpress(tissue_plot)

                    adv_pathology_id  <- glue("{spatial_plot_id}_pathology")
                    output[[adv_pathology_id]] <- renderCanvasXpress(pathology_plot)

                    adv_dotbox_id <- glue("{spatial_plot_id}_dot_box")
                    # render the dot plot differently depending if a canvasXpress object or error message is returned
                    if ("shiny.tag" %in% class(dot_box_plot)) {
                        output[[adv_dotbox_id]] <- renderUI(dot_box_plot)
                        dot_box_output          <- uiOutput(session$ns(adv_dotbox_id))
                    } else {
                        output[[adv_dotbox_id]] <- renderCanvasXpress(dot_box_plot)
                        dot_box_output          <- canvasXpressOutput(session$ns(adv_dotbox_id),
                                                                       height = "300px",
                                                                       width  = "1250px")
                    }


                    tagList(
                        imageOutput(session$ns(spatial_plot_adv_img_id), height = "0px"),
                        tags$div(style = glue("display:flex; position:relative; height:{250 * dim_factor}px;"),
                                 tags$div(style = "display:flex; position:absolute; justify-content:space-between; align-items:flex-start;
                                               top:0; left:0; right:0; bottom:0; overflow-x:auto; overflow-y:hidden;",
                                          tags$div(style = "flex:none;",
                                                   canvasXpressOutput(session$ns(adv_spatial_id),
                                                                      height = height,
                                                                      width  = width)),
                                          tags$div(style = "flex:none;",
                                                   canvasXpressOutput(session$ns(adv_cluster_id),
                                                                      height = height,
                                                                      width  = width)),
                                          tags$div(style = "flex:none;",
                                                   canvasXpressOutput(session$ns(adv_pathology_id),
                                                                      height = height,
                                                                      width  = width)),
                                          tags$div(style = "flex:none;",
                                                   canvasXpressOutput(session$ns(adv_tissue_id),
                                                                      height = height,
                                                                      width  = width))
                                 )
                        ),
                        tags$div(style = "display:flex; position:relative; height:330px; justify-content:space-around; align-items:center;",
                                 tags$div(style = "display:flex; top:0, left:0; right:0, bottom:0, overflow-x:auto; overflow-y:hidden;",
                                          tags$div(style = "flex:none;",
                                                   dot_box_output)
                                 )
                        )
                    )
                }
            })
        }) #plot button

        return(internal$valid)
    }) #moduleServer
}
