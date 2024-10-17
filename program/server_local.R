# ----------------------------------------
# --       PROGRAM server_local.R       --
# ----------------------------------------
# USE: Session-specific variables and
#      functions for the main reactive
#      shiny server functionality.  All
#      code in this file will be put into
#      the framework inside the call to
#      shinyServer(function(input, output, session)
#      in server.R
#
# NOTEs:
#   - All variables/functions here are
#     SESSION scoped and are ONLY
#     available to a single session and
#     not to the UI
#
#   - For globally scoped session items
#     put var/fxns in server_global.R
#
# FRAMEWORK VARIABLES
#     input, output, session - Shiny
#     ss_userAction.Log - Reactive Logger S4 object
# ----------------------------------------

# -- IMPORTS --


# -- VARIABLES --
userSel <- reactiveValues()
    userSel$logged_in              <- FALSE
    userSel$singleBoxesIds         <- list()
    userSel$pathologyBoxesIds      <- list()
    userSel$custom_signatures      <- list()
    userSel$matrixSignature        <- NULL
    userSel$matrixSample           <- NULL
    userSel$matrixGene             <- NULL
    userSel$matrixColor            <- NULL
    userSel$matrixGridSelected     <- FALSE
    userSel$matrixConfigSamples    <- NULL
    userSel$matrixConfigSignatures <- NULL
    userSel$matrixConfigGenes      <- NULL
    userSel$matrixConfigMetadata   <- NULL
    userSel$matrixConfigColoring   <- NULL
    userSel$matrixColoringScaling  <- NULL

selectedMatrixCol <- reactiveVal(NULL)
selectedMatrixRow <- reactiveVal(NULL)
matrix_samples    <- reactiveVal()
user_st_objects   <- reactiveVal(list())


# -- MODULES --
                                                                                                         
matrix_sample_row_sel <- downloadableTable(
    id               = "matrixSampleTbl",
    logger           = ss_userAction.Log,
    downloaddatafxns = NULL,
    filenameroot     = "",
    tabledata        = matrix_samples,
    class            = "nowrap display",
    formatStyle      = list(columns    = c("Sample"),
                            fontWeight = "bold"),
    filter           = list(position = "top",
                            clear    = FALSE),
    columnDefs       = list(list(orderable  = FALSE,
                                 searchable = FALSE,
                                 targets    = -1),
                            list(visible = FALSE, # hide the ST_ID column
                                 targets = 0)))



# -- REACTIVES & FUNCTIONS --

user_samples <- reactive({
    all_samples <- sg_all_samples

    if (identical(tolower(get_env_value("data_mode")), "file")) {
        if (NROW(all_samples) == 0) {
            logerror("There are no samples available.")
            createAlert(session,
                        "bodyAlert",
                        "noSamplesAlert",
                        style   = "warning",
                        content = get_message_text("no_samples"),
                        dismiss = FALSE,
                        append  = FALSE)
        } else {
            user_st_objects(build_seurat_st_data_objects(samples = all_samples))
        }
    } else {
        if (is.null(g_api_connection)) {
            createAlert(session,
                        "bodyAlert",
                        "databaseConnectionAlert",
                        style   = "error",
                        content = get_message_text("backend_unavailable"),
                        dismiss = FALSE,
                        append  = FALSE)
        } else {
            if (NROW(all_samples) > 0) {
                all_samples <- filter_api_samples_by_app_user(all_samples, session$user)
            }

            if (NROW(all_samples) == 0) {
                logerror(glue("All samples were filtered out for user: {session$user}"))
                createAlert(session,
                            "bodyAlert",
                            "noSamplesAlert",
                            style   = "warning",
                            content = get_message_text("no_user_samples", message_parameters = list(user_id = session$user)),
                            dismiss = FALSE,
                            append  = FALSE)
            } else {
                user_st_objects(build_api_st_data_objects(samples = all_samples))
            }
        }
    }

    all_samples
})

pathology_user_samples <- reactive({
    data <- NULL

    if ("PATHOLOGY_TAB_READY" %in% colnames(user_samples())) {
        data <- filter(user_samples(), PATHOLOGY_TAB_READY)
    }

    data
})

matrix_gene_choices <- reactive({
    read_organism_and_probes_genes(
        samples           = user_samples(),
        limit_to_organism = TRUE)
})

box_gene_choices <- reactive({
    result <- NULL

    if (is.null(g_api_connection)) {
        # The seurat file-based system needs to use file lists because otherwise
        # each file has to be read for the list of genes (not tenable in the ui)
        result <- read_organism_and_probes_genes(
            samples           = user_samples(),
            limit_to_organism = FALSE)
    }

    result
})

user_metadata_control_cols <- reactive({
    read_metadata_column_control(user_samples())
})


#' add_single_box
#'   Add a new or cloned box in the Single tab and return the ID of the created box
#'
#' @param cloned_box_source_id - ID of the box being cloned or NULL for a new box
#'
#' @return character
add_single_box <- function(cloned_box_source_id = NULL) {
    plotImmediately <- FALSE
    datasetSel      <- character(0)
    signatureSel    <- character(0)
    sampleSel       <- character(0)
    singleGeneSel   <- character(0)
    customSigSel    <- character(0)
    coloring        <- character(0)

    if (!is.null(cloned_box_source_id)) {
        datasetSel    <- input[[glue("{cloned_box_source_id}-datasetSel")]]
        signatureSel  <- input[[glue("{cloned_box_source_id}-signatureSel")]]
        sampleSel     <- input[[glue("{cloned_box_source_id}-sampleSel")]]
        singleGeneSel <- input[[glue("{cloned_box_source_id}-singleGeneSel")]]
        customSigSel  <- input[[glue("{cloned_box_source_id}-customSigSel")]]
    } else if (userSel$matrixGridSelected) {
        datasetSel      <- input$matrixOrganismTissueSel
        signatureSel    <- userSel$matrixSignature
        sampleSel       <- userSel$matrixSample
        singleGeneSel   <- userSel$matrixGene
        coloring        <- userSel$matrixColor
        customSigSel    <- character(0)
        plotImmediately <- TRUE
    }

    box_id                 <- glue("single-box-{UUIDgenerate()}")
    userSel$singleBoxesIds <- c(userSel$singleBoxesIds, box_id)

    insertUI(selector  = "#singleBoxList",
             where     = "afterBegin",
             ui        = stSingleBoxUI(box_id),
             immediate = TRUE)

    stSingleBoxServer(id                    = box_id,
                      genes_choices         = box_gene_choices(),
                      signature_choices     = sg_signature_choices,
                      signatures            = sg_signatures,
                      samples               = user_samples(),
                      st_objects            = user_st_objects,
                      custom_signatures     = userSel$custom_signatures,
                      selected_dataset      = datasetSel,
                      selected_signature    = signatureSel,
                      selected_sample_st_id = sampleSel,
                      selected_singleGene   = singleGeneSel,
                      selected_custom_sig   = customSigSel,
                      selected_coloring     = coloring,
                      plot_on_create        = plotImmediately,
                      logger                = ss_userAction.Log)

    observeEvent(input[[glue("{box_id}-closeBoxBtn")]], {
        removeUI(selector  = glue("#{box_id}-stBoxContainerDiv"),
                 immediate = TRUE)
        remove_shiny_inputs(box_id, input)
        gc()
        userSel$singleBoxesIds <- userSel$singleBoxesIds[userSel$singleBoxesIds != box_id]
    })

    observeEvent(input[[glue("{box_id}-cloneButton")]], {
        add_single_box(cloned_box_source_id = box_id)
    })

    observeEvent(input[[glue("{box_id}-customSigSel")]], {
        userSel$custom_signatures <- unique(c(userSel$custom_signatures, input[[glue("{box_id}-customSigSel")]]))

        for (boxid in c(userSel$singleBoxesIds[userSel$singleBoxesIds != box_id],
                        userSel$pathologyBoxesIds[userSel$pathologyBoxesIds != box_id])) {
            updateSelectInput(session,
                              inputId  = glue("{boxid}-customSigSel"),
                              selected = input[[glue("{boxid}-customSigSel")]],
                              choices  = userSel$custom_signatures)
        }
    })

    userSel$matrixGridSelected <- FALSE

    box_id
}


#' add_pathology_box
#'   Add a new or cloned box in the Pathology tab and return the ID of the created box
#'
#' @param cloned_box_source_id - ID of the box being cloned or NULL for a new box
#'
#' @return character
add_pathology_box <- function(cloned_box_source_id = NULL) {
    datasetSel    <- character(0)
    signatureSel  <- character(0)
    sampleSel     <- character(0)
    singleGeneSel <- character(0)
    customSigSel  <- character(0)
    metadataSel   <- character(0)

    if (!is.null(cloned_box_source_id)) {
        datasetSel    <- input[[glue("{cloned_box_source_id}-datasetSel")]]
        signatureSel  <- input[[glue("{cloned_box_source_id}-signatureSel")]]
        sampleSel     <- input[[glue("{cloned_box_source_id}-sampleSel")]]
        singleGeneSel <- input[[glue("{cloned_box_source_id}-singleGeneSel")]]
        customSigSel  <- input[[glue("{cloned_box_source_id}-customSigSel")]]
        metadataSel   <- input[[glue("{cloned_box_source_id}-metadataSel")]]
    }

    box_id                    <- glue("pathology-box-{UUIDgenerate()}")
    userSel$pathologyBoxesIds <- c(userSel$pathologyBoxesIds, box_id)

    insertUI(selector  = "#pathologyBoxList",
             where     = "afterBegin",
             ui        = stPathologyBoxUI(box_id),
             immediate = TRUE)

    stPathologyBoxServer(id                    = box_id,
                         genes_choices         = box_gene_choices(),
                         signature_choices     = sg_signature_choices,
                         signatures            = sg_signatures,
                         samples               = pathology_user_samples(),
                         st_objects            = user_st_objects,
                         custom_signatures     = userSel$custom_signatures,
                         selected_dataset      = datasetSel,
                         selected_signature    = signatureSel,
                         selected_sample_st_id = sampleSel,
                         selected_singleGene   = singleGeneSel,
                         selected_custom_sig   = customSigSel,
                         selected_metadata     = metadataSel,
                         logger                = ss_userAction.Log)

    observeEvent(input[[glue("{box_id}-closeBoxBtn")]], {
        removeUI(selector  = glue("#{box_id}-stBoxContainerDiv"),
                 immediate = TRUE)
        remove_shiny_inputs(box_id, input)
        gc()
        userSel$pathologyBoxesIds <- userSel$pathologyBoxesIds[userSel$pathologyBoxesIds != box_id]
    })

    observeEvent(input[[glue("{box_id}-cloneButton")]], {
        add_pathology_box(cloned_box_source_id = box_id)
    })

    observeEvent(input[[glue("{box_id}-customSigSel")]], {
        userSel$custom_signatures <- unique(c(userSel$custom_signatures, input[[glue("{box_id}-customSigSel")]]))

        for (boxid in c(userSel$singleBoxesIds[userSel$singleBoxesIds != box_id],
                        userSel$pathologyBoxesIds[userSel$pathologyBoxesIds != box_id])) {
            updateSelectInput(session,
                              inputId  = glue("{boxid}-customSigSel"),
                              selected = input[[glue("{boxid}-customSigSel")]],
                              choices  = userSel$custom_signatures)
        }
    })

    box_id
}

app_walkthrough <- reactive({
    wt_data <- read.csv(PATH_WALKTHROUGH_FILE,
                        comment.char     = "#",
                        stringsAsFactors = FALSE)

    steps <- wt_data %>%
        group_by(element) %>%
        mutate(intro = ifelse(is.na(title) || is.null(title) || (title == ""), intro, glue("{tags$h4(title)}{intro}"))) %>%
        select(element, intro, position)

    selected_steps <- steps %>%
        filter((element %in% c(".btn-box-tool", ".box-header")) ||
                   grepl(input$stTabs, tolower(element)))

    selected_steps
})

matrix_metadata_factors <- reactive({
    result          <- character(0)

    if (!is.null(g_api_connection)) {
        samples <- matrix_samples()

        if (NROW(samples) > 0) {
            samples <- user_samples() %>%
                filter(ST_ID %in% samples$ST_ID)

            tags <- get_api_tag_space(samples)
            if (!is.null(tags)) {
                result <- tags
            }
        }
    } else {
        # file based values
        organism_tissue <- get_organism_tissue(input$matrixOrganismTissueSel)

        if (all(organism_tissue != "")) {
            filtered_samples <- user_samples() %>%
                filter(Organism == organism_tissue[1],
                       Tissue == organism_tissue[2])

            meta_fields <- filtered_samples %>%
                pull(SEURAT_METADATA_FIELDS) %>%
                strsplit(REGEX_ALL_OTHER_META_FIELD_SEP)

            cluster_fields <- filtered_samples$CLUSTER_FIELD
            all_choices    <- METADATA_ITEM_TISSUE

            # for each sample, rename the clustering field to the desired display name
            for (i in seq_along(meta_fields)) {
                m_fields <- meta_fields[[i]]
                c_field  <- cluster_fields[[i]]

                if (!is.na(c_field) && (c_field %in% m_fields)) {
                    m_fields[m_fields == c_field] <- METADATA_ITEM_CLUSTERING
                }

                all_choices <- c(all_choices, m_fields)
            }

            if (any(filtered_samples$HAS_PATHOLOGY)) {
                all_choices <- c(all_choices, METADATA_ITEM_PATHOLOGY)
            }

            all_choices   <- unique(all_choices)
            first_choices <- c(METADATA_ITEM_TISSUE, METADATA_ITEM_CLUSTERING, METADATA_ITEM_PATHOLOGY)

            # Provide a list of all the possible options for the samples in the selected organism/tissue group
            result <- c(first_choices[first_choices %in% all_choices],
                        sort(all_choices[!(all_choices %in% first_choices)]))
        }
    }

    result
})


# ----------------------------------------
# --          SHINY SERVER CODE         --
# ----------------------------------------

if (g_configuration == "dev-online") {
    session$user <- get_env_value('shiny_session_user')
}

#### Initialization Only ####
observeEvent(TRUE, {
    show_pathology_tab      <- as.logical(get_env_value('show_pathology_tab', unset = FALSE))
    pathology_samples_exist <- (NROW(pathology_user_samples()) > 0)

    if (show_pathology_tab && !pathology_samples_exist) {
        logwarn("There are no available samples for pathology tab. Disabling pathology tab")
    }

    if (show_pathology_tab && pathology_samples_exist) {
        shinyjs::runjs("$('#stTabs li a[data-value=pathology]').css('display', 'block')");
        add_pathology_box()
    }

    add_single_box()

    loginfo(glue("Application started with {NROW(unique(user_samples()$Organism))} organisms, {NROW(sg_signatures)} signatures and {NROW(user_samples())} samples"),
            logger = ss_userAction.Log)

    updateSelectizeInput(
        session  = session,
        inputId  = "matrixOrganismTissueSel",
        server   = TRUE,
        choices  = get_unique_organism_tissue(user_samples()),
        selected = get_unique_organism_tissue(user_samples())[[1]],
        options  = list(placeholder   = PLACE_HOLDER_SELECTIZE_INPUTS ,
                        optgroupField = "Organism",
                        labelField    = "Organism_Tissue",
                        valueField    = "Organism_Tissue",
                        searchField   = c("Organism", "Tissue"),
                        sortField     = "Organism",
                        render        = I("{ option: function(item, escape) {
                                           return '<div>' + item.Tissue + '</div>'; }}")
        )
    )
}, once = TRUE)

observeEvent(input$loginButton, {
    valid_password <- !is.null(input$password) && (input$password == sg_password)

    feedbackDanger("password", !valid_password, "Incorrect password")

    if (valid_password) {
        userSel$logged_in <- TRUE
        toggleModal(session, "loginModal", "close")
    } else {
        updateTextInput(inputId = "password", value = "")
    }
})

output$logged_in <- reactive({
    userSel$logged_in
}); outputOptions(output, "logged_in", suspendWhenHidden = FALSE)

observeEvent(input$singleNewBoxBtn, {
    add_single_box()
})

observeEvent(input$pathologyNewBoxBtn, {
    add_pathology_box()
})


observeEvent(input$matrixOrganismTissueSel, {
    organism_tissue         <- input$matrixOrganismTissueSel
    matrix_metadata_options <- NULL

    if (length(organism_tissue) == 0) {
        organism_tissue <- ""
    }

    choices <- matrix_gene_choices()[[get_organism_tissue(organism_tissue)[1]]]

    updateSelectizeInput(inputId  = "matrixSingleGeneSel",
                         server   = TRUE,
                         selected = input$matrixSingleGeneSel,
                         choices  = choices,
                         options  = list(
                             placeholder = PLACE_HOLDER_SELECTIZE_INPUTS ,
                             render      = I("{ option: function(item, escape) {
                                           return '<div><b>' + item.value + '</b></div>'; }}")))

    primary_cols_df   <- NULL
    primary_col_names <- NULL

    if (!is.null(user_metadata_control_cols())) {
            primary_cols <- user_metadata_control_cols() %>%
                filter(Type == "primary")

            primary_cols_df <- primary_cols %>%
                column_to_rownames("Column_Name") %>%
                t() %>%
                as.data.frame()

            primary_col_names <- primary_cols$Column_Name
        }

    if (organism_tissue != "") {
        organism_tissue_list        <- get_organism_tissue(organism_tissue)
        selected_organism_tissue_df <- user_samples() %>%
            filter(Organism == organism_tissue_list[1], Tissue == organism_tissue_list[2]) %>%
            rename(Sample = "Sample_App_Name")

        if (identical(tolower(get_env_value("data_mode")), "file")) {
            all_other_metadata_col <- selected_organism_tissue_df$ALL_OTHER_METADATA
            all_other_metadata_df  <- data.frame()

            if (!all(is.na(all_other_metadata_col))) {
                all_other_metadata <- str_split(all_other_metadata_col, REGEX_ALL_OTHER_META_FIELD_SEP)
            }

            if (!is.null(primary_col_names)) {
                for (sample in seq_along(all_other_metadata)) {
                    primary_data <- data.frame(sample = all_other_metadata[[sample]]) %>%
                        separate(col = "sample", into = c("key", "value"), sep = REGEX_ALL_OTHER_META_KEYVAL_SEP) %>%
                        filter(key %in% primary_col_names,
                               !(key %in% c("Tissue", "Organism"))) %>% # avoid duplication when combining with the rest of the data after
                        pivot_wider(names_from = "key", values_from = "value")
                    all_other_metadata_df <- bind_rows(all_other_metadata_df, primary_data)
                }
            }

            if (NROW(all_other_metadata_df) > 0) {
                selected_organism_tissue_df <- bind_cols(selected_organism_tissue_df, all_other_metadata_df)
            }
        }

        all_cols <- selected_organism_tissue_df %>%
            arrange(match(ST_ID, str_sort(ST_ID, numeric = TRUE))) %>%
            select(ST_ID,
                   Sample,
                   any_of(primary_col_names)) %>% # rearrange the columns according to the prescribed order
            select(where(~!all(is.na(.)))) %>%
            as.data.frame()

        matrix_samples(all_cols)

        # provide a unique identifier to grab the additional metadata after
        sample_ids <- matrix_samples()$ST_ID
        matrix_samples(bind_cols(matrix_samples(),
                                 " " = get_ui_matrix_table_button(button_id  = "button",
                                                                  sample_ids = sample_ids)))
    } else {
        blank_df <- data.frame(ST_ID = "", Sample = "")

        if (!is.null(primary_cols_df)) {
            #add the sample column and the blank column for the button
            formatted_df <- bind_cols(ST_ID = "", Sample = "", primary_cols_df, " " = "")

            blank_df <- data.frame(formatted_df[FALSE, ],
                                   check.names = FALSE)
        }

        matrix_samples(blank_df)
    }

    shinyjs::runjs("hideSearchRow();")

    updateSelectizeInput(
        session  = session,
        inputId  = "matrixMetadataSel",
        server   = TRUE,
        choices  = matrix_metadata_factors())
})


observeEvent(input$stTabs == "matrix", {
    updateSelectizeInput(
        inputId  = "matrixSignatureSel",
        server   = TRUE,
        selected = input$matrixSignatureSel,
        choices  = sg_signatures[match(sg_signature_choices, sg_signatures$Name),],
        options  = list(
            placeholder = PLACE_HOLDER_SELECTIZE_INPUTS ,
            labelField  = "Name",
            searchField = c("Name", "Gene_List"),
            valueField  = "Name",
            render      = I("{ option: function(item, escape) {
                                           return '<div><b>' + item.Name +
                                                  '</b><br>' + '<i><span style=\"display:inline-block;font-size:x-small;line-height:1.3\">' +
                                                  item.Gene_List + '</span></i></div>'; }}")
        )
    )
})


observeEvent(input$extended_button, {
    sample_id   <- strsplit(input$extended_button, "extended_button_")[[1]][2]
    sample      <- user_samples() %>% filter(ST_ID == sample_id)
    modal_title <- div(id = "extendedModalTitle",
                       glue("Extended Attributes: {sample$Sample_App_Name}"),
                       div(style = "float:right;",
                           modalButton(label = NULL, icon = icon("xmark"))))

    meta_tbl_tags    <- tagList()
    additional_types <- NULL

    if (NROW(user_metadata_control_cols()) > 0) {
        additional_types <- user_metadata_control_cols() %>%
            filter(Type != "primary", Type != "secondary") %>%
            pull(Type) %>%
            unique() %>%
            str_to_title() %>%
            sort()

        # Make sure secondary is the first table of the additional types to show
        additional_types <- c("Secondary", additional_types)
    }

    samples <- user_samples() %>%
        select(-any_of("Sample_App_Name"))

    for (type in c(additional_types)) {
        meta_tbl_tags <- tagList(
            meta_tbl_tags,
            get_ui_addl_metadata_table(samples_data     = samples,
                                       metadata_cols    = user_metadata_control_cols(),
                                       smp_id           = sample_id,
                                       type             = type))
    }

    showModal(modalDialog(title     = modal_title,
                          size      = "l",
                          easyClose = TRUE,
                          footer    = NULL,
                          meta_tbl_tags))
})


observeEvent(input$plotMatrixBtn, {
    closeAlert(session, "plotsAlert")
    matrix_samples_sel <- matrix_sample_row_sel()
    signature_gene <- c(input$matrixSignatureSel, input$matrixSingleGeneSel, input$matrixMetadataSel)

    if (NROW(matrix_samples_sel) == 0) {
        createAlert(session,
                    "bodyAlert",
                    "plotsAlert",
                    style   = "warning",
                    content = get_message_text("chart_warning_samples"),
                    append  = FALSE)
    } else if (length(signature_gene) == 0) {
        createAlert(session,
                    "bodyAlert",
                    "plotsAlert",
                    style   = "warning",
                    content = get_message_text("chart_warning_gene_metadata"),
                    append  = FALSE)
    } else {
        selected_samples  <- matrix_samples_sel$Sample
        num_spatial_chart <- length(selected_samples) * length(signature_gene)

        if (num_spatial_chart == 1) {
            num_spatial_chart_msg <- get_message_text("one_spatial_chart")
        } else {
            num_spatial_chart_msg <- get_message_text("multi_spatial_chart",
                                                      message_parameters = list(num_spatial_chart = num_spatial_chart))
        }

        output$plotMatrixConfirmText  <- renderUI({tagList(
            fluidRow(style = "text-align:center; word-wrap:break-word;",
                     column(width = 5,
                            tags$b(tags$u("Tissue/Organism")),
                            tags$br(),
                            tags$br(),
                            HTML(glue_collapse(selected_samples, sep = "<br>")),
                     ),
                     column(width = 2,
                            tags$b("X")),
                     column(width = 5,
                            tags$b(tags$u("Signature/Genes/Metadata")),
                            tags$br(),
                            tags$br(),
                            tags$div(HTML(glue_collapse(signature_gene, sep = "<br>")))
                     )),
            tags$br(),
            tags$p(align = "center", num_spatial_chart_msg),
            tags$p(align = "center", get_message_text("calc_time_warning"))
        )})

        toggleModal(session = session, modalId = "plotMatrixConfirmationModal", toggle = "open")
        shinyjs::runjs("$('#plotMatrixConfirmationModal').on('shown.bs.modal', function(e) {$('#continueMatrixPlot').focus();});")
    }

})


observeEvent(input$cancelMatrixPlot, {
    toggleModal(session, "plotMatrixConfirmationModal", toggle = "close")
})


observeEvent(input$continueMatrixPlot, {
    toggleModal(session, "plotMatrixConfirmationModal", toggle = "close")
    output$matrixPlotOutput <- renderUI(div())
    shinyjs::runjs("Shiny.setInputValue('renderMatrix', 'TRUE', {priority: 'event'});")
})


observeEvent(input$renderMatrix, {
    show_modal_spinner(spin = "circle",
                       text = "Building Matrix Plots")
    userSel$matrixConfigSamples    <- matrix_sample_row_sel()
    userSel$matrixConfigSignatures <- input$matrixSignatureSel
    userSel$matrixConfigGenes      <- input$matrixSingleGeneSel
    userSel$matrixConfigMetadata   <- input$matrixMetadataSel
    userSel$matrixConfigColoring   <- input$coloringSelectize
    userSel$matrixColoringScaling  <- input$scalingSelectize

    matrix_plot_result <- NULL
    samples_signals    <- NULL
    signatures_sel     <- sg_signatures[(sg_signatures$Name %in% userSel$matrixConfigSignatures), , drop = F]
    samples_objects    <- get_full_sample_data_objects(selected_samples = userSel$matrixConfigSamples$ST_ID,
                                                       all_data_objects = user_st_objects)
    sample_names_by_id <- setNames(sapply(samples_objects, function(x) { x@name }),
                                   nm = sapply(samples_objects, function(x) { x@st_id }))
    samples_keep       <- c()

    # remove any samples where critical data failed to be retrieved from API
    samples_keep <- sapply(samples_objects, function(x) {
        if (identical(tolower(get_env_value("data_mode")), "api")) {
            if (all(!is.null(x@api_data$features),
                    !is.null(x@api_data$cells),
                    !is.null(x@api_data$coordinates),
                    !is.null(x@slice_image),
                    !is.null(x@top_genes))) {
                x@st_id
            } else {
                NA
            }
        } else {
            x@st_id
        }
    })

    if ((length(samples_keep) > 0) && all(!is.na(samples_keep))) {
        signals_result  <- get_sample_signals(st_ids           = samples_keep,
                                              all_data_objects = user_st_objects,
                                              signatures       = signatures_sel,
                                              genes            = userSel$matrixConfigGenes)
        samples_signals <- signals_result$signals

        # remove any samples where there was an error retrieving the signals
        samples_objects <- samples_objects[!(names(samples_objects) %in% signals_result$error_samples)]
    }

    excluded_samples <- unname(sample_names_by_id[setdiff(names(sample_names_by_id), names(samples_objects))])

    if (length(excluded_samples) > 0) {
        excluded_samples <- glue_collapse(glue("'{excluded_samples}'"), sep = ", ")
        dataset          <- userSel$matrixConfigSamples %>%
            mutate(dataset = glue("{Organism} > {Tissue}")) %>%
            pull(dataset) %>%
            unique()

        logwarn(glue("Unable to retrieve critical components for '{dataset}' sample(s): {excluded_samples}"))
        createAlert(session,
                    "bodyAlert",
                    "plotsAlert",
                    style   = "warning",
                    content = get_message_text("matrix_samples_excluded",
                                               message_parameters = list(dataset          = dataset,
                                                                         excluded_samples = excluded_samples)),
                    append  = FALSE)
    }

    if (length(samples_objects) == 0) {
        matrix_plot_result <- tags$div(style = "text-align:center;",
                                       tags$h4(class = "text-danger",
                                               "Plots Unavailable"))
        shinyjs::runjs("Shiny.setInputValue('matrixRenderIsComplete', 'TRUE', {priority: 'event'})")
    } else {
        matrix_plot_result <- get_ui_matrix_plot_table(
            selected_samples     = samples_objects,
            selected_signatures  = signatures_sel,
            selected_genes       = userSel$matrixConfigGenes,
            selected_metadata    = userSel$matrixConfigMetadata,
            coloring             = userSel$matrixConfigColoring,
            coloring_scaling     = userSel$matrixColoringScaling,
            samples_signals      = samples_signals,
            selectedMatrixCol    = selectedMatrixCol,
            selectedMatrixRow    = selectedMatrixRow,
            samples              = user_samples(),
            logger               = ss_userAction.Log)
    }

    output$matrixPlotOutput <- renderUI({
        matrix_plot_result
    })
}, ignoreNULL = TRUE)


observeEvent(selectedMatrixRow(), {
    show_modal_spinner(spin = "circle",
                       text = "Switching to Single Tab")
    selectedRow <- selectedMatrixRow()

    # switch to the single tab
    updateTabsetPanel(inputId  = "stTabs",
                      selected = "single")

    # clear out the current Single Tab contents
    lapply(userSel$singleBoxesIds, function(x) {
        removeUI(selector  = glue("#{x}-stBoxContainerDiv"),
                 immediate = TRUE)
        gc()
    })

    lapply(userSel$singleBoxesIds, remove_shiny_inputs, input = input)
    userSel$singleBoxesIds <- NULL

    # create the necessary boxes in the single tab
    for (col in rev(c(selectedRow$signatures, selectedRow$genes))) {
        userSel$matrixSample <- selectedRow$sample_ids

        if (col %in% selectedRow$genes) {
            userSel$matrixSignature <- character(0)
            userSel$matrixGene      <- col
        }

        if (col %in% selectedRow$signatures) {
            userSel$matrixSignature <- col
            userSel$matrixGene      <- character(0)
        }

        userSel$matrixColor        <- selectedRow$coloring
        userSel$matrixGridSelected <- TRUE
        add_single_box()
    }
})


observeEvent(selectedMatrixCol(), {
    show_modal_spinner(spin = "circle",
                       text = "Switching to Single Tab")
    selectedCol <- selectedMatrixCol()

    # switch to the single tab
    updateTabsetPanel(inputId  = "stTabs",
                      selected = "single")

    # clear out the current Single Tab contents
    lapply(userSel$singleBoxesIds, function(x) {
        removeUI(selector  = glue("#{x}-stBoxContainerDiv"),
                 immediate = TRUE)
        gc()
    })

    lapply(userSel$singleBoxesIds, remove_shiny_inputs, input = input)
    userSel$singleBoxesIds <- NULL

    # create the necessary boxes in the single tab
    if (is.null(selectedCol$genes)) {
        userSel$matrixSignature <- selectedCol$signatures
        userSel$matrixGene      <- character(0)
    }

    if (is.null(selectedCol$signatures)) {
        userSel$matrixSignature <- character(0)
        userSel$matrixGene      <- selectedCol$genes
    }

    for (row in rev(selectedCol$sample_ids)) {
        userSel$matrixSample       <- row
        userSel$matrixColor        <- selectedCol$coloring
        userSel$matrixGridSelected <- TRUE
        add_single_box()
    }
})


observeEvent(c(input$stTabs == "single", input$matrixRenderIsComplete), {
    remove_modal_spinner()
})
