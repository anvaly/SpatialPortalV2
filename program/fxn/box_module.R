# Supporting functions that are used across the box creating modules (ie in the Pathology and Single tabs)

#' clear_chart_box_inputs_feedback
#'   Clear warning feedback on inputs in chart boxes (Single and Pathology tabs)
#'
#' @param additional_input_ids - optional vector of additional inputs to clear feedback
#'
#' @return NULL
clear_chart_box_inputs_feedback <- function(additional_input_ids = character(0)) {
    input_ids <- c("datasetSel", "sampleSel", "signatureSel", "singleGeneSel", "customSigSel", additional_input_ids)

    for (input_id in input_ids) {
        hideFeedback(input_id)
    }

    shinyjs::hide(id = "selectionFeedbackDiv")
}


#' get_overall_top_genes
#'   Get the overall top genes for the selected sample and selected gene signature or single gene.
#'   * If a gene signature is selected, get the top N genes within this signature.
#'   * If a single gene is selected, this gene is used.
#'
#' @param selected_signature    - currently selected gene or gene signature
#' @param top_represented_genes - data.frame of top represented genes for each spot in the sample
#'
#' @return character
get_overall_top_genes <- function(selected_signature, top_represented_genes) {
    output <- NULL

    if (length(selected_signature) == 1) {
        output <- selected_signature
    } else if (!is.null(top_represented_genes)) {
        output <- top_represented_genes %>%
            filter(Symbol %in% selected_signature) %>%
            group_by(Symbol) %>%
            summarise(score_mean = mean(Value)) %>%
            arrange(desc(score_mean)) %>%
            select(Symbol) %>%
            top_n(g_display_top_n_sig_genes) %>%
            pull(Symbol)
    }

    if (length(output) > 0) {
        output <- glue_collapse(output, sep = "<br/>")
    }
    output
}


#' is_valid_selection
#'   Determine if selected value is valid
#'
#' @param value
#'
#' @return logical
is_valid_selection <- function(value){
    # selection is valid if it has a real value
    !is.null(value) && length(value) > 0 && value != ""
}


#' update_dataset_selectized_input
#'   Update dataset input choices with the datasets from the available samples
#'
#' @param session          - module session object
#' @param samples          - data.frame of metadata for all available samples
#' @param selected_dataset - currently selected dataset
update_dataset_selectized_input <- function(session, samples, selected_dataset = character(0)) {
    choices <- get_unique_organism_tissue(samples)
    updateSelectizeInput(
        session  = session,
        inputId  = "datasetSel",
        server   = TRUE,
        selected = choices[[1]],
        choices  = choices,
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
}


#' update_signatures_selectized_input
#'   Update gene signature input choices with the applicable genes for the selected dataset
#'
#' @param session            - module session object
#' @param signatures         - data.frame of genes for each gene signature
#' @param signature_choices  - character vector of all gene signature names
#' @param selected_signature - currently selected gene signature
update_signatures_selectized_input <- function(session, signatures, signature_choices, selected_signature = character(0)) {
    updateSelectizeInput(
        session,
        "signatureSel",
        server   = TRUE,
        selected = selected_signature,
        choices  = signatures[match(signature_choices, signatures$Name),],
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
}


#' update_customSigSel_selectized_input
#'   Update custom signature input choices with the so far entered signature during current session
#'
#' @param session               - module session object
#' @param custom_signatures     - character vector of all custom signatures
#' @param selected_customSigSel - currently selected custom signature
update_customSigSel_selectized_input <- function(session, custom_signatures, selected_customSigSel = character(0)) {
    JS_render_genes <- I("{
                       option_create: function(data, escape) {
                                          label = data.input.toString().trim().toUpperCase();
                                          return '<div class=\"create\">' + label + '</div>';
                                      },
                      option: function(item, escape) {
                                   var gene_symbol = item.value.toString().trim().toUpperCase();
                                   return '<div>' + gene_symbol + '</div>';
                              }
                         }")

    JS_create_custom_genes <- I("function(input) {
                                genes_list = input.trim().toUpperCase().split(/[ ,;]+/);
                                genes_list = [...new Set(genes_list)].toString();
                                custom_gene = {'value': genes_list, 'text': genes_list};
                                return custom_gene;
                            }")

    updateSelectizeInput(session, "customSigSel",
                         selected = selected_customSigSel,
                         choices  = custom_signatures,
                         options  = list(placeholder  = PLACE_HOLDER_CUSTOM_GENE_INPUT,
                                         labelField   = "text",
                                         searchField  = "value",
                                         valueField   = "value",
                                         render       = JS_render_genes,
                                         delimiter    = ":",
                                         create       = JS_create_custom_genes))
}


#' update_samples_selectized_input
#'   Update sample input choices with the applicable samples for the selected dataset
#'
#' @param session         - module session object
#' @param samples         - data.frame of metadata for all available samples
#' @param datasetSel      - currently selected dataset
#' @param selected_sample - currently selected sample
update_samples_selectized_input <- function(session, samples, datasetSel, selected_sample = character(0)) {
    datasetSel <- ifelse(length(datasetSel) == 0, "", datasetSel)

    sample_choices  <- character(0)
    organism_tissue <- get_organism_tissue(datasetSel)

    if (NROW(samples) > 0) {
        samples <- samples %>%
            filter(Organism %in% organism_tissue[1], Tissue %in% organism_tissue[2])

        if (NROW(samples) > 0) {
            sample_choices <- samples %>%
                arrange(match(ST_ID, str_sort(ST_ID, numeric = TRUE))) %>%
                select(ST_ID, Sample_App_Name)
        }
    }

    updateSelectizeInput(session,
                         "sampleSel",
                         server   = TRUE,
                         selected = selected_sample,
                         choices  = sample_choices,
                         options  = list(
                             placeholder = PLACE_HOLDER_SELECTIZE_INPUTS ,
                             labelField  = "Sample_App_Name",
                             searchField = "Sample_App_Name",
                             valueField  = "ST_ID",
                             render      = I("{ option: function(item, escape) {
                                                            return '<div><b>' + item.Sample_App_Name + '</b></div>'; }}")
                         )
    )
}
