#' get_message_text
#'     Build user message given message name and list of message parameters if any
#'
#' @param message_name       - string of message name
#' @param message_parameters - optional list of message parameters,
#'                             note that the number of message parameters should match the message selected
#'                             and must have the same names (default is empty list)
#'
#' @return string of formatted message
get_message_text <- function(message_name, message_parameters = list()) {
    message_text   <- NULL
    support_text   <- ""

    if (is.null(message_name) || message_name == "") {
        logwarn("Message cannot be NULL or empty")
    } else {
        message_text <- switch(
            message_name,
            # NOTE: if gluing on one of the text messages (support, etc) you must escape any other parameters with double {{braces}}
            backend_unavailable                    = glue("<b>The backend server/database is not available.</b><br/>{support_text}"),
            calc_time_warning                      = "This may take several minutes to calculate. Would you like to proceed?",
            chart_warning_gene_metadata            = "Charts cannot be created without choosing one or more Gene Signatures, Genes or Additional Metadata.",
            chart_warning_samples                  = "Charts cannot be created without choosing one or more samples. Please choose an organism/tissue and then select one or more lines in the resulting samples table.",
            contact_app_author                     = "Contact the app author or data owner for more information if needed.",
            dataset_required                       = "Dataset must be selected",
            enter_one_gene                         = "Please enter at least one gene",
            genes_not_represented                  = "This may occur if the selected gene(s) are not well represented in the sample.",
            genes_not_represented_info_unavailable = "This may occur if the selected gene(s) are not well represented, or extended information (clustering, pathology) is not available in the sample.",
            multi_spatial_chart                    = "{num_spatial_chart} Spatial charts have been requested.",
            no_samples                             = "There are no samples available to the application at this time. Please contact the application author for assistance.",
            no_user_samples                        = "There are no samples that can be visualized for userID: '{user_id}' . Please contact the application author for assistance.",
            one_spatial_chart                      = "1 Spatial chart has been requested.",
            only_one_selection_single_box          = "Only one of a Gene Signature, a Single Gene, or a Custom value can be selected",
            only_one_selection_pathology_box       = "Only one of a Gene Signature, Single Gene, Custom, or Metadata value can be selected",
            sample_required                        = "Sample must be selected",
            small_dynamic_range                    = "The signal dynamic range is too small for the chosen gene.",
            unable_create_spatial_extended         = "Unable to create the requested spatial plots",
            unable_create_spatial_overview         = "Unable to create the requested spatial plot",
            data_item_missing                      = "Unable to retrieve {item} for '{dataset}' sample '{sample}'.<br/>Contact the app owner or data owner for more information if this persists.",
            unexpected_data_retrieval_error        = "Data could not be retrieved due to an unexpected error.",
            matrix_samples_excluded                = "Unable to retrieve critical components for '{dataset}' sample(s): {excluded_samples}.<br/>Contact the app owner or data owner for more information if this persists.",
            top_genes_unexpected_error             = "Unable to retrieve top represented genes due to an unexpected error",
            extended_dot_plot_zero_variance        = "There is zero variance in the data so the plot cannot be created",
            pathology_select_spatial_points        = "Select one or more points on the spatial chart to the left<br/> then press the &quot;Get Pathology Image&quot; button below",
            pathology_tissue_area                  = "The tissue area shown will be a rectangular area<br/> around the outer bounds of all selected points",
            pathology_image_size                   = "Images are reduced from maximum resolution for screen viewing. Download the image for a larger version.",
            pathology_maximum_spots                = "The maximum number of spots allowed in a selection is {value} for performance reasons. Larger selections take significantly longer to render and download.",
            pathology_hints                        = "<div>Helpful Hints:<br/>
                                                      <ul>
                                                          <li>Click to select a single point</li>
                                                          <li>Hold ctrl while selecting to add points to a selection</li>
                                                          <li>Hold shift while dragging to select an area of points</li>
                                                      </ul></div>",
            pathology_no_points_selected           = "There are no points selected in the chart to the left.<br/>Select one or more points to define the pathology area.",
            pathology_sel_max_points_exceeded      = "The selected area contains more than {value} points. Please select a smaller area of points.",
            pathology_image_retrieval_error        = "Unable to retrieve the requested image."
            )

        if (!is.null(message_text)) {
            message_encoded_parameters <- regmatches(message_text,
                                                     gregexpr(REGEX_MESSAGE_ENCODED_PARAMETERS, message_text, perl = TRUE))[[1]]
            missed_params              <- setdiff(message_encoded_parameters, names(message_parameters))
            missed_params_names        <-  ""

            if ((length(missed_params) > 0) || (length(message_parameters) != length(message_encoded_parameters))) {

                if (length(missed_params) > 0) {
                    few_or_many         <- "few"
                    missed_params_names <- glue("and missing parameters '{glue_collapse(missed_params, sep = \", \")}'")

                    # loop over missed params to initialize and add them to main param list
                    for (missed_param in missed_params) {
                        message_parameters[[missed_param]] <- ""
                    }
                } else if (length(message_parameters) > 0) {
                    few_or_many <- "many"
                }

                logwarn(glue("Message '{message_name}' requested with too {few_or_many} parameters. This message requires {length(message_encoded_parameters)} value(s) {missed_params_names}"))
            }
            message_text <- do.call(glue, c(message_text, message_parameters, .null = "", .na = ""))
        } else {
            logwarn(glue("Message '{message_name}' could not be found"))
        }
    }

    message_text
}
