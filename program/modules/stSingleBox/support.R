### SINGLE BOX SPECIFIC FUNCTIONS####
#' get_unique_organism_tissue
#'   Construct a data.frame of unique organism-tissue combinations
#'
#' @param samples - data.frame of samples metadata
#'
#' @return data.frame
get_unique_organism_tissue <- function(samples) {
    results <- NULL
    if (!is.null(samples) && NROW(samples) > 0) {
        results <- samples %>%
            select(Organism, Tissue) %>%
            arrange(Tissue) %>%
            unique() %>%
            unite(col = "Organism_Tissue", sep = REGEX_ORGANISM_TISSUE_SEP, remove = FALSE)
    }
    results
}


#' switch_view_type
#'   Switch chart view in the UI box to standard ("Regular") or advanced ("Extended") view
#'
#' @param session  - module session object
#' @param standard - logical - TRUE for standard view, FALSE for advanced view
#'
#' @return
switch_view_type <- function(session, standard = TRUE) {
    if (standard) {
        updateButton(session, session$ns("viewToggleBtn_std"), value = TRUE,  style = "warning")
        updateButton(session, session$ns("viewToggleBtn_adv"), value = FALSE, style = "default")

        updateTabsetPanel(inputId  = "viewTypeTabs",
                          selected = "Regular")
        updateTabsetPanel(inputId  = "bottomExtension",
                          selected = "Regular_bottom")
    } else {
        updateButton(session, session$ns("viewToggleBtn_std"), value = FALSE, style = "default")
        updateButton(session, session$ns("viewToggleBtn_adv"), value = TRUE,  style = "warning")

        updateTabsetPanel(inputId  = "viewTypeTabs",
                          selected = "Advanced")
        updateTabsetPanel(inputId  = "bottomExtension",
                          selected = "Advanced_bottom")
    }
}
