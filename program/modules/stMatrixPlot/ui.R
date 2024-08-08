#' stMatrixPlotUI
#'   Create visualization spot for a spatial transcriptomics sample
#'
#' @param id - character - id to use in the module's namespace
#'
#' @return shiny.tag
stMatrixPlotUI <- function(id) {
    ns <- NS(id)
    uiOutput(ns("st_plot"))
}
