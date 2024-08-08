#' matrixGridtoSingleUI
#'   Create the action link in the UI
#'
#' @param id           - a unique id to use in the module's namespace
#' @param display_name - text to show in the UI
#'
#' @return taglist
matrixGridtoSingleUI <- function(id, display_name) {
    ns <- NS(id)
    tagList(
        shinyjs::useShinyjs(),
        tags$a(actionLink(ns(id), display_name), href = "#"),
        bsModal(
            tags$head(tags$style(
                paste(glue("#{ns(\"gridToSingleConfirmationModal\")}"), ".modal-footer{display:none}"))),
            id      = ns("gridToSingleConfirmationModal"),
            title   = tags$h3(align = "center", "Load Single Samples"),
            trigger = NULL,
            size    = "large",
            uiOutput(ns("gridToSingleConfirmText")),
            tags$br(),
            footer = tagList(
                fluidRow(
                    column(width = 4),
                    column(width = 2, actionButton(ns("cancelGridToSingle"), "Cancel")),
                    column(width = 2, actionButton(ns("continueGridToSingle"), "Continue")),
                    column(width = 4)
                )
            ))
    )

}


#' matrixGridtoSingleServer
#'   Server functionality for reacting to the matrix row/column clicks
#'
#' @param id                 - plot id
#' @param sample_names_by_id - selected samples names as a named vector with ST_IDs as the vector names
#' @param signatures         - name of the signature selected
#' @param genes              - name of the gene selected
#' @param coloring           - selected coloring
#' @param selectedMatrix     - reactive variable to store the selected values in
#'
#' @return Nothing
matrixGridtoSingleServer <- function(id,
                                     sample_names_by_id,
                                     signatures,
                                     genes,
                                     coloring,
                                     selectedMatrix) {
    moduleServer(
        id,
        function(input, output, session) {
            # setup internal variables so it can be used later
            internal_columns    <- c(genes, signatures)
            internal_samples    <- sample_names_by_id
            internal_genes      <- genes
            internal_coloring   <- coloring
            internal_signatures <- signatures

            # action link clicked
            observeEvent(input[[id]], {
                # add some content based on which row/column was clicked
                # paste is used right now because it can be easily vectorized compared to glue
                charts <- paste(unname(internal_samples), "-", internal_columns)

                output$gridToSingleConfirmText  <- renderUI({tagList(
                    fluidRow(style = "text-align:center; word-wrap:break-word;",
                             tags$p(align = "center", "The following charts have been requested:"),
                             tags$p(align = "center", tags$div(HTML(glue_collapse(charts, sep = "<br>")))),
                             tags$br(),
                             tags$p(align = "center",
                                    "All existing charts will be cleared out of the Single tab and replaced with the above charts. Would you like to proceed?")
                    ))})

                toggleModal(session, "gridToSingleConfirmationModal", toggle = "open")
                shinyjs::runjs("$('#gridToSingleConfirmationModal').on('shown.bs.modal', function(e) {$('#continueGridToSingle').focus();});")

            }, ignoreInit = TRUE)

            observeEvent(input$cancelGridToSingle, {
                toggleModal(session, "gridToSingleConfirmationModal", toggle = "close")
            })

            observeEvent(input$continueGridToSingle, {
                toggleModal(session, "gridToSingleConfirmationModal", toggle = "close")
                shinyjs::runjs(glue("Shiny.setInputValue('{{session$ns('switchToSingle')}}', 'TRUE', {priority: 'event'});",
                           .open = "{{",
                           .close = "}}"))
            })

            observeEvent(input$switchToSingle, {
                selectedMatrix(list(sample_ids = names(internal_samples),
                                    genes      = internal_genes,
                                    signatures = internal_signatures,
                                    coloring   = internal_coloring,
                                    guid       = UUIDgenerate()))
            })
        })
}
