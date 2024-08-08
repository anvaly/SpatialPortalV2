# ----------------------------------------
# --          PROGRAM ui_body.R         --
# ----------------------------------------
# USE: Create UI elements for the
#      application body (right side on the
#      desktop; contains output) and
#      ATTACH them to the UI by calling
#      add_ui_body()
#
# NOTEs:
#   - All variables/functions here are
#     not available to the UI or Server
#     scopes - this is isolated
# ----------------------------------------

# -- IMPORTS --


# ----------------------------------------
# --      BODY ELEMENT CREATION         --
# ----------------------------------------

# -- Create Elements

headerChanges <- tags$head(
    shinyjs::useShinyjs(),
    introjsUI(),
    useShinyFeedback(),
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$script(src = "custom.js")
)

topElements <- column(
    width = 12,
    tags$table(
        id    = "searchWidget",
        width = "100%",
        tags$tr(width = "100%",
                tags$td(width = "100%",
                        box(id          = "aboutBox",
                            title       = "About This Application",
                            width       = "100%",
                            collapsible = TRUE,
                            collapsed   = FALSE,
                            read_about())))
))

matrix_tab <- tabPanel(
    title = "Matrix",
    value = "matrix",
    tags$div(
        id = "matrixTab",
        fluidRow(
            column(
                width = 8,
                fluidRow(
                    style = list("padding-left: 20px;"),
                    tags$table(
                        tags$tr(width = "100%",
                                tags$td(width = "20%",
                                        style = "padding-bottom: 15px",
                                        ui_tooltip(id    = "organismTissueTxt",
                                                   label = "Organism/Tissue:",
                                                   text  = "Choose from the available Organisms and Tissues")),
                                tags$td(width = "80%",
                                        selectizeInput(inputId = "matrixOrganismTissueSel",
                                                       label   = NULL,
                                                       choices = NULL)))),
                    downloadableTableUI(id            = "matrixSampleTbl",
                                        downloadtypes = NULL,
                                        hovertext     = NULL,
                                        contentHeight = "300px",
                                        singleSelect  = FALSE))
            ),
            column(
                width = 4,
                style = "padding-left: 45px",
                tags$h3("Matrix Plot Configuration", class = "plotConfigurationHeader"),
                tags$hr(class = "plotConfigurationHeaderSeparator"),
                tags$br(),
                tags$div(id = "matrixOptions",
                         selectizeInput(inputId  = "matrixSignatureSel",
                                        label    =  ui_tooltip("signatureTooltip",
                                                               "Gene Signature",
                                                               "Choose custom gene signature(s)"),
                                        width    = "100%",
                                        choices  = NULL,
                                        multiple = TRUE,
                                        options  = list(placeholder = PLACE_HOLDER_SELECTIZE_INPUTS ,
                                                        searchField = "value",
                                                        plugins     = list("remove_button"))),
                         selectizeInput(
                             inputId  = "matrixSingleGeneSel",
                             label    = ui_tooltip("geneTooltip",
                                                   "Single Gene",
                                                   "Choose gene(s) available in the sample(s)"),
                             width    = "100%",
                             choices  = NULL,
                             multiple = TRUE,
                             options  = list(placeholder = PLACE_HOLDER_SELECTIZE_INPUTS ,
                                             plugins     = list("remove_button"))),
                         selectizeInput(
                             inputId  = "matrixMetadataSel",
                             label    = ui_tooltip("metadataTooltip",
                                                   "Additional Metadata",
                                                   "Choose additional data fields to plot in the matrix"),
                             width    = "100%",
                             choices  = NULL,
                             multiple = TRUE,
                             options  = list(placeholder = PLACE_HOLDER_SELECTIZE_INPUTS ,
                                             plugins     = list("remove_button")))
                ),
                fluidRow(
                    id = "matrixColoring",
                    column(
                        width = 6,
                        selectizeInput(
                            inputId   = "coloringSelectize",
                            label     =  ui_tooltip(id    = "ColoringTxt",
                                                    label = "Spatial Plot Coloring",
                                                    text  =  glue("Plots display pre-calculated normalized expression values (SCT). ",
                                                                  "The options modify the range covered by the colormap:<br/>",
                                                                  "<b>Native</b>: Full data range &nbsp;&nbsp;&nbsp;<br/>",
                                                                  "<b>Trimmed</b>: 1st to 99th percentile<br/>",
                                                                  "<b>Fixed</b>: Fixed range of [0-1]<br/>")),
                            choices  = OPTIONS_SPATIAL_PLOT_COLORING,
                            selected = "trim"
                        )
                    ),
                    column(
                        width = 6,
                        selectizeInput(
                            inputId   = "scalingSelectize",
                            label     =  ui_tooltip(id    = "ScalingTxt",
                                                    label = "Spatial Plot Scaling",
                                                    text  =  "Choose how the scaling should be applied"),
                            choices  = OPTIONS_SPATIAL_PLOT_SCALING,
                            selected = "individually"
                        )
                    )
                ),
                tags$span(style = "float:right;",
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          bsButton(inputId = "plotMatrixBtn",
                                   label   = tags$div("Create Plot Matrix", icon("arrow-right")),
                                   style   = "primary"))

            ) # column
        ), # fluidRow
        tags$hr(),
        fluidRow(uiOutput("matrixPlotOutput"))
    )
)


single_tab <- tabPanel(title = "Single",
                       value = "single",
                       tags$div(
                           style = "margin-bottom:10px; margin-left:15px",
                           bsButton(inputId = "singleNewBoxBtn",
                                    label   = "",
                                    size    = "large",
                                    icon    = icon("plus")),
                           bsTooltip(id        = "singleNewBoxBtn",
                                     title     = "Create a new UI Box Element",
                                     placement = "bottom"),
                           tags$span(style = "margin-left:5px; color:#3C8DBC; font-weight:bold",
                                     "Click to create a new chart widget below."),
                           tags$span(style = "color:gray; font-style:italic",
                                     "Chart widgets can be collapsed using the icon on the top right,
                            or closed/cloned using the buttons in the left-side configuration pane.")),
                       tagList(
                           fluidRow(
                               column(width = 12,
                                      tags$div(id = "singleBoxList"))))
)


pathology_tab <- tabPanel(title = "Pathology",
                          value = "pathology",
                          tags$div(
                              style = "margin-bottom:10px; margin-left:15px",
                              bsButton(inputId = "pathologyNewBoxBtn",
                                       label   = "",
                                       size    = "large",
                                       icon    = icon("plus")),
                              bsTooltip(id        = "pathologyNewBoxBtn",
                                        title     = "Create a new UI Box Element",
                                        placement = "bottom"),
                              tags$span(style = "margin-left:5px; color:#3C8DBC; font-weight:bold",
                                        "Click to create a new chart widget below."),
                              tags$span(style = "color:gray; font-style:italic",
                                        "Chart widgets can be collapsed using the icon on the top right,
                            or closed/cloned using the buttons in the left-side configuration pane.")),
                          fluidRow(
                              column(width = 12,
                                     tags$div(id = "pathologyBoxList")))
)


stTabs <- tabBox(
    id       = "stTabs",
    selected = "matrix",
    width    = 12,
    matrix_tab,
    single_tab,
    pathology_tab
)


# -- Modals
plotMatrixConfirmationModal <- bsModal(tags$head(tags$style("#plotMatrixConfirmationModal .modal-footer{display:none}")),
                                       id      = "plotMatrixConfirmationModal",
                                       title   = tags$h3(align = "center", "Selection Confirmation"),
                                       trigger = NULL,
                                       size    = "large",
                                       uiOutput("plotMatrixConfirmText"),
                                       tags$br(),
                                       footer = tagList(
                                           fluidRow(
                                               column(width = 4),
                                               column(width = 2, actionButton("cancelMatrixPlot", "Cancel")),
                                               column(width = 2, actionButton("continueMatrixPlot", "Continue")),
                                               column(width = 4)
                                           )
                                       ))

loginModal <- bsModal(id      = "loginModal",
                      title   = "Spatial Portal",
                      trigger = NULL,
                      size    = "small",
                      tags$head(tags$script("$(document).ready(function(){$('#loginModal').modal({backdrop: 'static', keyboard: false});});"),
                                tags$style("#loginModal .modal-header button.close{display:none;}"),
                                tags$style("#loginModal .modal-title {font-weight:bold;}"),
                                tags$style("#loginModal .modal-footer {display:none;}")),
                      tags$p(style = "font-weight: bold; text-align: center;",
                             "Pre-Release Password Required"),
                      tags$p("A password is required to access this application prior to the publication date for the paper containing the data.  Please contact the collaborators you are working with for the password if you do not know it."),
                      passwordInput(inputId     = "password",
                                    label       = NULL,
                                    value       = "",
                                    width       = "100%",
                                    placeholder = "Password"),
                      bsButton(inputId = "loginButton",
                               label   = "Continue",
                               width   = "100%"))

appBody <- if (get_env_value("require_password")) {
    tagList(loginModal,
            conditionalPanel(condition = "output.logged_in",
                             topElements,
                             stTabs,
                             plotMatrixConfirmationModal))
} else {
    tagList(topElements,
            stTabs,
            plotMatrixConfirmationModal)
}

# -- Register Elements in the ORDER SHOWN in the UI
add_ui_body(list(headerChanges, appBody))
