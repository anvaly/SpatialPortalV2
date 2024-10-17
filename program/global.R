# ----------------------------------------
# --          PROGRAM global.R          --
# ----------------------------------------
# USE: Global variables and functions
#
# NOTEs:
#   - All variables/functions here are
#     globally scoped and will be available
#     to server, UI and session scopes
# ----------------------------------------

suppressPackageStartupMessages({
    library(assertthat)
    library(dplyr)
    library(tidyr)
    library(glue)
    library(readr)
    library(canvasXpress)
    library(grDevices)
    library(stringr)
    library(rintrojs)
    library(DT)
    library(shinyFeedback)
    library(tictoc)
    library(config)
    library(shinybusy)
    library(uuid)
    library(tibble)
    library(future)
})

options(dplyr.summarise.inform = FALSE)
plan("future::multisession")


set_app_parameters(title       = "BMS Spatial Portal for Treatment Naïve PDAC",
                   titleinfo   = NULL,
                   loglevel    = "DEBUG",
                   showlog     = FALSE,
                   app_version = "5.2")

support_files <- list.files(path = "program/fxn", pattern = "^.*.R$", full.names = T, recursive = T)
for (sf in support_files) {
    source(sf, local = T)
}

source("program/modules/_control.R")

g_display_top_n_sig_genes   <- NULL
g_display_top_n_spot_genes  <- NULL
g_display_dynamic_range_min <- NULL
g_allowed_organisms         <- NULL
g_api_connection            <- NULL
g_configuration             <- Sys.getenv("PUBLIC_ST_APP_CONFIG", unset = "dev-public")
g_config_check              <- application_config_check()
g_cx_license                <- get_cx_license()

if (is.logical(g_config_check) && !g_config_check) {
    stop(glue("The application is misconfigured for configuration: {g_configuration}"))
} else {
    mock_lib <- get_env_value('concentriq_function_file')

    if (!is.null(mock_lib)) {
        # override sourced functions from data_concentriq.R
        source(mock_lib)
    }
}

