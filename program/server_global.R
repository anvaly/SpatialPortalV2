# ----------------------------------------
# --      PROGRAM server_global.R       --
# ----------------------------------------
# USE: Server-specific variables and
#      functions for the main reactive
#      shiny server functionality.  All
#      code in this file will be put into
#      the framework outside the call to
#      shinyServer(function(input, output, session)
#      in server.R
#
# NOTEs:
#   - All variables/functions here are
#     SERVER scoped and are available
#     across all user sessions, but not to
#     the UI.
#
#   - For user session-scoped items
#     put var/fxns in server_local.R
#
# FRAMEWORK VARIABLES
#     none
# ----------------------------------------

# -- IMPORTS --


# -- VARIABLES --

# app can be password-protected using a file-based password in any mode
sg_password <- NULL
if (get_env_value("require_password")) {
    sg_password <- get_app_password()
}

sg_signatures        <- read_signatures()
sg_signature_choices <- sg_signatures$Name[order(sg_signatures$Name)] %>% as.character()

# get all sample metadata - this is done at a server-global level for efficiency and will
# be filtered later (when applicable) to the user's allowed samples for local sessions
sg_all_samples       <- {
    if (identical(tolower(get_env_value("data_mode")), "file")) {
        get_cached_metadata()
    } else {
        get_api_all_sample_metadata()
    }
}

# -- FUNCTIONS --

