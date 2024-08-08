#' application_config_check
#'    This is called at application startup and it will ensure the
#'    application setup is correct and the necessary environment
#'    variables are present and either can be defaulted or appear to be correct.
#'
#' @return logical, TRUE if the config is valid or FALSE if invalid
application_config_check <- function() {
    result                    <- FALSE
    seurat_pkg_check_result   <- FALSE
    revealsc_pkg_check_result <- FALSE
    scidb_pkg_check_result    <- FALSE
    ucell_pkg_check_result    <- FALSE
    ebimage_pkg_check_result  <- FALSE

    if (exists("g_configuration") && !is.null(g_configuration)) {
        config_file <- get_config_file()

        # check if the config files contain the required fields
        if (grepl("api", g_configuration)) {
            result <- check_required_fields(
                required_fields = c(CONFIG_COMMON_REQUIRED_FIELDS, CONFIG_API_REQUIRED_FIELDS))
        } else {
            result <- check_required_fields(
                required_fields = c(CONFIG_COMMON_REQUIRED_FIELDS, CONFIG_SEURAT_REQUIRED_FIELDS))
        }
    }

    # if all required fields exist, perform custom setup for different configurations if needed
    if (result) {
        # Common configurations check
        ## use global environment variable .GlobalEnv for assigning values to global variables only
        ## as it is safer from using "<<-" and shorter than using "assign" function
        ## no need to use it upon global variables values retrieval
        .GlobalEnv$g_display_top_n_sig_genes   <- get_top_genes(top_gene_config_var = "display_top_n_sig_genes", default_value = 10)
        .GlobalEnv$g_display_top_n_spot_genes  <- get_top_genes(top_gene_config_var = "display_top_n_spot_genes", default_value = 3)
        .GlobalEnv$g_display_dynamic_range_min <- get_numeric_config_vars(config_var = "display_dynamic_range_min", default_value = 0)

        # perform other section specific configurations check if all common configurations are correct
        if (result) {
            #ALL Configurations
            tryCatch({
                ucell_pkg_check_result <- check_package_version(package_name = "UCell",
                                                                min_version  = get_env_value("pkg_min_ucell_version"))
                if (ucell_pkg_check_result) {
                    suppressPackageStartupMessages({
                        library(UCell)
                    })
                }

                ebimage_pkg_check_result <- all(requireNamespace("EBImage", quietly = TRUE),
                                                requireNamespace("httr", quietly = TRUE))

                if (ebimage_pkg_check_result) {
                    suppressPackageStartupMessages({
                        library(EBImage)
                        library(httr)
                    })
                } else {
                    warning("Required package EBImage is not installed")
                }
            },
            error = function(e) {
                warning(e)
            })

            # FILE-BASED Application
            if (identical(tolower(get_env_value("data_mode")), "file")) {
                # check that we have the Seurat package needed
                 tryCatch({
                     seurat_pkg_check_result <- check_package_version(package_name = "Seurat",
                                                                      min_version  = get_env_value("pkg_min_Seurat_version"),
                                                                      max_version  = get_env_value("pkg_max_Seurat_version"))

                    if (seurat_pkg_check_result) {
                        suppressPackageStartupMessages({
                            library(Seurat)
                        })
                    }
                },
                error = function(e) {
                    warning(e)
                })

            # API-BASED APPLICATION
            } else if (identical(tolower(get_env_value("data_mode")), "api")) {
                .GlobalEnv$g_allowed_organisms <- get_env_value("allowed_organisms") %>%
                    strsplit(split = ",") %>%
                    unlist() %>%
                    trimws()

                if ((length(g_allowed_organisms) == 0) || any(g_allowed_organisms == "")) {
                    logerror("Allowed organisms list is not configured correctly or contains empty values")
                    result <- FALSE
                }

                # check that we have the API packages needed
                tryCatch({
                    scidb_pkg_check_result <- check_package_version(package_name = "scidb",
                                                                    min_version  = get_env_value("pkg_min_scidb_version"))
                },
                error = function(e) {
                    warning(e)
                })

                tryCatch({
                    revealsc_pkg_check_result <- check_package_version(package_name = "revealsc",
                                                                       min_version  = get_env_value("pkg_min_revealsc_version"))
                    if (revealsc_pkg_check_result) {
                        suppressPackageStartupMessages({
                            library(revealsc)
                        })
                    }
                },
                error = function(e) {
                    warning(e)
                })

                test_connection <- NULL
                tryCatch({
                    host <- get_env_value('revealsc_host_name')

                    tic(glue("REVEALSC::connect to {host}"))
                    test_connection <- revealsc::revealsc_connect(
                        host              = host,
                        username          = get_env_value('revealsc_user_key',  TRUE),
                        token             = get_env_value('revealsc_token_key', TRUE),
                        result_size_limit = get_numeric_config_vars(config_var = "revealsc_api_size_limit", default_value = NULL),
                        local             = FALSE)
                },
                warning = function(w) {
                    logwarn(glue("Unable to connect to database due to {w$message}"))
                },
                error = function(e) {
                    logwarn(glue("Unable to connect to database due to {e$message}"))
                })
                toc()

                if (!is.null(test_connection)) {
                    g_api_connection <<- test_connection
                } else {
                    warning('Unable to connect to API as required in this configuration')
                }
            }

        }
    }

    if (result) {
        if (identical(tolower(get_env_value("data_mode")), "file")) {
            result <- seurat_pkg_check_result && ucell_pkg_check_result
            if (as.logical(get_env_value('show_pathology_tab', unset = FALSE))) {
                result <- result && ebimage_pkg_check_result
            }
        } else if (identical(tolower(get_env_value("data_mode")), "api")) {
            if (as.logical(get_env_value('show_pathology_tab', unset = FALSE))) {
                logwarn(glue("API configurations are not compatible with the Pathology Tab"))
                result <- FALSE
            } else {
                result <- scidb_pkg_check_result && revealsc_pkg_check_result && ucell_pkg_check_result
            }
        }
    }

    result
}


#' get_config_file
#'    A convenience function to wrap around the config::get call to have one
#'    place to change the hard coded file path
#'
#' @return contents of the configuration file as a list
get_config_file <- function() {
    config::get(config     = g_configuration,
                file       = PATH_CONFIGURATION_FILE,
                use_parent = FALSE)
}


#' get_env_value
#'      takes an environment variable name and returns the value to the caller if available
#'
#' @param envvar_name - the name of the config file value to retrieve
#' @param external    - execute a Sys.getenv() on the config file value (default = FALSE)
#' @param unset       - value to use if a variable is unset (default = NULL, which returns the config value)
#'
#' @return the environment variable value as a character or NULL if it does not exist
get_env_value <- function(envvar_name, external = FALSE, unset = NULL) {
    result <- NULL
    config_file <- get_config_file()

    if (envvar_name %in% names(config_file)) {
        result <- config_file[envvar_name]
    }

    if (!is.null(result)) {
        result <- as.character(result)

        if (external) {
            result.sys <- Sys.getenv(result, unset = "INVALID")
            if (result.sys != "INVALID") {
                result <- result.sys
            }
        }
    }

    if (is.null(result) && !is.null(unset)) {
        result <- unset
    }

    result
}


#' check_required_fields
#'      Helper method to check passed required fields existence st_app_config.yaml
#'
#' @param required_fields - list of required fields (common required fields, Seurat based fields or API based fields)
#'
#' @return boolean
check_required_fields <- function(required_fields) {
    valid_required_fields <- TRUE
    missed_fields         <- c()

    config_file <- get_config_file()

    for (field in required_fields) {
        if (!(field %in% names(config_file)) ||
            is.na(config_file[[field]]) ||
            (length(config_file[[field]]) == 0) ||
            (config_file[[field]] == "")) {
            missed_fields <- c(missed_fields, field)
        }
    }

    if (length(missed_fields) > 0) {
        logerror(glue("The following required fields are missing: {glue_collapse(missed_fields, sep = ', ')}"))
        valid_required_fields <- FALSE
    }

    valid_required_fields
}


#' check_package_version
#'   Helper method to check if the passed package version is valid.
#'   Package version is valid if:
#'     - "min_version" is NULL or it is less than or equal package version
#'     - and "max_version" is NULL or it is greater than package version
#'   The function throws warnings if:
#'     - the library is not found
#'     - or min_version and max_version are both null
#'
#' @param package_name - package name to check
#' @param min_version  - package required minimum version (default = NULL)
#' @param max_version  - package required maximum version (default = NULL)
#'
#' @return boolean
check_package_version <- function(package_name,
                                  min_version = NULL,
                                  max_version = NULL) {
    valid_package_version <- FALSE
    valid_min             <- FALSE
    valid_max             <- FALSE

    if (!requireNamespace(package_name, quietly = TRUE)) {
        warning(glue("Package '{package_name}' is not installed"))
    } else if (is.null(min_version) && is.null(max_version)) {
        warning(glue("Package '{package_name}' 'min_version' and 'max_version' configurations are both NULL"))
    } else {
        package_version <- as.character(packageVersion(package_name))

        if (is.null(min_version) ||
            (compareVersion(min_version, package_version) != 1)) {
            valid_min <- TRUE
        } else {
            warning(glue("{package_name} package version '{package_version}' is below the minimum required '{min_version}'"))
        }

        if (valid_min) {
            if (is.null(max_version) ||
                (compareVersion(max_version, package_version) == 1)) {
                valid_max <- TRUE
            } else {
                warning(glue("{package_name} package version '{package_version}' is greater than or equal the maximum required '{max_version}'"))
            }
        }

        valid_package_version <- valid_min && valid_max
    }

    valid_package_version
}


#' get_top_genes
#'   Parses top gene config character value to proper integer value
#'
#' @param top_gene_config_var - either 'display_top_n_spot_genes' or 'display_top_n_sig_genes'
#' @param default_value       - default value to use if the config variable cannot be parsed as a positive integer or zero
#'
#' @return integer
get_top_genes <- function(top_gene_config_var, default_value) {
    suppressWarnings({ top_n_genes <- as.integer(get_env_value(top_gene_config_var)) })
    if ((length(top_n_genes) == 0) ||
        is.na(top_n_genes) ||
        (top_n_genes < 0)) {
        logwarn(glue("Could not parse {top_gene_config_var} with value: {get_env_value(top_gene_config_var)}.",
                     " It must be a valid positive numeric value. Setting it to default value '{default_value}'",
                     .na = "", .null = ""))
        top_n_genes <- default_value
    }

    top_n_genes
}

#' get_numeric_config_vars
#'   Parses config character value to proper numeric value
#'
#' @param config_var    - either 'display_dynamic_range_min' or 'revealsc_api_size_limit'
#' @param default_value - default value to use if the config variable cannot be parsed as a numeric value
#'
#' @return numeric
get_numeric_config_vars <- function(config_var, default_value) {
    suppressWarnings({ numeric_config_var <- as.numeric(get_env_value(config_var)) })
    if ((length(numeric_config_var) == 0) || is.na(numeric_config_var)) {
        logwarn(glue("Could not parse {config_var} with value: {get_env_value(config_var)}.",
                     " It must be a valid numeric value. Setting it to default value '{default_value}'",
                     .na = "", .null = ""))
        numeric_config_var <- default_value
    }

    numeric_config_var
}
