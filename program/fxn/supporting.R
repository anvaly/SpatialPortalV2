# Please only put functions that really don't have any other place to go here

st_data_class <- setClass("st_data",
                          slots = list(is_api_based   = "logical",
                                       name           = "character",
                                       st_id          = "character",
                                       # common data objects
                                       top_genes      = "ANY",
                                       slice_image    = "ANY",
                                       # api-based information
                                       api_meta       = "ANY",
                                       api_data       = "ANY",
                                       # file-based information
                                       file_path      = "character",
                                       file_f_seurat  = "ANY"
                          ))

#' get_app_password
#'   Get app password from password file. If file cannot be found or cannot be read, stop the app.
#'
#' @return char
get_app_password <- function() {
    app_password <- NULL
    tryCatch({
        filename <- get_env_value("password_file")
        if (is.null(filename) ||
            filename == ""    ||
            !file.exists(filename)) {
            stopApp("Current 'password_file' setup is invalid or file does not exist. Please check your setup")
        } else {
            app_password <- readLines(filename, 1)
            if (any(is.null(app_password),
                    length(app_password) == 0,
                    app_password == "")) {
                stopApp("Current 'password_file' points to an empty or corrupted file. Please check your setup")
            }
        }

    },
    error = function(e) {
        stopApp(glue("Unable to read application password from file: ", get_env_value("password_file"), " due to: ", e$message, .null = "NULL"))
    })

    app_password
}


#' get_full_sample_data_objects
#'   Filter loaded samples objects to return selected sample objects
#'
#' @param selected_samples - ST_IDs of selected samples
#' @param all_data_objects - all current loaded data objects reactive variable
#'
#' @return selected sample data object
get_full_sample_data_objects <- function(selected_samples, all_data_objects) {
    result          <- list()
    updated_st_list <- all_data_objects()

    if (length(selected_samples) > 0) {
        # test the kind of samples (api/file) using the first item (mixed use is NOT supported!)
        if (updated_st_list[[selected_samples[1]]]@is_api_based) {
            # api-based: vectorized

            samples_to_update <- updated_st_list[which(names(updated_st_list) %in% selected_samples)]
            updated_samples   <- get_api_sample_data(samples_to_update)

            for (smp in updated_samples) {
                updated_st_list[[smp@st_id]] <- smp
            }
            result <- updated_samples
        } else {
            # file-based: must process 1 by 1

            for (st_id in selected_samples) {
                selected_sample <- updated_st_list[[st_id]]

                if (is.null(selected_sample@top_genes)) {
                    selected_sample@top_genes <- get_seurat_top_represented_genes(value(selected_sample@file_f_seurat))
                }

                if (is.null(selected_sample@slice_image)) {
                    selected_sample@slice_image <- get_seurat_slice_image(value(selected_sample@file_f_seurat))
                }

                updated_st_list[[st_id]] <- selected_sample
                result[[st_id]]          <- selected_sample
            }
        }

    }

    # update the cache with the new data
    all_data_objects(updated_st_list)

    result
}


#' get_organism_tissue
#'   Get a character vector of organism, tissue for the selected dataset
#'
#' @param datasetSel - currently selected dataset
#'
#' @return character
get_organism_tissue <- function(datasetSel) {
    str_split(datasetSel, pattern = REGEX_ORGANISM_TISSUE_SEP) %>% unlist()
}


#' get_sample_genes
#'  Get genes list for the selected sample (API or Seurat-based)
#'
#' @param dataset_sel        - name of selected dataset
#' @param sample_sel         - ST_ID of selected sample
#' @param st_objects         - reactive variable containing list of st_data_class objects for all available samples
#' @param initial_genes_list - list of genes for all available organisms and references (Seurat-based app) or NULL (API-based app)
#' @param sample_reference   - reference column value of the selected sample
#'
#' @return character vector
get_sample_genes <- function(dataset_sel,
                             sample_sel,
                             st_objects,
                             initial_genes_list = NULL,
                             sample_reference   = NULL) {
    sample_genes <- character(0)

    if (is.null(initial_genes_list)) {
        # API-based app - genes list for this sample comes from the API
        s_obj <- st_objects()[[sample_sel]]

        if (!is.null(s_obj) && !is.null(s_obj@api_data) && !is.null(s_obj@api_data$features)) {
            # feature data already retrieved
            sample_genes <- unique(s_obj@api_data$features$gene_symbol)
        } else if (!is.null(s_obj)) {
            if (is.null(s_obj@api_data$features)) {
                #NOTE - this is always singular, thus the unlist
                s_obj <- add_api_features(samples = list(s_obj))[[1]]
            }

            if (!is.null(s_obj@api_data$features)) {
                sample_genes <- s_obj@api_data$features %>%
                    pull(gene_symbol) %>%
                    unique()

                # update sample object in reactive stored list with api updated features
                updated_st_objects <- st_objects()
                updated_st_objects[[which(names(updated_st_objects) == sample_sel)]] <- s_obj
                st_objects(updated_st_objects)
            }
        }
    } else {
        # Seurat-based app - extract the applicable subset for this sample from the initial genes list
        if ((length(sample_reference) > 0) &&
            !is.na(sample_reference) &&
            (sample_reference %in% names(initial_genes_list))) {
            # sample has reference gene list
            sample_genes <- initial_genes_list[[sample_reference]]
        } else {
            s_dataset <- ""

            # sample doesn't have reference gene list, use organism genes
            if ((length(sample_sel) > 0) && (sample_sel != "")) {
                s_dataset <- dataset_sel
            }

            sample_genes <- initial_genes_list[[get_organism_tissue(s_dataset)[1]]]
        }
    }

    sample_genes
}


#' get_sample_signals
#'   Get given sample signals for genes and/or signatures
#'
#' @param st_ids           - sample st_ids for selected samples
#' @param all_data_objects - all current loaded data objects reactive variable
#' @param signatures       - selected signature objects (default = NULL)
#' @param genes            - selected genes (default = NULL)
#'
#' @return list of signals data for each sample and ST_IDs of any samples whose data could not be retrieved due to API error
get_sample_signals <- function(st_ids, all_data_objects, signatures = NULL, genes = NULL) {
    signals       <- list()
    error_samples <- character(0) # any samples that have an API error retrieving scores
    all_samples   <- all_data_objects()

    if (length(st_ids) > 0) {
        for (sample_index in 1:length(st_ids)) {
            sample        <- all_samples[[st_ids[sample_index]]]
            sample_signal <- data.frame()

            if (NROW(signatures) > 0) {
                #TODO - rewrite to get multiple signatures at once from the api
                #       after bugs are resolved w/r/t invalid genes - currently will not score
                for (signature_index in 1:NROW(signatures)) {
                    signature <- signatures$Name[signature_index]
                    features  <- list(signatures$Genes[signature_index][[1]])

                    if (sample@is_api_based) {
                        #TODO - remove validity check once bugs fixed, invalid genes cause scoring issues
                        valid_genes <- intersect(unlist(features), unique(sample@api_data$features$gene_symbol))

                        if (length(valid_genes) == 0) {
                            signal <- NULL
                            logwarn(glue("Selected gene(s) {glue_collapse(unlist(features), sep = ', ')} not found in sample {sample@name}"))
                        } else {
                            signal <- get_api_ucell_score(st_id = sample@st_id, features = list(valid_genes), all_data_objects)

                            if (is.null(signal)) {
                                # if the returned data is null even though we only requested valid genes for this sample,
                                # then the API retrieval failed
                                error_samples <- c(error_samples, sample@st_id)
                                logwarn(glue("Unable to retrieve scores for sample {sample@name}"))
                            } else {
                                signal <- signal %>% rename_with(~signature, "score")
                            }
                        }
                    } else {
                        signal_name   <- glue("sig.{signature}")
                        seurat.scored <- add_ucell_module_score(seurat   = value(sample@file_f_seurat),
                                                                features = features,
                                                                name     = signal_name)
                        signal <- get_seurat_data(seurat = seurat.scored,
                                                  items  = glue("{signal_name}"),
                                                  slot   = SEURAT_DATA_SLOT_NAME)
                    }

                    if (NROW(sample_signal) == 0) {
                        sample_signal <- signal
                    } else {
                        sample_signal           <- base::merge(sample_signal, signal, by = 0, all = TRUE)
                        rownames(sample_signal) <- sample_signal$Row.names
                        sample_signal$Row.names <- NULL
                    }
                }
            }

            if (!(sample@st_id %in% error_samples) && (NROW(genes) > 0)) {
                #TODO - rewrite to get multiple genes at once from the api
                #       after bugs are resolved w/r/t invalid genes
                for (gene_index in 1:length(genes)) {
                    gene          <- str_trim(genes[gene_index])
                    features      <- list(gene)

                    if (sample@is_api_based) {
                        #TODO - remove validity check once bugs fixed, invalid genes cause issues
                        if (!(gene %in% sample@api_data$features$gene_symbol)) {
                            signal <- NULL
                            logwarn(glue("Selected gene {gene} not found in sample {sample@name}"))
                        } else {
                            signal <- get_api_ucell_score(st_id = sample@st_id, features = features, all_data_objects)

                            if (is.null(signal)) {
                                # if the returned data is null even though we requested a valid gene for this sample,
                                # then the API retrieval failed
                                error_samples <- c(error_samples, sample@st_id)
                                logwarn(glue("Unable to retrieve module scores for sample {sample@name}"))
                            } else {
                                signal <- signal %>% rename_with(~gene, "score")
                            }
                        }
                    } else {
                        signal <- get_seurat_data(seurat = value(sample@file_f_seurat),
                                                  items  = unlist(features),
                                                  assay  = SEURAT_ASSAY)
                    }

                    if (NROW(sample_signal) == 0) {
                        sample_signal <- signal
                    } else {
                        sample_signal           <- base::merge(sample_signal, signal, by = 0, all = TRUE)
                        rownames(sample_signal) <- sample_signal$Row.names
                        sample_signal$Row.names <- NULL
                    }
                }
            }

            if (!(sample@st_id %in% error_samples)) {
                signals[[sample@st_id]] <- sample_signal
            }
        }
    }

    list(signals       = signals,
         # provide the ST_IDs of samples that were excluded due to API error (to distinguish from samples
         # that are missing from the signals data because none of the requested gene(s) were found)
         error_samples = unique(error_samples))
}


#' get_all_selected_spots
#'     Helper function to get all selected spots given a list of one or more spots
#'
#' @param selected_point_ids  - IDs for selected spots (one or more)
#' @param full_coordinates    - the coordinates data frame for the sample (NOT FILTERED)
#'                            - rownames: IDs
#'                            - columns: imagerow, imagecol
#'
#' @return list of spot IDs or NULL
get_all_selected_spots <- function(selected_point_ids, full_coordinates) {
    result <- NULL

    if (NROW(selected_point_ids) == 0 || (NROW(full_coordinates) == 0)) {
        # error
        #TODO
    } else if (NROW(selected_point_ids) == 1) {
        # expand size to a default minimum
        #TODO
        warning('single spot not yet implemented, returning NULL')
    } else {
        sel_spots <- full_coordinates %>%
            rownames_to_column('spot_id') %>%
            filter(spot_id %in% selected_point_ids)

        all_sel_spots <- full_coordinates %>%
            filter(imagerow >= min(sel_spots$imagerow) & imagerow <= max(sel_spots$imagerow),
                   imagecol >= min(sel_spots$imagecol) & imagecol <= max(sel_spots$imagecol)) %>%
            rownames()

        if (NROW(all_sel_spots) > 0) {
            result <- all_sel_spots
        }
    }

    result
}


#' read_about
#'   Read HTML content for "About This Application" box in the UI
#'
#' @return html
read_about <- function() {
    assert_that(file.exists(PATH_ABOUT_FILE),
                msg = "About file does not exist")
    HTML(readLines(PATH_ABOUT_FILE))
}


#' read_metadata_column_control
#'   Get the column categories from the designated file
#'
#' @param samples_df samples dataframe to be used to compare the metadata columns
#'
#' @return data.frame or NULL
read_metadata_column_control <- function(samples_df) {
    metadata_columns_file <- glue("{get_env_value('dir_metadata_files')}/{get_env_value('metadata_col_type_file')}")
    data                  <- NULL

    if (NROW(samples_df) > 0) {
        if (file.exists(metadata_columns_file)) {
            data <- suppressMessages(read_tsv(metadata_columns_file,
                                                col_types = "cc",
                                                trim_ws   = TRUE,
                                                comment   = get_env_value("files_comment_character"))) %>%
                mutate(Column_Name_Lower = tolower(Column_Name),
                       Type              = tolower(Type)) %>%
                filter(Column_Name_Lower != "sample_app_name") # avoid showing the Sample_App_Name

        } else {
            message(glue("Metadata columns file {metadata_columns_file} does not exist"))
        }
    }
    data
}


#' read_organism_and_probes_genes
#'   Construct a list of genes for each available organism and samples reference gene list
#'
#' @param samples           - data.frame of samples metadata
#' @param limit_to_organism - limit read to organism gene files only
#'
#' @return list
read_organism_and_probes_genes <- function(samples, limit_to_organism) {
    genes     <- list()
    organisms <- samples$Organism %>% unique()

    for (organism in organisms) {
        gene_file <- glue("{get_env_value('dir_metadata_files')}/{get_env_value('metadata_gene_file_prefix')}{tolower(organism)}.txt")
        assert_that(file.exists(gene_file),
                    msg = glue("Genes file {gene_file} does not exist for organism: {organism}"))
        genes[organism] <- suppressMessages(read_tsv(gene_file, comment = get_env_value("files_comment_character")))
    }

    if (!limit_to_organism && ("Reference" %in% names(samples))) {
        samples_with_probes <- samples %>%
            filter(!is.na(Reference), Reference !=  "")

        if (NROW(samples_with_probes) > 0) {
            probe_files <- unique(samples_with_probes$Reference)
            probes      <- list()

            for (pf in probe_files) {
                p_file_name <- glue("{get_env_value('dir_metadata_files')}/{get_env_value('metadata_gene_file_prefix')}{pf}.txt")

                if (file.exists(p_file_name)) {
                    probes[pf] <- suppressMessages(read_tsv(p_file_name, comment = get_env_value("files_comment_character")))
                }
            }

            for (idx in 1:NROW(samples_with_probes)) {
                sample_reference <- samples_with_probes[[idx, "Reference"]]

                if (sample_reference %in% names(probes)) {
                    genes[sample_reference] <- probes[sample_reference]
                } else {
                    message(glue("Reference gene list file {sample_reference} does not exist"))
                }
            }
        }
    }

    genes
}


#' read_signatures
#'   Construct a data.frame of genes for each gene signature
#'
#' @return data.frame
read_signatures <- function() {
    signatures_file <- glue("{get_env_value('dir_metadata_files')}/{get_env_value('metadata_signatures_file')}")
    assert_that(file.exists(signatures_file),
                msg = glue("Signatures file {signatures_file} does not exist"))

    data <- suppressMessages(read_tsv(signatures_file,
                                      col_types = "cc",
                                      trim_ws   = TRUE,
                                      comment   = get_env_value("files_comment_character")))
    data %>%
        rename(Name      = Signature_name,
               Gene_List = Gene_list) %>%
        mutate(Name      = make.names(Name, unique = T),
               Genes     = strsplit(Gene_List, ',', fixed = T),
               Gene_List = gsub(',', ', ', Gene_List, fixed = T))
}


#' remove_shiny_inputs
#'   Remove spatial visualization box from session's input object. This function is called when the user closes a box.
#'
#' @param id    - id of box to remove
#' @param input - shiny session's input object
#'
#' @return list
remove_shiny_inputs <- function(id, input) {
    lapply(grep(id, names(input), value = TRUE), function(i) {
        .subset2(input, "impl")$.values$remove(i)
    })
}


#' write_png_bkg_image
#'   Write passed tissue background image for selected sample into a temp png file.
#'   The image will be used later using renderImage for visualization.
#'
#' @param raw_img     - Tissue background image for selected sample
#' @param target_file - Temp file to be used for storing image raw data
#'
write_png_bkg_image <- function(raw_img, target_file) {
    write_img <- raw_img

    if (dim(write_img)[1] < 600) {
        add_rows <- array(1, dim = c(600 - dim(write_img)[1], 600, 3))
        write_img  <- abind::abind(write_img, add_rows, along = 1)
    }

    if (dim(write_img)[2] < 600) {
        add_cols <- array(1, dim = c(600, 600 - dim(write_img)[2], 3))
        write_img <- abind::abind(write_img, add_cols, along = 2)
    }

    png::writePNG(image  = write_img,
                  target = target_file)
}
