#' get_cached_metadata
#'     - read the cached metadata file if available
#'
#' @return NULL or data.frame
get_cached_metadata <- function() {
    result    <- NULL
    file_path <- get_env_value("metadata_cache_file")

    tryCatch({
        result <- readRDS(file_path) %>%
            filter(APPLICATION_READY)
    },
    error = function(e) {
        warning(glue("Could not read {file_path} "), e$message)
    })

    result
}


#' get_file_metadata_table
#'     build the metadata table using the cached metadata
#'
#' @param samples_data     - dataframe containing all the sample data
#' @param metadata_cols    - metadata column data frame or NULL if not found
#' @param smp_id           - character sample identifier
#' @param type             - type of table to show (ie. Secondary, Technical, Other, ...)
#'
#' @return dataframe
get_file_metadata_table <- function(samples_data, metadata_cols, smp_id, type) {
    metadata_table <- data.frame()

    parsed_all_other_metadata <- get_all_other_metadata(samples_data = samples_data,
                                                        sample_id    = smp_id)

    if (NROW(samples_data) > 0) {
        #Determine which fields we need
        if (tolower(type) == "other") {
            # if the metadata_cols is not available, everything goes into the Other dataframe
            # anything in ALL_OTHER_METADATA that is not included in metadata col
            metadata_table <- parsed_all_other_metadata %>%
                filter(!(key %in% metadata_cols$Column_Name))

            # pull the values that may be in the main dataframe
            main_df <- samples_data %>%
                filter(smp_id == ST_ID) %>%
                select("App Metadata File" = ASSOCIATED_METADATA_FILE,
                       "Sample File"       = ASSOCIATED_SAMPLE_FILE) %>%
                t() %>%
                as.data.frame() %>%
                rownames_to_column("key") %>%
                rename(value = V1)
        } else {
            metadata_col_type <- metadata_cols %>%
                filter(Type == tolower(type)) %>%
                pull(Column_Name)

            metadata_table <- parsed_all_other_metadata %>%
                filter(key %in% metadata_col_type)

            # pull the values that may be in the main dataframe
            main_df <- samples_data %>%
                filter(smp_id == ST_ID) %>%
                select(any_of(metadata_col_type)) %>%
                t() %>%
                as.data.frame() %>%
                rownames_to_column("key") %>%
                rename(value = V1)
        }

        if (NROW(main_df) > 0) {
            # Read values in extended module so metadata fields such as Clustering_color_scheme, etc.
            # take precedence over any corresponding values from the Seurat object that exist in
            # 'all_other_metadata' field
            metadata_table <- metadata_table %>%
                filter(!(key %in% main_df$key))
            metadata_table <- bind_rows(metadata_table, main_df)
        }

        metadata_table <- metadata_table %>%
            filter(!is.na(value),
                   value != "") %>%
            arrange(key)
    }

    metadata_table
}


#' get_all_other_metadata
#'    - Parses the all other metadata column into a data frame
#'
#' @param samples_data - the sample metadata cache table
#' @param sample_id    - the sample id
#'
#' @return the ALL_OTHER_METADATA column as a table
get_all_other_metadata <- function(samples_data, sample_id) {
    result <- NULL
    if ((NROW(samples_data) > 0)  &&
        !is.null(sample_id) &&
        (sample_id != "")) {
        all_other_metadata <- samples_data %>%
            filter(ST_ID == sample_id) %>%
            pull(ALL_OTHER_METADATA) %>%
            str_split(REGEX_ALL_OTHER_META_FIELD_SEP) %>%
            unlist()

        result <- data.frame(all_other_metadata) %>%
            separate(col  = "all_other_metadata",
                     into = c("key", "value"),
                     sep  = REGEX_ALL_OTHER_META_KEYVAL_SEP)
    }

    result
}
