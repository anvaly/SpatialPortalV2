#' add_ucell_module_score
#'   Return a seurat object with the added scored name for the given features,
#'   NOTE: score only added if possible, it is not always possible to add a score
#'
#' @param seurat   - seurat object
#' @param features - feature list for scoring
#' @param name     - score name (NOTE: make.names will be called on this!)
#' @param assay    - the assay to use (default = "SCT")
#'
#' @return seurat object
add_ucell_module_score <- function(seurat, features, name, assay = SEURAT_ASSAY) {
    result <- seurat

    try({
        suppressWarnings({
            result <- UCell::AddModuleScore_UCell(obj      = seurat,
                                                  features = features,
                                                  assay    = assay)
        })
    }, silent = T)

    # due to naming differences between UCell and Seurat's Module Score output
    # we need to rename the slot from UCell regardless so are using the base name

    result@meta.data[[make.names(name)]] <- result@meta.data$signature_1_UCell
    suppressWarnings({
        result@meta.data$signature_1_UCell <- NULL
    })

    result
}


#' build_seurat_st_data_objects
#'   Construct a list of st_data_class objects from available sample files
#'
#' @param samples - data.frame of samples metadata
#'
#' @return list
build_seurat_st_data_objects <- function(samples) {
    results <- list()

    if (NROW(samples) > 0) {
        # returns st_data_class objects constructed from available sample files
        results <- apply(samples, 1, function(x) {
            st_data_class(
                st_id          = x[["ST_ID"]],
                is_api_based   = FALSE,
                name           = x[["Sample_App_Name"]],
                top_genes      = NULL,
                slice_image    = NULL,
                api_meta       = NULL,
                api_data       = NULL,
                file_path      = x[["APPLICATION_SAMPLE_FILE"]],
                file_f_seurat  = future({
                    sobj <- NULL
                    tryCatch({
                        sobj <- readRDS(x[["APPLICATION_SAMPLE_FILE"]])
                    },
                    error = function(e) {
                        message(glue("Could not read rds for sample {x[['APPLICATION_SAMPLE_FILE']]} at ",
                                     "path {x[['APPLICATION_SAMPLE_FILE']]} due to error: {e}"))
                    })
                    sobj
                }, lazy = TRUE, seed = NULL)
            )
        })
        names(results) <- samples$ST_ID
    }

    results
}


#' get_seurat_clusters
#'   Get clusters - if a cluster field is provided, those values are used.
#'   If not, then the values in `seurat@active.ident` are used, if they exist and have more than one level.
#'
#' @param seurat        - seurat object for selected sample
#' @param cluster_field - name of clustering field for selected sample
#'
#' @return data.frame or NULL
get_seurat_clusters <- function(seurat, cluster_field) {
    result <- NULL

    if (!is.null(seurat)) {
        if (!is.null(cluster_field) && !is.na(cluster_field)) {
            result <- get_seurat_data(seurat = seurat ,
                                      items  = cluster_field,
                                      slot   = SEURAT_DATA_SLOT_NAME)

            if (!is.null(result)) {
                result <- result %>%
                    rename(cluster = 1)
            }
        }

        # backup
        if (all(is.null(result),
                is.factor(seurat@active.ident),
                length(levels(seurat@active.ident)) > 1)) {
            result <- get_seurat_data(seurat = seurat,
                                      items  = "ident",
                                      slot   = SEURAT_DATA_SLOT_NAME) %>%
                rename(cluster = 1)
        }
    }
    result
}


#' get_seurat_coordinates
#'   Get the tissue coordinates
#'
#' @param seurat      - seurat object for selected sample

#' @return data.frame
get_seurat_coordinates <- function(seurat) {
    Seurat::GetTissueCoordinates(seurat)
}


#' get_seurat_data
#'   Get variable(s) from Seurat object
#'
#' @param seurat - seurat object
#' @param items  - variable(s) to fetch
#' @param slot   - slot to pull data from
#'
#' @return data.frame or NULL
get_seurat_data <- function(seurat, items, slot = NULL, assay = NULL) {
    result <- NULL

    try({
        if (is.null(slot)) {
            result <- Seurat::FetchData(object = seurat,
                                        vars   = items,
                                        assay  = assay)
        } else if (is.null(assay)) {
            result <- Seurat::FetchData(object = seurat,
                                        vars   = items,
                                        slot   = slot)
        } else {
            result <- Seurat::FetchData(object = seurat,
                                        vars   = items,
                                        slot   = slot,
                                        assay  = assay)
        }
    }, silent = T)

    result
}


#' get_seurat_pathology
#'   Get pathology data from Seurat object
#'
#' @param seurat          - seurat object for selected sample
#' @param pathology_field - name of pathology field for selected sample
#'
#' @return data.frame or NULL
get_seurat_pathology <- function(seurat, pathology_field) {
    result <- NULL

    if (!is.null(seurat) && !is.null(pathology_field) && !is.na(pathology_field)) {
        result <- get_seurat_data(seurat = seurat,
                                  items  = pathology_field,
                                  slot   = SEURAT_DATA_SLOT_NAME)

        if (!is.null(result)) {
            result <- result %>%
                rename(pathology = 1)
        }
    }

    result
}


#' get_seurat_slice_image
#'   Get the background tissue image
#'
#' @param seurat      - seurat object for selected sample
#'
#' @return raw png value
get_seurat_slice_image <- function(seurat) {
    Seurat::GetImage(seurat, mode = "raw")
}


#' get_seurat_top_represented_genes
#'   Calculate the top represented genes for each spot in a given sample
#'
#' @param data          - sample seurat object
#' @param assay         - name of assay variable in seurat object to use for gene expression
#'
#' @return data.frame
get_seurat_top_represented_genes <- function(data, assay = SEURAT_ASSAY){
    top_genes <- NULL
    tryCatch({
        assay_data <- Seurat::GetAssayData(data, assay = assay)
        if (NROW(assay_data) > 0) {
            top_genes <- assay_data %>%
                as_tibble(rownames = "Symbol") %>%
                gather(key = "Coord", val = "Value", -Symbol) %>%
                filter(Value > 1)
        }
    },
    error = function(e) {
        warning("Could not calculate top represented genes due to:", e$message)
    })

    top_genes
}


#' get_seurat_spot_diameter
#'   Get spot diameter from Seurat object
#'
#' @param seurat      - seurat object for selected sample
#' @param sample_name - name of selected sample
#'
#' @return numeric or NULL
get_seurat_spot_diameter <- function(seurat, sample_name) {
    spot_diameter <- NULL

    tryCatch({
        spot_diameter <- seurat@images[[sample_name]]@scale.factors$spot_diameter_fullres
    },
    warning = function(w) {
        message(glue("Unable to get spot diameter for sample '{sample_name}' due to: {w$message}"))
    },
    error = function(e) {
        message(glue("Unable to get spot diameter for sample '{sample_name}' due to: {e$message}"))
    })

    spot_diameter
}


#' get_seurat_image_full_coordinates
#'   Get the image full coordinates from Seurat object
#'
#' @param seurat      - seurat object for selected sample
#' @param sample_name - name of selected sample
#'
#' @return data.frame or NULL
get_seurat_image_full_coordinates <- function(seurat, sample_name) {
    full_coords <- NULL

    tryCatch({
        full_coords <- seurat@images[[sample_name]]@coordinates
    },
    warning = function(w) {
        message(glue("Unable to get full coordinates for sample '{sample_name}' due to: {w$message}"))
    },
    error = function(e) {
        message(glue("Unable to get full coordinates for sample '{sample_name}' due to: {e$message}"))
    })

    full_coords
}
