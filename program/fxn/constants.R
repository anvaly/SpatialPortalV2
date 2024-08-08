#### Application Constant Variables File ####
# Put any constant definition used in the app in that file
# Constant variable format is: all caps, separated by underscores
# START the variables with prefixes indicating usage


# UI related constants
PLACE_HOLDER_SELECTIZE_INPUTS <- "Type/Click then Select"


## Matrix Tab UI Constants
OPTIONS_SPATIAL_PLOT_SCALING  <- c("Individually"      = "individually",
                                   "By Sample"         = "bySample",
                                   "By Signature/Gene" = "bySignatureGene",
                                   "All Together"      = "allTogether")
OPTIONS_SPATIAL_PLOT_COLORING <- c("Native"  = "native",
                                   "Trimmed" = "trim",
                                   "Fixed"   = "center")


## Single Tab UI Constants
OPTIONS_OVERVIEW_SPATIAL_PLOT_COLORING <- c(OPTIONS_SPATIAL_PLOT_COLORING, "Custom"  = "custom")
PLACE_HOLDER_CUSTOM_GENE_INPUT         <- "Select or Add a Custom Gene List"


## Pathology Tab UI Constants
OPTIONS_PATHOLOGY_METADATA_API_EXCLUDE <- c("Coord") # choices excluded from Metadata dropdown input in API-based app
HEIGHT_SCREEN_PATHOLOGY_IMAGE          <- "500"      # in pixels
HEIGHT_DOWNLOAD_PATHOLOGY_IMAGE        <- "1440"     # in pixels
MAX_ALLOWED_SELECTED_PATHOLOGY_SPOTS   <- 50


# IDs prefix/suffix literals
SUFFIX_SPATIAL_PLOT_ID          <- "bkg"
SUFFIX_ADVANCED_SPATIAL_PLOT_ID <- "adv_bkg"


# Files paths
PATH_WALKTHROUGH_FILE   <- "program/data/walkthrough_data.csv"
PATH_CONFIGURATION_FILE <- "program/data/st_app_config.yaml"
PATH_ABOUT_FILE         <- "program/data/about.html"


# Regular expressions
REGEX_USER_COLUMNS               <- "^user\\."
REGEX_MESSAGE_ENCODED_PARAMETERS <- "(?<=\\{)[^{}]+(?=\\})"
REGEX_SPACE                      <- "[[:space:]]"
REGEX_ORGANISM_TISSUE_SEP        <- " > "
REGEX_DOWNLOAD_NAME_FIRST_PART   <- ", |,| "
REGEX_ALL_OTHER_META_FIELD_SEP   <- "\\|\\|\\|"
REGEX_ALL_OTHER_META_KEYVAL_SEP  <- ":::"


# Files extensions
CX_BKG_IMG_FILE_EXT    <- ".png"
PATHOLOGY_DOWNLOAD_EXT <- ".tiff"


# Color schemes
COLOR_SCHEME_DEFAULT          <- c("navy", "blue", "cyan", "darkgreen", "green",  "yellow", "orange", "brown", "red", "magenta", "purple", "gray")
COLOR_SCHEME_SPECTRUM_SCALING <- list("#5E4FA2", "#3E96B7", "#88CFA4", "#D7EF9B", "#FFFFBF", "#FDD380", "#F88D52", "#DC494C", "#9E0142")
MISSING_DATA_COLOR            <- "#F0F0F0"


# meta data items
METADATA_ITEM_PATHOLOGY       <- "Pathology.Group"
METADATA_ITEM_CLUSTERING      <- "Clustering"
METADATA_ITEM_TISSUE          <- "Tissue Image"
SEURAT_DATA_SLOT_NAME         <- "data"
SEURAT_ASSAY                  <- "SCT"


# Required fields
CONFIG_COMMON_REQUIRED_FIELDS <- c(
    "dir_metadata_files"         # contains application required metadata files
)


CONFIG_SEURAT_REQUIRED_FIELDS <- c(
    "metadata_signatures_file",  # signatures metadata file
    "pkg_min_Seurat_version",    # Seurat package minimum required version
    "metadata_cache_file"        # metadata cache
)


CONFIG_API_REQUIRED_FIELDS    <- c(
    "allowed_organisms",         # contains application allowed organisms (default = "human, mouse, rat")
    "revealsc_host_name",
    "revealsc_user_key",
    "revealsc_token_key",
    "pkg_min_scidb_version",
    "pkg_min_revealsc_version",
    "api_counts_type"
)

CONFIG_API_SAMPLES_REQUIRED_COLUMNS <- c("Sample_ID", "Sample_Name", "Organism", "Tissue", "Protocol")
