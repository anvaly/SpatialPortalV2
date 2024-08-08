# Default color scheme to use when custom color scheme is not specified for a sample
generated_colors <- colorRampPalette(COLOR_SCHEME_DEFAULT)

#' assign_metadata_to_colors
#'     - sort pathology/cluster names based on color scheme
#'
#' @param colors   - named color vector of  the associated metadata
#' @param metadata - named vector of clusters or pathology
#'
#' @return named sorted vector
assign_metadata_to_colors <- function(colors, metadata) {
    # sort metadata values based on color scheme
    colors_order   <- names(colors)
    metadata_order <- names(metadata)
    metadata       <- metadata[metadata_order[order(match(metadata_order, colors_order))]]

    # assign colors to the available metadata values
    color_values                <- which(names(colors) %in% names(metadata))
    names(colors)[color_values] <- metadata[which(names(metadata) %in% names(colors))]
    colors
}


#' get_color_scale_breaks
#'   Calculate color spectrum scaling only
#'
#' @param signal       - dataframe        - used for color scaling based on color scaling options
#' @param scaling      - character/number - scaling should be "center", "trim" or a pair of values -- if anything else
#'                                          it is "native"
#' @param color_count  - number           - colors vector length
#'
#' @return list - color scaling values
get_color_scale_breaks <- function(signal,
                                   scaling,
                                   color_count) {
    breaks <- list()
    col.pr <- 2

    if ((is.list(scaling)) && (length(scaling) == 2)) {
        breaks <- seq(from       = scaling[[1]],
                      to         = scaling[[2]],
                      length.out = color_count)
    } else if (scaling == "trim") { # q1-q99
        q <- quantile(signal, probs = c(0, 0.01, 0.99, 1), na.rm = TRUE, names = FALSE)
        breaks <- seq(from       = q[2],
                      to         = q[3],
                      length.out = color_count)
    } else if (scaling == "center") {
        breaks <- seq(from       = 0,
                      to         = 1,
                      length.out = color_count)
    } else {
        breaks <- seq(from       = round(min(signal, na.rm = TRUE) - 5 * 10^(-2), 1),
                      to         = round(max(signal, na.rm = TRUE) + 5 * 10^(-2), 1),
                      length.out = color_count)
    }

    if (length(breaks) > 0) {
        testround <- round(breaks, digits = col.pr)
        if (length(unique(testround)) == length(testround)) {
            breaks <- testround
        }
    }

    breaks
}


#' get_download_filename
#'   Construct a file name for downloading a canvasXpress chart
#'
#' @param text - vector of identifiers to use in the file name, such as sample name and gene label
#'
#' @return character
get_download_filename <- function(text) {
    # only take the first part of anything that is multi-line
    lines <- sapply(text, function(x) { str_split(x, '\n', simplify = TRUE)[1] })

    glue(format(Sys.time(),'%Y_%m_%d_%H_%M'), ".", gsub(REGEX_DOWNLOAD_NAME_FIRST_PART, "_", paste(lines, collapse = ".")))
}


#' sort_plot_data_by_colors
#'   Re-order data.frame for plotting, to ensure the legend displays according to the specified color order
#'
#' @param data                 - data.frame with 1 column to reorder
#' @param colors               - named list of colors in desired legend order (names must match the values in the
#'                               data.frame column)
#' @param convert_to_character - logical - when TRUE, forces the column to character in the returned data.frame
#'
#' @return data.frame
sort_plot_data_by_colors <- function(data, colors, convert_to_character = TRUE) {

    # we expect that the data frame only has one column, otherwise return the original data frame
    if (NCOL(data) == 1) {
        col_name <- names(data)
    } else {
        warning("More than 1 column in the dataframe detected, using original dataframe")
        return(data)
    }

    data <- data %>%
        mutate({{ col_name }} := factor(.data[[col_name]], levels = names(colors), ordered = TRUE)) %>%
        arrange(.data[[col_name]])

    if (convert_to_character) {
        data <- data %>%
            mutate({{ col_name }} := as.character(.data[[col_name]]))
    }
    data
}


#' add_cx_large_chart_config
#'     - large plot format used in the single overview and pathology tabs
#'
#' @param cx_plot - canvasXpress plot object
#'
#' @return canvasXpress plot object
add_cx_large_chart_config <- function(cx_plot) {
    cx_plot %>%
        canvasXpress(disableToolbar            = FALSE,
                     toolbarItems              = c("Save", "History", "Table", "Explore", "Lasso", "Customize", "Maximize"),
                     showLegend                = TRUE,
                     colorByShowLegend         = TRUE,
                     showLegendTitle           = FALSE,
                     titleScaleFontFactor      = 0.7,
                     subtitleScaleFontFactor   = 0.6,
                     legendTextScaleFontFactor = 1)
}


#' add_cx_small_chart_config
#'     - small plot format used in the single extended tab
#'
#' @param cx_plot - canvasXpress plot object
#'
#' @return canvasXpress plot object
add_cx_small_chart_config <- function(cx_plot) {
    cx_plot %>%
        canvasXpress(toolbarItems              = c("Save", "Lasso", "Customize", "Maximize"),
                     showLegend                = TRUE,
                     colorByShowLegend         = TRUE,
                     showLegendTitle           = FALSE,
                     titleScaleFontFactor      = 0.8,
                     legendTextScaleFontFactor = 0.6)
}


#' add_cx_xsmall_chart_config
#'     - extra small plot format used in the matrix tab
#'
#' @param cx_plot - canvasXpress plot object
#'
#' @return canvasXpress plot object
add_cx_xsmall_chart_config <- function(cx_plot) {
    cx_plot %>%
        canvasXpress(disableToolbar            = TRUE,
                     showLegend                = FALSE,
                     colorByShowLegend         = FALSE)
}


#' add_cx_common_config
#'   Add different plots common configurations to passed CX plot
#'
#' @param cxObject - CX plot object
#'
#' @return CX Object
add_cx_common_config <- function(cxObject) {
    result <- cxObject

    if (!is.null(result) && is.list(result) && !is.null(result$x) && !is.null(result$x$config)) {

        common_config <- list(
            graphType                = "Scatter2D",
            scatterType              = "visium",
            backgroundType           = "windowImage",

            xAxisMajorTicks          = FALSE,
            xAxisMinorTicks          = FALSE,
            xAxisShow                = FALSE,
            yAxisMajorTicks          = FALSE,
            yAxisMinorTicks          = FALSE,
            yAxisShow                = FALSE,
            noValidate               = TRUE,
            printMagnification       = 3,
            selectionColor           = "#000000",

            zoomDisable              = TRUE,
            disableWheel             = TRUE,

            cX                       = g_cx_license
        )

        for (item in names(common_config)) {
            # do not override a value set in the chart configuration
            if (is.null(result$x$config[[item]])) {
                result$x$config[item] <- common_config[[item]]
            }
        }
    }
    result
}


#' get_cx_license
#'   Get configured CX package license
#'
#' @return character
get_cx_license <- function() {
    result <- ""

    try({
        result <- Sys.getenv("CX_CONFIG_LICENSE")
    }, silent = T)

    result
}


#' build_plot_custom_names
#'   Helper function to format passed slot value pairs (separated by equals and comma) into a named vector
#'
#' @param naming_string - slot value pairs (separated by equals and comma)
#' @param lower_case    - boolean to convert case of names values to lower case
#'
#' @return named vector or NULL
build_plot_custom_names <- function(naming_string, lower_case) {
    user_defined_names <- NULL

    if (!is.null(naming_string) && !is.na(naming_string) && (nchar(naming_string) > 3)) {
        # separate the user_defined_names and values into two lists
        name_order <- unlist(str_split(naming_string, ",")) %>%
            sapply(function(x) {
                user_defined_name <- unlist(str_split(x, "="))
                if (length(user_defined_name) != 2) {
                    user_defined_name[2] <- ""
                }
                user_defined_name
            })

        # clean up the values and create a named list
        # note all validation should be done in the caching script
        user_defined_names <- name_order[2, ] %>%
            str_remove_all(pattern = "'|\"")
        vector_name        <- as.character(trimws(name_order[1, ]))

        if (lower_case) {
            user_defined_names <- tolower(user_defined_names)
            vector_name        <- tolower(vector_name)
        }

        names(user_defined_names) <- vector_name
    }

    user_defined_names
}


#' build_plot_color_scheme
#'   Get the color specification object for plots
#'
#'   If a valid color string is given then this scheme is used to specify custom colors and group order
#'   * A color scheme is valid if it contains all valid hex color codes and all the data groups
#'   If color string is missing or invalid, default colors are used and no custom order is used for groups
#'
#' @param color_string - pairwise string specification for colors (equals, commas)
#' @param data_groups  - char or factor vector with all the unique groups in the selected sample
#'
#' @return list of two values:
#'           colors - char - hex color codes, optionally with names (ordered)
#'           object_colors - boolean - object_colors indicates the returned colors are from seurat object or not
build_plot_color_scheme <- function(color_string, data_groups) {
    colors        <- NULL
    object_colors <- TRUE

    # convert from factor to character so that group names are assigned correctly in the color scheme
    data_groups <- as.character(data_groups)

    try({
        # extract the order
        color_order <- unlist(str_split(color_string, ",")) %>%
            sapply(function(x) {unlist(str_split(x, "="))})

        # extract the colors
        # note all validation should be done in the caching script
        color_scheme_colors        <- color_order[2, ] %>% as.character()
        color_scheme_groups        <- color_order[1, ] %>% as.character()

        groups_without_spaces      <- gsub(REGEX_SPACE, "", color_scheme_groups) %>% tolower()
        data_groups_without_spaces <- gsub(REGEX_SPACE, "", data_groups) %>% tolower()

        if (all(length(data_groups) <= length(color_scheme_groups),
                data_groups_without_spaces %in% groups_without_spaces)) {
            # construct final groups colors names to be exactly as data groups names in case it exists
            final_groups <- sapply(color_scheme_groups, function(group) {
                group_index <- match(gsub(REGEX_SPACE, "", group) %>% tolower(), data_groups_without_spaces)
                if (is.na(group_index)) {
                    group
                } else {
                    data_groups[group_index]
                }
            })
            colors        <- color_scheme_colors
            names(colors) <- final_groups
        }
    }, silent = T)

    if (is.null(colors)) {
        data_groups   <- str_sort(data_groups, numeric = TRUE)
        colors        <- generated_colors(length(data_groups))
        names(colors) <- unique(data_groups)
        object_colors <- FALSE
    }

    list(colors = colors, object_colors = object_colors)
}


#' get_spot_scaling
#'     Helper function to calculate spot scaling factor given tissue image dimensions
#'
#' @param image_dim - dimensions of the background image plot (result of dim())
#'
#' @return numeric
get_spot_scaling <- function(image_dim) {
    result <- 1.4 * min(image_dim[1]/image_dim[2],
                        image_dim[2]/image_dim[1])

    max(1, floor(10*result)/10)
}


#' get_all_spots_colors
#'   Get the colors of all spots from the plotted data
#'
#' @param coordinates           - data.frame of spot coordinates
#' @param plot_data             - data.frame of plotted data
#' @param column_name           - name of column in plot_data that will be mapped to colors
#' @param color_key             - color key, for categorical data column
#' @param color_spectrum_breaks - color spectrum breaks, for numeric data column
#'
#' @return data.frame or NULL
get_all_spots_colors <- function(coordinates,
                                 plot_data,
                                 column_name,
                                 color_key             = NULL,
                                 color_spectrum_breaks = NULL) {

    spot_colors <- NULL

    if ((NROW(coordinates) > 0) && (NROW(plot_data) > 0)) {
        if (!is.null(color_key)) {
            # use color key for categorical data
            color_key_df <- data.frame(value = names(color_key),
                                       color = color_key)

            plot_data <- plot_data %>%
                rownames_to_column(var = "spot_id") %>%
                mutate(value = as.character(.data[[column_name]])) %>%
                select(spot_id, value)

            spot_colors <- coordinates %>%
                rownames_to_column(var = "spot_id") %>%
                left_join(plot_data, by = "spot_id") %>%
                left_join(color_key_df, by = "value")
        } else if (is.numeric(plot_data[[column_name]]) && !is.null(color_spectrum_breaks)) {
            # use color spectrum for numeric data
            plot_data <- plot_data %>%
                rownames_to_column(var = "spot_id") %>%
                mutate(color = NA_character_) %>%
                select(spot_id,
                       setNames(column_name, nm = "signal"),
                       color)

            min_break <- min(color_spectrum_breaks)
            max_break <- max(color_spectrum_breaks)

            if (!is.na(min_break) && !is.na(max_break) && (max_break > min_break)) {
                plot_data <- plot_data %>%
                    mutate(trimmed_signal = case_when(signal < min_break ~ min_break,
                                                      signal > max_break ~ max_break,
                                                      TRUE               ~ signal),
                           scaled_signal = (trimmed_signal - min(trimmed_signal, na.rm = TRUE)) / (max_break - min_break))
            } else {
                plot_data <- plot_data %>%
                    mutate(scaled_signal = signal)
            }

            colors_df1 <- plot_data %>%
                filter(!is.na(scaled_signal))

            tryCatch({
                colors_df1 <- colors_df1 %>%
                    mutate(color = rgb(colorRamp(COLOR_SCHEME_SPECTRUM_SCALING)(scaled_signal), maxColorValue = 255))
            },
            warning = function(w) {
                warning(w$message)
            },
            error = function(e) {
                warning(e$message)
            })

            colors_df2 <- plot_data %>%
                filter(is.na(scaled_signal))

            colors_df_merged <- bind_rows(colors_df1, colors_df2)

            spot_colors <- coordinates %>%
                rownames_to_column(var = "spot_id") %>%
                left_join(colors_df_merged, by = "spot_id")
        }
    }

    if (NROW(spot_colors) > 0) {
        spot_colors <- spot_colors %>%
            mutate(color = replace_na(color, MISSING_DATA_COLOR)) %>%
            select(spot_id, imagecol, imagerow, color)
    }

    spot_colors
}
