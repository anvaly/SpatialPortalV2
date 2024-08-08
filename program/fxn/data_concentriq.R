#' get_concentriq_visualization
#'
#' Get the high-resolution pathology image from the concentriq APIs and (optionally) draw spots on it
#'
#' @param concentriq_image_id  - the concentriq image ID
#' @param sample_spot_diameter - the spot diameter (in pixels)
#'                               NOTE: comes from @images[name]@scale.factors$spot_diameter_fullres in Seurat objects
#' @param selected_points_df   - data.frame with columns: spot_id, imagecol, imagerow, color
#'                               NOTE: if a color value is NA the outline will be black
#'                               NOTE: to remove or not draw a spot, do not include it in the data.frame
#' @param show_spot_outlines   - boolean of whether to draw spots on the image
#' @param image_height         - integer maximum height of the output image (this is the fixed dimension on plots)
#'                               returned plot will be less than or equal to this height.  Width varies accordingly
#'
#' @return NULL or an EBImage::Image object
get_concentriq_visualization <- function(concentriq_image_id,
                                         sample_spot_diameter,
                                         selected_points_df,
                                         show_spot_outlines,
                                         image_height) {
    result <- NULL

    if (!all(c("spot_id", "imagecol", "imagerow", "color") %in% colnames(selected_points_df))) {
        warning('selected_points_df must contain the following columns: "spot_id", "imagecol", "imagerow", "color"')
    }

    if (NROW(selected_points_df) < 1) {
        warning('selected_points_df has no rows!')
    }

    if ((NROW(selected_points_df) > 0) &&
        all(c("spot_id", "imagecol", "imagerow", "color") %in% colnames(selected_points_df))) {

        spot_images <- list()
        radius      <- sample_spot_diameter / 2

        if (show_spot_outlines) {
            # ensure every spot has a color specified
            selected_points_df <- selected_points_df %>%
                mutate(color = ifelse(is.na(color), "#000000", color))
        }

        result <- fetch_spot_patch(
            image_id         = concentriq_image_id,
            spot_df          = selected_points_df,
            radius_px        = radius,
            show_outline     = show_spot_outlines,
            image_height_max = as.integer(image_height))
    }

    result
}


# INTERNAL functionality - should not be called directly, may be removed in the future
fetch_spot_patch <- function(image_id,
                             spot_df,
                             radius_px,
                             show_outline,
                             image_height_max) {
    result <- NULL

    tryCatch({
        httr::set_config(httr::config(http_version = 0))

        # read in Concentriq credentials
        api_user <- get_env_value('concentriq_api_user', TRUE)
        api_pass <- get_env_value('concentriq_api_pass', TRUE)

        url_call1      <- paste0(get_env_value('concentriq_url'), "images/", image_id)
        get_image_info <- httr::GET(url_call1,
                                    httr::authenticate(api_user, api_pass, type = 'basic'),
                                    config = httr::config(ssl_verifypeer = FALSE))

        api_return1 <- httr::content(get_image_info)

        if (httr::http_error(get_image_info)) {
            message(glue("Concentriq GET image info failed [{httr::status_code(get_image_info)}]"))
        } else {
            iip_url_base <- gsub("info.json", "",
                                 api_return1$data$imageData$imageSources[[1]]$imageServerUrl)

            expand_val  <- round(radius_px*2)
            SW          <- c(min(spot_df$imagecol) - expand_val,
                             min(spot_df$imagerow) - expand_val)

            rect_width  <- max(spot_df$imagecol) - min(spot_df$imagecol) + expand_val*2+1
            rect_height <- max(spot_df$imagerow) - min(spot_df$imagerow) + expand_val*2+1

            output_height <- min(rect_height, image_height_max)

            # format: origin-x, origin-y, width, height/size-width, size-height/rotation/quality.format
            frame_coord <- paste0(SW[1], ",", SW[2], ",", rect_width, ",", rect_height,
                                  "/,", output_height,
                                  "/0",
                                  "/default.tiff")

            url_call2 <- paste0(iip_url_base, frame_coord)

            get_patch <- httr::GET(url_call2,
                                   httr::authenticate(api_user, api_pass, type = 'basic'),
                                   config = httr::config(ssl_verifypeer = FALSE))

            if (httr::http_error(get_patch)) {
                message(glue("Concentriq GET image patch failed [{httr::status_code(get_patch)}]"))
            } else {
                image <- tiff::readTIFF(httr::content(get_patch)) %>%
                    EBImage::Image(colormode = "Color")

                if (show_outline) {
                    for(i in 1:NROW(spot_df)){
                        scalar    <- output_height/rect_height
                        linewidth <- max(ceiling(7*scalar), 1)
                        for(a in c(-linewidth:linewidth)){
                            image <- EBImage::drawCircle(image,
                                                         (spot_df$imagerow[i]-SW[2])*scalar,
                                                         (spot_df$imagecol[i]-SW[1])*scalar,
                                                         radius = ceiling(scalar*radius_px) + a,
                                                         col    = spot_df$color[i])
                        }
                    }
                }
                result <- EBImage::transpose(image)
            }
        }
    }, warning = function(w) {
        warning(w)
    }, error = function(e) {
        warning(e)
    })

    result
}

