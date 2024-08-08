suppressPackageStartupMessages({
    library(shinyFeedback)
    library(RColorBrewer)
    library(dplyr)
    library(tibble)
})

source('program/fxn/plots.R', local = TRUE)

source('program/modules/stSingleBox/ui.R',      local = TRUE)
source('program/modules/stSingleBox/server.R',  local = TRUE)
source('program/modules/stSingleBox/support.R', local = TRUE)

source('program/modules/stMatrixPlot/ui.R',     local = TRUE)
source('program/modules/stMatrixPlot/server.R', local = TRUE)

source('program/modules/matrixGridToSingle/matrixGridToSingle.R', local = TRUE)
