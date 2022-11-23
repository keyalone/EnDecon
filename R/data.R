#' raw count matrix of reference scRNA-seq dataset
#'
#' We obtain the filtered human breast cancers scRNA-seq dataset from 
#' [Zenodo data repository] (<https://doi.org/10.5281/zenodo.4739739>) 
#' The dataset contains the expression levels of 11920 genes and 3024 cells. 
#'
#' @name breast.sc.ref
#' @docType data
#' @usage data(breast.sc.ref)
#' @keywords datasets
#' @format a large matrix
#' @examples
#' data(breast.sc.ref)
#'
"breast.sc.ref"

#' cell type labels of scRNA-seq data
#'
#' @name breast.sc.cell.label
#' @docType data
#' @usage data(breast.sc.cell.label)
#' @keywords datasets
#' @format a vector
#' @examples
#' data(breast.sc.cell.label)
#'
"breast.sc.cell.label"


#' raw count matrix of spatial transcriptomics data
#' 
#' We obtain the filtered human breast cancers scRNA-seq dataset from 
#' [Zenodo data repository] (<https://doi.org/10.5281/zenodo.4739739>). 
#' The dataset contains the expression levels of 11920 genes and 306 spots. 
#'
#' @name breast.st
#' @docType data
#' @usage data(breast.st)
#' @keywords datasets
#' @format a large matrix
#' @examples
#' data(breast.st)
#'
"breast.st"

#' coordinate of spots for the spatially resolved transcriptomics data
#'
#' coordinate of spots for the breast.st. The center of the grids are served as the coordinates of the corresponding generated spots.
#'
#' @name breast.st.loc
#'
#' @docType data
#' @usage data(breast.st.loc)
#' @keywords datasets
#' @format a list
#' @examples
#' data(breast.st.loc)
#'
"breast.st.loc"