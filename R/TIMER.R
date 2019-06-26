#' TIMER Data
#' Downloaded May 14th 2018
#' TIMER data was published by  Li B, et al. Comprehensive analyses of tumor immunity:
#' implications for cancer immunotherapy. Genome Biology. 2016;17(1):174. [PMID: 27549193]
#'
#' @docType data
#'
#' @usage data(TIMER)
#'
#' @format An object of class \code{"DataFrame"};  with 11509 rows and 7 columns
#'
#' @keywords datasets
#'
#' @references Li B, et al. Comprehensive analyses of tumor immunity:
#' implications for cancer immunotherapy. Genome Biology. 2016;17(1):174.
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/27549193}{PubMed})
#'
#'
#' @source \href{https://cistrome.shinyapps.io/timer/_w_b38e86ce/immuneEstimation.txt"}{Download TIMER data}
#'
#' @examples
#' data(TIMER)
#' summary(TIMER)
#' cor(TIMER$B_cell, TIMER$CD4_Tcell)
#' \donttest{cor(TIMER$B_cell, TIMER$CD4_Tcell)}
"TIMER"



