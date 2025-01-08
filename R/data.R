#' Kimmey_5hr_stim
#'
#' This data was originally published by Kimmey et al. in 2019.
#' It is human PBMC mass cytometry data from a single donor.
#' The sample underwent a 5-hour phorbol 12-myristate 13-acetate (PMA) ionomycin stimulation prior to being inputted into the cytometer.
#'
#' @format A data frame with 95492 rows and 12 variables:
#'
#'   \describe{
#'   \item{Cluster}{Cluster the cell was grouped into (1-9)}
#'   \item{CD16}{}
#'   \item{CD56}{}
#'   \item{CD123}{}
#'   \item{CD4}{}
#'   \item{CD20}{}
#'   \item{CD14}{}
#'   \item{CD19}{}
#'   \item{CD38}{}
#'   \item{CD3}{}
#'   \item{CD11c}{}
#'   \item{HLA.DR}{}
#'   }
#'
#' @source <https://doi.org/10.1038/s41467-019-09128-7a>
"Kimmey_5hr_stim"

#' Dusoswa_OMIP_54_markers
#'
#' This data was originally published by Dusoswa et al. in 2019.
#' Marker definitions are extracted from the published Optimized Multicolor Immunofluorescence Panel (OMIP), which was designed for studies using mice brain, spleen, and bone marrow samples on the mass cytometer.
#'
#' @format A data frame with 8 rows and 2 variables:
#'
#'   \describe{
#'   \item{Cluster}{Cluster the cell was grouped into (1-9)}
#'   \item{Markers}{Markers and marker signs, seperated by commas}
#'   }
#'
#' @source <https://doi.org/10.1002/cyto.a.23725>
"Dusoswa_OMIP_54_markers"

#' Lee_AML_cell_types_markers
#'
#' This cell type-marker table was originally published by Lee et al. in 2017 to test the Automated Cell-type Discovery and Classification (ACDC) algorithm.
#' That table was generated based on a gating hierarchy on Cytobank.
#' The original data was published by Levine et al. in 2015 and was collected from healthy human bone marrow.
#'
#' @format A data frame with 14 rows and 2 columns:
#' \describe{
#'   \item{Name}{Cell type name}
#'   \item{Markers}{Markers and marker signs, seperated by commas}
#' }
#'
#' @source <https://doi.org/10.1093/bioinformatics/btx054>
"Lee_AML_cell_types_markers"
