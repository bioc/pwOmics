#' Omics example dataset.
#'
#' A dataset as input example for readOmics: A list containing two input lists,
#' one for protein data, one for gene data, both including a vector of all 
#' measured IDs as first element and a list as second element including for 
#' each tp a dataframe with IDs and log foldchanges per timepoint.
#'
#' @docType data
#' @keywords datasets
#' @name OmicsExampleData
#' @usage data(OmicsExampleData)
#' @return List with 2 sublists.
#' @format A list with a 'P' sublist containing protein data and a 'G' sublist
#' containing gene/transcript data. Each of the sublists has a first element 
#' with all measured protein/gene IDs and a second element with a list of the 
#' length of the number of measured time points. Each of these lists contains 
#' a dataframe with a first column of significant protein/gene IDs at that time 
#' point and a second column with the corresponding logFCs.
NULL