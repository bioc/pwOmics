#' Get Omics timepoints.
#'
#' This function returns the timepoints of the OmicsData.
#'  
#' @param data_omics OmicsData object.
#' @return list of protein time points and gene time points; in case of single 
#' time point measurement experiment number entered in OmicsData object.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' getOmicsTimepoints(data_omics)
getOmicsTimepoints <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    timepoints = data_omics[[1]][[1]][[1]]
    return(timepoints)
}

#' Get all protein IDs.
#'
#' This function returns the protein IDs of all proteins measured.
#'  
#' @param data_omics OmicsData object.
#' @return all protein IDs.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' getOmicsallProteinIDs(data_omics)
getOmicsallProteinIDs <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    allprotIDs = unique(data_omics[[1]][[1]][[2]])
    colnames(allprotIDs) = "all.protein.IDs"
    return(allprotIDs)
}


#' Get all gene IDs.
#'
#' This function returns the gene IDs of all genes (transcripts) measured.
#'  
#' @param data_omics OmicsData object.
#' @return all gene IDs.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' getOmicsallGeneIDs(data_omics)
getOmicsallGeneIDs <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    allgeneIDs = unique(data_omics[[1]][[1]][[3]])
    colnames(allgeneIDs) = "all.gene.IDs"
    return(allgeneIDs)
}


#' Get Omics dataset.
#'
#' This function returns the omics datasets of the experiment.
#'  
#' @param data_omics OmicsData object.
#' @param writeData boolean value indicating if datasets should be written into 
#' csv file.
#' @return list with protein data set as first element and gene data set as
#' second element; both elements contain matrices with significant proteins/
#' genes/transcripts at the given time points.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' getOmicsDataset(data_omics)
getOmicsDataset <- function(data_omics, writeData = FALSE) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    lenpdata = vector()
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    {lenpdata[plen] = dim(data_omics[[1]][[2]][[1]][[plen]])[1]}
    matrix_prot = matrix(ncol = 2*length(data_omics[[1]][[1]][[1]][[1]]), 
                         nrow = max(lenpdata))
    lengdata = vector()
    for(glen in 1: length(data_omics[[1]][[1]][[1]][[2]]))
    {lengdata[glen] = dim(data_omics[[1]][[2]][[2]][[glen]])[1]}
    matrix_gen = matrix(ncol = 2*length(data_omics[[1]][[1]][[1]][[2]]), 
                        nrow = max(lengdata))
    
    for(k in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    {matrix_prot[1: dim(data_omics[[1]][[2]][[1]][[k]])[1],(k*2-1):(k*2)] = 
         as.matrix(data_omics[[1]][[2]][[1]][[k]][,1:2])}
    for(j in 1: length(data_omics[[1]][[1]][[1]][[2]]))
    {matrix_gen[1: dim(data_omics[[1]][[2]][[2]][[j]])[1],(j*2-1):(j*2)] = 
         as.matrix(data_omics[[1]][[2]][[2]][[j]][,1:2])}
    
    colnames(matrix_prot) = rep(names(data_omics[[1]][[2]][[1]]),each = 2)
    colnames(matrix_gen) = rep(names(data_omics[[1]][[2]][[2]]),each = 2)
    if(writeData == TRUE)
    {write.csv(matrix_prot, "protein_dataset.csv")
     write.csv(matrix_gen, "gene_dataset.csv")}
    
    return(list(ProteinDataset = matrix_prot, GeneDataset = matrix_gen))
}

#' Get downstream analysis pathways.
#'
#' This function returns pathways identified in the downstream analysis 
#' containing the significantly abundant proteins.
#'  
#' @param data_omics OmicsData object.
#' @return list of length = number of protein time points, each element
#' containing a data frame with the pathway IDs in the generated biopax model,
#' corresponding pathway names and flags if those pathways are enriched
#' (1 = enriched, NA = not enriched).
#' @keywords manip
#' @export
#' @examples
#' #please run with whole database files (prepared according to vignette)
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("chea", "pazar"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWTFTGs(data_omics, only_enriched = FALSE)
#' getDS_PWs(data_omics)
#' }
getDS_PWs <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    DS_PWs = list()
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    {
        no_de_prots = dim(data_omics[[1]][[2]][[1]][plen][[1]])[1]
        PW_tp_NA = vector()
        for(g in 1: no_de_prots)
        {PW_tp_NA[g] = is.na(data_omics[[1]][[3]][[1]][[plen+1]][[g]][[1]][1])}
        temp_omics = data_omics[[1]][[3]][[1]][[plen+1]][1:no_de_prots]
        tps_PWs = rbindlist(temp_omics[!PW_tp_NA])
        DS_PWs[[plen]] = as.data.frame(tps_PWs)[,1:3]
    }
    names(DS_PWs) = names(data_omics[[1]][[2]][[1]])
    
    return(DS_PWs)
}

#' Get downstream analysis transcription factors in pathways.
#'
#' This function returns the genes identified in the downstream analysis and 
#' a column indicating if the genes are transcription factors. 
#'  
#' @param data_omics OmicsData object.
#' @return list of length = number of protein time points, each element
#' containing a character vector with identified transcription factors.
#' @keywords manip
#' @export
#' @examples
#' #please run with whole database files (prepared according to vignette)
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("chea", "pazar"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWTFTGs(data_omics, only_enriched = FALSE)
#' getDS_TFs(data_omics)
#' }
getDS_TFs <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    DS_TFs = list()
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    {
        len_omics = length(data_omics[[1]][[3]][[1]][[plen+1]])
        DS_TFs[[plen]] = data_omics[[1]][[3]][[1]][[plen+1]][[len_omics-1]]
        if(is.na(DS_TFs[[plen]][[1]][1]))
        {message("No transcription factors could be identified for time point ",
                 data_omics[[1]][[1]][[1]][[1]][plen],
                 " in downstream analyis. \n", sep = "")
        }else{
            if(length(as.character(DS_TFs[[plen]][,1][which(DS_TFs[[plen]][,2] == 1)])) == 0)
            {message("No transcription factors could be identified for time point ",
                     data_omics[[1]][[1]][[1]][[1]][plen],
                     " in downstream analyis. \n", sep = "")
            }else{
                DS_TFs[[plen]] = 
                    as.character(DS_TFs[[plen]][,1][which(DS_TFs[[plen]][,2] == 1)])}
        }
    }
    names(DS_TFs) = names(data_omics[[1]][[2]][[1]])
    
    return(DS_TFs)
}

#' Get downstream analysis target genes of TFs found in pathways.
#'  
#' @param data_omics OmicsData object.
#' @return list of length = number of protein time points, each element
#' containing a character vector with identified target genes.
#' @keywords manip
#' @export
#' @examples
#' #please run with whole database files (prepared according to vignette)
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("chea", "pazar"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWTFTGs(data_omics, only_enriched = FALSE)
#' getDS_TGs(data_omics)
#' }
getDS_TGs <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    DS_TGs = list()
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    {
        len_omics = length(data_omics[[1]][[3]][[1]][[plen+1]])
        DS_TGs[[plen]] = 
            as.character(data_omics[[1]][[3]][[1]][[plen+1]][[len_omics]])
    }
    names(DS_TGs) = names(data_omics[[1]][[2]][[1]])
    
    return(DS_TGs)
}


#' Get upstream TFs.
#'  
#' @param data_omics OmicsData object.
#' @return list of length = number of gene/transcript time points, each element
#' containing a data frame with upstream transcription factors and flag if
#' these transcription factors are enriched (1 = enriched, NA = not enriched).
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyRsofTFs(data_omics, only_enriched = FALSE, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' getUS_TFs(data_omics)
#' }
getUS_TFs <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    US_TFs = list()
    for(glen in 1: length(data_omics[[1]][[1]][[1]][[2]]))
    {
        no_de_genes = dim(data_omics[[1]][[2]][[2]][glen][[1]])[1]
        TF_tp_NA = vector()
        for(g in 1: no_de_genes)
        {TF_tp_NA[g] = is.na(data_omics[[1]][[3]][[2]][[glen+1]][[g]][[1]][1])}
        temp_omics = data_omics[[1]][[3]][[2]][[glen+1]][1:no_de_genes]
        tps_TFs = rbindlist(temp_omics[!TF_tp_NA])
        US_TFs[[glen]] = as.data.frame(tps_TFs)
    }
    names(US_TFs) = names(data_omics[[1]][[2]][[2]])
    
    return(US_TFs)
}

#' Get upstream pathways of identified transcription factors.
#'  
#' @param data_omics OmicsData object.
#' @return list of length = number of gene/transcript time points, each element
#' containing a list of transcription factors; these transcription factor
#' elements contain data frame with pathway IDs and pathway names.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyRsofTFs(data_omics, only_enriched = FALSE, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' getUS_PWs(data_omics)
#' }
getUS_PWs <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    US_PWs = list()
    for(glen in 1: length(data_omics[[1]][[1]][[1]][[2]]))
    {
        length_list_TFs = 
            length(data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])-1]])
        TF_tp_NA = 
            is.na(data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])-1]])
        US_PWs[[glen]] = 
            data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])-1]][!TF_tp_NA]
    }
    names(US_PWs) = names(data_omics[[1]][[2]][[2]])
    
    return(US_PWs)
}

#' Get upstream regulators of identified transcription factors.
#'  
#' @param data_omics OmicsData object.
#' @return list of length = number of gene/transcript time points, each element
#' containing a vector of protein regulator IDs.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyRsofTFs(data_omics, only_enriched = FALSE, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' getUS_regulators(data_omics)
#' }
getUS_regulators <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    US_regulators = list()
    for(glen in 1: length(data_omics[[1]][[1]][[1]][[2]]))
    {
        US_regulators[[glen]] = 
            data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]]
    }
    names(US_regulators) = names(data_omics[[1]][[2]][[2]])
    
    return(US_regulators)
}

#' Get upstream regulators of identified transcription factors.
#'  
#' @param data_omics OmicsData object.
#' @return biopax model generated as consensus biopax models from all 
#' pathway databases selected for analysis.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' getBiopaxModel(data_omics)
#' }
getBiopaxModel <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    BiopaxModel = data_omics[[2]][[2]]
    return(BiopaxModel)
}


#' Get protein intersection for the omics data on the different time points.
#'  
#' The timepoints or measurement names for comparison have to be defined in 
#' tp_prot and tp_genes as given in the readOmics function.
#'  
#' @param data_omics OmicsData object.
#' @param tp_prot numeric integer defining protein timepoint measurement chosen 
#' for comparison.
#' @param tp_genes numeric integer defining gene/transcript timepoint 
#' measurement chosen for comparison.
#' @return list with three elements: 1) character vector of protein IDs
#' identified in both upstream and downstream analysis 2) protein time point
#' 3) gene/transcript time point.
#' @keywords manip
#' @export
#' @examples
#' #please run with whole database files (prepared according to vignette)
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("chea", "pazar"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, only_enriched = FALSE, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' data_omics = identifyPWTFTGs(data_omics, only_enriched = FALSE)
#' getProteinIntersection(data_omics, tp_prot = 4, tp_genes = 4)
#' }
getProteinIntersection <- function(data_omics, tp_prot, tp_genes) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    if(!tp_prot %in% data_omics[[1]][[1]][[1]][[1]])
    {stop("tp_prot is not found in protein time points of 
          OmicsData object.")}
    if(!tp_genes %in% data_omics[[1]][[1]][[1]][[2]])
    {stop("tp_genes is not found in gene/transcript time points 
          of OmicsData object.")}
    
    plen = which(data_omics[[1]][[1]][[1]][[1]] == tp_prot)
    glen = which(data_omics[[1]][[1]][[1]][[2]] == tp_genes)
    prot_inters = vector()
    prot_data_prot = data_omics[[1]][[2]][[1]][[plen]][,1]
    prot_data_genes = data_omics[[1]][[3]][[2]][[glen+1]]$regulatorsPW
    prot_inters = prot_data_prot[which(prot_data_prot %in%
                                           as.character(prot_data_genes))]
    return(list(Protein_Intersection = prot_inters, Protein_Timepoint = tp_prot, 
                Gene_Timepoint = tp_genes))
    }


#' Get TF intersection for the omics data on the different time points.
#'  
#' @param data_omics OmicsData object.
#' @param tp_prot numeric integer defining protein timepoint measurement chosen 
#' for comparison.
#' @param tp_genes numeric integer defining gene/transcript timepoint
#' measurement chosen for comparison.
#' @param only_enriched Boolean value defining if transcription factors should 
#' be identified only for enriched pathways (TRUE); or for all identified 
#' pathways (FALSE); default is TRUE.
#' @return list with three elements: 1) character vector of transcription factor
#' IDs identified in both upstream and downstream analysis 2) protein time point
#' 3) gene/transcript time point.
#' @keywords manip
#' @export
#' @examples
#' #please run with whole database files (prepared according to vignette)
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("chea", "pazar"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, only_enriched = FALSE, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' data_omics = identifyPWTFTGs(data_omics, only_enriched = FALSE)
#' getTFIntersection(data_omics, 4,4,only_enriched = TRUE)
#' }
getTFIntersection <- function(data_omics, tp_prot, tp_genes,
                              only_enriched = TRUE) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    if(!tp_prot %in% data_omics[[1]][[1]][[1]][[1]])
    {stop("tp_prot is not found in protein time points of 
          OmicsData object.")}
    if(!tp_genes %in% data_omics[[1]][[1]][[1]][[2]])
    {stop("tp_genes is not found in gene/transcript time points 
          of OmicsData object.")}
    
    plen = which(data_omics[[1]][[1]][[1]][[1]] == tp_prot)
    glen = which(data_omics[[1]][[1]][[1]][[2]] == tp_genes)
    TF_inters = vector()
    prot_TFs = unique(getDS_TFs(data_omics)[[plen]])
    if(only_enriched == TRUE)
    {gene_TFs = unique(as.character(getUS_TFs(data_omics)[[glen]][,1]))
    }else{
        gene_TFs = 
            unique(as.character(getUS_TFs(data_omics)[[glen]][which(getUS_TFs(data_omics)[[glen]][,2]==1),1])) 
    }
    TF_inters = prot_TFs[which(prot_TFs %in% gene_TFs)]
    return(list(TF_Intersection =TF_inters, Protein_Timepoint = tp_prot, 
                Gene_Timepoint = tp_genes))
    }


#' Get genes intersection for the omics data on the different time points.
#'  
#' @param data_omics OmicsData object.
#' @param tp_prot numeric integer defining protein timepoint measurement chosen 
#' for comparison.
#' @param tp_genes numeric integer defining gene/transcript timepoint 
#' measurement chosen for comparison.
#' @return list with three elements: 1) character vector of gene
#' IDs identified in both upstream and downstream analysis 2) protein time point
#' 3) gene/transcript time point.
#' @keywords manip
#' @export
#' @examples
#' #please run with whole database files (prepared according to vignette)
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("chea", "pazar"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, only_enriched = FALSE, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' data_omics = identifyPWTFTGs(data_omics, only_enriched = FALSE)
#' getGenesIntersection(data_omics, tp_prot = 4, tp_genes = 4)
#' }
getGenesIntersection <- function(data_omics, tp_prot, tp_genes) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    if(!tp_prot %in% data_omics[[1]][[1]][[1]][[1]])
    {stop("tp_prot is not found in protein time points of 
          OmicsData object.")}
    if(!tp_genes %in% data_omics[[1]][[1]][[1]][[2]])
    {stop("tp_genes is not found in gene/transcript time points 
          of OmicsData object.")}
    
    plen = which(data_omics[[1]][[1]][[1]][[1]] == tp_prot)
    glen = which(data_omics[[1]][[1]][[1]][[2]] == tp_genes)
    genes_inters = vector()
    prot_tgenes = unique(getDS_TGs(data_omics)[[plen]])
    genes = unique(data_omics[[1]][[2]][[2]][[glen]][,1])
    genes_inters = prot_tgenes[which(prot_tgenes %in% genes)]
    
    return(list(Genes_Intersection = genes_inters, Protein_Timepoint = tp_prot, 
                Gene_Timepoint = tp_genes))
    }


#' Get omics data intersection on the three levels.
#'  
#' Get intersection for the omics data on all three levels (proteins, TFs,
#' genes) on corresponding time points.
#'  
#' @param data_omics OmicsData object.
#' @return list with three elements: 
#' 1) protein intersection
#' 2) transcription factor intersection
#' 3) gene intersection
#' each element contains a list with overlapping time points of both upstream
#' and downstream analyses.
#' @keywords manip
#' @export
#' @examples
#' #please run with whole database files (prepared according to vignette)
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("chea", "pazar"))
#' \dontrun{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, only_enriched = FALSE, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' data_omics = identifyPWTFTGs(data_omics, only_enriched = FALSE)
#' gettpIntersection(data_omics)
#' }
gettpIntersection <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    same_tps = data_omics[[1]][[1]][[1]][[1]][which(data_omics[[1]][[1]][[1]][[1]] %in% 
                                                        data_omics[[1]][[1]][[1]][[2]])]
    prot = list()
    TF = list()
    genes = list()
    for(g in same_tps)
    { temp_ind = which(same_tps==g)
      prot[[temp_ind]] = getProteinIntersection(data_omics, g,g)[[1]]
      TF[[temp_ind]] = getTFIntersection(data_omics, g,g)[[1]]
      genes[[temp_ind]] = getGenesIntersection(data_omics, g,g)[[1]]
    }
    names(prot) = paste("tp", same_tps, sep = "")
    names(TF) = paste("tp", same_tps, sep = "")
    names(genes) = paste("tp", same_tps, sep = "")
    return(list(Intersection = list(Protein_Intersection = prot, 
                                    TF_Intersection = TF, 
                                    Genes_Intersection = genes)))
}