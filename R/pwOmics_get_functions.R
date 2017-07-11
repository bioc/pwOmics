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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
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
#' containing a data frame with the pathway IDs in the generated biopax model and
#' corresponding pathway names.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
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
        DS_PWs[[plen]] = as.data.frame(tps_PWs)[,1:4]
        colnames(DS_PWs[[plen]])[3] = "upreg"
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
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
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
            if(length(as.character(DS_TFs[[plen]][,1][which(DS_TFs[[plen]][,4] == 1)])) == 0)
            {message("No transcription factors could be identified for time point ",
                     data_omics[[1]][[1]][[1]][[1]][plen],
                     " in downstream analyis. \n", sep = "")
            }else{
                DS_TFs[[plen]] = 
                    DS_TFs[[plen]][which(DS_TFs[[plen]][,4] == 1),1:4]}
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
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
#' getDS_TGs(data_omics)
#' }
getDS_TGs <- function(data_omics) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    DS_TGs = list()
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    {
        len_omics = length(data_omics[[1]][[3]][[1]][[plen+1]])
        DS_TGs[[plen]] = data_omics[[1]][[3]][[1]][[plen+1]][[len_omics]]
    }
    names(DS_TGs) = names(data_omics[[1]][[2]][[1]])
    
    return(DS_TGs)
}


#' Get upstream TFs.
#'  
#' @param data_omics OmicsData object.
#' @return list of length = number of gene/transcript time points, each element
#' containing a data frame with upstream transcription factors.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics_plus = identifyPR(data_omics_plus)
#' \dontrun{
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, 
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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics_plus = identifyPR(data_omics_plus)
#' \dontrun{
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, 
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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics_plus = identifyPR(data_omics_plus)
#' \dontrun{
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, 
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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
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
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection. Type of checking is defined with
#' parameter 'phospho'.
#' @param phospho boolean value; TRUE in case up- and downregulation should be
#' checked based on provided downstream phosphoprotein influence from 
#' identifyPR function; FALSE in case up- and downregulation should be checked
#' for without phosphoprotein database knowledge. Default is TRUE.
#' @return list with three elements: 1) character vector of protein IDs
#' identified in both upstream and downstream analysis 2) protein time point
#' 3) gene/transcript time point.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, noTFs_inPW = 1, order_neighbors = 10)
#' getProteinIntersection(data_omics, tp_prot = 4, tp_genes = 4, 
#' updown = FALSE, phospho = TRUE)
#' }
getProteinIntersection <- function(data_omics, tp_prot, tp_genes, updown = FALSE, phospho = TRUE) {
    
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
    if(updown == FALSE)
    {   prot_data_prot = data_omics[[1]][[2]][[1]][[plen]][,1]
        prot_data_genes = as.character(data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]]$regulators)    
        prot_inters = prot_data_prot[which(prot_data_prot %in% as.character(prot_data_genes))]
    }else{
        if(phospho == TRUE)
        {prot_data_prot = data_omics[[1]][[2]][[1]][[plen]]
        prot_data_prot[,4] = (prot_data_prot[,2] * prot_data_prot[,3])>0 
        
        prot_data_genes = data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]]
        temp_reg = unique(merge(data_omics[[5]], prot_data_genes, by.x = "Phosphoprotein", by.y = "regulators"))
        temp_reg[,4] = (temp_reg[,2] * (temp_reg[,3]=="upreg"))>0
        prot_data_genes[,3] = prot_data_genes[,2]
        prot_data_genes = unique(prot_data_genes)
        prot_data_genes[,3] = as.character(prot_data_genes[,3])
        
        for(j in 1: dim(prot_data_genes)[1])
        {if(as.character(prot_data_genes[j,1]) %in% as.character(temp_reg[,1]))
            {
              if(length(unique(temp_reg[which(temp_reg[,1] == as.character(prot_data_genes[j,1])),4])) == 1)  ##if there is no ambiguity about regulation: TRUE/FALSE
              {prot_data_genes[j,3] = temp_reg[which(temp_reg[,1] == as.character(prot_data_genes[j,1])),4]
              }else{
               prot_data_genes[j,3] = NA ##matches everything, if there is ambiguity
              }
            }    
        }
        #in case phosphorylation data is ambiguous (NA) both options should be ok
        for(k in 1:dim(prot_data_prot)[1])
        {  if(length(which(as.character(prot_data_genes$regulators) %in% as.character(prot_data_prot[k,1])))> 0)
           {  if(is.na(prot_data_genes[which(as.character(prot_data_genes$regulators) %in% prot_data_prot[k,1]),3]) | is.na(prot_data_prot[k,4]))
                { 
                prot_inters[k] = as.character(prot_data_genes[which(as.character(prot_data_genes$regulators) %in% as.character(prot_data_prot[k,1])),]$regulators[1])
                }else if(prot_data_genes[which(as.character(prot_data_genes$regulators) %in% prot_data_prot[k,1]),3] == prot_data_prot[k,4]){
                prot_inters[k] = as.character(prot_data_genes[which(as.character(prot_data_genes$regulators) %in% as.character(prot_data_prot[k,1])),]$regulators[1])    
                }
           }
        }    
        prot_inters = unique(na.omit(prot_inters))  
        }else{
            prot_data_prot = data_omics[[1]][[2]][[1]][[plen]]
            prot_data_prot[,4] = (prot_data_prot[,2])>0 
            
            prot_data_genes = data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]]
            prot_data_genes = unique(prot_data_genes)
            prot_data_genes[,3] = as.character(prot_data_genes[,2]) == "upreg"
            
            #in case phosphorylation data is ambiguous (NA) both options should be ok
            for(k in 1:dim(prot_data_prot)[1])
            {  if(length(which(as.character(prot_data_genes$regulators) %in% as.character(prot_data_prot[k,1])))> 0)
                {  if(is.na(prot_data_genes[which(as.character(prot_data_genes$regulators) %in% prot_data_prot[k,1]),3]) | is.na(prot_data_prot[k,4]))
                    {   prot_inters[k] = as.character(prot_data_genes[which(as.character(prot_data_genes$regulators) %in% as.character(prot_data_prot[k,1])),]$regulators[1])
                    }else if(prot_data_genes[which(as.character(prot_data_genes$regulators) %in% prot_data_prot[k,1]),3] == prot_data_prot[k,4]){
                        prot_inters[k] = as.character(prot_data_genes[which(as.character(prot_data_genes$regulators) %in% as.character(prot_data_prot[k,1])),]$regulators[1])    
                    }
                }
            }    
            prot_inters = unique(na.omit(prot_inters))  
        }
    }
    return(list(Protein_Intersection = as.character(prot_inters), Protein_Timepoint = tp_prot, 
                Gene_Timepoint = tp_genes))
}


#' Get TF intersection for the omics data on the different time points.
#'  
#' @param data_omics OmicsData object.
#' @param tp_prot numeric integer defining protein timepoint measurement chosen 
#' for comparison.
#' @param tp_genes numeric integer defining gene/transcript timepoint
#' measurement chosen for comparison.
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection. Type of checking is defined with
#' parameter 'phospho'.
#' @param phospho boolean value; TRUE in case up- and downregulation should be
#' checked based on provided downstream phosphoprotein influence from 
#' identifyPR function; FALSE in case up- and downregulation should be checked
#' for without phosphoprotein database knowledge. Default is TRUE.
#' @return list with three elements: 1) character vector of transcription factor
#' IDs identified in both upstream and downstream analysis 2) protein time point
#' 3) gene/transcript time point.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, noTFs_inPW = 1, order_neighbors = 10)
#' getTFIntersection(data_omics, tp_prot = 4, tp_genes = 4, updown = FALSE, 
#' phospho = TRUE)
#' }
getTFIntersection <- function(data_omics, tp_prot, tp_genes,
                              updown = FALSE, phospho = TRUE) {
    
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
    if(updown == FALSE)
    {   prot_TFs = as.character(unique(getDS_TFs(data_omics)[[plen]])$genes_PW)
        if(dim(getUS_TFs(data_omics)[[glen]])[1]!= 0)
           {gene_TFs = unique(as.character(getUS_TFs(data_omics)[[glen]][,1]))
            }else{
            gene_TFs = NA   
           }
        TF_inters = prot_TFs[which(prot_TFs %in% gene_TFs)]
    }else{
        if(phospho == TRUE)
        {prot_TFs = getDS_TFs(data_omics)[[plen]]
        
        temp_reg_updown = data.frame(TF = rownames(data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]]),
                                     regulators = data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][,1],
                                     regulation = data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][,2],
                                     new_regulation = data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][,2])
        temp_reg_updown[,4] = as.character(temp_reg_updown[,4])
        ##gehe durch alle upstream protein regulators
        for(k in 1: dim(data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]])[1])
        {
            #falls einer davon in phospholiste
            if(sum(data_omics[[5]][,1] %in% data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][k,1])>0)
            { 
                #und phospholiste eindeutig
                if(length(unique(data_omics[[5]][which(data_omics[[5]][,1] %in% data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][k,1]),2]))==1)
                {temp_reg_updown[k,4] = unique(data_omics[[5]][which(data_omics[[5]][,1] %in% data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][k,1]),2] * 
                                                   (data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][k,2]=="upreg"))
                }else{ #sonst
                    temp_reg_updown[k,4] = NA   
                }
            }
        }
        temp_reg_updown[,1] = gsub("\\.[0-9]*", "", as.character(temp_reg_updown[,1]))
        Tfs = data.frame(TFs = unique(temp_reg_updown[,1]), finalreg = rep(NA, times = length(unique(temp_reg_updown[,1]))))
        Tfs = Tfs[which(Tfs[,1] %in% getUS_TFs(data_omics)[[glen]][,1]),]
        for(s in 1: length(Tfs[,1]))
        {#falls für diesen TF ambiguous regulation --> keep in
            if(NA %in% subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4])  ##either NA
            {Tfs[s,2] = NA
            }else if(1 %in% subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4] &  ##or up AND down
                    -1 %in% subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4])
            {Tfs[s,2] = NA
            }else if(length(unique(subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4])) == 1) ##if unambigous take it
             {Tfs[s,2] =   unique(subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4])
            }
        }
        
        prot_TFs = cbind(prot_TFs, rep(NA, times = dim(prot_TFs)[1]))
        colnames(prot_TFs)[5] = "final_reg"
        for(s in 1: dim(prot_TFs)[1])
        {   if(length(which(as.character(Tfs[,1]) %in% prot_TFs[s,1]))> 0)
               {if(!is.na(Tfs[which(as.character(Tfs[,1]) %in% prot_TFs[s,1]),]$finalreg == prot_TFs[s,2] |
                   is.na(Tfs[which(as.character(Tfs[,1]) %in% prot_TFs[s,1]),]$finalreg)))
                  {prot_TFs$final_reg[s] = TRUE
                     }else{
                   prot_TFs$final_reg[s] = FALSE 
                     }
               }
        }
        
        TF_inters = as.character(unique(prot_TFs[which(prot_TFs$final_reg == TRUE),1]))
        }else{
            prot_TFs = data_omics[[1]][[3]][[1]][[plen+1]][[length(data_omics[[1]][[3]][[1]][[plen+1]])-1]]
            prot_TFs = prot_TFs[which(prot_TFs[,4] == 1),]
            prot_TFs = unique(prot_TFs)
            
            Tfs = unique(getUS_TFs(data_omics)[[glen]])
            for(s in 1: dim(prot_TFs)[1])
            {   if(as.character(prot_TFs[s,1]) %in% as.character(Tfs[,1])) 
                {   if( TRUE %in% Tfs[which(as.character(Tfs[,1]) %in% as.character(prot_TFs[s,1])),2] &
                        prot_TFs[which(as.character(prot_TFs[s,1]) %in% as.character(Tfs[,1])),2] == TRUE)  ##if unambigous take it
                       {TF_inters[s] = as.character(prot_TFs[s,1])}
                }
            }
            TF_inters = unique(na.omit(TF_inters))    
        }
    }
    
    return(list(TF_Intersection = as.character(TF_inters), Protein_Timepoint = tp_prot, 
                Gene_Timepoint = tp_genes))
}


#' Get genes intersection for the omics data on the different time points.
#'  
#' @param data_omics OmicsData object.
#' @param tp_prot numeric integer defining protein timepoint measurement chosen 
#' for comparison.
#' @param tp_genes numeric integer defining gene/transcript timepoint 
#' measurement chosen for comparison.
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection. Type of checking is defined with
#' parameter 'phospho'.
#' @param phospho boolean value; TRUE in case up- and downregulation should be
#' checked based on provided downstream phosphoprotein influence from 
#' identifyPR function; FALSE in case up- and downregulation should be checked
#' for without phosphoprotein database knowledge. Default is TRUE.
#' @return list with three elements: 1) character vector of gene
#' IDs identified in both upstream and downstream analysis 2) protein time point
#' 3) gene/transcript time point.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, noTFs_inPW = 1, order_neighbors = 10)
#' getGenesIntersection(data_omics, tp_prot = 4, tp_genes = 4, updown = FALSE, 
#' phospho = TRUE)
#' }
getGenesIntersection <- function(data_omics, tp_prot, tp_genes,
                                 updown = FALSE, phospho = TRUE) {
    
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
    if(updown == FALSE)
    {  prot_tgenes = unique(getDS_TGs(data_omics)[[plen]][,1])
       genes = unique(data_omics[[1]][[2]][[2]][[glen]][,1])
       genes_inters = prot_tgenes[which(prot_tgenes %in% genes)]
    }else{
       if(phospho == TRUE)
       {prot_tgenes = getDS_TGs(data_omics)[[plen]]
       prot_tgenes[,2] = gsub(" ", "",prot_tgenes[,2])
       
       temp_reg_updown = data.frame(TF = rownames(data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]]),
                                    regulators = data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][,1],
                                    regulation = data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][,2],
                                    new_regulation = data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][,2])
       temp_reg_updown[,4] = as.character(temp_reg_updown[,4])
       ##search through all upstream protein regulators
       for(k in 1: dim(data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]])[1])
       {
           #in case one of them is in phosphoprotein list
           if(sum(data_omics[[5]][,1] %in% data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][k,1])>0)
           { 
               #and phosphoprotein list is unambiguous calculate new influence (phospho * regulation)
               if(length(unique(data_omics[[5]][which(data_omics[[5]][,1] %in% data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][k,1]),2]))==1)
               {temp_reg_updown[k,4] = unique(data_omics[[5]][which(data_omics[[5]][,1] %in% data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][k,1]),2] 
                                              * (data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])]][k,2]=="upreg"))
               }else{
                   temp_reg_updown[k,4] = NA   
               }
           }
       }
       temp_reg_updown[,1] = gsub("\\.[0-9]*", "", as.character(temp_reg_updown[,1]))
       Tfs = data.frame(TFs = unique(temp_reg_updown[,1]), finalreg = rep(NA, times = length(unique(temp_reg_updown[,1]))))
       Tfs = Tfs[which(Tfs[,1] %in% getUS_TFs(data_omics)[[glen]][,1]),]
       for(s in 1: length(Tfs[,1]))
       {#falls für diesen TF ambiguous regulation --> keep in
           if(NA %in% subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4])
           {Tfs[s,2] = NA
           }else if(1 %in% subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4] & 
                    -1 %in% subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4])
           {Tfs[s,2] = NA
           }else if(length(unique(subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4])) == 1)
           {Tfs[s,2] =   unique(subset(temp_reg_updown, temp_reg_updown[,1] == as.character(Tfs[s,1]))[,4])
           }
       }
       
       ##für die TFs, für die klare regulation von oben ist:
       check_ds_t_list = as.character(Tfs[which(!is.na(Tfs[,2])),1])
       
       #identify for this time point all upstream TFs
       temp = data.frame()
       for(j in 1: (length(data_omics[[1]][[3]][[2]][[glen+1]])-2))
       {   if(!names(data_omics[[1]][[3]][[2]][[glen+1]][[j]])[1] == "NA.")
           temp = rbind(temp, data_omics[[1]][[3]][[2]][[glen+1]][[j]]) 
       }
       #get regulation: up --> TRUE, down --> FALSE
       for(s in 1: dim(temp)[1])
       {
           if(temp[s,2] == TRUE)
           {temp[s,2] = 1
           }else if(temp[s,2] == FALSE){
            temp[s,2] = -1   
           } #when identified in downstream analysis
           if(temp[s,1] %in% Tfs[,1])
            {  if(!is.na(Tfs[which(Tfs[,1] %in% temp[s,1]),2])) #when finalreg defined in DS analysis calculate resulting regulation
                 {temp[s,2] = temp[s,2]*(Tfs[which(Tfs[,1] %in% temp[s,1]),2]=="upreg") }
            }
       }
       ##all TFs with corrected downstream regulation
       genes = data.frame(genes = data_omics[[1]][[2]][[2]][[glen]][,1], 
                          final_reg = rep(NA, times = dim(data_omics[[1]][[2]][[2]][[glen]])[1]))  ##genes at that time point
       for(j in 1: (length(data_omics[[1]][[3]][[2]][[glen+1]])-2))
       {  if(length(check_ds_t_list)> 0)
          {for(k in 1: length(check_ds_t_list)) 
            {
              if(!names(data_omics[[1]][[3]][[2]][[glen+1]][[j]])[1] == "NA.")
              {##falls TFs in klarer regulation von oben nehme aus temp liste, sonst NA
              if(check_ds_t_list[k] %in% as.character(data_omics[[1]][[3]][[2]][[glen+1]][[j]]$upstreamTFs))
              {genes[j, 2] = temp[which(temp[,1] == check_ds_t_list[k]),2]}
              }
            }
         }else{
         genes[j, 2] = NA}
       }

       genes_inters = vector()
       for(k in 1:dim(genes)[1])
       {  if(length(which(as.character(prot_tgenes[,1]) %in% genes[k,1]))> 0)
           {if(is.na(genes[k,2])){   #in case of ambiguity
               genes_inters[k] = as.character(genes[k,1]) 
             }else if(as.character(genes[k,2]) %in% unique(prot_tgenes[which(as.character(prot_tgenes[,1]) %in% genes[k,1]),2]) ){
                 genes_inters[k] = as.character(genes[k,1])
             }else if(length(unique(prot_tgenes[which(as.character(prot_tgenes[,1]) %in% genes[k,1]),2])) == 1){ #no ambiguity
                if(unique(prot_tgenes[which(as.character(prot_tgenes[,1]) %in% genes[k,1]),2]) == genes[k,2]) 
                {genes_inters[k] = as.character(genes[k,1])}
             }
           }else{
             genes_inters[k] =  NA   
           }  
       }
       genes_inters = unique(na.omit(genes_inters)) 
     }else{
         prot_tgenes = getDS_TGs(data_omics)[[plen]]
         prot_tgenes = unique(prot_tgenes)
         prot_tgenes[,2] = gsub(" ", "", prot_tgenes[,2])
         
         gene_reg = data_omics[[1]][[2]][[2]][[glen]]
         gene_reg[,3] = gene_reg[,2]>0
         
         for(k in 1: dim(prot_tgenes)[1])
         {
             if(as.character(prot_tgenes[k,1]) %in% gene_reg[,1])
             {   if(gene_reg[which(gene_reg[,1] == as.character(prot_tgenes[k,1])),3] == TRUE & as.character(prot_tgenes[k,2]) == "TRUE")
                 {genes_inters[k] = as.character(prot_tgenes[k,1])} 
             }
         }
         genes_inters = unique(na.omit(genes_inters))      
     }
    }
    return(list(Genes_Intersection = as.character(genes_inters), Protein_Timepoint = tp_prot, 
                Gene_Timepoint = tp_genes))
    }


#' Get omics data intersection on the three levels.
#'  
#' Get intersection for the omics data on all three levels (proteins, TFs,
#' genes) on corresponding time points.
#'  
#' @param data_omics OmicsData object.
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection. Type of checking is defined with
#' parameter 'phospho'.
#' @param phospho boolean value; TRUE in case up- and downregulation should be
#' checked based on provided downstream phosphoprotein influence from 
#' identifyPR function; FALSE in case up- and downregulation should be checked
#' for without phosphoprotein database knowledge. Default is TRUE.
#' @return list with three elements: 
#' 1) protein intersection
#' 2) transcription factor intersection
#' 3) gene intersection
#' each element contains a list with overlapping time points of both upstream
#' and downstream analyses.
#' @keywords manip
#' @export
#' @examples
#' data(OmicsExampleData)
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24), 
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics.newupdown")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics.newupdown"))
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, noTFs_inPW = 1, order_neighbors = 10)
#' gettpIntersection(data_omics, updown = FALSE, phospho = TRUE)
#' }
gettpIntersection <- function(data_omics, updown = FALSE, phospho = TRUE) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    same_tps = data_omics[[1]][[1]][[1]][[1]][which(data_omics[[1]][[1]][[1]][[1]] %in% 
                                                        data_omics[[1]][[1]][[1]][[2]])]
    prot = list()
    TF = list()
    genes = list()
    for(g in same_tps)
    { temp_ind = which(same_tps==g)
      prot[[temp_ind]] = getProteinIntersection(data_omics, g,g, updown, phospho)[[1]]
      TF[[temp_ind]] = getTFIntersection(data_omics, g,g, updown, phospho)[[1]]
      genes[[temp_ind]] = getGenesIntersection(data_omics, g,g, updown, phospho)[[1]]
    }
    names(prot) = paste("tp", same_tps, sep = "")
    names(TF) = paste("tp", same_tps, sep = "")
    names(genes) = paste("tp", same_tps, sep = "")
    return(list(Intersection = list(Protein_Intersection = prot, 
                                    TF_Intersection = TF, 
                                    Genes_Intersection = genes)))
}