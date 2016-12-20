#' Identify phosphorylation regulation influence downstream
#'
#' This function identifies the downstream regulation influence
#' of phosphoprotein regulation for further downstream analysis steps.
#'  
#' @param data_omics_plus output list of readPWdata function; first element
#' contains an OmicsData object, secons element the genelist data corresponding 
#' to the selected pathway database.
#' @return OmicsData object: list of 4 elements (OmicsD, PathwayD, TFtargetsD, 
#' Status); OmicsD containing omics data set + results (after analysis);
#' PathwayD containing selected pathway databases + biopax model;
#' TFtargetsD containing selected TF target gene databases + TF target gene data.
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
#' package = "pwOmics")) 
#' \dontrun{
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics_plus = identifyPR(data_omics_plus)
#' }
identifyPR <- function(data_omics_plus){
    
    updown = NULL
    for(s in 1:length(data_omics_plus[[1]][[1]][[1]][[1]][[1]]))
    {
        updown = data_omics_plus[[1]][[5]][match(as.character(data_omics_plus[[1]][[1]][[2]][[1]][[s]][,1]), as.character(data_omics_plus[[1]][[5]][,1])),2]
        data_omics_plus[[1]][[1]][[2]][[1]][[s]] = cbind(data_omics_plus[[1]][[1]][[2]][[1]][[s]], updown)
    }
    message("Phosphoprotein downstream regulation information is
            considered in downstream analysis. \n")
    return(data_omics_plus)
}



#' Identify pathway IDs and pathway names of differentially abundant proteins
#'
#' This function identifies the pathways of the differentially abundant 
#' phosphoproteins dependent on the chosen database. 
#' Requires rBiopaxParser package. Takes a 
#' lot of time for a high number of proteins and/or if all databases are chosen.
#' First, chosen databases are loaded, then new internal pathway IDs are 
#' generated. 
#' Afterwards the genelists of the different databases are loaded or generated, 
#' depending on the loadgenelists option. After pathway identification for the 
#' reference time point, also pathway identification for different time points 
#' is performed. Pathway ID mapping takes some time, especially for such big 
#' databases as reactome, so use savegenelists and loadgenelists for easier and
#' faster usage...
#'  
#' @param data_omics_plus output list of readPWdata function; first element
#' contains an OmicsData object, secons element the genelist data corresponding 
#' to the selected pathway database.
#' @return OmicsData object: list of 4 elements (OmicsD, PathwayD, TFtargetsD, 
#' Status); OmicsD containing omics data set + results (after analysis);
#' PathwayD containing selected pathway databases + biopax model;
#' TFtargetsD containing selected TF target gene databases + TF target gene data.
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
#' package = "pwOmics")) 
#' \dontrun{
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWs(data_omics_plus)
#' }
identifyPWs <- function(data_omics_plus){
    
    if(data_omics_plus[[1]][[4]] == 1)
    {stop("Please read in the omics data set and both pathway
          database information and TF-target gene information first with
          readOmics, readTFdata and readPWdata functions.")}
    
    if(class(data_omics_plus[[1]]) != "OmicsData")
    {stop("'data_omics_plus[[1]]' is not an OmicsData
          object.")}
    
    if(length(data_omics_plus[[2]]) < length(data_omics[[2]][[1]]))
    {stop("'data_omics_plus[[2]]' does not contain all genelists of the
          selected pathway databases. Please check if all genelists are
          present in the working directory and if necessary run readPWdata
          again with loadgenelists set to FALSE.")}
    
    data_omics = PWidentallprots(data_omics_plus[[1]], data_omics_plus[[2]])
    message("Pathways are identified for all proteins measured. \n")
    
    data_omics = PWidenttps(data_omics)
    message("Pathways are identified for the different timepoints. \n")
    
    if(length(data_omics[[1]][[3]][[1]])!= 0 & 
           length(data_omics[[1]][[3]][[2]])!= 0)
    {data_omics[[4]] = data_omics[[4]] +1
    }else{
        data_omics[[4]] = data_omics[[4]] 
    }
    
    return(data_omics)
}





#' Identify TFs in pathways and their target genes - downstream 
#' analysis.
#'
#' This function identifies the transcription factors being part of the
#' pathways of downstream analysis. Subsequently it finds the target genes of
#' these transcription factors from the selected TF-target gene database.
#'  
#' @param data_omics OmicsData object.
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection; FALSE = default, if only 
#' deregulation should be checked for.
#' @return OmicsData object: list of 4 elements (OmicsD, PathwayD, TFtargetsD,
#' Status); OmicsD containing omics data set + results (after analysis);
#' PathwayD containing selected pathway databases + biopax model;
#' TFtargetsD containing selected TF target gene databases + TF target gene data.
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
#' package = "pwOmics")) 
#' \dontrun{
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
#' }
identifyPWTFTGs <- function(data_omics, updown = FALSE) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    genelists = loadGenelists()
    if(length(genelists) < length(data_omics[[2]][[1]]))
    {stop("Current working directory does not contain all genelists of 
          the selected pathway databases. Please check if all genelists are
          present in the working directory and if necessary run readPWdata
          again with loadgenelists set to FALSE.")}    
    message("Genelists of databases are loaded/generated. \n\n")  
    pathwayIDs = pathwayNames = NULL
    
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    { 
        PWinfo = preparePWinfo(data_omics, plen) 
        data_omics = PWinfo[[1]]
        tps_PWs = unique(PWinfo[[2]])
        PWofinterest = PWinfo[[3]]

        if(length(tps_PWs) != 0)
        {genelist_n = apply(rbindlist(genelists),2,as.character)
         genes_PW = list()
         for(k in 1: dim(tps_PWs)[1])
         {   ind_genesPW = vector()
             ind_genesPW = which(as.character(genelist_n[,2]) == 
                                     as.character(tps_PWs[,pathwayIDs])[k])
             genes_PW[[k]] = unique(genelist_n[ind_genesPW,1])
             genes_PW[[k]] = data.frame(genes_PW = genes_PW[[k]], "upreg" = rep(tps_PWs[,upreg][k], times = length(genes_PW[[k]])),
                                        "phosphoeffect" = rep(tps_PWs[,phosphoeffect][k], times = length(genes_PW[[k]])))  #new
         }
         names(genes_PW) = tps_PWs[,pathwayNames]
         genes_PW_ov = genes_PW     
         message("Gene sets of pathways are identified: ", 
                 names(data_omics[[1]][[2]][[1]][plen]), "\n")   
         
         temp_genelist = do.call("rbind", genes_PW)  #new
         temp_lists = identTFTGsinPWs(data_omics, temp_genelist)  
         temp_Genelist = temp_lists[[1]]
         temp_targetlist = apply(temp_lists[[2]], 2, as.character)
         colnames(temp_Genelist) = c("genes_PW","upreg", "phosphoeffect", "TF_PW")
         data_omics[[1]][[3]][[1]][[plen+1]][[length(data_omics[[1]][[3]][[1]][[plen+1]])+1]] =
             temp_Genelist
         names( data_omics[[1]][[3]][[1]][[plen+1]][length(data_omics[[1]][[3]][[1]][[plen+1]])]) = 
             "Genelist_PW"
         data_omics[[1]][[3]][[1]][[plen+1]][[length(data_omics[[1]][[3]][[1]][[plen+1]])+1]] = 
             temp_targetlist 
         names( data_omics[[1]][[3]][[1]][[plen+1]][length(data_omics[[1]][[3]][[1]][[plen+1]])]) = 
             "Target_Genelist_PW"
         message("Transcription factors in pathways of downstream analysis are identified: ", 
                 names(data_omics[[1]][[2]][[1]][plen]), "\n")
         message("Target genes of TFs in downstream analysis are identified: ", 
                 names(data_omics[[1]][[2]][[1]][plen]), "\n\n")
         
        }else{
            message(paste("No downstream signaling could be found for time point ", 
                          data_omics[[1]][[1]][[1]][[1]][[plen]],"\n", sep = ""))
        }
    }
    return(data_omics)
}


