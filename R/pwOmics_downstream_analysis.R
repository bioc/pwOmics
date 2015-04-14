#' Identify pathway IDs and pathway names of differentially abundant proteins
#'
#' This function identifies the pathways of the differentially abundant proteins
#' dependent on the chosen database. Requires rBiopaxParser package. Takes a 
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
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' \donttest{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = "Genelists")
#' }
#' \dontrun{
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



#' Pathway enrichment - downstream analysis.
#'
#' This function does pathway enrichment for pathways determined via identifyPWs
#' function. 
#'  
#' @param data_omics OmicsData object.
#' @param method correction method for multiple testing correction as specified
#' in p.adjust documentation; default is Benjamini & Hochberg correction.
#' @param alpha significance level for pathway enrichment; 
#' default is alpha = 0.05.
#' @param ... further input parameters for multiple comparison adjustment.
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
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' \donttest{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = "Genelists")
#' }
#' \dontrun{
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' }
enrichPWs <- function(data_omics, method = "BH", alpha = 0.05, ...) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    if(data_omics[[4]][[1]] < 3)
    {stop("Transcription factors in upstream analysis or pathways 
          in downstream analysis were not yet identified. Please use the 
          identifyTFs and the identifyPWs function to do so.")}
    pathwayIDs = NULL
    
    no_all_proteins = dim(unique(data_omics[[1]][[1]][[2]]))[1]
    PW_NA = vector()
    for(j in 1: no_all_proteins)
    {PW_NA[j] = is.na(data_omics[[1]][[3]][[1]][[1]][[j]][[1]][1])}
    allprots_PWs = rbindlist(data_omics[[1]][[3]][[1]][[1]][!PW_NA])
    
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    { 
        no_de_prots = dim(data_omics[[1]][[2]][[1]][plen][[1]])[1]
        PW_tp_NA = vector()
        for(g in 1: no_de_prots)
        {PW_tp_NA[g] = is.na(data_omics[[1]][[3]][[1]][[plen+1]][[g]][[1]][1])}
        
        tps_PWs = rbindlist(data_omics[[1]][[3]][[1]][[plen+1]][!PW_tp_NA])
        if(length(tps_PWs) != 0)
        {tps_PWs = as.character(tps_PWs[,pathwayIDs])
         p_vals = vector()
         count_all = vector()
         count_de = vector()
         for(PW in 1: length(tps_PWs))
         {   count_all[PW] = length(which(allprots_PWs[,pathwayIDs] == tps_PWs[PW]))
             count_de[PW] = length(which(tps_PWs == tps_PWs[PW]))
             p_vals[PW] = fisher.test(matrix(c(no_de_prots, no_all_proteins - no_de_prots, 
                                               count_de[PW], count_all[PW] - count_de[PW]), nrow = 2, 
                                             dimnames = list(pa = c("de(tp)", "all-de(tp)"), 
                                                             counts = c("prots", "pws"))))$p.value
         }
         q_vals = p.adjust(p_vals, method = "BH")
         enr_PWs = unique(tps_PWs[which(q_vals < alpha)])
         for(k in 1: no_de_prots)
         {if(!is.na(data_omics[[1]][[3]][[1]][[plen+1]][[k]][[1]][1]))
         { data_omics[[1]][[3]][[1]][[plen+1]][[k]][,3] = NA
           colnames(data_omics[[1]][[3]][[1]][[plen+1]][[k]])= 
               c("pathwayIDs", "pathwayNames", "enrichedPWs")
           for(s in 1: length(enr_PWs))
           { ind_match = which(data_omics[[1]][[3]][[1]][[plen+1]][[k]][,1] == enr_PWs[s])
             if(length(ind_match)!= 0)
             {data_omics[[1]][[3]][[1]][[plen+1]][[k]][ind_match,3] = 1}
           }
         }
         }
        }else{
            message(paste("No enriched PWs found for time point ", 
                          data_omics[[1]][[1]][[1]][[1]][[plen]], ".\n", sep = ""))
        }
    }
    message("Pathway enrichment is completed.\n")
    
    return(data_omics)  
}



#' Identify TFs in enriched pathways and their target genes - downstream 
#' analysis.
#'
#' This function identifies the transcription factors being part of the enriched 
#' pathways of downstream analysis. Subsequently it finds the target genes of
#' these transcription factors from the selected TF-target gene database.
#'  
#' @param data_omics OmicsData object.
#' @param only_enriched boolean value defining if transcription factors should
#' be identified only for enriched pathways (TRUE); or for all identified
#' pathways (FALSE); default is TRUE.
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
#' PWdatabase = c("biocarta"), 
#' TFtargetdatabase = c("chea"))
#' \donttest{
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = FALSE)
#' }
#' \dontrun{
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichPWs(data_omics)
#' identifyPWTFTGs(data_omics, only_enriched = FALSE)
#' }
identifyPWTFTGs <- function(data_omics, only_enriched = TRUE) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    genelists = loadGenelists()
    if(length(genelists) < length(data_omics[[2]][[1]]))
    {stop("Current working directory does not contain all genelists of 
          the selected pathway databases. Please check if all genelists are
          present in the working directory and if necessary run readPWdata
          again with loadgenelists set to FALSE.")}    
    message("Genelists of databases are loaded/generated. \n\n")  
    enrichedPWs = pathwayIDs = pathwayNames = NULL
    
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    { 
        PWinfo = preparePWinfo(data_omics, plen) 
        data_omics = PWinfo[[1]]
        tps_PWs = PWinfo[[2]]
        PWofinterest = PWinfo[[3]]
        
        if(only_enriched == TRUE)
        {tps_PWs = tps_PWs[which(tps_PWs[,enrichedPWs] == 1),]}
        
        if(length(tps_PWs) != 0)
        {genelist_n = apply(rbindlist(genelists),2,as.character)
         genes_PW = list()
         for(k in 1: dim(tps_PWs)[1])
         {   ind_genesPW = vector()
             ind_genesPW = which(as.character(genelist_n[,2]) == 
                                     as.character(tps_PWs[,pathwayIDs])[k])
             genes_PW[[k]] = unique(genelist_n[ind_genesPW,1])
         }
         names(genes_PW) = tps_PWs[,pathwayNames]
         genes_PW_ov = genes_PW     
         message("Gene sets of pathways are identified: ", 
                 names(data_omics[[1]][[2]][[1]][plen]), "\n")   
         
         temp_genelist = as.data.frame(unique(unlist(genes_PW))) 
         temp_lists = identTFTGsinPWs(data_omics, temp_genelist)  
         temp_Genelist = temp_lists[[1]]
         temp_targetlist = as.character(temp_lists[[2]])
         colnames(temp_Genelist) = c("genes_PW","TF_PW")
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


