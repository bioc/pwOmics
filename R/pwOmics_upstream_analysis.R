#' Transcription factor identification.
#'
#' This function identifies the upstream transcription factors of the provided
#'  gene IDs.
#'
#' @param data_omics OmicsData object.
#' @return OmicsData object: list of 4 elements (OmicsD, PathwayD, TFtargetsD,
#'  Status); OmicsD containing omics data set + results (after analysis);
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
#' }
identifyTFs <- function(data_omics) {
    
    if(data_omics[[4]] == 1)
    {stop("Please read in the omics data set and both pathway
          database information and TF-target gene information first with
          readOmics, readTFdata and readPWdata functions.")}
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    data_omics = TFidentallgenes(data_omics)
    message("Upstream TFs are identified for all target genes measured.\n")
    
    data_omics = TFidenttps(data_omics)
    message("Upstream TFs are identified for the different timepoints.\n")
    
    if(length(data_omics[[1]][[3]][[1]])!= 0 & 
           length(data_omics[[1]][[3]][[2]])!= 0)
    {data_omics[[4]] = data_omics[[4]] +1
    }else{
        data_omics[[4]] = data_omics[[4]] 
    }
    
    return(data_omics)
}


#' Transcription factor enrichment - upstream analysis.
#'
#' This function does transcription factor enrichment for transcription factor 
#' identified to be upstream of the diff. expressed genes/transcripts in the
#' omics data set imported with readOmics function. In order to use this
#' function first read in the transcription factor target gene information via
#' readTFdata and identify the upstream TFs of the differentially expressed 
#' genes/transcripts with the identifyTFs function.
#'
#' @param data_omics OmicsData object.
#' @param method correction method for multiple testing correction as specified 
#' in p.adjust documentation; default is Benjamini & Hochberg correction.
#' @param alpha significance level for transcription factor enrichment; 
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
#' setwd(system.file("extdata", package = "pwOmics"))
#' data_omics = readTFdata(data_omics)
#' data_omics_plus = readPWdata(data_omics, 
#' loadgenelists = "Genelists")
#' }
#' \dontrun{
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = enrichTFs(data_omics)
#' }
enrichTFs <- function(data_omics, method = "BH", alpha = 0.05, ...) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    if(data_omics[[4]][[1]] < 3)
    {stop("Transcription factors in upstream analysis or pathways in 
           downstream analysis were not yet identified. Please use the 
           identifyTFs and the identifyPWs function to do so.")}
    
    no_all_genes = dim(data_omics[[1]][[1]][[3]])[1]
    allgenes_TFs = rbindlist(data_omics[[1]][[3]][[2]][[1]])
    allgenes_TFs = allgenes_TFs[which(!is.na(allgenes_TFs))]
    for(glen in 1: length(data_omics[[1]][[1]][[1]][[2]]))
    { 
        no_de_genes = dim(data_omics[[1]][[2]][[2]][glen][[1]])[1]
        tps_TFs = rbindlist(data_omics[[1]][[3]][[2]][[glen+1]])
        tps_TFs = tps_TFs[which(!is.na(tps_TFs))]
        tps_TFs = as.vector(apply(tps_TFs,2,as.character))
        p_vals = vector()
        count_all = vector()
        count_de = vector()
        for(TF in 1: length(tps_TFs))
        {   count_all[TF] = length(which(apply(allgenes_TFs,2,as.character) == tps_TFs[TF]))
            count_de[TF] = length(which(tps_TFs == tps_TFs[TF]))
            
            p_vals[TF] = fisher.test(matrix(c(count_all[TF]- count_de[TF], 
                                              no_all_genes - no_de_genes, count_de[TF], no_de_genes), 
                                            nrow = 2, dimnames = list(number = c("TF", "genes"), 
                                                                      genes = c("not d.e.", "d.e."))))$p.value
        }
        q_vals = p.adjust(p_vals, ...)
        enr_TFs = unique(tps_TFs[which(q_vals < alpha)])
        for(k in 1: no_de_genes)
        {if(!is.na(data_omics[[1]][[3]][[2]][[glen+1]][[k]][[1]][1]))
        { data_omics[[1]][[3]][[2]][[glen+1]][[k]][,2] = NA
          colnames(data_omics[[1]][[3]][[2]][[glen+1]][[k]])= c("upstreamTFs", "enrichedTFs")
          for(s in 1: length(enr_TFs))
          { ind_match = which(data_omics[[1]][[3]][[2]][[glen+1]][[k]][,1] == enr_TFs[s])
            if(length(ind_match)!= 0)
            {data_omics[[1]][[3]][[2]][[glen+1]][[k]][ind_match,2] = 1}
          }
        }
        }
    }
    message("Transcription factor enrichment is completed.\n")
    return(data_omics)
}




#' Identify regulators of enriched transcription factors - upstream analysis.
#'
#' This function identifies the regulators upstream of the enriched/identified 
#' transcription factors in upstream analysis.
#' Converting the pathway information to a regulatory graph needs some time...
#' Warnings regarding the skipping of edges in building the regulatory graph can 
#' be ignored.
#'  
#' @param data_omics OmicsData object.
#' @param only_enriched boolean value defining if transcription factors should
#' be identifies only for enriched pathways (TRUE); or for all identified
#' pathways (FALSE); default is TRUE.
#' @param noTFs_inPW integer; only regulators in upstream pathways with more
#' than this number of TFs are identified.
#' @param order_neighbors integer specifiying the order of the neighborhood:
#' order 1 is TF plus its immediate neighbors.
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
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' identifyRsofTFs(data_omics, only_enriched = FALSE, noTFs_inPW = 1, 
#' order_neighbors = 10)
#' }
identifyRsofTFs <- function(data_omics, only_enriched = TRUE, noTFs_inPW = 2, 
                            order_neighbors = 6) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    genelists = loadGenelists()
    if(length(genelists) < length(data_omics[[2]][[1]]))
    {stop("Current working directory does not contain all genelists of 
          the selected pathway databases. Please check if all genelists are
          present in the working directory and if necessary run readPWdata
          again with loadgenelists set to FALSE.")}  
    message("Genelists of databases are loaded/generated. \n")    
    enrichedTFs = NULL
    for(glen in 1: length(data_omics[[1]][[1]][[1]][[2]]))
    { 
        tps_TFs = identTFs(data_omics, glen)
        if(only_enriched == TRUE)
        {tps_TFs = tps_TFs[which(tps_TFs[,enrichedTFs] == 1),]}
        pathway_info = identPWsofTFs(genelists, tps_TFs)
        pws_morex_TFs = selectPWsofTFs(pathway_info[[1]], pathway_info[[2]],
                                       noTFs_inPW)
        up_regulators = identRegulators(pws_morex_TFs, data_omics,
                                        order_neighbors, noTFs_inPW)
        
        data_omics[[1]][[3]][[2]][[glen+1]] = 
            c(data_omics[[1]][[3]][[2]][[glen+1]], list(pathway_info[[1]]))
        names(data_omics[[1]][[3]][[2]][[glen+1]])[dim(data_omics[[1]][[2]][[2]][[glen]])[1] +1] = 
            "upstreamPW"
        
        data_omics[[1]][[3]][[2]][[glen+1]] = 
            c(data_omics[[1]][[3]][[2]][[glen+1]], 
              as.data.frame(na.omit(as.character(unlist(up_regulators)))))
        names(data_omics[[1]][[3]][[2]][[glen+1]])[dim(data_omics[[1]][[2]][[2]][[glen]])[1] +2] = 
            "regulatorsPW"
        message("Regulators for time point ", 
                data_omics[[1]][[1]][[1]][[2]][glen] ," were identified. \n")    
    }
    message("Regulatory pathway elements of transcription factors in upstream 
        analysis are identified.\n")
    return(data_omics)
}
