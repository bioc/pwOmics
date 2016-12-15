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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics")) 
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
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


#' Identify regulators of transcription factors - upstream analysis.
#'
#' This function identifies the regulators upstream of the identified 
#' transcription factors in upstream analysis.
#' Converting the pathway information to a regulatory graph needs some time...
#' Warnings regarding the skipping of edges in building the regulatory graph can 
#' be ignored.
#'  
#' @param data_omics OmicsData object.
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
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("userspec"))
#' data_omics = readPhosphodata(data_omics, 
#' phosphoreg = system.file("extdata", "phospho_reg_table.txt", 
#' package = "pwOmics")) 
#' data_omics = readTFdata(data_omics, 
#' TF_target_path = system.file("extdata", "TF_targets.txt", 
#' package = "pwOmics"))
#' data_omics_plus = readPWdata(data_omics,  
#' loadgenelists = system.file("extdata/Genelists", package = "pwOmics")) 
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyRsofTFs(data_omics, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' }
identifyRsofTFs <- function(data_omics, noTFs_inPW = 2, 
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
        if(!length(tps_TFs) == 0)
        {
        pathway_info = identPWsofTFs(genelists, tps_TFs)                        
        pws_morex_TFs = selectPWsofTFs(pathway_info[[1]], pathway_info[[2]],
                                       noTFs_inPW)
        up_regulators = identRegulators(pws_morex_TFs, data_omics,   
                                        order_neighbors, noTFs_inPW) 
        up_regulators_bound = list()
        for(s in 1: length(pws_morex_TFs))
        {   if(!is.na(names(up_regulators[[s]]))[1])
            {up_regulators_bound[[s]] = do.call("rbind", up_regulators[[s]])
            }else{
             up_regulators_bound[[s]] = up_regulators[[s]][[1]]   
            }
        }   
        up_regulators_bound = na.omit(do.call("rbind", up_regulators_bound))
        
        data_omics[[1]][[3]][[2]][[glen+1]] = 
            c(data_omics[[1]][[3]][[2]][[glen+1]], list(pathway_info[[1]]))
        names(data_omics[[1]][[3]][[2]][[glen+1]])[dim(data_omics[[1]][[2]][[2]][[glen]])[1] +1] = 
            "upstreamPW"

        data_omics[[1]][[3]][[2]][[glen+1]][[length(data_omics[[1]][[3]][[2]][[glen+1]])+1]] = as.data.frame(up_regulators_bound)
        names(data_omics[[1]][[3]][[2]][[glen+1]])[dim(data_omics[[1]][[2]][[2]][[glen]])[1] +1] = 
            "regulatorsPW"
        message("Regulators for time point ", 
                data_omics[[1]][[1]][[1]][[2]][glen] ," were identified. \n")  
        }else{
        message("No upstream TFs and upstream regulators identified for time point ", 
                    data_omics[[1]][[1]][[1]][[2]][glen] ,". \n")      
        }
    }
    message("Regulatory pathway elements of transcription factors in upstream 
        analysis are identified.\n")
    return(data_omics)
}
