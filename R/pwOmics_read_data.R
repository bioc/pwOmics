#' Read in omics data.
#'
#' This function reads omics data: differentially expressed genes and relatively
#' differentially abundant proteins with corresponding fold changes for each 
#' time point. 
#'
#' @param tp_prots numeric vector of protein timepoints used in experiment; 
#' in case of single time point experiments simply assign an experiment number
#' (e.g. 1).
#' @param tp_genes numeric vector of gene/transcript timepoints used in 
#' experiment; in case of single time point experiments simply assign an 
#' experiment number (e.g. 1).
#' @param omics list containing protein and gene IDs and fold changes: 
#' Input are two lists, one for protein data, one for gene data,
#' both including a vector of all measured IDs as first element and a list as 
#' second element including for each tp a dataframe with IDs and foldchanges 
#' per timepoint.
#' @param PWdatabase character vector of pathway database names which should be
#' used for pathway identification, possible choices are "biocarta", "kegg", 
#' "nci", "reactome".
#' @param TFtargetdatabase character vector of TF target database names which 
#' should be used for transcription factor/target gene identification,
#' possible choices are "chea", "pazar", "userspec". In case the user is able 
#' to provide an own list of transcription factor target gene matches, he can
#' indicate this via "userspec".
#' @return OmicsData object:
#' list of 4 elements (OmicsD, PathwayD, TFtargetsD, 
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
#' \dontrun{
#' data_omics = readOmics(tp_prots = c(0.25, 1, 4, 8, 13, 18, 24),
#' tp_genes = c(1, 4, 8, 13, 18, 24), OmicsExampleData,
#' PWdatabase = c("biocarta", "kegg", "nci", "reactome"), 
#' TFtargetdatabase = c("chea", "pazar", "userspec"))
#' }
readOmics <- function(tp_prots, tp_genes, omics, PWdatabase, TFtargetdatabase) {
    
    OmicsData = list(OmicsD = list(OmicsDescription = list(
        timepoints = list(Proteins = tp_prots, Genes = tp_genes), 
        allProteinIDs = unique(omics$P[[1]]),
        allGeneIDs = unique(omics$G[[1]])),
                                   OmicsDataSet = list(
                                       ProteinData = list(),
                                       GeneData = list()),
                                   OmicsResults = list(
                                       ProteinData = list(),
                                       GeneData = list())),
                     PathwayD = list(PWdatabaseNames = PWdatabase, 
                                     BiopaxModel = list()),
                     TFtargetsD = list(TFtargetdatabaseNames = TFtargetdatabase,
                                       TFtargetData = list()),
                     Status = 1)
    uni_protdata = data.frame()
    for(plen in 1: length(tp_prots))
    {uni_protdata = omics$P[[2]][[plen]][isUnique(omics$P[[2]][[plen]][,1]),]
     OmicsData$OmicsD$OmicsDataSet$ProteinData[[plen]] = uni_protdata}
    names(OmicsData$OmicsD$OmicsDataSet$ProteinData) = paste("tp", tp_prots, 
                                                             sep = "")
    uni_genedata = data.frame()
    for(glen in 1: length(tp_genes))
    {uni_genedata = omics$G[[2]][[glen]][isUnique(omics$G[[2]][[glen]][,1]),]
     OmicsData$OmicsD$OmicsDataSet$GeneData[[glen]] = uni_genedata}
    names(OmicsData$OmicsD$OmicsDataSet$GeneData) = paste("tp", tp_genes, 
                                                          sep = "")
    
    class(OmicsData) = "OmicsData"
    return(OmicsData)
}


#' Reads in chosen transcription factor target database information.
#'
#' This function reads in transcription factor information given the selected 
#' transcription factor target gene database. The information is downloaded
#' via the AnnotationHub package and merged, if necessary.
#'
#' @param data_omics OmicsData object.
#' @param TF_target_path character vector indicating path of the txt file of 
#' matching transcription factors and target genes; the file should be a txt 
#' file with first column transcription factors and second column target gene 
#' symbols without a header.
#' @param cell_match character indicating the cell line/cells for which the TF target 
#' gene data should be extracted from the database; this is only possible for 
#' chea database. Available cell-specific data from chea for matching are 
#' "Hs578T", "Raji B cells and iDC", "MCF7", "THP-1", "Hela cells", "STHdh",
#' "H3396 breast cancer cells","HL60","HESC","T-ALL", "HPC-7", 
#' "ovarian surface epithelium", "HaCaT", "HCT116","U2OS", "Wilms tumor-derived CCG99-9611",
#' "HepG2","HUMAN INTESTINAL CELL LINE CACO-2", "HEK293T","K562", "AK7",
#' "NEUROBLASTOMA","JURKAT","T-47D","LS174T", "MULTIPLE HUMAN CANCER CELL TYPES",
#' "501MEL", "PC3", "CACO-2", "FETAL_BRAIN", "HELA", "U937_AND_SAOS2", 
#' "CD4_POS_T", "ERYTHROLEUKEMIA", "RHABDOMYOSARCOMA", "293T", "SW620",
#' "LYMPHOBLASTOID", "VCAP", "SK-N-MC", "CADO-ES1", "MEDULLOBLASTOMA", "M12",
#' "K562_HELA_HEPG2_GM12878", "NT2", "SHEP-21N", "LN229_GBM", "MCF-7", "MELANOMA",
#' "MYOFIBROBLAST", "NTERA2", "MEGAKARYOCYTES", "HMVEC", "ZR75-1", "TREG",
#' "TLL", "A2780", "MONOCYTES", "BEAS2B", "LNCAP PROSTATE CANCER CELL LINES", 
#' "MCF10A", "GC-B", "BL", "IMR90", "EOC", "PCA", "PROSTATE_CANCER", "OVCAR3",
#' "MALME-3M", "HFKS", "HEK293", "HELA-AND-SCP4", "CD34+", "IB4-LCL", "MDA-MB-231",
#' "U87", "T47D", "Z138-A519-JVM2", "DLD1", "ATHEROSCLEROTIC-FOAM", "LCL-AND-THP1",
#' "NB4", "PFSK-1 AND SK-N-MC", "EP156T","GBM1-GSC","CD4+", "FIBROSARCOMA",
#' "LGR5+ INTESTINAL STEM CELL","NEUROBLASTOMA BE2-C". 
#' If no tissue is given the data from all cells/cell lines are merged.
#' @param TF_filter_threshold integer defining a threshold number to 
#' filter out those transcription factors having higher numbers of target genes
#' than 'TF_filter_threshold' from the further analysis
#' 
#' @return OmicsData object - a list containing information about the user data 
#' (timepoints, IDs, fold changes) and the selected databases chosen for the 
#' analysis.
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
#' TFtargetdatabase = c("pazar"))
#' data_omics = readTFdata(data_omics)
#' data_omics[[3]]
readTFdata <- function(data_omics, TF_target_path, cell_match = 0, TF_filter_threshold = 0) {
    
    if(class(data_omics) != "OmicsData")
    { stop("Parameter 'data_omics' is not an OmicsData object.\n")} 
    
    if(!"chea" %in% data_omics[[3]][[1]] & cell_match != 0)
    { stop("Matching of cell line/ cells option is only available for the chea
           database.\n")} 
    
    if("chea" %in% data_omics[[3]][[1]] | "pazar" %in% data_omics[[3]][[1]])
    {ah = AnnotationHub()}
    
    TF_data_comb = data.frame()
    if("chea" %in% data_omics[[3]][[1]])
    { 
        chea = query(ah, "ChEA")
        TF_data_chea = chea[[1]]
        TF_data_chea = TF_data_chea[which(TF_data_chea[,"Species"]=="human" | 
                                              TF_data_chea[,"Species"]=="HUMAN"),]
        if(cell_match != 0 && is.character(cell_match))
        {TF_data_chea = TF_data_chea[which(TF_data_chea[,"CellType"] == cell_match),]}
        TF_data_chea = data.frame(TF = TF_data_chea[,2], 
                                  target = TF_data_chea[,4])
        TF_data_comb = rbind(TF_data_comb,TF_data_chea)
        message("ChEA database information was read.\n")
    }
    if("pazar" %in% data_omics[[3]][[1]])
    { 
        pazar = query(ah, "Pazar")
        TF_data_pazar = list()
        ret_list = c(1:2, 4:35, 37:46, 48:49, 51:62, 64:69, 71, 73:76, 78:84, 88: length(pazar))
        for(k in ret_list)
        {   if(class(pazar[[k]]) == "GRanges")
        { TF_data_pazar[[k]] = as.data.frame(pazar[[k]])
        }else{
            TF_data_pazar[[k]] = pazar[[k]]
        }
        }
        TF_data_pazar = rbindlist(TF_data_pazar,use.names = TRUE, fill = TRUE)
        TF_data_pazar = data.frame(TF_data_pazar)
        TF_data_pazar = TF_data_pazar[which(TF_data_pazar[,"Species"]=="Homo sapiens"),]
        target = getAlias_Ensemble(TF_data_pazar[,"EnsemblGeneAccession"])
        target = target[match(TF_data_pazar[,"EnsemblGeneAccession"], target[,1] ),2]
        TF_data_pazar = data.frame(TF = gsub("_HUMAN", "", TF_data_pazar[,"TFName"]), 
                                   target)
        TF_data_pazar = unique(TF_data_pazar)
        mouse_ind = grep("_MOUSE", 
                         apply(as.data.frame(TF_data_pazar),2,as.character)[,1])
        rat_ind = grep("_RAT", 
                       apply(as.data.frame(TF_data_pazar),2,as.character)[,1])
        TF_data_pazar = TF_data_pazar[-c(mouse_ind, rat_ind),]
        TF_data_comb = rbind(TF_data_comb, TF_data_pazar)
        message("Pazar database information was read.\n")
    }  
    if("userspec" %in% data_omics[[3]][[1]])
    { if(!is.character(TF_target_path))
    {stop("No path for user specified TF target gene list was given. 
          Either provide a pathname and the file or chose
          another transcription factor target database.")
    }else{
        TF_data = readTFtargets(data_omics, TF_target_path)
        colnames(TF_data)=c("TF", "target")
        TF_data_comb = rbind(TF_data_comb, TF_data)}
      message("User specified database information was read.\n")
    }
    
    if(is.numeric(TF_filter_threshold) & TF_filter_threshold != 0)
    {TF_matrix = matrix(nrow = length(as.character(unique(TF_data_comb[,1]))), ncol = 2)
    TF_matrix[,1] = as.character(unique(TF_data_comb[,1]))
    for(k in 1: length(unique(TF_data_comb[,1])))
    { TF_matrix[k,2] = length(which(TF_data_comb[,1] == TF_matrix[k,1]))}
    TF_matrix = TF_matrix[order(as.numeric(TF_matrix[,2])),]
    TF_matrix = TF_matrix[-which(as.numeric(TF_matrix[,2])> TF_filter_threshold),]
    TF_data_comb = TF_data_comb[which(TF_data_comb[,1] %in% TF_matrix[,1]),]}

    data_omics[[3]][[2]] = TF_data_comb
    #increment status
    if(length(data_omics[[2]][[2]])!=0 & length(data_omics[[3]][[2]]) != 0)
    {data_omics[[4]] = data_omics[[4]] +1
    }else{
        data_omics[[4]] = data_omics[[4]] 
    }
    return(data_omics)
}

#' Read in pathway database data needed for pathway identification.
#'
#' This function reads pathway data of the chosen database(s) via the 
#' AnnotationHub [1] package and rBiopaxParser [2] package. 
#' Takes a lot of time for a high number of proteins and/or if all databases 
#' are chosen.
#' First, chosen databases are retrieved, then new internal pathway IDs are 
#' generated. Afterwards the genelists of the different databases are loaded or
#' generated, depending on the loadgenelists option. Pathway ID mapping takes 
#' some time, especially for such big databases as reactome, so the genelists 
#' are automatically stored in the current working folder and can be used 
#' via loadgenelists in case you use this function again for easier and faster 
#' usage...
#' Biopax level of retrieved databases is 2 by default.
#'   
#' @param data_omics OmicsData object.
#' @param loadgenelists path of genelist RData files stored previously; all 
#' genelists stored in this path are read in and used automatically if path is 
#' given; if loadgenelists = FALSE, then genelists from pathway databases have 
#' to be generated first.
#' @param biopax_level integer indicating biopax level of pathway database
#' information. default level is 2.
#'  
#' @references 1. Morgan M, Carlson M, Tenenbaum D and Arora S. AnnotationHub: 
#' Client to access AnnotationHub resources. R package version 1.99.75. 
#' @references 2. Kramer F, Bayerlova M, Klemm F, Bleckmann A and Beissbarth T. 
#' rBiopaxParser - an R package to parse, modify and visualize BioPAX data.
#' 
#' @return list of OmicsData object and genelists for selected pathway databases.
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
#' data_omics_plus = readPWdata(data_omics, loadgenelists = "Genelists")
#' data_omics_plus[[2]][[1]]
#' }
readPWdata <- function(data_omics, loadgenelists, biopax_level = 2) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    db_info = loadPWs(data_omics[[2]][[1]], biopax_level)
    message("Pathway database information is loaded. \n")
    
    intIDs = createIntIDs(data_omics, db_info)
    message("New internal IDs for pathways were generated. \n")
    
    data_omics[[2]][[2]] = createBiopaxnew(intIDs, data_omics[[2]][[1]])   
    message("New Biopax model was created with new internal pathway IDs. \n")
    
    if(loadgenelists == FALSE)
    { genelists = genGenelists(intIDs, data_omics[[2]][[1]])
    }else{
        files_genelists = list.files(path = loadgenelists, pattern = "*.RData")
        if(sum(file.exists(paste(loadgenelists,"/Genelist_", 
                                 data_omics[[2]][[1]], ".RData", sep = ""))) == 0)  
        {stop("No genelist files corresponding to your pathway database
              choice have been found in the working directory. Check your 
              working folder or run readPWdata with loadgenelists = FALSE
              to generate these files.")}
        genelists = list()
        for(s in 1: length(files_genelists))
        {genelists[[s]] = get(load(paste(loadgenelists, files_genelists[s], 
                                         sep = "/")))}
    }  
    message("Genelists of databases are loaded/generated. \n")    
    
    if(length(data_omics[[2]][[2]])!=0 & length(data_omics[[3]][[2]]) != 0)
    {data_omics[[4]] = data_omics[[4]] +1
    }else{
        data_omics[[4]] = data_omics[[4]] 
    }
    return(list(data_omics, genelists))
}

