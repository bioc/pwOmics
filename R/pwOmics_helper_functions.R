#' Identification of upstream transcription factors for all genes.
#'
#' @param data_omics OmicsData object. 
#' @return data_omics OmicsData object with upstream TFs identified for all 
#' genes.
#' 
#' @keywords manip
TFidentallgenes <- function(data_omics){
    
    all_genes = unique(data_omics[[1]][[1]][[3]])
    for(glen in 1: dim(all_genes)[1])
    {data_omics[[1]][[3]][[2]]$all_gen[[glen]] = list(NA)}
    names(data_omics[[1]][[3]][[2]]$all_gen) = 
        as.character(all_genes[[1]])
    
    ind_allgenes = which(all_genes[,1] %in% data_omics[[3]][[2]][,2])
    for(f in ind_allgenes)
    { ind_match = which(data_omics[[3]][[2]][,2] == 
                            as.character(all_genes[f,1]))
      df_tf = data.frame(data_omics[[3]][[2]][ind_match,1])
      colnames(df_tf) = c("upstreamTFs")
      if(is.na(data_omics[[1]][[3]][[2]]$all_gen[[f]][[1]][1]))
      {data_omics[[1]][[3]][[2]]$all_gen[[f]] =  df_tf
      }else{
          data_omics[[1]][[3]][[2]]$all_gen[[f]] = 
              rbind(data_omics[[1]][[3]][[2]]$all_gen[[f]], df_tf)
      }
      message(paste("Upstream TFs identified for gene ",
                    as.character(all_genes[f,1]),
                    ", no. ", f, " from ", dim(all_genes)[1], "\n", sep = ""))
    }
    return(data_omics)
}


#' Identification of upstream transcription factors.
#' 
#' Identification of upstream transcription factors for the differentially 
#' expressed genes of the different timepoints.
#'
#' @param data_omics OmicsData object.
#' @return data_omics OmicsData object with upstream TFs of differentially 
#' expressed genes for separate timepoints identified.
#' 
#' @keywords manip
TFidenttps <- function(data_omics){

    data_omics[[1]][[3]][[2]] = 
        c(data_omics[[1]][[3]][[2]], 
          data_omics[[1]][[2]][[2]])
    for(glen in 1: length(data_omics[[1]][[1]][[1]][[2]]))
    {data_omics[[1]][[3]][[2]][[glen+1]] = 
     as.list(as.character(data_omics[[1]][[3]][[2]][[glen+1]][[1]]))
     names(data_omics[[1]][[3]][[2]][[glen+1]]) = 
         data_omics[[1]][[3]][[2]][[glen+1]] }
    
    for(k in 1: length(data_omics[[1]][[1]][[1]][[2]]))  
    {
        names_genes = names(data_omics[[1]][[3]][[2]][[k+1]])  
        for(no_pro in 1: length(names_genes))                 
        {if(names_genes[no_pro] %in% names(data_omics[[1]][[3]][[2]][[1]]))
        { data_omics[[1]][[3]][[2]][[k+1]][[no_pro]] = 
          data_omics[[1]][[3]][[2]][[1]][which(names(data_omics[[1]][[3]][[2]][[1]]) == 
                                                       names_genes[no_pro])][[1]]   
        }else{
            data_omics[[1]][[3]][[2]][[k+1]][[no_pro]] = NA
        }
        }
    }
    return(data_omics)
}

#' Identification of pathwayIDs and pathway names for proteins at individual
#' timepoints.
#'
#' Take the identified pathways from the list of all proteins and transfer this 
#' information for the individual timepoints.
#' 
#' @param data_omics OmicsData object. 
#' @return data_omics OmicsData object with all pathways identified for the 
#' individual timepoints. 
#' 
#' @keywords manip
PWidenttps <- function(data_omics)
{
    data_omics[[1]][[3]][[1]] = 
        c(data_omics[[1]][[3]][[1]], 
          data_omics[[1]][[2]][[1]])
    for(plen in 1: length(data_omics[[1]][[1]][[1]][[1]]))
    {data_omics[[1]][[3]][[1]][[plen+1]] = 
         as.list(as.character(data_omics[[1]][[3]][[1]][[plen+1]][,1]))
     names(data_omics[[1]][[3]][[1]][[plen+1]]) = 
         data_omics[[1]][[3]][[1]][[plen+1]] }

    for(k in 1: length(data_omics[[1]][[1]][[1]][[1]]))  
    {
        names_prot = names(data_omics[[1]][[3]][[1]][[k+1]])  
        for(no_pro in 1: length(names_prot))                  
        {data_omics[[1]][[3]][[1]][[k+1]][[no_pro]] = 
             data_omics[[1]][[3]][[1]][[1]][which(names(data_omics[[1]][[3]][[1]][[1]]) 
                                                  == names_prot[no_pro])][[1]]}   
    }
    return(data_omics)
}

#' Identification of pathwayIDs and pathway names for all proteins.
#' 
#' Identification of pathwayIDs and pathway names for all proteins.
#'
#' @param data_omics OmicsData object. 
#' @param genelists lists of genelists from chosen pathway databases.
#' @return OmicsData object with identified pathway IDs for list of all
#' proteins.
#' 
#' @keywords manip
PWidentallprots<- function(data_omics, genelists){
    
    all_prots = as.character(unique(data_omics[[1]][[1]][[2]][[1]]))
    for(plen in 1: length(unique(data_omics[[1]][[1]][[2]][[1]])))
    {data_omics[[1]][[3]][[1]]$all_prot[[plen]] = list(NA)}
    names(data_omics[[1]][[3]][[1]]$all_prot) = all_prots
    
    if("biocarta" %in% data_omics[[2]][[1]])
    { genelist_ind = 1
      datab = "biocarta"
      data_omics = PWidentallprotssub(data_omics,
                                      genelists, genelist_ind, datab)
    }
    if("kegg" %in% data_omics[[2]][[1]])
    { genelist_ind = 2
      datab = "kegg"
      data_omics = PWidentallprotssub(data_omics, 
                                      genelists, genelist_ind, datab)
    } 
    if("nci" %in% data_omics[[2]][[1]])
    { genelist_ind = 3
      datab = "nci"
      data_omics = PWidentallprotssub(data_omics, 
                                      genelists, genelist_ind, datab)
    } 
    if("reactome" %in% data_omics[[2]][[1]])
    { genelist_ind = 4
      datab = "reactome"
      data_omics = PWidentallprotssub(data_omics, 
                                      genelists, genelist_ind, datab)
    }  
    for(plen in 1: length(unique(data_omics[[1]][[1]][[2]][[1]])))
    {data_omics[[1]][[3]][[1]]$all_prot[[plen]] = 
         unique(data_omics[[1]][[3]][[1]]$all_prot[[plen]])}
    return(data_omics)
}

#' Internal subfunction for all protein pathway identification.
#'
#' Internal subfunction for all protein pathway identification.
#'
#' @param data_omics OmicsData object.
#' @param genelists lists of genelists from chosen pathway databases
#' @param genelist_ind integer specifying pathway database genelist matching;
#' 1 = biocarta, 2 = kegg, 3 = nci, 4 = reactome.
#' @param datab character vector indicating database name for message.
#' @return OmicsData object with identified pathways for each protein.
#' 
#' @keywords manip
PWidentallprotssub <- function(data_omics, genelists, genelist_ind, datab){
    
    V2 = V3 = NULL
    all_prots = unique(data_omics[[1]][[1]][[2]])
    genes_cha = apply(genelists[[genelist_ind]][,1, with = FALSE], 2, 
                      as.character)
    ind_allprots = which(all_prots[,1] %in% genes_cha)
    ind_allprotslist = which(genes_cha %in% as.character(all_prots[,1]))
    for (f in ind_allprots)
    {for(s in ind_allprotslist)
    { if(as.character(all_prots[f,1]) == genes_cha[s] )
    {df_pw = data.frame(genelists[[genelist_ind]][s,V2], 
                        genelists[[genelist_ind]][s,V3], NA)
     colnames(df_pw) = c("pathwayIDs", "pathwayNames", "enrichedPWs")
     if(is.na(data_omics[[1]][[3]][[1]]$all_prot[[f]][[1]][1]))
     {data_omics[[1]][[3]][[1]]$all_prot[[f]] =  df_pw
     }else{
         data_omics[[1]][[3]][[1]]$all_prot[[f]] = 
             rbind(data_omics[[1]][[3]][[1]]$all_prot[[f]],
                   df_pw)
     }}
    }
     message(datab, " pathways were identified for protein no. ",
             f," - ", all_prots[f,1], "\n")
    }
    return(data_omics)
}



#' Create internal IDs.
#'
#' Create new internal ids in biopax$df:
#' 
#' @param data_omics OmicsData object. 
#' @param PWinfo pathway database information from chosen pathway databases as
#' from loadPWs.
#' @return list of internal IDs from specified pathway databases.
#' 
#' @keywords manip
createIntIDs <- function(data_omics, PWinfo) {
    
    bp_biocarta = "NA"
    bp_kegg = "NA"
    bp_nci = "NA"
    bp_react_homos = "NA"
    if("biocarta" %in% data_omics[[2]][[1]])
    {PWinfo_ind = 1
     PWDBname = "biocarta"
     bp_biocarta = genIntIDs(data_omics, PWinfo, PWinfo_ind, PWDBname)}
    if("kegg" %in% data_omics[[2]][[1]])
    {PWinfo_ind = 2
     PWDBname = "kegg"
     bp_kegg = genIntIDs(data_omics, PWinfo, PWinfo_ind, PWDBname)}
    if("nci" %in% data_omics[[2]][[1]])
    {PWinfo_ind = 3
     PWDBname = "nci"
     bp_nci = genIntIDs(data_omics, PWinfo, PWinfo_ind, PWDBname)}
    if("reactome" %in% data_omics[[2]][[1]])
    {PWinfo_ind = 4
     PWDBname = "reactome"
     bp_react_homos = genIntIDs(data_omics, PWinfo, PWinfo_ind, PWDBname)}
    
    intIDs = list(bp_biocarta, bp_kegg, bp_nci, bp_react_homos)
    return(intIDs)
}


#' Internal function for generation of pathway database specific internal IDs.
#'
#' Generates new internal ids (database-specific) in biopax$df.
#' 
#' @param data_omics OmicsData object. 
#' @param PWinfo pathway database information from chosen pathway databases as
#'  from loadPWs.
#' @param PWinfo_ind integer specifying element of loadPWs output matching the
#'  chosen pathway database.
#' @param PWDBname character; pathway database name.
#' @return data.table with newly generated internal IDs of biopax model.
#' 
#' @keywords manip
genIntIDs <- function(data_omics, PWinfo, PWinfo_ind, PWDBname) {
    
    bp_mod = PWinfo[[PWinfo_ind]] 
    bp_mod_dt = bp_mod$dt
    col_names = colnames(bp_mod$dt)
    bp_mod_dt = sapply(bp_mod_dt, as.character)
    for (i in 1:dim(bp_mod_dt)[1]) 
    {
        bp_mod_dt[i,"id"] = paste("bp3_", PWDBname, "_",bp_mod_dt[i,"id"],sep="")
        if (bp_mod_dt[i,"property_attr"] == "rdf:resource") 
        {
            bp_mod_dt[i,"property_attr_value"] = paste("#bp3_", PWDBname,
                "_",sub("#","",bp_mod_dt[i,"property_attr_value"]), sep = "")
        }
    }
    bp_mod_dt = data.table(bp_mod_dt)
    for (col in col_names) 
        {set(bp_mod_dt, j=col, value=as.character(bp_mod_dt[[col]]))}
    bp_mod$dt = bp_mod_dt
    message("Internal IDs for ", PWDBname," database were generated.\n")
    
    return(bp_mod)
}



#' Read in matching transcription factor target gene pairs.
#'
#' In case the user is able to provide a file with transcription factor - 
#' target gene matches (e.g. from TRANSFAC database) this function can read in 
#' the information.
#' The file needs to be a txt file with first column transcription factors and 
#' second column target gene symbols without a header.
#' 
#' @param data_omics OmicsData object. 
#' @param TF_target_path path of txt file containing the transcription factor
#' target gene information as specified above.
#' @return data frame with user-specified TF target gene information.
#' 
#' @keywords manip
readTFtargets <- function(data_omics, TF_target_path) {
    TF_database_p = read.delim(TF_target_path, header = FALSE)
    
    if(dim(TF_database_p)[2] != 2)
    {stop("User-specified TF-target gene information file has wrong format.")}
    
    return(TF_database_p)
}


#' Get Gene Symbols from Ensemble Gene Accession IDs.
#'
#' @param ids vector of Ensemble Gene Accession IDs.
#' @return ids character vector of gene symbols.
#' 
#' @keywords manip
getAlias_Ensemble <- function(ids)
{
  ensembl = useDataset("hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "www.ensembl.org"))
  hgnc = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
               filters = "ensembl_gene_id", values = ids, mart = ensembl)
  return(hgnc)
}


#' Load pathway database information.
#'
#' This function loads the pathway information from pathway databases. 
#' Needed in the identifyPWs function.
#'  
#' @param pwdatabases vector indicating with which pathway database the 
#' pathways should be determined; possible choices
#' are "biocarta", "kegg", "nci", "reactome".
#' @param biopax_level integer indicating biopax level of pathway database 
#' information to be retrieved.
#' @return list of biopax model corresponding to specified pathway databases.
#' 
#' @keywords manip
loadPWs <- function(pwdatabases, biopax_level) {
    
    ah = AnnotationHub()
    biopax = query(ah, "NIH Pathway Interaction Database")
    if(biopax_level == 2)
    {ind_pw = which(grepl("bp2", biopax$title)== TRUE)
    }else{
     stop("Biopax level is not supported momentarily.")
    }
    
    if("biocarta" %in% pwdatabases)
    {bp_biocarta = biopax[[grep("BioCarta", biopax$title)[which(grep("BioCarta", biopax$title) %in% ind_pw)][[1]] ]]
     }else{
     bp_biocarta = NA}
    if("kegg" %in% pwdatabases)
    {bp_kegg = biopax[[grep("KEGG", biopax$title)[which(grep("KEGG", biopax$title) %in% ind_pw)][[1]] ]]
     }else{
     bp_kegg = NA}
    if("nci" %in% pwdatabases)
    {bp_nci = biopax[[grep("NCI-Nature_Curated", biopax$title)[which(grep("NCI-Nature_Curated", biopax$title) %in% ind_pw)][[1]] ]]
     }else{
     bp_nci = NA}
    if("reactome" %in% pwdatabases)
    {bp_react_homos = biopax[[grep("Reactome", biopax$title)[which(grep("Reactome", biopax$title) %in% ind_pw)][[1]] ]]
     }else{
     bp_react_homos = NA}
    return(list(bp_biocarta, bp_kegg, bp_nci, bp_react_homos))
}

#' Generate genelists from pathway databases.
#'
#' This function generates genelists from the chosen pathway databases for 
#' further processing in detPathways function.
#'  
#' @param intIDs list containing Biopax model with newly generated internal IDs
#'  as processed with genintIDs function.
#' The components of the list are biopax models for "biocarta", "kegg", "nci",
#'  "reactome". In case a database was not chosen the list 
#' entry contains a NA.
#' @param pwdatabases vector indicating with which pathway database the pathways
#'  should be determined; possible choices
#' are "biocarta", "kegg", "nci", "reactome".
#' @return list of genelists of specified pathway databases.
#' 
#' @keywords manip
genGenelists <- function(intIDs, pwdatabases) {
    
    Genelist_biocarta_all = "NA"
    Genelist_kegg_all = "NA"
    Genelist_nci_all = "NA"
    Genelist_react_homos_all = "NA"
    
    if("biocarta" %in% pwdatabases){
        database_int = 1
        PWDB_name = "biocarta"
        Genelist_biocarta_all = genGenelistssub(intIDs, database_int, PWDB_name)
    }   
    if("kegg" %in% pwdatabases){
        database_int = 2
        PWDB_name = "kegg"
        Genelist_kegg_all = genGenelistssub(intIDs, database_int, PWDB_name)
    }   
    if("nci" %in% pwdatabases){
        database_int = 3
        PWDB_name = "nci"
        Genelist_nci_all = genGenelistssub(intIDs, database_int, PWDB_name)
    }   
    if("reactome" %in% pwdatabases){
        database_int = 4
        PWDB_name = "reactome"
        Genelist_react_homos_all = genGenelistssub(intIDs, database_int,
                                                   PWDB_name)
    }   
    return(list(Genelist_biocarta_all, Genelist_kegg_all, Genelist_nci_all,
                Genelist_react_homos_all))
}

#' Generate internally genelists from pathway databases.
#'
#' This function generates genelists for a particular pathway database for
#'  further processing in detPathways function.
#'  
#' @param intIDs list containing Biopax model with newly generated internal IDs
#'  as processed with genintIDs function.
#' The components of the list are biopax models for "biocarta", "kegg", "nci",
#'  "reactome". In case a database was not chosen the list 
#' entry contains a NA.
#' @param database_int integer indicating database entry in indIDs (output of
#' genintIDs); biocarta = 1, kegg = 2, nci = 4, reactome = 4.
#' @param PWDB_name character; pathway database name.
#' @return data.table of genelist of particular pathway database.
#' 
#' @keywords manip
genGenelistssub <- function(intIDs, database_int, PWDB_name) {
    
    pathways = listPathways(intIDs[[database_int]])
    Genelist = data.table()
    Genelist_all = data.table()
    for(d in 1:dim(pathways)[1])
    { 
        Genelist = pathway2Geneset(intIDs[[database_int]], 
                                   pwid = pathways[d,1])$name
        Genelist = Genelist[which(Genelist!="")]
        Genelist = cbind(Genelist, pathways[d,1],
                         listInstances(intIDs[[database_int]],
                                       pathways[d,1])$name)
        
        if(dim(Genelist)[1] != 1)
        {Genelist_all = rbind(Genelist_all, as.data.table(Genelist))}
        message( PWDB_name, " Genelist is generated, pathway no. ", d, "\n")
    }
    save(Genelist_all, file = paste(getwd(), "/Genelist_",
                                    PWDB_name, ".RData", sep = "") ) 
    return(Genelist_all)
}



#' Create a new Biopax model containing all database information.
#'
#' This function creates a new biopax model depending on which 
#' pathway databases are chosen for analysis.
#'  
#' @param intIDs output list of genintIDs function.
#' @param pwdatabases vector indicating with which pathway database the
#' pathways should be determined; possible choices are "biocarta", "kegg",
#' "nci", "reactome".
#' @return biopax object generated from the specified pathway databases.
#' 
#' @keywords manip
createBiopaxnew <- function(intIDs, pwdatabases) {
    all_databases = createBiopax(level = 2)
    if("biocarta" %in% pwdatabases)
    {all_databases$dt = rbind(all_databases$dt, intIDs[[1]]$dt)}
    if("kegg" %in% pwdatabases)
    {all_databases$dt = rbind(all_databases$dt, intIDs[[2]]$dt)}
    if("nci" %in% pwdatabases)
    {all_databases$dt = rbind(all_databases$dt, intIDs[[3]]$dt)}
    if("reactome" %in% pwdatabases)
    {all_databases$dt = rbind(all_databases$dt, intIDs[[4]]$dt)}  
    return(all_databases)
}


#' Loading of genelists
#'
#' This function automatically loads the genelists corresponding to the selected
#' pathway databases stored as RData file in the current working directory.
#' 
#' @return genelist of specified pathway database.
#' 
#' @keywords manip
loadGenelists <- function()
{   
    files_genelists = list.files(path = getwd(), pattern = "*.RData")
    genelists = list()
    for(s in 1: length(files_genelists))
    { genelists[[s]] = base::get(load(paste(getwd(), files_genelists[s], sep = "/")))}
    return(genelists)
}

#' Identify overlapping upstream regulators of x transcription factors
#'
#' @param pws_morex_TFs list of transcription factors in identified pathways.
#' @param data_omics OmicsData object.
#' @param order_neighbors integer specifiying the order of the neighborhood: 
#' order 1 is TF plus its immediate neighbors.
#' @param noTFs_inPW integer; only regulators in upstream pathways with more 
#' than this number of TFs are identified.
#' @return list of possible proteomic regulators. 
#' 
#' @keywords manip
identRegulators <- function(pws_morex_TFs, data_omics, order_neighbors,
                            noTFs_inPW) {
    regulators = list()
    regulators[1: length(pws_morex_TFs)] = NA
    names(regulators) = names(pws_morex_TFs)
    for(pwxTFs in 1: length(pws_morex_TFs))
    { 
        neighbors = findxnextneighbors(data_omics, pws_morex_TFs, pwxTFs,
                                       order_neighbors)
        if(length(neighbors) > 0)
        {
            if(noTFs_inPW > 1)
            {regulators[[pwxTFs]] = findxneighborsoverlap(neighbors, noTFs_inPW,
                                                     regulators[[pwxTFs]])
            }else{
             regulators[[pwxTFs]] = neighbors
            }
        }else{
         regulators[[pwxTFs]] = NA    
        }
        if(is.na(regulators[[pwxTFs]])[[1]])
        {message("Pathway '", names(pws_morex_TFs[pwxTFs]), 
                 "' is checked for regulators...\n The pathway seems to have no inhibitory/activating pathway components.\n")
        }else{
         message("Pathway '", names(pws_morex_TFs[pwxTFs]), 
                 "' is checked for regulators...\n Regulators identified: ",
                 regulators[[pwxTFs]], "\n")   
        }
    } 
    return(regulators)
}


#' Find overlap of next neighbors of transcription factors in identified
#' pathways. 
#' 
#' Find the overlap of x next neighbors of transcription factors in identified
#' pathways. Writes the overlap into a given list called 'regulators'.
#'
#' @param neighbors list of x next neighbors for each transcription factor in
#' the pathway as provided by findxnextneighbors function.
#' @param regul list element of regulators list for current pathway.
#' @param noTFs_inPW numeric value specifying number of TFs being at least part 
#' of the pathway.
#' @return list of regulators identified in x next neighbors of TFs.
#' 
#' @keywords manip
findxneighborsoverlap <- function(neighbors, noTFs_inPW, regul) {
    if(length(unlist(neighbors))>0)
    {for(nover in 1: length(unlist(neighbors)))
    { neighbors_found = 0
      for(no_TFs in 1: length(neighbors))
      {if(unlist(neighbors)[nover] %in% neighbors[[no_TFs]] )
      {neighbors_found = neighbors_found + 1
       if(neighbors_found > (noTFs_inPW-1) )
       {regul = c(regul, neighbors[[no_TFs]][which(neighbors[[no_TFs]] == 
                                                unlist(neighbors)[nover])])}
      }
      }
    }
    }
    regul = unique(regul[2: length(regul)])
    return(regul)
}

#' Find next neighbors of transcription factors in identified pathways.
#'  
#' Produces a list of x next neighbors for each transcription factor in the
#' pathway.
#'
#' @param data_omics OmicsData object.
#' @param pws_morex_TFs list of transcription factors in identified pathways.
#' @param pwxTFs numeric variable of pathway currently investigated
#'  (from pws_morexTFs).
#' @param order_neighbors integer specifiying the order of the neighborhood:
#' order 1 is TF plus its immediate neighbors.
#' @return list of x next neighbors for each TF in the pathway. 
#' 
#' @keywords manip
findxnextneighbors <- function(data_omics, pws_morex_TFs, pwxTFs,
                               order_neighbors) {
    
    mygraph = pathway2RegulatoryGraph(data_omics[[2]][[2]],
                                      names(pws_morex_TFs)[pwxTFs], 
                                      expandSubpathways = TRUE,
                                      splitComplexMolecules = TRUE, 
                                      useIDasNodenames = FALSE, verbose = TRUE)
    
    if(class(mygraph)== "graphNEL")
    {   mygraph_i = igraph.from.graphNEL(mygraph, name = TRUE)
        neighbors = list()
        for(neigh in 1: length(pws_morex_TFs[[pwxTFs]]))
        { if(pws_morex_TFs[[pwxTFs]][neigh] %in%
                 get.vertex.attribute(mygraph_i,"name"))
        {neighbors[[neigh]] = get.vertex.attribute(mygraph_i,"name",index = 
                              unlist(neighborhood(mygraph_i,order_neighbors,
                              pws_morex_TFs[[pwxTFs]][neigh],  mode="in")))
         neighbors[[neigh]] = neighbors[[neigh]][which(neighbors[[neigh]] != 
                              pws_morex_TFs[[pwxTFs]][neigh]) ]
        }else{
            neighbors[[neigh]] = NA
        }
          names(neighbors)[neigh] = pws_morex_TFs[[pwxTFs]][neigh]
        }
    }else{
    neighbors = list()    
    }
    return(neighbors)
}

#' Select pathways with more than x TFs
#'
#' @param pathway_list first element of list returned from identPWsofTFs
#' function; contains a list of pathways.
#' @param pathway_frame second element of list returned from identPWsofTFs 
#' function; contains a data.frame of pathways.
#' @param noTFs_inPW numeric value specifying number of TFs being at least
#' part of the pathway.
#' @return list of pathways with more than x TFs.
#' 
#' @keywords manip
selectPWsofTFs <- function(pathway_list, pathway_frame, noTFs_inPW) {
    countTFsinPW = 0
    pws_more_TFs = list()
    ind_no_TFs = vector()
    for(TF_found in 1: dim(unique(pathway_frame))[1])
    { pws_more_TFs[[TF_found]] = NA
      names(pws_more_TFs)[TF_found] = 
          as.character(unique(pathway_frame[TF_found,1]))
      for(pw_found in 1: length(pathway_list))
      {if(as.character(unique(pathway_frame[TF_found,1])) %in% 
              pathway_list[[pw_found]][,1])
      {pws_more_TFs[[TF_found]] = append(pws_more_TFs[[TF_found]], 
                                          names(pathway_list)[pw_found])}
      }
      pws_more_TFs[[TF_found]] = 
          pws_more_TFs[[TF_found]][2: length(pws_more_TFs[[TF_found]])]
      if(length(pws_more_TFs[[TF_found]]) >= noTFs_inPW )
      {ind_no_TFs[TF_found] = TF_found}
    }
    ind_no_TFs = na.omit(ind_no_TFs)
    pws_more_TFs = pws_more_TFs[ind_no_TFs]
    return(pws_more_TFs)
}

#' Identification of pathways containing the transcription factors identified in 
#' upstream analysis
#'
#' @param genelists data.table as read/loaded by loadGenelist function.
#' @param tps_TFs data.table of upstream transcription factors and the flag for 
#' enrichment as returned from identTFs function.
#' @return list with first element being a pathway list and second being a 
#' pathway dataframe of pathways including the TFs of the specified timepoint.
#' 
#' @keywords manip
identPWsofTFs <- function(genelists, tps_TFs) {
    upstreamTFs = NULL
    genelist_n = apply(rbindlist(genelists),2,as.character)
    pathway_list = list()
    pathway_frame = data.frame()
    PWTF_tp_NA = vector()
    for(k in 1: dim(tps_TFs)[1])
    { ind_match = which(as.character(genelist_n[,1]) == 
                            as.character(tps_TFs[,upstreamTFs])[k]) 
      if(length(ind_match)> 0)
      {if(length(ind_match) == 1)
      { df_to_bind = data.frame(matrix(unique(genelist_n[ind_match,2:3]),
                                       ncol = 2))
        colnames(df_to_bind) = c("V2","V3")
        pathway_frame = rbind(pathway_frame, df_to_bind)
        pathway_list[[k]] = df_to_bind
      }else
      {pathway_frame = rbind(pathway_frame, unique(genelist_n[ind_match,2:3]))
       pathway_list[[k]] = unique(genelist_n[ind_match,2:3])}
      }else{
          pathway_list[[k]] = NA
      }
      PWTF_tp_NA[k] = is.na(pathway_list[[k]][[1]][1])
    }
    names(pathway_list) = as.character(tps_TFs[,upstreamTFs])
    pathway_list  = pathway_list[!PWTF_tp_NA]
    return(list(pathway_list, pathway_frame))
}

#' This function provides a data.table of upstream transcription factors and the 
#' flag for enrichment.
#'  
#' @param data_omics OmicsData object.
#' @param glen numeric value; identifier for current timepoint.
#' @return data.table of upstream TFs and an enrichment flag.
#' 
#' @keywords manip
identTFs <- function(data_omics, glen) {
    no_de_genes = dim(data_omics[[1]][[2]][[2]][glen][[1]])[1]
    genes_tp_NA = vector()
    for(g in 1: no_de_genes)
    {genes_tp_NA[g] = is.na(data_omics[[1]][[3]][[2]][[glen+1]][[g]][[1]][1])}
    tps_TFs = 
        unique(rbindlist(data_omics[[1]][[3]][[2]][[glen+1]][!genes_tp_NA]))
    return(tps_TFs)
}

#' Prepare OmicsData object for pathway information.
#'
#' This function identifies the TFs in the pathway genes and determines their
#' target genes on basis of the given (chosen) TF-target database(s).
#'  
#' @param data_omics OmicsData object.
#' @param temp_genelist dataframe of unique gene IDs in enriched/not enriched
#' PWs.
#' @return list with first element being a genelist of the pathways and second
#' being a target gene list of TFs.
#' 
#' @keywords manip
identTFTGsinPWs <- function(data_omics, temp_genelist) {
    
    temp_genelist$TFs_PW = NA
    temp_targetlist = vector()
    for(gene_no in 1: dim(temp_genelist)[1])
    { if(as.character(temp_genelist[gene_no,1]) %in% 
             as.character(data_omics[[3]][[2]][,1]))
      {temp_genelist[gene_no,2] = 1  #identification TF
       match_gene = which(as.character(data_omics[[3]][[2]][,1]) == 
                            as.character(temp_genelist[gene_no,1]))
       if(length(match_gene)>0)
       {temp_targetlist = 
          c(temp_targetlist, as.character(data_omics[[3]][[2]][match_gene,2]))}
      }else{
       temp_genelist$TFs_PW[gene_no] = NA
    } 
    }
    return(list(temp_genelist, temp_targetlist))
}

#' Prepare OmicsData object for pathway information.
#'
#' This function prepares the OmicsData object for the identified pathway
#' information.
#'  
#' @param data_omics OmicsData object.
#' @param plen integer indicating the timepoint currently investigated.
#' @return list of OmicsData object, current timepoint and pathways of interest.
#' 
#' @keywords manip
preparePWinfo <- function(data_omics, plen) {
    no_de_prots = dim(data_omics[[1]][[2]][[1]][plen][[1]])[1]
    PW_tp_NA = vector()
    for(g in 1: no_de_prots)
    {PW_tp_NA[g] = is.na(data_omics[[1]][[3]][[1]][[plen+1]][[g]][[1]][1])}
    temp_prot_list = data_omics[[1]][[3]][[1]][[plen+1]][1:no_de_prots]
    tps_PWs = rbindlist(temp_prot_list[!PW_tp_NA])
    PWofinterest = data_omics[[1]][[3]][[1]][[plen+1]][!PW_tp_NA] 
    return(list(data_omics,tps_PWs, PWofinterest))
}
