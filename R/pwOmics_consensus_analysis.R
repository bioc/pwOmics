#' Static analysis.
#' 
#' Identify for each corresponding timepoint of the two datasets the consensus
#' network. Protein intersection of the omics data and TF intersection are 
#' linked via SteinerTree algorithm applied on STRING protein-protein
#' interaction database. 
#' The Steiner tree algorithm refers to the shortest path heuristic algorithm
#' of [1,2].
#' Target genes of this consensus network are identified via the chosen
#' TF-target gene database(s). Please note that the consensus graphs can be
#' different as in the Steiner Tree algorithm the start terminal node is 
#' picked arbitrarily and there are always several shortest path distances. 
#' By default the same time points of both data sets are considered.
#'  
#' @references 1. Path heuristic and Original path heuristic, Section 4.1.3 of 
#' the book "The Steiner tree Problem", Peter L. Hammer
#' @references 2. "An approximate solution for the Steiner problem in graphs", 
#' H Takahashi, A Matsuyama
#'  
#'  
#' @param data_omics OmicsData object.
#' @param run_times integer specifying number of times to run SP Steiner tree 
#' algorithm to find minimal graph, default is 3.
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection. Type of checking is defined with
#' parameter 'phospho'.
#' @param phospho boolean value; TRUE in case up- and downregulation should be
#' checked based on provided downstream phosphoprotein influence from 
#' identifyPR function; FALSE in case up- and downregulation should be checked
#' for without phosphoprotein database knowledge. Default is TRUE.
#' @param tp_prot integer specifying the time point that should be included 
#' into the static consensus net for the phosphoprotein data
#' @param tp_gene integer specifying the time point that should be included 
#' into the static consensus net for the transcriptome data
#' 
#' @return list of igraph objects; length corresponds to number of overlapping 
#' time points from upstream and downstream analysis.
#' @keywords manip
#' @export
#' @examples
#' \dontrun{
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
#'
#' data_omics_plus = identifyPR(data_omics_plus)    
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' data_omics = identifyPWTFTGs(data_omics)
#' statConsNet = staticConsensusNet(data_omics)
#' }
staticConsensusNet <- function(data_omics, run_times = 3, updown = FALSE,
                               tp_prot = NULL, tp_gene = NULL, phospho = TRUE) {
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    requireNamespace("STRINGdb", quietly = TRUE)
    string_db = STRINGdb$new(version = "10", species = 9606,
                             score_threshold = 0, input_directory = "") 
    PPI_graph = getSTRING_graph(string_db) 
    
    if(is.null(tp_prot) & is.null(tp_gene))
        {same_tps = data_omics[[1]][[1]][[1]][[1]][which(data_omics[[1]][[1]][[1]][[1]] 
                                                        %in% data_omics[[1]][[1]][[1]][[2]])]
        if(length(same_tps) == 0)
        {stop("No matching time points were found in the two 
              corresponding data sets.")}
        
        ST_proteins = vector()
        consensus_graph = list()
        for(tps in 1: length(same_tps))
        {
            STRINGIDs = getConsensusSTRINGIDs(data_omics, tps, string_db, updown, phospho)
            if(dim(STRINGIDs)[1]>0)
            {message("Consensus graph for time point ", same_tps[tps], 
                    " was generated.\n")   
    
            ST_net = SteinerTree_cons(STRINGIDs$STRING_id, 
                                      PPI_graph, run_times) 
            ST_net = getAliasfromSTRINGIDs(data_omics, ST_net, updown, phospho,
                                           STRINGIDs, tps, string_db)
            
            V(ST_net)$label.cex = 1
            plot(ST_net, main = paste("Steiner tree of consensus graph\n time ", 
                                      same_tps[tps], sep = ""))
            legend(x = 0, y = -1.2 , legend = c("consensus proteins",
                                                "steiner node proteins", "consensus TFs"), 
                   fill = c("red", "yellow", "lightblue"), cex = 0.7)
            message("Steiner tree of consensus graph for time point ", same_tps[tps],
                    " was build.\n")
            
            ST_TFTG = getbipartitegraphInfo(data_omics, tps, updown, phospho)
            ST_net_targets = genfullConsensusGraph(ST_net, ST_TFTG)
            V(ST_net_targets)$label.cex = 0.6
            plot.igraph(ST_net_targets, main = paste("Consensus graph\n time ", 
                                                     same_tps[tps], sep = ""), vertex.size = 18)
            if("yellow" %in% V(ST_net_targets)$color)
            {legend(x = 0, y = -1.2 , legend = c("consensus proteins", 
                "steiner node proteins", "consensus TFs", "consensus target genes"), 
                   fill = c("red", "yellow", "lightblue", "green"), cex = 0.7)
            }else{
            legend(x = 0, y = -1.2 , legend = c("consensus proteins", 
                     "consensus TFs", "consensus target genes"), 
                     fill = c("red", "lightblue", "green"), cex = 0.7)   
            }
            
            ST_net_targets = addFeedbackLoops(ST_net_targets)
            consensus_graph[[tps]] = ST_net_targets
            }else{
            consensus_graph[[tps]] = NA
                message("Not enough intersecting molecules to generate
                        a consensus graph for time point ", same_tps[tps], "\n")   
            }
        }
        names(consensus_graph) = as.character(same_tps)
    }else{
        ST_proteins = vector()
        consensus_graph = list()
        TF_proteins = getTFIntersection(data_omics, tp_prot, tp_gene, updown, phospho)$TF_Intersection
        proteins = c(getProteinIntersection(data_omics, tp_prot, tp_gene, updown, phospho)$Protein_Intersection,
                     TF_proteins)
        ST_proteins_STRINGid = string_db$map(as.data.frame(proteins), "proteins", takeFirst = TRUE)
        doubID = which(duplicated(ST_proteins_STRINGid$ST_proteins))
        if(length(doubID)>0)
        {STRINGIDs = ST_proteins_STRINGid[-doubID,]  
        }else{
         STRINGIDs = na.omit(ST_proteins_STRINGid)}
        if(dim(STRINGIDs)[1]>0)
        {message("Consensus graph for protein time point ", tp_prot, " and transcript/gene time point ", tp_gene,
                 " was generated.\n")   
            
            ST_net = SteinerTree_cons(STRINGIDs$STRING_id, 
                                      PPI_graph, run_times) 
            ###
            ST_netnames = V(ST_net)$name
            ind_name = which(ST_netnames %in% STRINGIDs$STRING_id)
            ST_netnames[ind_name] = STRINGIDs$proteins[na.omit(match(ST_netnames, STRINGIDs$STRING_id))]
            temp_STR = as.data.frame(ST_netnames[-ind_name])
            colnames(temp_STR) = "STRING_id"
            ST_netnames[-ind_name] = 
                string_db$add_proteins_description(temp_STR)$preferred_name
            V(ST_net)$name = ST_netnames
            V(ST_net)$color[which(V(ST_net)$name %in% TF_proteins)] = "lightblue"
            ###
            
            V(ST_net)$label.cex = 1
            plot(ST_net, main = paste("Steiner tree of phosphoprotein time\n ", 
                                      tp_prot, " and transcript time ", tp_gene, sep = ""))
            legend(x = 0, y = -1.2 , legend = c("consensus proteins",
                                                "steiner node proteins", "consensus TFs"), 
                   fill = c("red", "yellow", "lightblue"), cex = 0.7)
            message("Steiner tree of graph phosphoprotein time point", 
                    tp_prot, " and transcript time point ", tp_gene," was build.\n")
            
            ###
            ST_TFs = TF_proteins
            ST_targets = getGenesIntersection(data_omics, tp_prot, tp_gene, updown, phospho)$Genes_Intersection
            ST_TFTG = list()
            for(k in 1: length(ST_TFs))
            {ST_TFTG[[k]] = as.character(data_omics[[3]][[2]][which(as.character(data_omics[[3]][[2]][,1]) == ST_TFs[k]),2])
            ST_TFTG[[k]] = ST_TFTG[[k]][which(ST_TFTG[[k]] %in% ST_targets)]
            }
            names(ST_TFTG) = ST_TFs
            ###
            ST_net_targets = genfullConsensusGraph(ST_net, ST_TFTG)
            V(ST_net_targets)$label.cex = 0.6
            plot.igraph(ST_net_targets, main = paste("Consensus graph phosphoprotein time\n ", 
                                                     tp_prot, " and transcript time ", tp_gene, 
                                                     sep = ""), vertex.size = 18)
            if("yellow" %in% V(ST_net_targets)$color)
            {legend(x = 0, y = -1.2 , legend = c("consensus proteins", 
                                                 "steiner node proteins", "consensus TFs", "consensus target genes"), 
                    fill = c("red", "yellow", "lightblue", "green"), cex = 0.7)
            }else{
                legend(x = 0, y = -1.2 , legend = c("consensus proteins", 
                                                    "consensus TFs", "consensus target genes"), 
                       fill = c("red", "lightblue", "green"), cex = 0.7)   
            }
            
            ST_net_targets = addFeedbackLoops(ST_net_targets)
            consensus_graph = ST_net_targets
        }else{
            consensus_graph[[tps]] = NA
            message("Not enough intersecting molecules to generate
                        a consensus graph for time point ", same_tps[tps], "\n")   
        }
    }   
    return(consensus_graph)
}




#' Add Consensus Graph information.
#' 
#' Adds phosphoprotein information based on phosphoprotein data table and 
#' redraws Consensus Graph edges
#'    
#' @param data_omics OmicsData object.
#' @param consensusGraph result from static analysis: consensus graph generated 
#' by staticConsensusNet function.
#' @param phosphotab dataframe with phosphoprotein information annotated in 
#' columns 'Gene.names', 'Amino acid', 'Position' (of phosphosite).
#' @return graph of igraph class containing complemented consensus graph
#' information
#' @keywords manip
#' @export
#' @examples
#' \dontrun{
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
#'
#' data_omics_plus = identifyPR(data_omics_plus) 
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' data_omics = identifyPWTFTGs(data_omics)
#' statConsNet = staticConsensusNet(data_omics)
#' }
infoConsensusGraph <- function(data_omics, consensusGraph, phosphotab)
{
    if(length(consensusGraph[[1]]) == 1)
    {consensusGraphs = list(consensusGraph)
    temp = 1
    }else{
    consensusGraphs = consensusGraph
    temp = 0     
    }
    for(k in 1: length(consensusGraphs))
    {   p_names_graph = V(consensusGraphs[[k]])$name[which(V(consensusGraphs[[k]])$color == "red")]
        p_newnames_graph = paste(phosphotab$Gene.names[match(p_names_graph, phosphotab$Gene.names)], "_",
                             phosphotab$Amino.acid[match(p_names_graph, phosphotab$Gene.names)],
                             phosphotab$Position[match(p_names_graph, phosphotab$Gene.names)], 
                             sep = "")
    
    V(consensusGraphs[[k]])$name[which(V(consensusGraphs[[k]])$color == "red")] = p_newnames_graph
    }
    for(j in 1: length(consensusGraphs))
    {   
        ind_green = which(V(consensusGraphs[[j]])$color == "green")
        E(consensusGraphs[[j]])$lty = 1
        E(consensusGraphs[[j]])[from(ind_green)]$lty = 3
    }
    if(temp==1)
    { consensusGraphs = consensusGraphs[[1]]
    temp = 0}
    return(consensusGraphs) 
}



#' Generate STRING PPI graph.
#' 
#' Generates connected graph with undirected edges from STRING PPI-database.
#' 
#' @param string_db STRING_db object generated by getConsensusSTRINGIDs 
#' function.
#' @return igraph object connected graph from STRING PPI database.
#' 
#' @keywords manip
getSTRING_graph <- function(string_db){
    string_homo_new = string_db$load()
    conn_intgraph = decompose.graph(string_homo_new)[[1]]
    graph_STRING = conn_intgraph
    
    return(graph_STRING)
} 


#' Add feedback loops from target genes to proteins/TFs if present.
#' 
#' @param ST_net_targets full consensus graph.
#' @return igraph object with feedback loops added.
#' 
#' @keywords manip
addFeedbackLoops <- function(ST_net_targets){
    
    prot_nodes = V(ST_net_targets)$name[which(V(ST_net_targets)$color == "red" | 
                                                  V(ST_net_targets)$color == "yellow" |
                                                  V(ST_net_targets)$color == "lightblue")]
    gene_nodes = 
        V(ST_net_targets)$name[which(V(ST_net_targets)$color == "green")]
    for(h in 1: length(prot_nodes))
    { feedback_ind = vector()
      if(prot_nodes[h] %in% gene_nodes)
      {
          feedback_ind = which(V(ST_net_targets)$name == prot_nodes[h])
          ST_net_targets[feedback_ind[1],feedback_ind[2] ] = 1
      }
    }    
    return(ST_net_targets)
}

#' Combine SteinerNet with bipartite graph to get full consensus network.
#' 
#' @param ST_net steiner tree graph generated by SteinerTree_cons function.
#' @param ST_TFTG steiner tree graph extended with consensus target genes and 
#' the edges between TFs and target genes.
#' @return igraph object of network comprising steiner tree graph and 
#' TF - target gene interactions.
#' 
#' @keywords manip
genfullConsensusGraph <- function(ST_net, ST_TFTG){
    
    ST_net_targets = ST_net
    match_vec = which(names(ST_TFTG) %in% V(ST_net)$name)
    for(s in match_vec)
    {
        if(s>1)
        {targets_in_graph = V(ST_net_targets)$name[which(V(ST_net_targets)$color == "green")]
         targets_to_add = unique(ST_TFTG[[s]])[which(!unique(ST_TFTG[[s]]) %in% targets_in_graph)]
        }else{
            targets_to_add = unique(ST_TFTG[[s]])    
        }
        temp = ST_net_targets
        ST_net_targets = ST_net_targets + targets_to_add
        len_targets = length(targets_to_add)
        if(len_targets >=1)
        {V(ST_net_targets)$color[(length(V(temp))+1):(length(V(temp)) + len_targets)] = "green"}
        
        ind = which(V(ST_net_targets)$name == names(ST_TFTG[s]) &
                        V(ST_net_targets)$color == "lightblue")
        ind_green = which(V(ST_net_targets)$color == "green" &
                              V(ST_net_targets)$name %in% unique(ST_TFTG[[s]]) )
        if(length(ind)>0 & length(ind_green)>0)
        {ST_net_targets = add.edges(ST_net_targets, rbind(ind_green, ind))}
        
    }
    return(ST_net_targets)
}



#' Get TF-target gene information for the consensus graph.
#' 
#' @param data_omics OmicsData object.
#' @param tps integer specifying current timepoint under consideration.
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection. Type of checking is defined with
#' parameter 'phospho'.
#' @param phospho boolean value; TRUE in case up- and downregulation should be
#' checked based on provided downstream phosphoprotein influence from 
#' identifyPR function; FALSE in case up- and downregulation should be checked
#' for without phosphoprotein database knowledge. Default is TRUE.
#' @return list of transcription factor target gene interactions.
#' 
#' @keywords manip
getbipartitegraphInfo <- function(data_omics, tps, updown = FALSE, phospho = TRUE){
    same_tps = data_omics[[1]][[1]][[1]][[1]][which(data_omics[[1]][[1]][[1]][[1]] 
                                                    %in% data_omics[[1]][[1]][[1]][[2]])]
    ST_TFs = 
        getTFIntersection(data_omics, same_tps[tps], same_tps[tps], updown, phospho)$TF_Intersection
    ST_targets = 
        getGenesIntersection(data_omics, same_tps[tps], same_tps[tps], updown, phospho)$Genes_Intersection
    ST_TFTG = list()
    for(k in 1: length(ST_TFs))
    {ST_TFTG[[k]] = as.character(data_omics[[3]][[2]][which(as.character
                                                            (data_omics[[3]][[2]][,1]) == ST_TFs[k]),2])
     ST_TFTG[[k]] = ST_TFTG[[k]][which(ST_TFTG[[k]] %in% ST_targets)]
    }
    names(ST_TFTG) = ST_TFs
    return(ST_TFTG)  
}

#' Map alias names to STRING IDs of consensus graph.
#' 
#' @param data_omics OmicsData object.
#' @param ST_net steiner tree graph generated by SteinerTree_cons function.
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection. Type of checking is defined with
#' parameter 'phospho'.
#' @param phospho boolean value; TRUE in case up- and downregulation should be
#' checked based on provided downstream phosphoprotein influence from 
#' identifyPR function; FALSE in case up- and downregulation should be checked
#' for without phosphoprotein database knowledge. Default is TRUE.
#' @param consSTRINGIDs first element of list generated by getConsensusSTRINGIDs 
#' function; a data.frame including the proteins to be considered as terminal 
#' nodes in Steiner tree with colnames ST_proteins and the corresponding STRING 
#' IDs in column 'STRING_id'.
#' @param tps integer specifying current timepoint under consideration.
#' @param string_db second element of list generated by getConsensusSTRINGIDs 
#' function; species table (for human) of STRING database.
#' @return igraph object with alias name annotation.
#' 
#' @keywords manip
getAliasfromSTRINGIDs <- function(data_omics, ST_net, updown = FALSE, phospho = TRUE,
                                  consSTRINGIDs, tps, string_db){
    same_tps = data_omics[[1]][[1]][[1]][[1]][which(data_omics[[1]][[1]][[1]][[1]] 
                                                    %in% data_omics[[1]][[1]][[1]][[2]])]
    ST_netnames = V(ST_net)$name
    ind_name = which(ST_netnames %in% consSTRINGIDs$STRING_id)
    ST_netnames[ind_name] = consSTRINGIDs$ST_proteins[na.omit(match(ST_netnames, 
                                                                    consSTRINGIDs$STRING_id))]
    temp_STR = as.data.frame(ST_netnames[-ind_name])
    colnames(temp_STR) = "STRING_id"
    ST_netnames[-ind_name] = 
        string_db$add_proteins_description(temp_STR)$preferred_name
    V(ST_net)$name = ST_netnames
    TF_proteins = 
        getTFIntersection(data_omics, same_tps[tps], same_tps[tps], updown, phospho)$TF_Intersection
    V(ST_net)$color[which(V(ST_net)$name %in% TF_proteins)] = "lightblue"
    
    return(ST_net)
}

#' Get consensus graph in STRING IDs.
#' 
#' @param data_omics OmicsData object.
#' @param tps integer specifying current timepoint under consideration.
#' @param string_db STRING_db object.
#' @param updown boolean value; TRUE in case up- and downregulation should be
#' checked individually for intersection. Type of checking is defined with
#' parameter 'phospho'.
#' @param phospho boolean value; TRUE in case up- and downregulation should be
#' checked based on provided downstream phosphoprotein influence from 
#' identifyPR function; FALSE in case up- and downregulation should be checked
#' for without phosphoprotein database knowledge. Default is TRUE.
#' @return igraph object consensus graph with STRING IDs (only including 
#' proteins and transcription factors).
#' 
#' @keywords manip
getConsensusSTRINGIDs <- function(data_omics, tps, string_db, updown = FALSE, phospho = TRUE){
    same_tps = data_omics[[1]][[1]][[1]][[1]][which(data_omics[[1]][[1]][[1]][[1]] 
                                                    %in% data_omics[[1]][[1]][[1]][[2]])]
    ST_proteins = 
        as.character(getProteinIntersection(data_omics, same_tps[tps], same_tps[tps], updown, phospho)$Protein_Intersection)
    ST_proteins = unique(c(ST_proteins,
                           getTFIntersection(data_omics, same_tps[tps], same_tps[tps], updown, phospho)$TF_Intersection))
    ST_proteins_STRINGid = string_db$map(as.data.frame(ST_proteins),
                                         "ST_proteins", takeFirst = TRUE)
    doubID = which(duplicated(ST_proteins_STRINGid$ST_proteins))
    if(length(doubID)>0)
    {consSTRINGIDs = ST_proteins_STRINGid[-doubID,]  
    }else{
        consSTRINGIDs = na.omit(ST_proteins_STRINGid)}
    return(consSTRINGIDs)
}

#' Steiner tree algorithm.
#' 
#' Use this function to get the Steiner tree based on the STRING protein-protein 
#' interaction database.
#'  
#' @param terminal_nodes character vector of final nodes used for generation
#' of Steiner tree.
#' @param PPI_graph igraph object; graph should be connected and have
#' undirected edges.
#' @param run_times integer specifying number of times to run SP Steiner tree
#' algorithm to find minimal graph.
#' @return igraph object including Steiner tree.
#' 
#' @keywords manip
SteinerTree_cons <- function(terminal_nodes, PPI_graph, run_times) {
    
    color = NULL
    terminal_nodes = na.omit(terminal_nodes)
    V(PPI_graph)$color = "yellow"
    V(PPI_graph)[terminal_nodes]$color = "red"
    terminals = V(PPI_graph)[color == "red"]
    steinertmin = vector()
    steinertrees = list()
    for(runs in 1: run_times)
    {edges = c()
     prob = sample(1:length(terminals), 1, replace = FALSE)
     subtree = terminals$name[[prob]]
     nsubtree = setdiff(terminals$name, subtree)
     tparam = 1
     while(tparam <= length(terminals))
     {   
         paths = get.all.shortest.paths(PPI_graph,subtree, nsubtree)
         if(length(paths$res)>1)
         {
             paths_length = sapply(paths$res, length)
             sp = paths$res[which(paths_length == min(paths_length))][[1]]
             subtree = igraph::union(subtree, V(PPI_graph)$name[sp])
             nsubtree = setdiff(nsubtree, V(PPI_graph)$name[sp])
         }else{
             subtree = subtree
             nsubtree = nsubtree
         }
         tparam = tparam+1
     }
     steinert = minimum.spanning.tree(induced.subgraph(PPI_graph, subtree))
     for(i in length(which(V(steinert)$color == "yellow"))+1)
     { degr = degree(steinert, v = V(steinert), mode = c("all"))
       todel = names(which(degr == 1))
       todel = todel[which(!todel %in% terminals$name)]
       if(length(todel) > 0)
       {steinert = delete.vertices(steinert, todel)}
     } 
     steinertrees[[runs]] = steinert
     steinertmin[runs]  = length(V(steinert)$name)
    }
    return(steinertrees[[which(steinertmin == min(steinertmin))[1]]])  
} 


#' Dynamic analysis.
#' 
#' Generates continous data for dynamic analysis of protein, TF and gene data
#' via smoothing splines. 50 time points are generated this way.  
#' The following nodes are considered: 
#' Nodes which are part of the static consensus graphs from corresponding time 
#' points of the two different measurement types. In case a node is not
#' significantly changed at a certain point in time its FC is assumed to remain
#' constant at this time point.
#' Calculation of the consensus-based dynamic net parameters are based on the ebdbNet 
#' R package [1]. The number of time points generated via smoothing splines (50) 
#' is based on their results for median AUCs of ROC curves.
#' The number of forward time units a node is assumed to influence other nodes
#' can be specified via the laghankel parameter. The cutoff determining the
#' percent of total variance explained by the singular values generated by
#' singular value decomposition (SVD) of the block-Hankel matrix H in order to
#' specify the hidden state dimension K (for further details see [1]).
#'    
#' @references 1. A. Rau, F. Jaffrezic, J.-L. Foulley, R. W. Doerge (2010). 
#' An empirical Bayesian method for estimating biological networks from temporal
#' microarray data. Statistical Applications in Genetics and Molecular Biology, 
#' vol. 9, iss. 1, article 9.     
#'    
#' @param data_omics OmicsData object.
#' @param consensusGraphs result from static analysis: consensus graph generated 
#' by staticConsensusNet function.
#' @param laghankel integer specifying the maximum relevant time lag to be used 
#' in constructing the block-Hankel matrix.
#' @param cutoffhankel cutoff to determine desired percent of total variance 
#' explained; default = 0.9 as in [1].
#' @param conv.1 value of convergence criterion 1; default value is 0.15 
#' (for further details see [1]).
#' @param conv.2 value of convergence criterion 2; default value is 0.05 
#' (for further details see [1]).
#' @param conv.3 value of convergence criterion 3; default value is 0.05 
#' (for further details see [1]).
#' @param verbose boolean value, verbose output TRUE or FALSE
#' @param max.iter maximum overall iterations; default value is 100 
#' (for further details see [1]).
#' @param max.subiter maximum iterations for hyperparameter updates; 
#' default value is 200 (for further details see [1]).
#' @return list of 2 elements: 
#' 1) output parameters of dynamic network inference with ebdbNet package
#' 2) splines data generated.
#' @keywords manip
#' @export
#' @examples
#' \dontrun{
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
#' 
#' data_omics_plus = identifyPR(data_omics_plus) 
#' setwd(system.file("extdata/Genelists", package = "pwOmics"))
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, 
#' noTFs_inPW = 1, order_neighbors = 10)
#' data_omics = identifyPWTFTGs(data_omics)
#' statConsNet = staticConsensusNet(data_omics)
#' dynConsNet = consDynamicNet(data_omics, statConsNet)
#' }
consDynamicNet <- function(data_omics, consensusGraphs, laghankel = 3, 
                           cutoffhankel = 0.9, conv.1 = 0.15, 
                           conv.2 = 0.05, 
                           conv.3 = 0.05, verbose = TRUE, max.iter = 100, 
                           max.subiter = 200){
    
    requireNamespace("ebdbNet", quietly = TRUE)
    requireNamespace("longitudinal", quietly = TRUE)
    as.longitudinal = dataFormat = hankel = ebdbn = NULL
    
    if(class(data_omics) != "OmicsData")
    {stop("Parameter 'data_omics' is not an OmicsData object.")}
    
    if(TRUE %in% !(names(consensusGraphs) %in% 
                       as.character(data_omics[[1]][[1]][[1]][[1]]) &
                       names(consensusGraphs) %in% 
                       as.character(data_omics[[1]][[1]][[1]][[2]])) ) 
    {stop("ConsensusGraphs do not match OmicsData object.")}
    
    prot_nodes = vector()
    gene_nodes = vector()
    for(j in 1: length(consensusGraphs))
    {   prot_nodes = c(prot_nodes, 
                       V(consensusGraphs[[j]])$name[which(V(consensusGraphs[[j]])$color == "red" | 
                                                              V(consensusGraphs[[j]])$color == "yellow" |
                                                              V(consensusGraphs[[j]])$color == "lightblue")])
        gene_nodes = c(gene_nodes, 
                       V(consensusGraphs[[j]])$name[which(V(consensusGraphs[[j]])$color == "green")])
    }
    spl_prot = getFCsplines(data_omics, unique(prot_nodes), "proteins")
    spl_genes = getFCsplines(data_omics, unique(gene_nodes), "genes")
    
    mat_splines_prot = predictFCvals(data_omics, nopredpoints = 50, spl_prot,
                                     "Protein splines")
    mat_splines_genes = predictFCvals(data_omics, nopredpoints = 50, spl_genes,
                                      "Gene splines")
    timevals = mat_splines_prot[[2]]
    rownames(mat_splines_prot[[1]]) = paste(rownames(mat_splines_prot[[1]]),
                                            "_p", sep = "")
    rownames(mat_splines_genes[[1]]) = paste(rownames(mat_splines_genes[[1]]),
                                             "_g", sep = "")
    
    mat_splines = rbind(mat_splines_prot[[1]], mat_splines_genes[[1]])
    
    ts_data = ts(t(mat_splines), start = 0, end = max(timevals),
                 frequency = 50/max(timevals))
    ts_data = as.longitudinal(ts_data)
    ts_data = dataFormat(ts_data)
    K_hankel <- hankel(ts_data, laghankel, cutoffhankel)$dim
    dynConsensusNet = ebdbn(ts_data, K_hankel, input = "feedback", conv.1,
                            conv.2, conv.3, verbose, max.iter, max.subiter)
    return(list(dynConsensusNet, mat_splines))
}

#' Prediction of continous data points via smoothing splines.
#'
#' @param data_omics OmicsData object.
#' @param nopredpoints integer number; how many timpoints should be predicted?
#' @param splineslist list of protein or gene nodes for which splines were 
#' generated: output of getFCsplines function.
#' @param title character vector specifying name of title.
#' @return list of splines matrix values and calculated times. 
#' 
#' @keywords manip
predictFCvals <- function(data_omics, nopredpoints, splineslist, title) {
    
    uni_timepoints = BiocGenerics::union(data_omics[[1]][[1]][[1]][[1]], 
                           data_omics[[1]][[1]][[1]][[2]])
    timevals = seq(0, max(uni_timepoints), length.out = nopredpoints)
    
    splineslist =  splineslist[which(!is.na(splineslist))]
    splines_mat = matrix(ncol = nopredpoints, nrow = length(splineslist))
    rownames(splines_mat) = names(splineslist)
    colnames(splines_mat) = as.character(round(timevals, digits = 2))
    for(k in 1: length(which(!is.na(splineslist))))
    { splines_mat[k,1:50] = predict(splineslist[k][[1]], timevals)$y}
    
    return(list(splines_mat, timevals))
}

#' Get fold change splines.
#' 
#' Calculate the splines used for the dynamic analysis.
#'  
#' @param data_omics OmicsData object.
#' @param nodes character vector of nodes the fold change splines should be 
#' calculated for.
#' @param nodetype character indicating to calculate splines for "proteins" or 
#' "genes".
#' @return splines values used in dynamic analysis.
#' 
#' @keywords manip
getFCsplines <- function(data_omics, nodes, nodetype) {
    
    if(length(nodes)<1)
    {message("Warning: Node number equals 0, no splines can be calculated.\n")
     return(NA)}
    
    if(nodetype == "proteins")
    {timepoints = data_omics[[1]][[1]][[1]][[1]]
     used_list = 1
    }else if(nodetype == "genes"){
        timepoints = data_omics[[1]][[1]][[1]][[2]]
        used_list = 2
    }
    uni_timepoints = BiocGenerics::union(data_omics[[1]][[1]][[1]][[1]], 
                           data_omics[[1]][[1]][[1]][[2]])
    fc = list()
    for(s in 1: length(nodes))
    {
        fc[[s]] = vector()
        fc[[s]][1] = 0
        for(tp in uni_timepoints)
        {
            if(tp %in% timepoints)
            { ind_current_tp = gsub("tp", "", 
                                    names(data_omics[[1]][[2]][[used_list]]))
              if(nodes[s] %in% 
                     data_omics[[1]][[2]][[used_list]][[which(tp == ind_current_tp)]][,1])
              {fc[[s]][which(tp == uni_timepoints)+1] = 
                   data_omics[[1]][[2]][[used_list]][[which(tp == ind_current_tp)]][which(data_omics[[1]][[2]][[used_list]][[which(tp == ind_current_tp)]][,1] == nodes[s]),2]
              }else{
                  fc[[s]][which(tp == uni_timepoints)+1] = fc[[s]][which(tp == uni_timepoints)]
              }
            }else{
                fc[[s]][which(tp == uni_timepoints)+1] = fc[[s]][which(tp == uni_timepoints)]
            }
        }
    }
    names(fc) = nodes
    spl = list()
    for(k in 1: length(fc))
    {
        if(!all(fc[[k]] == 0))
        {spl[[k]] = smooth.spline(c(0,uni_timepoints), fc[[k]], df =4)
        }else{
            spl[[k]] = NA
        }
    }
    names(spl) = nodes
    return(spl)
}


