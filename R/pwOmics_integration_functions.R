#' Find downstream signaling axis.
#'
#' This function determines the regulated downstream structures of a selected
#' phosphoprotein at a specified time point.
#'  
#' @param data_omics OmicsData object.
#' @param phosphoprot character specifying the name of the phosphoprotein that is
#' selected as the starting point for downstream analysis.
#' @param tpDS integer specifying the time point considered for downstream analysis
#' of phosphoprotein data .
#' @return list of downstream pathways identified for this time point of 
#' phosphoprotein measurement with sublists containing information about the 
#' transcription factor, the regulation, the target genes of these transcription 
#' factors and the matching transcripts.
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
#' setwd(system.file("extdata/Genelists", package = "pwOmics.newupdown"))
#' \dontrun{
#' data_omics_plus = identifyPR(data_omics_plus)
#' data_omics = identifyPWs(data_omics_plus)
#' data_omics = identifyTFs(data_omics)
#' data_omics = identifyPWTFTGs(data_omics)
#' data_omics = identifyRsofTFs(data_omics, noTFs_inPW = 1, order_neighbors = 10)
#' SYK_axis = findSignalingAxes(data_omics,phosphoprot = "SYK", tpDS = 2)
#' }
findSignalingAxes <- function(data_omics, phosphoprot, tpDS)
{
    ##identify DS pathways of this phosphoprot
    ind_DS = which(data_omics[[1]][[1]][[1]][[1]] %in% tpDS)
    DS_length = dim(data_omics[[1]][[2]][[1]][[ind_DS]])[1]
    if(TRUE %in% grepl(paste(phosphoprot,"$", sep = ""), names(data_omics[[1]][[3]][[1]][[ind_DS+1]])))
    {phospho_ind = grep(paste(phosphoprot,"$", sep = ""), names(data_omics[[1]][[3]][[1]][[ind_DS+1]]))
    DS_pathways = data_omics[[1]][[3]][[1]][[ind_DS+1]][[phospho_ind]]

    res = as.list(DS_pathways$pathwayIDs)
    names(res) = (DS_pathways$pathwayIDs)
    if(length(res)!= 0)
    {
        ##identify downstream transcription factors and target genes 
        ##for each pathway
        DS_pathway_nodes = list()
        for(j in 1: dim(DS_pathways)[1])
        {
            mygraph = pathway2RegulatoryGraph(data_omics[[2]][[2]],
                                    as.character(DS_pathways$pathwayIDs[j]), 
                                    expandSubpathways = TRUE,
                                    splitComplexMolecules = TRUE, 
                                    useIDasNodenames = FALSE, verbose = TRUE)
            if(length(mygraph)>0 )
            { if(length(graph::nodes(mygraph))>0)
               { 
                ##in case of reactome pathway: check individually to protein ID 
                ##Annotation: "CSNK1A1" instead of "UniProt:P48729 CSNK1A1" 
                if(TRUE %in% grepl("UniProt:", graph::nodes(mygraph)))
                {temp_nodes = vector()
                    for(k in 1: length(graph::nodes(mygraph)))
                    {
                        if(grepl("UniProt:", graph::nodes(mygraph)[k]))
                        {temp_nodes[k] = strsplit(graph::nodes(mygraph)[k], " ")[[1]][2] 
                        }else{
                            temp_nodes[k] = graph::nodes(mygraph)[k] 
                        }
                    }
                    DS_pathway_nodes[[j]] = temp_nodes
                }else{
                    DS_pathway_nodes[[j]] = graph::nodes(mygraph)}
                DS_pathway_nodes[[j]] = na.omit(DS_pathway_nodes[[j]])
                
                res[[j]] = as.list(DS_pathway_nodes[[j]])
                names(res[[j]]) = DS_pathway_nodes[[j]]
                
                
                ##identify TFs
                ind_grep = vector()
                for(s in 1: length(DS_pathway_nodes[[j]]))
                {   
                    
                    if(DS_pathway_nodes[[j]][s] == "Ca++")
                    {DS_pathway_nodes[[j]][s] = "Ca2+"}
                    if(DS_pathway_nodes[[j]][s] == "Ca ++")
                    {DS_pathway_nodes[[j]][s] = "Ca2+"}
                    if(DS_pathway_nodes[[j]][s] == "Fe++")
                    {DS_pathway_nodes[[j]][s] = "Fe2+"}
                    if(DS_pathway_nodes[[j]][s] == "Zn++")
                    {DS_pathway_nodes[[j]][s] = "Zn2+"}
                    if(grepl("\\[",DS_pathway_nodes[[j]][s]) | grepl("\\]",DS_pathway_nodes[[j]][s]))
                    {DS_pathway_nodes[[j]][s] = gsub("\\]", "", gsub("\\[", "", DS_pathway_nodes[[j]][s]))}
                    
                    if(length( grep(DS_pathway_nodes[[j]][s],
                            as.character(data_omics[[3]][[2]][,1]), 
                            ignore.case = TRUE)) <1 )
                    {ind_grep[s] = NA
                    }else{
                        ind_grep[s] = paste(grep(paste("^",DS_pathway_nodes[[j]][s], sep = ""),
                                            as.character(data_omics[[3]][[2]][,1]), 
                                            ignore.case = TRUE), collapse = ",")}
                    if(!is.na(ind_grep[s]) & ind_grep[s] == "")
                    {ind_grep[s] = NA}
                    res[[j]][s][[1]] = list(TF = data.frame(is_TF = !is.na(ind_grep[s]), 
                                    is_upreg = DS_pathways$upreg[j]), 
                                    target_genes = data.frame())
                    
                    target_genes = list()
                    if(!is.na(ind_grep[s]))
                    {
                        ##identify target genes
                        target_genes = unique(data_omics[[3]][[2]][strsplit(ind_grep[s],"," )[[1]],2])
                        target_genes = data.frame(target_genes, 
                                                  upreg = DS_pathways$upreg[j], 
                                                  t(rep(NA, times = 2* length(data_omics[[1]][[1]][[1]][[2]]))))
                        
                        ###target genes that are found diff. expressed on transcript level...
                        for(k in 1: dim(target_genes)[1])
                        {
                            for(d in 1: length(data_omics[[1]][[1]][[1]][[2]]))
                            {target_genes[k,d+2] = as.character(target_genes[k,1]) %in% data_omics[[1]][[2]][[2]][[d]][,1]
                            if( as.character(target_genes[k,1]) %in% data_omics[[1]][[2]][[2]][[d]][,1])
                            {target_genes[k,d+2+length(data_omics[[1]][[1]][[1]][[2]])] = 
                                data_omics[[1]][[2]][[2]][[d]][which(data_omics[[1]][[2]][[2]][[d]][,1] %in% as.character(target_genes[k,1])),2] >0
                            }else{
                                target_genes[k,d+2+length(data_omics[[1]][[1]][[1]][[2]])] = NA   
                            }
                            }
                        }
                        colnames(target_genes)[1:(2+(2*length(data_omics[[1]][[1]][[1]][[2]])))] =  
                            c(as.character(DS_pathway_nodes[[j]][s]), 
                            "upreg" ,paste("tp_", data_omics[[1]][[1]][[1]][[2]], "_transcript", sep = ""), 
                            paste("tp_", data_omics[[1]][[1]][[1]][[2]], "_transcript_upreg", sep = ""))
                    }else{
                        target_genes = NA
                    }
                    res[[j]][[s]][[2]] = data.frame(target_genes)
                    
                    ##if matching, then check for each time point individually, 
                    ##which transcript is later on found in which pathway
                    if(is.data.frame(target_genes) )
                    {matching_transcripts = vector()
                    for(d in 1: length(data_omics[[1]][[1]][[1]][[2]]))
                    {
                        matching_transcripts = c(matching_transcripts,
                                which(grepl("TRUE", res[[j]][[s]][[2]][,d+2])))
                    }
                    res[[j]][[s]][[2]][unique(matching_transcripts),]
                    res[[j]][[s]][[3]] = res[[j]][[s]][[2]][unique(matching_transcripts),]
                    }else{
                        res[[j]][[s]][[3]] = NA    
                    }
                    names(res[[j]][[s]])[3] = "matching_transcripts"
                }
            }
           }
        }
    }
    return(res)
    }else{
        return(res = list())
    }
}



#' Summarize table of matching downstream transcripts.
#'
#' This function returns the summarized table of matching transcripts information
#' based on the result list from the findSignalingAxes function.
#'  
#' @param data_omics OmicsData object.
#' @param axis list output of findSignalingAxes.
#' @return dataframe containing the target genes of the axis matching the
#' transcript data in at least one of the time points in the first column, the
#' direction of their regulation inferred from the phosphoproteome data 
#' in the second column, the transcript regulation in the next columns (per time
#' point) and the summarized matched information in the following columns (per
#' time point). 
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
#' SYK_axis = findSignalingAxes(data_omics,phosphoprot = "SYK", tpDS = 2)
#' SYK_transcripts = get_matching_transcripts(data_omics, SYK_axis)
#' }
get_matching_transcripts <- function(data_omics, axis)
{
    if(length(axis)>0)
    {dat = data.frame(t(rep(NA, times = 2+2*length(data_omics[[1]][[1]][[1]][[2]]))))
    colnames(dat) = c("", "upreg", 
                      paste("tp_", data_omics[[1]][[1]][[1]][[2]], "_transcript", sep=""),
                      paste("tp_", data_omics[[1]][[1]][[1]][[2]], "_transcript_upreg", sep=""))
    for(p in 1: length(axis))
    {   
        temp = vector()
        for(i in 1: length(axis[[p]]))
        {
            if(length(axis[[p]][[i]])==3)
            {temp[i] = axis[[p]][[i]]$TF[[1]]}
        }
        if(TRUE %in% temp)
        {for(s in 1: length(grep("TRUE",temp)))
            {
                colnames(axis[[p]][[which(temp == TRUE)[s]]]$matching_transcripts)[1] = ""
                dat = rbind(dat,axis[[p]][[which(temp == TRUE)[s]]]$matching_transcripts)
            }
        }else{
            colnames(dat) = c("", "upreg", 
                 paste("tp_", data_omics[[1]][[1]][[1]][[2]], "_transcript", sep=""),
                 paste("tp_", data_omics[[1]][[1]][[1]][[2]], "_transcript_upreg", sep=""))
        }
    }
    dat = dat[-1,]
    dat = unique(dat)
    return(dat)
    }else{
        return(axis)
    }
}

#' Generate a folder with downstream information about all phosphoproteins.
#'
#' This function generates a folder structure with an RData object and csv
#' tables for each timepoint sepcified for each phosphoprotein considered in
#' the data_omics object.
#'  
#' @param data_omics OmicsData object.
#' @param timepoints integer vector specifying the timepoints of interest for
#' downstream analysis.
#' @return Folder structure in working directory, containing phosphoprotein
#' downstream information in individual folders.
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
#' setwd("~/Signaling_axes/")
#' generate_DSSignalingBase(data_omics, timepoints = c(0.25, 1, 4, 8, 13, 18, 24))
#' }
generate_DSSignalingBase <- function(data_omics, timepoints = c(0.25, 1, 4, 8, 13, 18, 24))
{
    all_phosphos = as.character(getOmicsallProteinIDs(data_omics)[,1])
    for(k in 1:length(all_phosphos))
    {
        all = all_phosphos[k]                                                      
        dir.create(paste(getwd(),"/", all, "_downstream/",sep = ""),               
                   showWarnings = TRUE, recursive = FALSE, mode = "0777")
        for(tp in timepoints)
        {   axis = findSignalingAxes(data_omics, phosphoprot = all, tpDS = tp)
        axis_ds = get_matching_transcripts(data_omics, axis)
        write.csv(axis_ds, paste(getwd(),"/", all, "_downstream/", all,"_tp", 
                                 tp, "_matching_transcripts.csv", sep = ""))
        save(axis, file = paste(getwd(),"/", all, "_downstream/", all,"_tp", 
                                tp, "_matching_transcripts.RData", sep = ""))
        rm(axis)
        }
    }
}


















