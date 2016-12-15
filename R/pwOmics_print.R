#' Print an OmicsData object.
#'
#' @param x an OmicsData object to print.
#' @param ... further arguments to be passed to print.
#' @return prints OmicsData object.
#'
#' @keywords manip
#' @export
print.OmicsData = function(x, ...) {
    
    cat("OmicsData class:\n")
    cat("******************************************************************\n")
    cat("Time points: \n")
    print(head(x[[1]][[1]][[1]]))
    cat("******************************************************************\n")
    cat("all Protein IDs: \n")
    print(head(x[[1]][[1]][[2]]))
    cat("\n")
    cat("all Gene IDs: \n")
    print(head(x[[1]][[1]][[3]]))
    cat("\n")
    cat("******************************************************************\n")
    #OmicsDataSet
    cat("Time point data proteins: \n")
    for(i in 1: length(x[[1]][[1]][[1]][[1]]))
    {print(head(x[[1]][[2]][[1]][[i]]))
     cat("\n")}
    cat("\n")
    cat("Time point data genes: \n")
    for(j in 1: length(x[[1]][[1]][[1]][[2]]))
    {print(head(x[[1]][[2]][[2]][[j]]))
     cat("\n")}
    cat("\n")
    cat("******************************************************************\n")
    if(length(x[[1]][[3]][[1]]) != 0)
    {cat("Time point downstream analysis results of protein data: \n")
     cat("All proteins: \n")
     print(head(as.matrix(names(x[[1]][[3]][[1]][[1]][1:3]))))
     
     cat("Time points: \n")
     for(k in 1: length(x[[1]][[1]][[1]][[1]]))
     { cat("Protein data pathways",  names(x[[1]][[3]][[1]])[k+1], "\n" )
       cat("---------------------\n")
       print(head(x[[1]][[3]][[1]][[k+1]][1]))
       cat("...")
       cat("\n\n\n")}
       if(names(head(x[[1]][[3]][[1]][[ length(x[[1]][[1]][[1]][[1]])+1]][[length(x[[1]][[3]][[1]][[ length(x[[1]][[1]][[1]][[1]])+1]])-1]])[1]) == "genes_PW")
       {cat("Protein data TFs: \n")
       print(head(x[[1]][[3]][[1]][[ length(x[[1]][[1]][[1]][[1]])+1]][[length(x[[1]][[3]][[1]][[ length(x[[1]][[1]][[1]][[1]])+1]])-1]])[1:3])
       cat("Protein data downstream target genes: \n")
       print(head(x[[1]][[3]][[1]][[ length(x[[1]][[1]][[1]][[1]])+1]][[length(x[[1]][[3]][[1]][[ length(x[[1]][[1]][[1]][[1]])+1]])]]))
       }
    }else{
        cat("No results have been generated yet in downstream analysis.\n 
            Please use functions readPWdata and identifyPWs for further 
            analysis.\n")
    }
    cat("******************************************************************\n")
    if(length(x[[1]][[3]][[2]]) != 0)
    {
        cat("Time point upstream analysis results of gene data: \n")
        cat("All genes: \n")
        print(head(as.matrix(names(x[[1]][[3]][[2]][[1]][1:3]))))

        cat("Time points: \n")
        for(h in 1: length(x[[1]][[1]][[1]][[2]]))
        { cat("Gene data upstream TFs",  names(x[[1]][[3]][[2]])[h+1], "\n" )
          cat("---------------------\n")
          print(x[[1]][[3]][[2]][[h+1]][1:3])}
        if(names((x[[1]][[3]][[2]][[ length(x[[1]][[1]][[1]][[2]])+1]][[length(x[[1]][[3]][[2]][[ length(x[[1]][[1]][[1]][[2]])+1]])]]))[1] != "upstreamTFs")
        {cat("Gene data upstream pathways: \n")
        print(head((x[[1]][[3]][[2]][[ length(x[[1]][[1]][[1]][[2]])+1]][[length(x[[1]][[3]][[2]][[ length(x[[1]][[1]][[1]][[2]])+1]])-1]]))[1:3])
        cat("Gene data upstream regulators: \n")
        print(head((x[[1]][[3]][[2]][[ length(x[[1]][[1]][[1]][[2]])+1]][[length(x[[1]][[3]][[2]][[ length(x[[1]][[1]][[1]][[2]])+1]])]])))
        }
    }else{
        cat("No results have been generated yet in upstream analysis.\n
            Please use functions readTFdata and identifyTFs for further 
            analysis.\n")
    }
    cat("******************************************************************\n")
    cat("Chosen pathway databases: \n")
    print(head(x[[2]][[1]]))
    cat("\n")
    cat("Biopax model: \n")
    print(head(x[[2]][[2]]))
    cat("\n")
    cat("******************************************************************\n")
    cat("Chosen TF target databases: \n")
    print(head(x[[3]][[1]]))
    cat("\n")
    cat("TF target data: \n")
    print(head(x[[3]][[2]]))
    cat("\n")
    cat("******************************************************************\n")
    if(dim(x[[5]])[1] > 0)
    {cat("Phosphoprotein regulation data: \n")
    print(head(x[[5]]))
    cat("\n")
    cat("******************************************************************\n")
    }
    cat("Status: \n")
    print(head(x[[4]]))
}
