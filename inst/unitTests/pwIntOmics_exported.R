#*******************************************************************************
#
# unit tests for exported function of pwOmics.R
#
#*******************************************************************************


test.readOmics <- function(){
    ##check for right input format
    ex_data = data(omics_example)
    checkTrue(length(readOmics(tp_prots = c(2,5,10,20), 
                                    tp_genes = c(10,20,60,120), 
                                    ex_data,
                                    PWdatabase = c("biocarta", "kegg", "nci"), 
                                    TFtargetdatabase = c("chea"))) == 4)
    ##check for right reading inof data
    checkTrue("OmicsD" %in% names(readOmics(tp_prots = c(2,5,10,20), 
                       tp_genes = c(10,20,60,120), 
                       ex_data,
                       PWdatabase = c("biocarta", "kegg", "nci"), 
                       TFtargetdatabase = c("chea")))
}

test.readTFdata <- function(){
    ##check for right input format
    checkException(readTFdata(4))
}

test.readPWdata <- function(){
    ##check for right input format
    checkException(readPWdata(4))
}

test.enrichPWs <- function(){
    ##check for right input format
    checkException(enrichPWs(4))
}

test.enrichTFs <- function(){
    ##check for right input format
    checkException(enrichTFs(4))
}

test.identifyPWTFTGs <- function(){
    ##check for right input format
    checkException(identifyPWTFTGs(4))
}

test.identifyRsofTFs <- function(){
    ##check for right input format
    checkException(identifyRsofTFs(4))
}

test.getOmicsTimepoints <- function(){
    ##check for right output format
    checkTrue(length(getOmicsTimepoints(data_omics)) == 2)
}
test.getOmicsallProteinIDs <- function(){
    ##check for right output format
    checkTrue(dim(getOmicsallProteinIDs(data_omics))[[2]] == 1)
}
test.getOmicsallGeneIDs  <- function(){
    ##check for right output format
    checkTrue(dim(getOmicsallGeneIDs(data_omics))[[2]] == 1)
}
test.getOmicsDataset   <- function(){
    ##check for right output format
    checkTrue(length(getOmicsDataset(data_omics)) == 2)
}
test.getDS_PWs   <- function(){
    ##check for right number of time points in output
    checkTrue(length(getDS_PWs(data_omics)) == 
                  length(data_omics[[1]][[1]][[1]][[1]]))
}
test.getDS_TFs   <- function(){
    ##check for right number of time points in output
    checkTrue(length(getDS_TFs(data_omics)) == 
                  length(data_omics[[1]][[1]][[1]][[1]]))
}
test.getDS_TGs   <- function(){
    ##check for right number of time points in output
    checkTrue(length(getDS_TGs(data_omics)) == 
                  length(data_omics[[1]][[1]][[1]][[1]]))
}
test.getUS_TFs   <- function(){
    ##check for right number of time points in output
    checkTrue(length(getUS_TFs(data_omics)) == 
                  length(data_omics[[1]][[1]][[1]][[2]]))
}
test.getUS_PWs   <- function(){
    ##check for right number of time points in output
    checkTrue(length(getUS_PWs(data_omics)) == 
                  length(data_omics[[1]][[1]][[1]][[2]]))
}
test.getUS_regulators   <- function(){
    ##check for right number of time points in output
    checkTrue(length(getUS_regulators(data_omics)) == 
                  length(data_omics[[1]][[1]][[1]][[2]]))
}
test.getBiopaxModel   <- function(){
    ##check for right output class
    checkTrue("biopax" == class(getBiopaxModel(data_omics)))
}
test.getProteinIntersection <- function(){
    ##check for right output format
    checkTrue(length(getProteinIntersection(data_omics)) == 3)
}
test.getTFIntersection  <- function(){
    ##check for right output format
    checkTrue(length(getTFIntersection(data_omics)) == 3)
}
test.getGenesIntersection  <- function(){
    ##check for right output format
    checkTrue(length(getGenesIntersection(data_omics)) == 3)
}
test.gettpIntersection  <- function(){
    ##check for right output format
    checkTrue(length(gettpIntersection(data_omics)$Intersection) == 3)
}
test.staticConsensusNet  <- function(){
    ##check for right input format
    checkException(staticConsensusNet(3))
}
test.dynamicConsensusNet  <- function(){
    ##check for right input format
    checkException(dynamicConsensusNet(3))
}
test.clusterTimeProfiles <- function(){
    ##check for right input format
    checkException(clusterTimeProfiles(3))
}


