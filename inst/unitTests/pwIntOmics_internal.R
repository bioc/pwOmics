#*******************************************************************************
#
# unit tests for internal functions of pwOmics.R
#
#*******************************************************************************

##########internal functions
test.loadPWs <- function(){
    ##check that pathways are correctly loaded
    checkTrue(length(pwIntOmics:::loadPWs(c("kegg", "biocarta", "nci", 
                                            "reactome"),
        PWdb_path = "/home/awachte/pwIntOmics2/pwIntOmics1/rBiopaxParser_parsed_files/")) == 4)
}

test.readTFtargets <- function(){
    ##check that wrong format results in error
    checkException(pwIntOmics:::readTFtargets(data_omics,  
                    "/home/awachte/pwIntOmics2/Code/pwIntOmics/inst/unitTests"))
}

test.genIntIDs <- function(){
    ##check that output is of biopax class
    loaded_PWs = pwIntOmics:::loadPWs("kegg", 
                 PWdb_path = "/home/awachte/pwIntOmics2/pwIntOmics1/rBiopaxParser_parsed_files/")
    checkTrue("biopax" == 
              class(pwIntOmics:::genIntIDs(data_omics, loaded_PWs, 2, "kegg"))[1])
}

test.createIntIDs <- function(){
    ##check that right number of pathways are loaded
    loaded_PWs = pwIntOmics:::loadPWs("kegg", 
                 PWdb_path = "/home/awachte/pwIntOmics2/pwIntOmics1/rBiopaxParser_parsed_files/")
    checkTrue(length(pwIntOmics:::createIntIDs(data_omics, loaded_PWs) == 4))
}

test.getAlias_Ensemble <- function(){
    ##check that translation of IDs is right
    test_ids = c("ENSG00000143632", "ENSG00000143632", "ENSG00000143632", 
                 "ENSG00000143632", "ENSG00000143632", "ENSG00000143632", 
                 "ENSG00000109846", "ENSG00000120738")
    gene_symbol_ids = c("ACTA1", "ACTA1", "ACTA1", "ACTA1", "ACTA1", "ACTA1", 
                        "CRYAB", "EGR1")
    checkTrue(all(pwIntOmics:::getAlias_Ensemble(test_ids) == gene_symbol_ids))
}

test.genGenelistssub <- function(){
    ##check for right output class
    loaded_PWs = pwIntOmics:::loadPWs("kegg", 
                 PWdb_path = "/home/awachte/pwIntOmics2/pwIntOmics1/rBiopaxParser_parsed_files/")
    intIDs_gen = pwIntOmics:::genintIDs(data_omics, loaded_PWs, 2, "kegg")
    checkTrue("data.table" = class(pwIntOmics:::genGenelistssub(intIDs_gen, 
                                                                2, "kegg")))
}

test.genGenelists <- function(){
    ##check that output is a list of length 4
    loaded_PWs = pwIntOmics:::loadPWs("kegg", 
                 PWdb_path = "/home/awachte/pwIntOmics2/pwIntOmics1/rBiopaxParser_parsed_files/")
    intIDs_gen = pwIntOmics:::genintIDs(data_omics, loaded_PWs, 2, "kegg")
    checkTrue(length(pwIntOmics:::genGenelists(intIDs_gen, "kegg"))==4)
}

test.createBiopaxnew <- function(){
    ##check that output is a biopax model
    loaded_PWs = pwIntOmics:::loadPWs("kegg", 
                 PWdb_path = "/home/awachte/pwIntOmics2/pwIntOmics1/rBiopaxParser_parsed_files/")
    intIDs_gen = pwIntOmics:::genintIDs(data_omics, loaded_PWs, 2, "kegg")
    checkTrue("biopax" == class(pwIntOmics:::createBiopaxnew(data_omics, "kegg")))
}

test.identTFs <- function(){
    ##check that output is a data.table with right number of columns
    checkTrue("data.table" == class(pwIntOmics:::identTFs(data_omics, 1)))
    checkTrue(dim(pwIntOmics:::identTFs(data_omics, 1))[2] == 2)
}

test.loadGenelists <- function(){
    ##check that output is a list
    checkTrue("list" == class(pwIntOmics:::loadGenelists()))
}

test.identPWsofTFs <- function(){
    ##check that output is a list of length 2
    setwd("/home/awachte/pwIntOmics2/pwIntOmics1/Genelists")
    genelists = pwIntOmics:::loadGenelists()
    dt = pwIntOmics:::identTFs(data_omics, 1)
    checkTrue(length(pwIntOmics:::identPWsofTFs(genelists, dt)) == 2)
}

test.preparePWinfo <- function(){
    ##check that output is a list of length 2
    checkTrue(length(pwIntOmics:::preparePWinfo(data_omics, 1)) == 2)
}

test.getSTRING_graph <- function(){
    ##check that STRING graph contains right number of vertices
    string_db = STRINGdb$new(version = "9_05", species = 9606,
                             score_threshold = 0, input_directory = "") 
    checkTrue(length(V(pwIntOmics:::getSTRING_graph(string_db))) == 20128)    
}
