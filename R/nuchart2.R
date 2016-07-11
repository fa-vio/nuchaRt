source("R/version.R")
source("R/global.R")

################################
#
# Function BuildGraph(...).
#
################################
BuildGraph <- function( genesList, numWrks=1, searchLvl=1, grain=0, plotGviz=1, plotIgraph=0, only=0,
                        fileColls="",
                        fileExpression="",
                        fileSAM="data/human/uniques.sam",
                        fileGenes="data/human/human_genes.txt",
                        fileFrags="data/human/DIG/HindIII.txt" ) {

    ################################################################################
    # Build the neighbourhood graph of the given gene(s).
    # Starting from one or more input genes, a neighbourhood graph is built,
    # according to the datasets used. This function sets up all data and calls a
    # C++ function, where a ParallelFor from FastFlow library skeleton is used to
    # build the graph in a fast and memory-efficient way.
    #
    # Args:
    #   genesList:  a list of gene from which the graph is built. It can be a list
    #               of symbols or a list of EntrezIDs (space separated), or a cluster
    #               of genes (like HLA or KRAB).
    #   numWrks:    number of working threads to use during computation.
    #               Default is 1.
    #   searchLvl:  maximum height of the graph starting from the given initial
    #               genes. Default is 1, that is all genes adjacents to the root(s).
    #   grain:      grain for parallel loop. Default is 0, which equally partitions
    #               tasks among workers
    #   fileSAM:    path to the file of nucleotide sequence alignments to be used.
    #               Default is a small sample file named 'uniques.sam'. Note that
    #               SAM files can be several GigaBytes in size.
    #   plotIgraph: Plot using iGraph instead of Grapviz. Default is 0
    #
    # Returns:
    #   this function does not return an object. Instead, it saves in the Global
    #   Environment a data.table object named 'edges' containing the resulting
    #   neighbourhood graph in the form of a list of edges, with information for
    #   each edge, concering chromosome fragment coordinates, genes symbols,
    #   nucleotide sequences, a 'Weight' column and a 'Probs' column (to be filled
    #   by the 'NormaliseEdges' function) and two hashed values needed for internal
    #   operations.
    #
    ################################################################################

    genesVec <- NA
    intrs <- 0

    # parse fragments and genes and save them to global environment.
    # useful for pretty visualisation.
    # SAM file not automatically parsed, due to its size and the concrete risk of
    # compromising R's stability.
    LoadFragmentFile(fileFrags)
    LoadGenesFile(fileGenes, fileFrags, numWrks)

    if(!file.exists(fileSAM)) {
        stop( 'ERROR: the SAM file provided does not exists.\nGraph cannot be built.' )
    } else {
        assign("file.SAM", fileSAM, envir = .GlobalEnv)
    }

    # Check function arguments, before invoking C++ function
    if( missing(genesList) )
        stop( 'ERROR: at least one gene is required to start building the graph.' )

    # check for clusters
    if( !missing(genesList) ) {
        if( grepl('chr',genesList) ) genesVec <- genesList
        else genesVec <- strsplit(genesList,split=" ")[[1]]
    }

    if( !missing(fileColls) )
        intrs <- 1

    if(searchLvl < 1) {
        warning( 'WARNING: Search Level must be greater than ZERO. Setting to 1.', immediate.=T )
        searchLvl <- 1
    }

    if(numWrks < 1) {
        warning( 'WARNING: Number of worker threads must be greater than ZERO. Setting to 1.', immediate.=T )
        numWrks <- 1
    }

    # get user login name
    info <- Sys.info()
    info <- unname(info[7])

    listFromC <- NA
    listFromC <- runNuChart2(genesVec, fileSAM, fileGenes, fileFrags, fileColls, fileExpression,
                             info, numWrks, searchLvl, grain, plotGviz, intrs, only)

    dtEdges <- rbindlist(listFromC, fill = T)
    vert  <- unique(dtEdges[,c(Gene1,Gene2)])
    rm(listFromC)
    gc()
    mallinfo::malloc.trim() # return any unused memory to OS - needed in order to avoid memory overflows

    attr(dtEdges, "source") <- "NuChaRt_Edges"
    assign("g.Edges", dtEdges, envir = .GlobalEnv)
    attr(vert, "source") <- "NuChaRt_Vertices"
    assign("g.Vertices", vert, envir = .GlobalEnv)
    attr(genesVec, "source") <- "NuChaRt_Start"
    assign("startGenes", genesVec, envir = .GlobalEnv)
}


################################
#
# Function NormaliseEdges(...)
#
###############################
NormaliseEdges <- function(edgesList, numWrks=1, grain=0, edgesProb=0,
                           fileSAM="data/human/uniques.sam",
                           featFile="data/human/GF/gf_full.txt") {

    ################################################################################
    # calculate the weight of each edge using the IWLS C++ function, that computes
    # a Poisson regression using genomic features.
    #
    # Args:
    #   edgesList:  a graph expressed as a list of edges, as returned by the
    #               'BuildGraph' function.
    #   numWrks:   number of working threads to use during computation.
    #               Default is 1.
    #   grain:      grain for parallel loop.
    #               Default is 0, which equally partitions tasks among workers
    #   edgesProb:  threshold for edges probability.
    #               Default is 0, which accepts all edges. Allowed probability
    #               values range from 0.1 to 1).
    #   fileSAM:    path to the file of nucleotide sequence alignments to be used.
    #               Default is a small sample file (SAM files may be several GB).
    #
    # Returns:
    #   a data.table containing start coordinates of the two chromosome
    #   fragments belonging to the two linked genes, a hashed value that
    #   identifies the chromosome, and a value for the weight (a score) and the
    #   probability of existence of the edge.
    #
    ################################################################################

    if( missing(edgesList) || attr(edgesList, "source") != "NuChaRt_Edges" )
        stop( 'ERROR: the \'edgesList\' is not a valid NuChaRt object.' )

    if(edgesProb != 0)  {
        if(edgesProb < 0.0 || edgesProb > 1.0) {
            warning( 'WARNING: threshold for edge plotting must fall within [0.1, 1.0]. Setting to 0.0.', immediate.=T )
            edgesProb <- 0.0
        }
    }

    if(numWrks < 1) {
        warning( 'WARNING: Number of worker threads must be greater than ZERO. Setting to 1.', immediate.=T )
        numWrks <- 1
    }

    miniEdges <- data.table::copy(edgesList)
    miniEdges[,c(1L,2L,4L,5L,6L,7L,9L,10L)] <- NULL
    miniEdges <- as.data.frame(miniEdges)

    nnodes <- length(g.Vertices)

    retEdges <- runNormalisation(miniEdges, fileSAM, featFile, numWrks, edgesProb, grain, nnodes)

    weightedEdges <- rbindlist(retEdges, fill = TRUE)
    rm(retEdges)
    rm(miniEdges)

    # TODO(fabio): assign by column number is not really the best solution!
    data.table::set(edgesList, i=NULL, j = 11L, value = weightedEdges$Weight)
    data.table::set(edgesList, i=NULL, j = 12L, value = weightedEdges$Prob)


    #vert <- unique(edgesList[,c(Gene1,Gene2)])
    #assign("g.Vertices", vert, envir = .GlobalEnv)
    rm(weightedEdges)
    gc()
    mallinfo::malloc.trim() # return any unused memory to OS - needed in order to avoid memory overflows

    edgsStat <- subset(edgesList, select = c(Gene1,Gene2,Weight,Prob))
    edgsStat[,CS:=ifelse(edgesList$Chr1==edgesList$Chr2, yes="CIS", no = "TRANS")]

    name <- character()
    # print edges basic statistics:
    # Gene1 Gene2 Weight Prob Cis/Trans
    for(g in startGenes)
        name <- paste(name, g, sep = "_")

    now <- Sys.time()
    tstmp <- as.character(format(now, "%Y%m%dh%H%M%S"))

    file_path <- paste(getwd(), "/csv/edges_stat/", tstmp, name, sep="")
    dir.create(file_path, recursive = T, showWarnings = T)

    write.table(edgsStat, file=paste(file_path,"/edgsStat.csv", sep=""), quote = FALSE, na = "",
                row.names = FALSE, sep = "\t")
    cat("Edges statistics written in ", paste(file_path,"/edgsStat.csv", sep=""), sep="")
}

