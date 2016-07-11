###################################
## GraphStatistics
##
## Provide some functions to describe and confront graph
##
###################################

GraphStatistics <- function(g){

    #Global statistics
    nv <- length(V(g))
    cat(sprintf("Number of vertices: %s; ",nv))
    ne <- length(E(g))
    cat(sprintf("Number of edges: %s.\n",ne))

    #     dens <- graph.density(g)
    #     cat(sprintf("Graph density: %s\n",dens))

    cl <- clusters(g)$no
    cat(sprintf("Number of Graph clusters: %s; ",cl))
    diam <- diameter(g)
    cat(sprintf("Graph diameter: %s.\n",diam))

    deg.dist <- degree.distribution(g,cumulative = T)
    ddDF <- data.table(Degree=1:length(deg.dist), NumberOfNodes=deg.dist)
    #     cat(sprintf("Degree distribution: %s\n",deg.dist))

    #Centrality Measures
    cat("Centrality measures\n")

    deg <- degree(g)
    deg <- data.table(Node=names(deg), Degree=deg)
    #     cat(sprintf("Graph degree:\n"))
    #     Display(g,deg)

    clsns <- closeness(g)
    clsns <- data.table(Node=names(clsns), Closeness=clsns)
    #     cat(sprintf("Graph closeness:\n"))
    #     Display(g,clsns)

    btwns <- betweenness(g)
    btwns <- data.table(Node=names(btwns), Betweenness=btwns)
    #     cat(sprintf("Graph betweenness:\n"))
    #     Display(g,btwns)

    clCoef <- transitivity(g, vids=NULL, type="local")
    clCoef <- data.table(Node=V(g)$label, Transitivity=clCoef)
    #     cat(sprintf("Clustering Cofficent:\n"))
    #     Display(g,clCoef)

    eigVec <- evcent(g)$vector
    eigVec <- data.table(Node=names(eigVec), Eigenvector=eigVec)
    #     cat(sprintf("Eigenvector centrality:\n"))
    #     Display(g,eigVec)

    #Degree distribution
    cat("Degree distribution\n")
    plot(degree.distribution(g), xlab="node degree")
    lines(degree.distribution(g))

    #Centrality Measures
    cat("Betweeness vs. evcent\n")
    order(betweenness(g))
    order(evcent(g)$vector)
    plot(evcent(g)$vector, betweenness(g))
    text(evcent(g)$vector, betweenness(g), 0:100, cex=0.6, pos=4)

    return ( list(DegrDistr=ddDF,
                  Degrees=deg,
                  Closeness=clsns,
                  Betweenness=btwns,
                  ClustCoeff=clCoef,
                  Eigenvector=eigVec) )
}


###################################
## graph_ergm
##
## A function to compute ergm analysis on the graph
##
## Mandatory:
## v1 = vertices table of graph
## e1 = edges table of graph
## t1 = features table
## form = statistical formula
##
###################################

graphERGM <- function(graph, edgesList, BEDfile="src/PnuChart/extdata/BED/full.bed", form="nodecov('features')") {

    if(class(graph) != "igraph")
        stop('ERROR: an iGraph graph must be provided')

    names <- V(g)$name
    features <- createFeaturesTable(graph, edgesList, BEDfile)

    m1 <- get.adjacency(g, sparse = F)
    #Create networks
    n1 <- network(m1,directed=FALSE)

    f1 <- NULL
    for (i in 1:length(names)[1]) {
        f1 <- c(f1,sum(features==names[i]))
    }

    #Check formula
    var <- strsplit(form, "\'")[[1]][2]

    #Load features
    n1 %v% var <- f1

    #Create formula
    fmla <- as.formula(paste("n1 ~ edges + ", form))

    #Compute model
    model <- ergm(fmla)

    #return report
    return(list(m1=m1,n1=n1,model=model))

}


###################################
## graph_correlation
##
## A function to compare graph
##
## Mandatory:
## v1 = vertices table of graph 1
## v2 = vertices table of graph 2
## e1 = edges table of graph 1
## e2 = edges table of graph 2
##
###################################

graphCorrelation <- function(graph1, graph2) {
   if(class(graph1) != "igraph" || class(graph2) != "igraph")
        stop('ERROR: an iGraph graph must be provided')

    m1 <- get.adjacency(graph1, sparse = F)
    m2 <- get.adjacency(graph2, sparse = F)
    c=cor(c(m1), c(m2))
    return(list(c=c,m1=m1,m2=m2))
}

createFeaturesTable <- function(g, edgesList, bed) {
    if(class(g) != "igraph")
        stop('ERROR: an iGraph graph must be provided')

    if(missing(bed))
        stop('BED File not specified')

    if( missing(edgesList) )
        stop( 'ERROR: list of edges must be provided.' )

    f <- LoadMACSbed(bed)
    nm <- V(g)$name
    tableBED <- f[0]
    tmpSymbols <- NULL

    # First fragment - chr == chr1 & st > st1 & end < end1
    chrF <- f[Chr %in% edgesList$Chr1,]
    chrF <- subset(unique(chrF,by = c("Start","Stop")))
    for(i in 1:nrow(edgesList)) {
        bds <- chrF[Chr==edgesList[i]$Chr1 & Start>edgesList[i]$Start1 & Stop<edgesList[i]$End1,]
        if(nrow(bds) != 0) {
            tableBED <- rbindlist(list(tableBED,bds))
            tmpSymbols <- c(tmpSymbols, rep(edgesList[i]$Gene2,nrow(bds)))
        }
    }

    # Second fragment - chr == chr2 & st > st2 & end < end2
    chrF <- f[Chr %in% edgesList$Chr2,]
    chrF <- subset(unique(chrF,by = c("Start","Stop")))
    for(i in 1:nrow(edgesList)) {
        bds <- chrF[Chr==edgesList[i]$Chr2 & Start>edgesList[i]$Start2 & Stop<edgesList[i]$End2,]
        if(nrow(bds) != 0) {
            tableBED <- rbindlist(list(tableBED,bds))
            tmpSymbols <- c(tmpSymbols, rep(edgesList[i]$Gene1,nrow(bds)))
        }
    }

    #setnames(tableBED, old=colnames(tableBED), new=c("Chr", "Start", "Stop", "Sequence", "Symbol"))
    tableBED$Symbol <- tmpSymbols

    return (tableBED)
}

createBEDTable <- function(g, edgesList, bed) {
    if(class(g) != "igraph")
        stop('ERROR: an iGraph graph must be provided')

    if(missing(bed))
        stop('BED File not specified')

    if( missing(edgesList) )
        stop( 'ERROR: list of edges must be provided.' )

    f <- LoadBindingSiteFiles(bed)
    nm <- V(g)$name
    tableBED <- f[0]
    tmpSymbols <- NULL

    # First fragment - chr == chr1 & st > st1 & end < end1
    chrF <- f[Chr %in% edgesList$Chr1,]
    chrF <- subset(unique(chrF,by = c("Start","Stop")))
    for(i in 1:nrow(edgesList)) {
        bds <- chrF[Chr==edgesList[i]$Chr1 & Start>edgesList[i]$Start1 & Stop<edgesList[i]$End1,]
        if(nrow(bds) != 0) {
            tableBED <- rbindlist(list(tableBED,bds))
            tmpSymbols <- c(tmpSymbols, rep(edgesList[i]$Gene2,nrow(bds)))
        }
    }

    # Second fragment - chr == chr2 & st > st2 & end < end2
    chrF <- f[Chr %in% edgesList$Chr2,]
    chrF <- subset(unique(chrF,by = c("Start","Stop")))
    for(i in 1:nrow(edgesList)) {
        bds <- chrF[Chr==edgesList[i]$Chr2 & Start>edgesList[i]$Start2 & Stop<edgesList[i]$End2,]
        if(nrow(bds) > 0) {
            tableBED <- rbindlist(list(tableBED,bds))
            tmpSymbols <- c(tmpSymbols, rep(edgesList[i]$Gene1,nrow(bds)))
        }
    }

    setnames(tableBED, old=colnames(tableBED), new=c("Chr", "Start", "Stop", "Sequence", "Symbol"))
    tableBED$Symbol <- tmpSymbols

    return (tableBED)
}
