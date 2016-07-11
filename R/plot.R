################
#
# GRAPH
# Functions needed to create and plot graphs, either flat (using iGraph) or interactive
# (using networkD3).
# Once the graph is constructed, some statistical analysis over the graph can
# be performed: a topologycal analysis of the graph, centrality measures,
# exponential random graph modelling, graph correlation.
#
################

PlotFlatGraph <- function(edgesList, expr.file, plotit=1) {
    if( missing(edgesList) )
        stop( 'ERROR: list of edges must be provided.' )

    if( attr(edgesList, "source") != "NuChaRt_Edges" )
        stop( 'ERROR: \'edgesList\' is not a valid NuChaRt object.\n',
              'It must be generated using the \'BuildGraph\' function.' )

    vert <- unique(edgesList[,c(Gene1,Gene2)])
    vert <- VertexInfo(vert,expr.file)
    vert <- as.data.frame(vert)
    edg <- subset(edgesList, select = c(Gene1, Gene2, Seq1, Seq2, Weight, Prob))
    edg[Prob==0, Prob:=1]
    edg <- as.data.frame(edg)

    g <- graph.empty()
    g <- as.undirected(g)
    g <- g + vertices( vert$Symbol,label=vert$Symbol,start=vert$Start,end=vert$Stop,
                       exp=vert$Expr )
    name=paste(edg[,1],"-",edg[,2])
    g <- g + edges(t(edg[,1:2]), seq1=edg[,3], seq2=edg[,4], weight=edg[,5], prob=edg[,6])

    V(g)$deg <- degree(g, V(g))

    for (i in 1:length(V(g))) {
        # intergenic - kept for compatibility (?)
        if(length(grep("chr",V(g)[i]$name))==1) {
            V(g)[i]$weight=2
            V(g)[i]$color="red"
        }
        if ((V(g)[i]$exp==0)) {
                V(g)[i]$weight=4
                V(g)[i]$color="white"
        } else if (V(g)$exp[i]<0) {
            V(g)[i]$weight=abs(V(g)[i]$exp)
            V(g)[i]$color="green"
        } else {
            V(g)[i]$weight=2*V(g)[i]$exp
            V(g)[i]$color="orange"
        }
    }
    V(g)$shape="circle"

    if(plotit) {
        #Put genes in Skyblue
        for (i in 1:length(V(g))){if (V(g)[i]$color=="white") { V(g)[i]$color="SkyBlue2"}}
        #tiff ("/home/imerelli/f6_4.tiff", width=10, height=10, units="in", res=300)
        plot( g,layout=layout.fruchterman.reingold(g),vertex.color=V(g)$color,
          edge.width=2.6*E(g)$prob, vertex.label=V(g)$label, vertex.label.cex=0.9, vertex.size=1.5*V(g)$weight )
    }

    return (g)
}

PlotD3Graph <- function(edgesList, expr.file="src/PnuChart/extdata/expression.txt") {
    if( missing(edgesList) || attr(edgesList, "source") != "NuChaRt_Edges" )
        stop( 'ERROR: the \'edgesList\' is not a valid NuChaRt object.' )

    vert <- unique(edgesList[,c(Gene1,Gene2)])
    vert <- VertexInfo(vert,expr.file)

    edg <- subset(edgesList, select = c(Gene1, Gene2, Seq1, Seq2, Weight, Prob))
    setkey(edg, Gene2)

    # NodeIds must starts from zero.
    #vert$Id <- 0:(nrow(vert)-1)
    vert[,Id:=0:(nrow(vert)-1)]
    edg[Gene1 %in% vert$Symbol, Gene1 := with(vert, as.character(Id[na.omit(match(Gene1, Symbol))]))]
    edg[Gene2 %in% vert$Symbol, Gene2 := with(vert, as.character(Id[na.omit(match(Gene2, Symbol))]))]

    # source and target columns must be numeric
    edg[,Gene1:=as.numeric(Gene1)]
    edg[,Gene2:=as.numeric(Gene2)]

    edg$Prob <- abs(round(10*edg$Prob, 0))
    edg[Prob==0, Prob:=1]

    #vert$Expr <- abs(round(10*vert$Expr, 0))
    vert[,Expr:=abs(round(10*Expr, 0))]
    vert[Expr==0, Expr:=1]

    vert <- as.data.frame(vert)
    edg <- as.data.frame(edg)

    #nD3 <-
    forceNetwork(Links = edg, Nodes = vert, Source = "Gene2",
                 Target = "Gene1", Value = "Prob", NodeID = "Symbol", Group = "Chr",
                 opacity = 0.7, colourScale = "d3.scale.category20b()", zoom = T,
                 legend = T, Nodesize = "Expr")
    #saveNetwork(network = nD3, file = 'Net1.html')
}

PlotSimpleD3Graph <- function(edgesList) {
    if( missing(edgesList) )
        stop( 'ERROR: the \'edgesList\' is not a valid NuChaRt object.' )

    edg <- subset(edgesList, select = c(Gene1, Gene2, Seq1, Seq2, Weight, Prob))
    setkey(edg, Gene2)

    simpleNetwork(edg, fontSize = 5, zoom = T, opacity = 0.7)
}

PlotFlatGraphERBS <- function(edgesList, plotit=1) {
    if( missing(edgesList) )
        stop( 'ERROR: list of edges must be provided.' )

    vert <- unique(edgesList[,c(Gene1,Gene2)])
    #vert <- VertexInfo(vert)
    vert <- Genes[Symbol%in%vert]

    vert <- as.data.frame(vert)
    edg <- subset(edgesList, select = c(Gene1, Gene2, Seq1, Seq2, Prob))
    edg[Prob==0, Prob:=1]
    edg <- as.data.frame(edg)

    g <- graph.empty()
    g <- as.undirected(g)
    g <- g + vertices( vert$Symbol,label=vert$Symbol,start=vert$Start,end=vert$Stop)
    name <- paste(edg[,1],"-",edg[,2])
    g <- g + edges(t(edg[,1:2]), seq1=edg[,3], seq2=edg[,4], prob=edg[,5])

    V(g)$deg <- degree(g, V(g))

    for (i in 1:length(V(g))) {
        # intergenic - kept for compatibility (?)
        if(length(grep("ERBS",V(g)[i]$name))==1) {
            V(g)[i]$weight=14
            V(g)[i]$color="red"
            V(g)[i]$shape="square"
        } else {
            if(V(g)[i]$name%in%startGenes) {
                V(g)[i]$weight=14
                V(g)[i]$color="Yellow"
                V(g)[i]$shape="circle"
            } else {
                V(g)[i]$weight=14
                V(g)[i]$color="SkyBlue2"
                V(g)[i]$shape="circle"
            }
        }
    }

    if(plotit) {
        plot_file <- paste(startGenes,"_",basename(file.SAM),".pdf",sep="")
        #tiff (plot_file, width=10, height=10, units="in", res=300)
        #pdf(file = plot_file,  paper = "a4")
        plot( g,layout=layout.fruchterman.reingold(g),vertex.color=V(g)$color,
              edge.width=5*E(g)$prob, vertex.label=V(g)$label, vertex.label.cex=0.9,
              vertex.size=1.5*V(g)$weight )
        #dev.off()
    }

    return (g)
}


###################################
## VertexInfo
##
## internal function used to create a list of 'augmented' vertices for the graph:
## for each vertex, some information concerning the gene are retrieved, so that
## thay can be mapped on the graph
##
##################################
VertexInfo <- function(vertices,fileEXPR) {
    if(length(vertices) == 0)
        stop( 'ERROR: vertices must be an array of genes symbols' )

    if(!exists("Genes")) {
        #LoadGenesFile("src/PnuChart/extdata/human_genes.txt")
        stop( 'ERROR: no Genes table found' )
    }

    v <- Genes[Symbol%in%vertices,]
    #v$Expr <- 0
    v[,Expr:=0]
    ex <- NULL

    if( missing(fileEXPR)) {
        data.table::setorder(v, "Symbol")
        v[,Expr:=1]
    } else if(grepl(".diff",fileEXPR)) {
        ex <- LoadDiffFile(fileEXPR)
        vEx <- subset(ex, subset=Chr %in% v$Chr)
        tmpLfc <- NULL
        for(j in 1:nrow(v)) {
            redEx <- vEx[Chr==v[j]$Chr & Start>=v[j]$Start & Stop<=v[j]$Stop,]
            if(nrow(redEx) > 0) {
                redEx[is.na(LogFC), LogFC:=0]
                v[j]$Expr <- mean(redEx$LogFC)
            } else {
                v[j]$Expr <- 0
            }
        }
        data.table::setkeyv(v, c("Chr","Start"))
    } else if(grepl(".txt",fileEXPR)) {
        ex <- LoadExpressionFiles(fileEXPR)
        vEx <- unique(ex[Symbol%in%vertices,],by="Symbol")
        m <- merge(v, vEx, by="Symbol", all=T)
        m[is.na(LogFC), LogFC:=0]
        m[is.na(pValue), pValue:=0.06]
        m <- m[`pValue` <= 0.05]
        data.table::setorder(v, "Symbol")
        #v$Expr <- m$LogFC
        v[,Expr:=m$LogFC]
    }

    rm(ex)
    rm(vEx)
    rm(m)
    gc()

    assign("file.EXPR", fileEXPR, envir = .GlobalEnv)
    assign("g.Vertices", v, envir = .GlobalEnv)
    return(v)
}

###################################
## Dsplay
##
## internal function from printing genes characteristics 1D
##
##
##
##################################

Display<-function(g,v){
    for (i in 1:length(V(g))){
        cat(sprintf("%s : %s\n",V(g)$name[i],v[i]))
    }
}

###################################
## Display2
##
## internal function from printing genes characteristics 2D
##
##
##
##################################

Display2<-function(g,v){
    for (i in 1:length(V(g))){
        for (j in 1:length(V(g))){
            cat(sprintf("%s - %s : %s\n",V(g)$name[i],V(g)$name[j],v[i,j]))
        }
    }
}
