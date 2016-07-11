## 
## write NuchaRt output files for muxViz
##
writeHicToMuxviz <- function(edgesList, startG) {
    path <- getwd()
    folder <- paste(path,"/muxviz/",as.character(startG),"/HiC/",sep="")
    
    if(!dir.exists(folder))
        dir.create(folder, showWarnings = F, recursive = T)
    
    gr_vert <- unique(edgesList[,c(Gene1,Gene2)])
    gr_edg <- subset(edgesList, select = c(Gene1, Gene2, Prob))
    gr_edg[,Prob:=round(Prob, digits = 3)]
    
    startG <- as.character(startG)
    outfile_ed <- paste(folder,as.character(startG),"_HiC_edges.txt", sep="")
    write.table(gr_edg, outfile_ed, sep=" ", row.names=FALSE, col.names=FALSE, quote = F)
    
    outfile_lo <- paste(folder,as.character(startG),"_HiC_layout.txt", sep="")
    write.table(gr_vert, outfile_lo, sep=" ", row.names=TRUE,col.names = c("nodeID nodeLabel"), quote = F)
    
    outfile_cf <- paste(folder,as.character(startG),"_HiC_config.txt", sep="")
    config <- paste(outfile_ed,";",as.character(startG),"_interactions",";", outfile_lo, sep="")
    flConn <- file(outfile_cf)
    writeLines(config, con = flConn)
    close(flConn)
    
    ##    layers file, make by hand.
    ##    Es:
    ##    layerID layerLabel
    ##    1 marriage_alliances
    ##    2 business
}

##
## write Hi-C graph with mapped RNApol to Muxviz
##  - edges list
##  - edges per layer (edges extended)
##  - genes list (single layer layout)
##  - nodes features (external)
##  - config file
##
rnapolToMuxviz <- function(edgesList, startG, time="0h") {
    path <- getwd()
    folder <- paste(path,"/muxviz/",as.character(startG),"/RNApol/",sep="")
    
    if(!dir.exists(folder))
        dir.create(folder, showWarnings = F, recursive = T)
    
    nodes <- workWithWIGs(edgesList,time=time)
    
    gr_vert <- unique(edgesList[,c(Gene1,Gene2)])
    gr_edg <- subset(edgesList, select = c(Gene1, Gene2, Prob))
    gr_edg[,Prob:=round(Prob, digits = 3)]
    gr_edg[Prob==0,Prob:=0.1]
    
    hic_rna <- subset(nodes, select = c(Maximum))
    setnames(hic_rna, old = colnames(hic_rna), new=c("size"))
    hic_rna[,size:=ceiling(size)]
    avg <- ceiling(mean(hic_rna[,size]))
    hic_rna[,color:=ifelse(size>avg,"\"#FF0000\"","\"#2E8B57\"")]
    hic_rna[,layerID:=2]
    hic_rna[,nodeID:=1:nrow(hic_rna)]
    setcolorder(hic_rna, c("nodeID","layerID","color","size"))
    hic_rna[size==0, size:=1]
    
    startG <- as.character(startG)
    outfile_ed <- paste(folder,startG,"_HiC_edges_",time,".txt", sep="")
    write.table(gr_edg, outfile_ed, sep=" ", row.names=FALSE, col.names=FALSE, quote = F)
    
    gr_edg[,L1:=2]
    gr_edg[,L2:=2]
    setcolorder(gr_edg, c("Gene1","L1","Gene2","L2","Prob"))
    outfile_extd <- paste(folder,startG,"_HiC_edges_extended_",time,".txt", sep="")
    write.table(gr_edg, outfile_extd, sep=" ", row.names=FALSE, col.names=FALSE, quote = F)
    
    outfile_lo <- paste(folder,startG,"_HiC_layout_",time,".txt", sep="")
    write.table(gr_vert, outfile_lo, sep=" ", row.names=TRUE, col.names = c("nodeID nodeLabel"), quote = F)
    
    outfile_ext <- paste(folder,startG,"_HiC_external_",time,".txt", sep="")
    write.table(hic_rna, outfile_ext, sep=" ", row.names=FALSE, col.names=TRUE, quote = F)
    
    outfile_cf <- paste(folder,startG,"_HiC_config_",time,".txt", sep="")
    config <- paste(outfile_ed,";",startG,"_interactions",";", outfile_lo, sep="")
    flConn <- file(outfile_cf)
    writeLines(config, con = flConn)
    close(flConn)
    
    return(list(edges=outfile_ed, extended=outfile_extd, layout=outfile_lo))
    
    ##    layers file, make by hand?.
    ##    Es:
    ##    layerID layerLabel
    ##    1 marriage_alliances
    ##    2 business
}

##
## read graph file from InteGO2, write outputs for muxviz:
##  - edges per layer (edges extended)
##  - genes list (single layer layout)
##  - genes with mutation (external)
##
go2ToMuxviz <- function(go2_file, mut_file, startG, time="0h") {
    path <- getwd()
    folder <- paste(path,"/muxviz/",as.character(startG),"/GO2/",sep="")
    
    if(!dir.exists(folder))
        dir.create(folder, showWarnings = F, recursive = T)
    
    go2_graph <- fread(go2_file)
    go2_graph <- go2_graph[V3>0.199]
    
    outfile_go2_ed <- paste(folder,as.character(startG),"_GO2_edges_",time,".txt", sep="")
    write.table(go2_graph, outfile_go2_ed, sep=" ", row.names=FALSE, col.names=FALSE, quote = F)
    
    go2_vert <- unique(go2_graph[,c(V1,V2)])
    
    outfile_go2_lo <- paste(folder,as.character(startG),"_GO2_layout_",time,".txt", sep="")
    write.table(go2_vert, outfile_go2_lo, sep=" ", row.names=T, col.names = c("nodeID nodeLabel"), quote = F)
    
    go2_graph[,L1:=1]
    go2_graph[,L2:=1]
    setcolorder(go2_graph, c("V1","L1","V2","L2","V3"))
    
    outfile_go2_extd <- paste(folder,as.character(startG),"_GO2_edges_extended_",time,".txt", sep="")
    write.table(go2_graph, outfile_go2_extd, sep=" ", row.names=FALSE, col.names=FALSE, quote = F)
    
    mutations <- fread(mut_file, skip = 1)
    mut_dt <- data.table(genes=go2_vert)
    setkey(mutations, "V1")
    mut_dt <- mutations[mut_dt]
    mut_dt[,c(1,3):=NULL, with=F]
    setnames(mut_dt, "V2", "size")
    mut_dt[is.na(size),size:=1]
    
    mut_dt[,nodeID:=1:nrow(mut_dt)]
    avg <- ceiling( ifelse(is.na(mean(mut_dt[,size])), yes=1, no=mean(mut_dt[,size])) )
    mut_dt[,color:=ifelse(size>avg,"\"#FF0000\"","\"#2E8B57\"")]
    mut_dt[,layerID:=1]
    setcolorder(mut_dt, c("nodeID","layerID","color","size"))
    
    outfile_go2 <- paste(folder,startG,"_GO2_external_",time,".txt", sep="")
    write.table(mut_dt, outfile_go2, sep=" ", row.names=FALSE, col.names=TRUE, quote = F)
    
    return(list(edges=outfile_go2_ed, extended=outfile_go2_extd, layout=outfile_go2_lo))
}

##
## merge HiC, Go2 and string graphs for a multi-layer graph with Muxviz
##  - one single layout file for all three graphs (all nodes included)
##  - one nodes feature file (external) to rule them all
##
## INPUT DATA:
##   fileGO2 := edges list downloaded from InteGO2, with weights telling ontological similarity
##   fileString := edges list downloaded from STRING, with weights telling likelihood of interaction
##   fileMuts := file containing genes mutations (COSMIC DB)
##   startGene := gene of interest
##   time := time point (0h or 4h)
##   suffix := suffix for directories and file names (character string)
graphsForMuxviz <- function(fileGO2, fileString, fileMuts, startGene, time, suffix="") {
    if( !exists("Genes", envir = .GlobalEnv) || attr(Genes, "source") != "NuChaRt_Genes" ) 
        stop( 'ERROR: the \'Genes\' table is required, and it must be a valid NuChaRt object.' )
    if( !exists("g.Edges", envir = .GlobalEnv) ) 
        stop( 'ERROR: the \'Edges\' table is required, and it must be a valid NuChaRt object.' )
    if( missing(time) )
        stop( 'ERROR: \'time\' parameter is required: a string (e.g.: \"0h\") describing the time point' )
        
    
    rnapol <- rnapolToMuxviz(g.Edges, startGene, time)
    strEd <- processStringGraphFiles(fileString)
    go2 <- go2ToMuxviz(fileGO2, fileMuts, startGene, time)
    
    # read graph files - all tables have 3 columns: "V1", "V2", "V3"
    hicGr <- fread(rnapol[[1]]) 
    strGr <- fread(strEd[[1]])
    go2Gr <- fread(go2[[1]])
    
    path <- getwd()
    if(suffix != "")
        suffix <- paste("_",suffix,sep="")
    
    folder <- paste(path,"/muxviz/",startGene,"/",time,suffix,"/", sep="")
    # directory containing .wig files with RNApol data. Should be a paramenter
    wigDir <- paste("/space/datasets/multilayer/RNApol/",time,"/",sep="")
    
    if(!dir.exists(folder))
        dir.create(folder, showWarnings = F, recursive = T)
    
    # merge all graphs and obtain all involved genes
    mrg <- merge(go2Gr, hicGr, by = c("V1","V2"), all = T)
    setnames(mrg, old=colnames(mrg), new=c("gene1","gene2","similarity","proximity"))
    mrg[is.na(similarity),similarity:=0]
    mrg[is.na(proximity), proximity:=0]

    mrg_all <- merge(mrg, strGr, by.x = c("gene1","gene2"), by.y = c("V1","V2"), all=T)
    mrg_all[is.na(proximity), proximity:=-1]
    mrg_all[is.na(similarity), similarity:=-1]
    mrg_all[is.na(V3), V3:=-1]
    setnames(mrg_all, "V3", "interaction")

    fName <- paste(folder,startGene,"_merged_",time,suffix,".txt", sep="")
    write.table(mrg_all, fName, sep="\t", row.names = F, col.names = T, quote = F)
    
    gn_names <- unique(mrg_all[,c(gene1,gene2)])
    dt <- data.table(genes=gn_names)
    dt[,nodeID:=1:nrow(dt)]
    setcolorder(dt, c("nodeID", "genes"))
    setnames(dt, "genes", "nodeLabel")
    
    fName_lo <- paste(folder,startGene,"_allnodes_layout_",time,suffix,".txt", sep="")
    write.table(dt, fName_lo, sep=" ", row.names = F, col.names = T, quote = F)
    
    # ---------------------------------------------------------------------------------------------
    # find RNApol, write externals - applies to hic nodes
    hicGr_vert <- unique(hicGr[,c(V1,V2)])
    hicDT <- data.table(verts=hicGr_vert)
    setkey(dt, "nodeLabel")
    setkey(Genes, "Symbol")
    
    nids <- dt[hicDT]
    nids[,nodeLabel:=NULL]
    
    nodes <- Genes[hicDT]
    nodes[, c(5,6) := NULL, with=F]
    nodes[,Span2H := ceiling(abs(Start-Stop)/200)]
    nodes[, Avg := 0]
    nodes[,Maximum := 0]
    
    margin=10000
    wigVec <- vector(length = 1000000)
    for(i in 1:nrow(nodes)) {
        wig <- fread( findWigFile(wigDir,nodes[i]$Chr), skip = 1 )
        wigVec <- wig[,V1]
        st <- round( (nodes[i]$Start - margin)/200 )
        range <- round( ((2*margin)/200) + (st+nodes[i]$Span2H) )
        nodes[i]$Avg  <- ifelse(is.na(mean(wigVec[st:range])), yes=0.0, no=mean(wigVec[st:range]))
        nodes[i]$Maximum <- ifelse(is.na(max(wigVec[st:range])), yes=0.0, no=max(wigVec[st:range]))
    }
    rm(wigVec)
    
    hic_rna <- subset(nodes, select = c(Maximum))
    setnames(hic_rna, old = colnames(hic_rna), new=c("size"))
    hic_rna[,size:=ceiling(size)]
    avg <- ceiling(mean(hic_rna[,size]))
    hic_rna[,color:=ifelse(size>avg,"\"#FF00FF\"","\"#00FF00\"")]
    hic_rna[,layerID:=2]
    hic_rna[,nodeID:=nids]
    setcolorder(hic_rna, c("nodeID","layerID","color","size"))
    hic_rna[size==0, size:=1]
    # ---------------------------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------------------------
    # find expression, write externals - applies to string nodes
    
    # genes expression file containing data for times oh, 1h, 2h, 4h.
    # obtained from MCF-7 files. - should be a parameter
    fr_mcf7 <- fread("/space/datasets/multilayer/expression/my_MCF-7.tsv")
    setkey(fr_mcf7, "Symbol")
    
    strGr_vert  <- unique(strGr[,c(V1,V2)])
    strDT <- data.table(genes=strGr_vert)
    
    nids <- dt[strDT]
    nids[,nodeLabel:=NULL]
    
    h0 <- fr_mcf7[strDT]
    final <- merge(strDT, h0, by.x = c("genes"), by.y = c("Symbol"), all.x = T)
    
    expcol <- ifelse(time=="0h", 3, 6)
    
    exp_val <- subset(final, select = c(1,expcol), with=F)
    setnames(exp_val, old = colnames(exp_val), new = c("V1","V2"))
    exp_val[is.na(V2), V2:=100]
    exp_val[,V2:=ceiling((V2/10)/2)]
    avg <- ceiling(mean(exp_val[,V2]))
    exp_val[,color:=ifelse(V2>avg,"\"#FF0000\"","\"#00FFFF\"")]
    exp_val[,nodeID:=nids]
    exp_val[,layerID:=3]
    setcolorder(exp_val, c("nodeID","layerID","color","V2","V1"))
    setnames(exp_val, old=colnames(exp_val), new=c("nodeID","layerID","color","size","V1"))
    exp_val[,V1:=NULL]
    # ---------------------------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------------------------
    # find mutations, write externals - applies to GO2 nodes
    go2_vert <- unique(go2Gr[,c(V1,V2)])
    
    mutations <- fread(fileMuts, skip = 1)
    mutDT <- data.table(genes=go2_vert)
    
    nids <- dt[mutDT]
    nids[,nodeLabel:=NULL]
    
    setkey(mutations, "V1")
    mut_dt <- mutations[mutDT]
    mut_dt[,c(1,3):=NULL, with=F]
    setnames(mut_dt, "V2", "size")
    mut_dt[is.na(size),size:=1]
    
    mut_dt[,nodeID:=nids]
    avg <- ceiling( ifelse(is.na(mean(mut_dt[,size])), yes=1, no=mean(mut_dt[,size])) )
    mut_dt[,color:=ifelse(size>avg,"\"#0000FF\"","\"#FFFF00\"")]
    mut_dt[,layerID:=1]
    setcolorder(mut_dt, c("nodeID","layerID","color","size"))
    # ---------------------------------------------------------------------------------------------
    
    dtlist = list(mut_dt, hic_rna, exp_val)
    finalDT <- rbindlist(dtlist, use.names = T, fill = T, idcol = NULL)
    setkey(finalDT, "nodeID")
    finaDT <- unique(finalDT)
    
    fName <- paste(folder,startGene,"_allnodes_external_",time,suffix,".txt", sep="")
    write.table(finalDT, fName, sep=" ", row.names = F, col.names = T, quote = F)
    
    assign("allNodes", mrg_all, envir = .GlobalEnv)
    assign("hicRNA", hic_rna, envir = .GlobalEnv)
    assign("nodesExpr", exp_val, envir = .GlobalEnv)
    assign("mutations", mut_dt, envir = .GlobalEnv)
    
    # write extended edges list
    dt2 <- fread(rnapol[[2]], header = F)
    dt3 <- fread(strEd[[2]], header = F)
    dt1 <- fread(go2[[2]], header = F)

    # inter-links L1-L2
    vdt12 <- data.table(v1=unique(dt1[,c(V1,V3)]))
    vdt2 <- data.table(v1=unique(dt2[,c(V1,V3)]))
    setkey(vdt12, "v1")
    vdt12[vdt2]
    vdt12[,l1:=1]
    vdt12[,v2:=v1]
    vdt12[,l2:=2]
    vdt12[,w:=0.999]
    
    # inter-links L1-L3
    vdt13 <- data.table(v1=unique(dt1[,c(V1,V3)]))
    vdt3 <- data.table(v1=unique(dt3[,c(V1,V3)]))
    setkey(vdt13, "v1")
    vdt13[vdt3]
    vdt13[,l1:=1]
    vdt13[,v2:=v1]
    vdt13[,l2:=3]
    vdt13[,w:=0.999]
    
    # inter-links L2-L3
    vdt23 <- data.table(v1=unique(dt2[,c(V1,V3)]))
    vdt3 <- data.table(v1=unique(dt3[,c(V1,V3)]))
    setkey(vdt23, "v1")
    vdt23[vdt3]
    vdt23[,l1:=2]
    vdt23[,v2:=v1]
    vdt23[,l2:=3]
    vdt23[,w:=0.999]
    
    fl <- list(go2=dt1, rnapol=dt2, string=dt3, l12=vdt12, l13=vdt13, l23=vdt23)
    eel <- rbindlist(fl, use.names = F, fill = F, idcol = NULL)
    eel[V5==0, V5:=0.1]
    
    fName_ed <- paste(folder,startGene,"_alledges_extended_",time,suffix,".txt", sep="")
    write.table(eel, fName_ed, sep=" ", row.names = F, col.names = F, quote = F)
    
    # write config file for independent layers
    # FORMAT:
    #   path/to/edges;layer_label;path/to/nodes/layout
    #   path/to/edges;layer_label;path/to/nodes/layout
    #   ....
    outfile_cf <- paste(folder,as.character(startGene),"_config_independent_",time,suffix,".txt", sep="")
    flConn <- file(outfile_cf)
    cat(go2[[1]],";","GO2+mutations",";",fName_lo, "\n",
        rnapol[[1]],";","HiC+RNApol",";",fName_lo,"\n", 
        strEd[[1]],";","STRING+expression",";",fName_lo, sep = "", file = flConn)
    close(flConn)
    
    # write layout and config file for multilayer
    fdir <- paste(path,"/muxviz/",startGene,"/", sep="")
    outfile_ly <- paste(fdir,as.character(startGene),"_layers.txt",sep="")
    flConn <- file(outfile_ly)
    cat("layerID layerLabel", "1 GO2+mutations", "2 HiC+RNApol", "3 STRING+expression", sep = "\n", file=flConn)
    close(flConn)
    
    outfile_cf <- paste(folder,as.character(startGene),"_config_multilayer_",time,suffix,".txt", sep="")
    config <- paste(fName_ed,";",outfile_ly,";",fName_lo, sep="")
    flConn <- file(outfile_cf)
    writeLines(config, con = flConn)
    close(flConn)
}

##
## stringFile is the graph file for a gene, downloaded from string.
## write those graphs in muxviz format, with expression on nodes
##
processStringGraphFiles <- function(stringFile) {
    strExtd <- parseStringGraphs(stringFile)
    fr_mcf7 <- fread("/space/datasets/multilayer/expression/my_MCF-7.tsv")
    
    base <- basename(stringFile)
    path <- dirname(stringFile)
    gene <- strsplit(base, split = '_')[[1]][1]
    curPath <- getwd()
    folder <- paste(curPath,"/muxviz/",gene,"/STRING/",sep="")
    
    file_layout <- paste(folder,gene,"_layout.txt", sep="")
    gLayout <- fread(file_layout)
    lbls <- paste(paste("^",gLayout$nodeLabel,"$",sep=""), collapse = "|")
    
    h0 <- fr_mcf7[Symbol%like%lbls, ]
    wh0 <- subset(h0, select=c(2,3), with=F)
    
    final <- merge(gLayout, wh0, by.x = c("nodeLabel"), by.y = c("Symbol"), all.x = T)
    setnames(final, old = colnames(final), new = c("V1","nodeID","V2"))
    
    final[is.na(V2), V2:=100]
    final[,V2:=ceiling((V2/10)/2)]
    avg <- ceiling(mean(final[,V2]))
    final[,color:=ifelse(V2>avg,"\"#FF0000\"","\"#2E8B57\"")]
    final[,V1:=NULL]
    final[,layerID:=2]
    setcolorder(final, c("nodeID","layerID","color","V2"))
    setnames(final, old=colnames(final), new=c("nodeID","layerID","color","size"))
    setkey(final, "nodeID")
    
    outfile_ext <- paste(folder,gene,"_nodes_external_0h.txt", sep="")
    write.table(final, outfile_ext, sep=" ", row.names=FALSE, col.names=TRUE, quote = F)
    
    wh0 <- subset(h0, select=c(2,6), with=F)
    
    final <- merge(gLayout, wh0, by.x = c("nodeLabel"), by.y = c("Symbol"), all.x = T)
    setnames(final, old = colnames(final), new = c("V1","nodeID","V2"))
    
    final[is.na(V2), V2:=100]
    final[,V2:=ceiling((V2/10)/2)]
    avg <- ceiling(mean(final[,V2]))
    final[,color:=ifelse(V2>avg,"\"#FF0000\"","\"#2E8B57\"")]
    final[,V1:=NULL]
    final[,layerID:=3]
    setcolorder(final, c("nodeID","layerID","color","V2"))
    setnames(final, old=colnames(final), new=c("nodeID","layerID","color","size"))
    setkey(final, "nodeID")
    
    outfile_ext <- paste(path,"/graphs/",gene,"_nodes_external_4h.txt", sep="")
    write.table(final, outfile_ext, sep=" ", row.names=FALSE, col.names=TRUE, quote = F)
    
    return(strExtd)
}

##
## INTERNAL
## parse plain-text graphs from string and write 
## them back into muxviz format
##
parseStringGraphs <- function(stringFile) {
    gene <- basename(stringFile)
    path <- dirname(stringFile)
    gene <- strsplit(gene, split = '_')[[1]][1]
    curPath <- getwd()
    folder <- paste(curPath,"/muxviz/",gene,"/STRING/",sep="")
    if(!dir.exists(folder))
        dir.create(folder, showWarnings = F, recursive = T)
    
    txtGr <- fread(stringFile, header = T)
    setnames(txtGr, "#node1", "node1")
    vert_Gr  <- unique(txtGr[,c(node1,node2)])
    gr_edg <- subset(txtGr, select = c(node1,node2,combined_score))
    gr_edg[,combined_score:=round(combined_score, digits = 3)]
    
    outfile_ed <- paste(folder,gene,"_edges.txt", sep="")
    write.table(gr_edg, outfile_ed, sep=" ", row.names=FALSE, col.names=FALSE, quote = F)
    
    gr_edg[,L1:=3]
    gr_edg[,L2:=3]
    setcolorder(gr_edg, c("node1","L1","node2","L2","combined_score"))
    outfile_extd <- paste(folder,gene,"_edges_extended.txt", sep="")
    write.table(gr_edg, outfile_extd, sep=" ", row.names=FALSE, col.names=FALSE, quote = F)
    
    outfile_lo <- paste(folder,gene,"_layout.txt", sep="")
    write.table(vert_Gr, outfile_lo, sep=" ", row.names=TRUE, 
                col.names = c("nodeID nodeLabel"), quote = F)
    
    outfile_cf <- paste(folder,gene,"_config.txt", sep="")
    config <- paste(outfile_ed,";",gene,"_string",";", outfile_lo, sep="")
    flConn <- file(outfile_cf)
    writeLines(config, con = flConn)
    close(flConn)
    
    return(list(edges=outfile_ed, extended=outfile_extd,laout=outfile_lo))
    
    ##    layers file, make by hand
    ##    Es:
    ##    layerID layerLabel
    ##    1 marriage_alliances
    ##    2 business
}

##
## parse and work with wig files
## to be fixed and tested
##
workWithWIGs <- function(edgesList,time="0h",margin=10000,wigPath="/space/datasets/multilayer/RNApol/") {
    if( missing(edgesList) ) 
        stop( 'ERROR: the \'edgesList\' is not a valid NuChaRt object.' )
    
    wigDir <- paste(wigPath,time,sep="")
    verts <- paste(paste("^",unique(edgesList[,c(Gene1,Gene2)]),"$",sep=""), collapse = "|")
    
    nodes <- Genes[Symbol%like%verts,]
    nodes[, c(5,6) := NULL, with=F]
    nodes[,Span2H := ceiling(abs(Start-Stop)/200)]
    nodes[, Avg := 0]
    nodes[,Maximum := 0]
    
    wigVec <- vector(length = 1000000)
    for(i in 1:nrow(nodes)) {
        wig <- fread( findWigFile(wigDir,nodes[i]$Chr), skip = 1 )
        wigVec <- wig[,V1]
        st <- round( (nodes[i]$Start - margin)/200 )
        range <- round( ((2*margin)/200) + (st+nodes[i]$Span2H) )
        nodes[i]$Avg  <- ifelse(is.na(mean(wigVec[st:range])), yes=0.0, no=mean(wigVec[st:range]))
        nodes[i]$Maximum <- ifelse(is.na(max(wigVec[st:range])), yes=0.0, no=max(wigVec[st:range]))
    }
    rm(wigVec)
    
    return(nodes)
}

##
## INTERNAL
## find .wig file by 'chr'
##
findWigFile <- function(dir_path, chr) {
    pt <- paste("*",chr,".wig$",sep="")
    wFile <- list.files(path=dir_path, pattern=pt, full.names=T, recursive=FALSE)
    return (wFile)
}