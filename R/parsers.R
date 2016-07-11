################
##
## PARSERS
## Functions useful for parsing datasets files. Once parsed, data is
## stored in the Global Environment as data.tables, and can be conveniently 
## visualised from R.
## 
## NOTE: some functions do not return any object. Instead, loaded data is stored 
## in the Global Environment.
##
################

LoadFragmentFile<-function( frag.file ) {
    if(missing(frag.file) || !file.exists(frag.file))
        stop('Fragments file not found')
    
    DF <- data.table::fread(frag.file, skip=2, header=F, sep='\t', stringsAsFactors=FALSE, select=c(1:3))
    data.table::setnames(DF, old=colnames(DF), new=c("Chr", "Start", "Stop"))
    DF[,Id:=c(0:(nrow(DF)-1))]
    
    attr(DF, "source") <- "NuChaRt_Frags"
    assign("Frags", DF, envir = .GlobalEnv)
}

LoadGenesFile<-function( genes.file, frags.file, num_wrks) {
    if(missing(genes.file) || !file.exists(genes.file))
        stop('Genes File not found')
    
    if(missing(frags.file) || !file.exists(frags.file))
        stop('Fragments File not found')
    
    DF <- fread(genes.file, header=TRUE, sep='\t', stringsAsFactors=FALSE, select=c(1,3,4,5,6))
    data.table::setnames(DF, old=colnames(DF), new=c("Chr", "Start", "Stop", "Symbol", "EntrezID"))
    data.table::setkey(DF, EntrezID)
    DF <- na.omit(DF)
    DF <- subset(unique(DF))
    
    data.table::setorder(DF, Chr, EntrezID, Start, Stop)
    DF[,Id:=c(0:(nrow(DF)-1))]
    
    stsp <- extendGenesCoords(genes.file, frags.file, num_wrks)
    DF[,Start:=stsp$starts[1:nrow(DF)]]
    DF[,Stop:=stsp$stops[1:nrow(DF)]]
    
    attr(DF, "source") <- "NuChaRt_Genes"
    assign("Genes", DF, envir = .GlobalEnv)
}

LoadSamFile<-function( sam.file ) {
    if(missing(sam.file) || !file.exists(sam.file))
        stop('SAM File not found')
    
    if( grepl('ord_', sam.file) ) {
        DF <- fread(sam.file, header=FALSE, sep='\t', stringsAsFactors=FALSE, select=c(2,3,4,5,6))
        data.table::setnames(DF, old=colnames(DF), new=c("Chr1", "Start1", "Chr2", "Start2", "Seq"))
    } else {
        DF <- fread(sam.file, header=FALSE, skip=26, sep='\t', stringsAsFactors=FALSE, select=c(3,4,7,8,10))
        data.table::setnames(DF, old=colnames(DF), new=c("Chr1", "Start1", "Chr2", "Start2", "Seq"))
        set(DF, i=which(DF$Chr2=='='), j=3L, DF$Chr1)
        data.table::setorder(DF, Chr1, Start1)
    }
    DF[,Id:=c(0:(nrow(DF)-1))]
    
    attr(DF, "source") <- "NuChaRt_Sam"
    assign("Sam", DF, envir = .GlobalEnv)
}

LoadDiffFile <- function(diff.file) {
    if(missing(diff.file) )
        stop('DIFF File not specified')
    
    if(!file.exists(diff.file))
        stop('DIFF File does not exist.')
    
    DF <- fread(diff.file, header=TRUE, sep='\t', stringsAsFactors=FALSE, select=c(4,10,12,13))
    data.table::setnames(DF, old=colnames(DF), new=c("Locus", "LogFC", "pValue", "qValue"))
    
    locus <- DF$Locus
    elems <- as.data.table(matrix(unlist(strsplit(locus, split = "[:-]")), ncol = 3, byrow = TRUE))
    
    DF[,Locus:=NULL]
    DF[,Chr:=elems$V1]
    DF[,Start:=as.numeric(elems$V2)]
    DF[,Stop:= as.numeric(elems$V3)]
    DF[grepl("chr",Chr,ignore.case = TRUE)==FALSE,Chr:=paste("Chr",Chr, sep="")]
    data.table::setkeyv(DF, c("Chr","Start"))
    #DF <- DF[`p-value` <= 0.05]
    #setorder(DF, -LogFC)
    
    data.table::setcolorder(DF, c("Chr", "Start", "Stop", "LogFC", "pValue", "qValue"))
    attr(DF, "source") <- "NuChaRt_Diff"
    
    write.table(x = DF, file = "src/PnuChart/extdata/mouse/gene_exp.txt", 
                append=FALSE, sep = '\t', col.names = FALSE, quote = FALSE, row.names = FALSE)
    
    return (DF)
}

LoadExpressionFiles <- function(expr.file) {
    if(missing(expr.file) )
        stop('Expressions File not specified')
    
    if(!file.exists(expr.file))
        stop('Expressions File does not exist.')
    
    #DF <- fread(expr.file, header=TRUE, stringsAsFactors=FALSE, na.strings = " ", sep = " ", sep2 = " ")
    DF <- read.csv(expr.file, header=TRUE, sep=" ", quote = "", fill = TRUE)
    DF <- data.table::as.data.table(DF)
    data.table::setnames(DF, old=colnames(DF), new=c("EntrezID", "Symbol", "LogFC", "pValue"))
    data.table::setkey(DF, Symbol)
    DF <- na.omit(DF)
    
    attr(DF, "source") <- "NuChaRt_Expr"
    
    return (DF)
}

LoadBindingSiteFiles <- function(bed.file) {
    if(missing(bed.file) )
        stop('BED File not specified')
    
    if(!file.exists(bed.file))
        stop('BED File does not exist.')
    
    DF <- fread(bed.file, header=TRUE, sep='\t', stringsAsFactors=FALSE, select=c(1,2,3,4,6))
    setnames(DF, old=colnames(DF), new=c("Chr", "Start", "Stop", "Sequence", "Strand"))
    
    attr(DF, "source") <- "NuChaRt_Bed"
    
    return (DF)
}

LoadMACSbed <- function(MACS.bed.file) {
    if(missing(MACS.bed.file) )
        stop('BED File not specified')
    
    if(!file.exists(MACS.bed.file))
        stop('BED File does not exist.')
    
    DF <- fread(MACS.bed.file, header=FALSE, sep='\t', stringsAsFactors=FALSE, select=c(1,2,3,5))
    setnames(DF, old=colnames(DF), new=c("Chr", "Start", "Stop", "Ref.Value"))
    
    attr(DF, "source") <- "NuChaRt_MACSbed"
    
    return (DF)
}

LoadFeaturesFile <- function(feat.file) {
    if(missing(feat.file) )
        stop('Features File not specified')
    
    if(!file.exists(feat.file))
        stop('Features File does not exist.')
    
    DF <- fread(feat.file, header = F, sep='\t', stringsAsFactors = F)
    setnames(DF, old=colnames(DF), new=c("Chr", "Start", "Stop", "Len", "GCc", "Map"))
    
    attr(DF, "source") <- "NuChaRt_Feat"
    
    return (DF)
}

parseChiapetBed <- function(chiapetFile) {
    base <- basename(chiapetFile)
    base <- strsplit(base,"[.]")[[1]][1]
    path <- dirname(chiapetFile)
    
    bedT <- fread(chiapetFile, header = F, sep="\t")
    bed4 <- subset(bedT, select=c("V4"))
    
    bed4[, c("chr1","stsp1","chr2","stsp2") := tstrsplit(V4, "[:-]")]
    bed4[,c("start1","empty1","end1") := tstrsplit(stsp1, "[..]")]
    bed4[,c("empty1","stsp1"):=NULL]
    
    bed4[,c("start2","empty2","end2","score") := tstrsplit(stsp2, "[..]|,")]
    bed4[,empty2:=NULL]
    bed4[,c(1,4,6,7):=NULL,with=F]
    setcolorder(bed4, c("chr1","start1", "chr2", "end2", "score"))
    setnames(bed4, "end2", "start2")
    bed4 <- bed4[score>2,]
    setkeyv(bed4, c("chr1","start1"))
    
    rm(bedT)
    tabname <- paste(path,"/","tsv_",base,".tsv",sep="")
    write.table(bed4, file=tabname, sep = "\t", quote = F, row.names = F, col.names = T)
}
