##################################################
# Global variables
##################################################

###
# write global variables here
###


welcomeFunction <- function(){

    cat("\n")
    cat(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n")
    cat("::: Welcome to NuchaRt\n")
    #cat("==========================\n")
    cat(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n")
    cat("\n")
    cat(":: NuchaRt: a high-level parallel computing tool for augmented Hi-C data analysis\n")
    cat("\n")
    cat(":: This software is released under GNU GPL v3:\n")
    cat(":: http://www.gnu.org/copyleft/gpl.html\n")
    cat("\n")
    cat(paste(":: Version:",nuchaRtVersion,"\n"))
    cat(paste(":: Last update:",nuchaRtUpdate,"\n"))
    cat("\n")
    cat(paste(":: Your system is a",Sys.info()["sysname"],"\n"))
    cat(paste("::",Sys.info()["version"],"\n"))
    cat(paste("::",version["version.string"][[1]],"\n"))
    cat("\n")
}

buildPath <- function(folder,objname){
    if( Sys.info()["sysname"]=="Windows" ){
        return( paste(getwd(),folder,objname,sep="\\") )
    }else{
        return( paste(getwd(),folder,objname,sep="/") )
    }
}

buildTmpPath <- function(objname){
    return(  buildPath("tmp",objname) )
}

concatenatePath <- function(folder,objname){
    if( Sys.info()["sysname"]=="Windows" ){
        return( paste(folder,objname,sep="\\") )
    }else{
        return( paste(folder,objname,sep="/") )
    }
}

getExecutablePath <- function(exec_name){
    path <- ""
    if( Sys.info()["sysname"]=="Windows" ){
        path <- buildPath("bin",paste0(exec_name,"_windows"))
    }else if( Sys.info()["sysname"]=="Linux" ){
        path <- buildPath("bin",paste0(exec_name,"_linux"))
    }else{
        path <- buildPath("bin",paste0(exec_name,"_macosx"))
    }

    return(path)
}

mdebug <- function(message){
    cat(paste("DEBUG:",message,"\n"))
}

