args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}

argsL <- as.list(do.call("cbind", parseArgs(args))[c(F,T)])
names(argsL) <- as.list(do.call("cbind", parseArgs(args))[c(T,F)])
args <- argsL;rm(argsL)

if(! is.null(args$help)) {
  cat("
      Mandatory arguments:
      --df_windows                 - Input Rdata file that contains all the windows to be splitted
      Optional arguments:
      --help \n\n")
  q(save="no")
}

if(is.null(args$df_windows)) {stop("Option --df_windows should be provided")} else{df_windows=args$df_windows}
load(df_windows)

# create one empty file per window to run the scripts in parallel across windows
for (id in 1:nrow(df_window_merge)){
  w = paste(df_window_merge[id,"chr"], ":", df_window_merge[id,"start"], "-", df_window_merge[id,"end"], sep="")
  system(paste("touch ", w, sep=""))
}