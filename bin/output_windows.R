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
      --nwindow_list               - number of chunks containing all the windows
        OR --genome_overlap        - a dataframe with header CHR START END, windows overlapping each line (genomic pos) will be grouped for the skat test
      Optional arguments:
      --help \n\n")
  q(save="no")
}

if(is.null(args$df_windows)) {stop("Option --df_windows should be provided")} else{df_windows=args$df_windows}
if(is.null(args$nwindow_list)) {stop("Option --nwindow_list should be provided")} else{nwindow_list=as.numeric(args$nwindow_list)}
load(df_windows)

# create one empty file per window to run the scripts in parallel across windows
list_w = c()
for (id in 1:nrow(df_window_merge)){
  w = paste(df_window_merge[id,"chr"], ":", df_window_merge[id,"start"], "-", df_window_merge[id,"end"], sep="")
  #system(paste("touch ", w, sep=""))
  list_w = c(list_w , w)
}
chunk <- function(x, n) split(x, sort(rank(x) %% n))
res = chunk(list_w, n = nwindow_list)
for(idl in 1:length(res)){
  r = res[[idl]]
  for(rr in r){
    cat(rr, file=paste("windows_id_", idl, ".txt", sep=""), append=T, sep= "\n")
  }
}