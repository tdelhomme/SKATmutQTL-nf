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
        OR --bed_overlap        - a dataframe with header CHR START END, windows overlapping each line (genomic pos) will be grouped for the skat test
      Optional arguments:
      --help \n\n")
  q(save="no")
}

if(is.null(args$df_windows)) {stop("Option --df_windows should be provided")} else{df_windows=args$df_windows}
if(is.null(args$nwindow_list) & is.null(args$bed_overlap)) {stop("Option --nwindow_list OR --bed_overlap should be provided")} 
if(!is.null(args$nwindow_list)) {nwindow_list=as.numeric(args$nwindow_list)} 
if(!is.null(args$bed_overlap)) {
  library(GenomicRanges)
  bed = makeGRangesFromDataFrame(read.table(args$bed_overlap, h=T, stringsAsFactors = F))
} 
load(df_windows)


if(!is.null(args$nwindow_list)){
  # create one empty file per window to run the scripts in parallel across windows
  list_w = c()
  for (id in 1:nrow(df_window_merge)){
    w = paste(df_window_merge[id,"chr"], ":", df_window_merge[id,"start"], "-", df_window_merge[id,"end"], sep="")
    #system(paste("touch ", w, sep=""))
    list_w = c(list_w , w)
  }
  chunk <- function(x, n) split(x, sort(rank(x) %% n))
  res = chunk(list_w, n = nwindow_list)
}

if(!is.null(args$bed_overlap)){
  res = list()
  df_window_merge_gr = makeGRangesFromDataFrame(df_window_merge[,c("chr", "start", "end")])
  tmp = as.data.frame(findOverlaps(bed, df_window_merge_gr))
  for(bed_wind_id in unique(tmp$queryHits) ){
    tmpw = df_window_merge[ tmp[which(tmp$queryHits == bed_wind_id), "subjectHits"], ] # for wind id in the bed, return overlap windows in df_window_merge
    ll = unlist(lapply(1:nrow(tmpw), function(id){
      paste(tmpw[id,"chr"], ":", tmpw[id,"start"], "-", tmpw[id,"end"], sep="")
    }))
    res[[bed_wind_id]] = ll
  }
  save(res, file="list_windows_per_signSNP.Rdata")
}

for(idl in 1:length(res)){
  r = res[[idl]]
  for(rr in r){
    cat(rr, file=paste("windows_id_", idl, ".txt", sep=""), append=T, sep= "\n")
  }
}