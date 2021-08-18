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
      --manhattan_plot_function                - Manhattan plot function

      Optional arguments:
      --help \n\n")
  q(save="no")
}
if(is.null(args$manhattan_plot_function)) {stop("Option --manhattan_plot_function should be provided")} else{manhattan_plot_function=args$manhattan_plot_function}

source(manhattan_plot_function)

all_p_files = list.files(path=".", pattern="pvalue")

all_p = c()
for(f in all_p_files){
  load(f)
  all_p = c(all_p, pval)
  rm(pval)
}

w =  gsub("_pvalue", "", all_p_files)
dd = data.frame(pvalues = all_p,
                chr=unlist(lapply(w, function(w0) unlist(strsplit(w0, ":"))[1])),
                pos=unlist(lapply(w, function(w0) { ( as.numeric(unlist(strsplit(unlist(strsplit(w0, ":"))[2],"-"))[1]) + 5000 ) })))

manhattan.plot(dd$chr, dd$pos, dd$pvalue)

# second plot
library(qqman)
dd$snp = "snp"
dd$chr_name = as.numeric(gsub("chr", "", dd$chr))
manhattan(dd, chr="chr_name", bp="pos", p="pvalues", snp="snp")

