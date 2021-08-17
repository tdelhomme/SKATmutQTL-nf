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
      --input_file                - .Rdata input file

      Optional arguments:
      --help \n\n")
  q(save="no")
}

library(SKAT)

if(is.null(args$input_file)) {stop("Option --input_file should be provided")} else{input_file=args$input_file}
load(input_file)

if(is.na(input_data$mat_geno)){ 
  pval = 1
} else {
  Z = input_data$mat_geno
  y = input_data$mat_pheno
  X = matrix(rep(0, nrow(Z)), ncol=1)
  
  obj <- SKAT_Null_Model(y ~ X, out_type="C")
  pval = SKAT(Z, obj)$p.value
}

save(pval, file = paste(gsub("_input_skat.Rdata", "", input_file), "_pvalue", sep=""))
