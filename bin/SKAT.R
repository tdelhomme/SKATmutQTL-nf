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

all_input_files = system(paste("ls *.Rdata"), intern=T)

for(input_file in all_input_files){
  
  load(input_file)
  
  if(length(table(input_data$mat_pheno))==1){ 
    pval = 1
  } else {
    Z = input_data$mat_geno
    y = input_data$mat_pheno
    X = matrix(rep(0, nrow(Z)), ncol=1)
    
    obj <- SKAT_Null_Model(y ~ X, out_type="C")
    pval = SKAT(Z, obj, method = "optimal.adj")$p.value
  }
  
  save(pval, file = paste(gsub("_input_skat.Rdata", "", input_file), "_pvalue", sep=""))
}
