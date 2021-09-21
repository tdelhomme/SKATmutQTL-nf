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
      --somatic_files             - Text file containing in one column the list of the files to consider in input somatic folder
      --somatic_folder            - Input folder of somatic files

      Optional arguments:
      --help \n\n")
  q(save="no")
}
library(data.table)
library(GenomicRanges)
library(readr)

# load somatic data 
if(is.null(args$somatic_files)) {stop("Option --somatic_files should be provided")} else{somatic_files=args$somatic_files}
if(is.null(args$somatic_folder)) {stop("Option --somatic_folder should be provided")} else{somatic_folder=args$somatic_folder}

#############################
##### SOMATIC PHENOTYPE #####

t2in<-fread(somatic_files, head=FALSE)
files <- as.character(t2in$V1)

#load data
all_gr_mut = lapply(files, function(f){
  ff = system(paste("ls ", somatic_folder, "/*", f, ".csv", sep=""), intern=T)
  mut = read_csv(ff)
  
  # change colnames for consistency with Hartwig
  if("reference_allele" %in% colnames(mut)) colnames(mut)[which(colnames(mut) == "reference_allele")] = "REF"
  if("mutated_to_allele" %in% colnames(mut)) colnames(mut)[which(colnames(mut) == "mutated_to_allele")] = "ALT"
  
  mut = mut[(mut$REF %in% c('A', 'T', 'C', 'G') & mut$ALT %in% c('A', 'T', 'C', 'G')),]
  mut = mut[!mut$chr %in% c('chrY'),]
  
  # make a VRange obj
  gr_mut = GRanges(seqnames = mut$chr, ranges = IRanges(start = mut$start, end = mut$end))
})
names(all_gr_mut) = files

save(all_gr_mut, file="gr_mut_somatic.Rdata")