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
      --w                         - Window (e.g. chr1_1_10000)
      --somatic_files             - Text file containing in one column the list of the files to consider in input somatic folder
      --somatic_Rdata             - Rdata object containing the pre-loaded somatic files (gr_mut)
      --somatic_folder            - Input folder of somatic files
      --germline_VCF              - Germline VCF containing genotypes
      
      Optional arguments:
      --help \n\n")
  q(save="no")
}
library(data.table)
library(GenomicRanges)
library(readr)

if(is.null(args$somatic_files)) {stop("Option --somatic_files should be provided")} else{somatic_files=args$somatic_files}
if(is.null(args$somatic_Rdata)) {stop("Option --somatic_Rdata should be provided")} else{ load(args$somatic_Rdata) }
if(is.null(args$somatic_folder)) {stop("Option --somatic_folder should be provided")} else{somatic_folder=args$somatic_folder}
if(is.null(args$germline_VCF)) {stop("Option --germline_VCF should be provided")} else{germline_VCF=args$germline_VCF}
if(is.null(args$wlist)) {stop("Option --wlist should be provided")} else{all_w = as.character(read.table(args$wlist)[,1])}

#############################
##### SOMATIC PHENOTYPE #####
t2in<-fread(somatic_files, head=FALSE)
files <- as.character(t2in$V1)

for(w in all_w){
  
  print(paste(date(), " INFO: working on the window ", w, sep=""))
  
  # extract window from germline VCF file
  system(paste(" bcftools view -r ", w, " ", germline_VCF, " | bgzip -c > window.vcf.gz", sep=""))
  system("tabix -p vcf window.vcf.gz")
  
    # check if empty VCF
  if( length(system("zcat window.vcf.gz | grep \"^chr\" | head -n1", intern=T)) == 1 ) {
    
    gr_wind = GRanges(seqnames = unlist(strsplit(w,":"))[1], 
                      ranges = IRanges(start = as.numeric(unlist(strsplit(unlist(strsplit(w,":"))[2], "-"))[1]),
                                       end = as.numeric(unlist(strsplit(unlist(strsplit(w,":"))[2], "-"))[2])))
    
    res = lapply(files, function(f){
      gr_mut = all_gr_mut[[f]] 
      counts <- sum(countOverlaps(gr_wind, gr_mut)) # gr_mut is computed in extern
    })
    
    mat_pheno = matrix(unlist(res), ncol = 1)
    rownames(mat_pheno) = files
    
    ####################
    ##### GENOTYPE #####
    library(VariantAnnotation)
    
    vcf <- open(VcfFile("window.vcf.gz",  yieldSize=1000000))
    vcf_chunk = readVcf(vcf, "hg19")
    
    # if we have no somatic mutation in the window OR no germline variants return NAs
    if(length(table(as.numeric(mat_pheno[,1]))) == 0 | nrow(vcf_chunk) == 0 ){
      input_data = list("mat_geno" = NA, "mat_pheno" = NA)
    } else {
      # continue with the genotype
      nbc = 1
      seen_rs = c()
      
      #and continue
      while(dim(vcf_chunk)[1] != 0) {
        print(paste(date(), "INFO processing chunk: ", nbc, sep=""))
        # remove duplicated RS IDs -- take the first one
        rs = rownames(vcf_chunk)
        vcf_chunk = vcf_chunk[which(!duplicated(rs) & !(rs %in% seen_rs)),]
        rs = rownames(vcf_chunk)
        seen_rs = c(seen_rs, rs)
        
        chrom = as.character(seqnames(rowRanges(vcf_chunk,"seqnames")))
        pos = start(ranges(rowRanges(vcf_chunk,"seqnames")))
        gts = geno(vcf_chunk)$GT
        dos = apply(gts, 2, function(c) {
          c[which(c %in% c("0/0", "0|0"))] = 0
          c[which(c %in% c("1/0", "0/1", "0|1", "1|0"))] = 1
          c[which(c %in% c("1/1", "1|1"))] = 2
          as.numeric(c)
        })
        if(class(dos) == "numeric") { dos = t(as.data.frame(dos)) } # when only one variant, it is not a dataframe
        rownames(dos) = rownames(gts)
        res = cbind( data.frame(rs = rownames(dos), chrom=chrom, pos=pos), dos)
        
        if(exists("all_res")){ all_res = rbind(all_res, res) } else { all_res = res }
        
        nbc = nbc + 1
        vcf_chunk = readVcf(vcf, "hg19")
      }
      
      mat_geno = t(all_res[, 4:ncol(all_res)])
      rownames(mat_geno) = colnames(all_res)[4:ncol(all_res)]
      colnames(mat_geno) = all_res$rs
      # match with somatic names
      rownames(mat_geno) = unlist(lapply(rownames(mat_geno), function(x){
        if(!(grepl("DNA", x))) {res = unlist(strsplit(x,"_"))[1]} else {res=x}
        if(grepl("DNA", x)) return(gsub("-DNA_", "_DNA_", res))
        if(grepl("SM", x)) return(gsub("SM-", "SM_", res))
        if(!(grepl("DNA", x)) & !(grepl("SM", x))) return(res)
      }))
      kept_sm = intersect(rownames(mat_pheno), rownames(mat_geno))
      mat_geno = mat_geno[kept_sm,] # re-order similarly to the somatic samples
      
      input_data = list("mat_geno" = mat_geno, "mat_pheno" = mat_pheno)
      
    }
    outputf = paste(w, "_input_skat.Rdata", sep="")
    save(input_data, file=outputf)
  } else { print("WARNING: empty germline VCF") }
}
