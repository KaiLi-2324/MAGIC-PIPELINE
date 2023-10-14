# AB filter
# The input is a vcf file. According to the position corresponding to the AD field in the format, extract the ad of each mutation sample (starting from 0/1) and calculate ab = ad[2] / sum(ad)
# If ab meets the condition, the output will not be changed. If ab does not meet the condition, 0/1 will be changed to ./. for output.

library(stringr)
library(magrittr)
library(readr)
library(getopt)

# args
spec <- matrix(
  # There are five in each row, and the fifth one is optional, which means that the fifth column does not need to be written.
  # byrow Fill the elements of a matrix row by row
  # ncol  Fill each row with five elements
  c("input_vcf_file",  "v", 1, "character", "vcf file after vep",
    "output_file", "o", 1, "character",  "output_file_prefix",
    "minAB",  "i", 1, "double",  "min AB",
    "maxAB",  "I", 1, "double",  "max AB",
    "help",   "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
print(opt)
if( !is.null(opt$help) || is.null(opt$input_vcf_file) || is.null(opt$output_file)  || 
    is.null(opt$minAB)){
  # ... Here you can also customize some things to put in it
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

date()
con <- file(opt$input_vcf_file, "r")
out.file <- opt$output_file
minAB <- opt$minAB
maxAB <- opt$maxAB

ln_num <- 1
ln=readLines(con,n=1)
while(length(ln) != 0) {
  if(exists("vcf.col.nam")){
    vcf.col <- ln %>% str_split(., pattern = "\t") %>% unlist # Get the columns of vcf file
    format.col <- vcf.col[format.idx] %>% str_split(., pattern = ":") %>% unlist # Get the columns in the FORMAT field in the vcf file
    ad.idx <- str_which(format.col, pattern = "AD")
    mtt.idx <- grep("^0/1", vcf.col)
    vcf.col[mtt.idx] <- 
      sapply(mtt.idx, simplify = T, FUN = function(x){
        dat <- vcf.col[x]
        ad <- str_split(dat, pattern = ":", simplify = T) %>% 
          extract(., 1, ad.idx) %>% 
          str_split(., pattern = ",") %>% unlist %>%
          as.numeric()
        ab <- ad[2]/sum(ad)
        if((ab < minAB) | is.na(ab) | (ab > maxAB) ){ # If AB does not meet the condition of 0.25, replace it with ./.; if 0/1:.:0,0 will report an error, which is a very special case.
          str_sub(dat, 1, 3) <- "./."
        }
        return(dat)
      })
    ln.new <- str_c(vcf.col, collapse = "\t")
    write_lines(ln.new, path = out.file, append = T)
  }else{
    write_lines(ln, path = out.file, append = T)
    if(str_detect(ln, "^#CHROM")){
      vcf.col.nam <- ln %>% str_split(., pattern = "#") %>% unlist %>% extract(.,2) %>%
        str_split(., pattern = "\t") %>% unlist
      info.idx <- which(vcf.col.nam=="INFO")
      format.idx <- which(vcf.col.nam=="FORMAT")
      sam.stt.idx <- format.idx + 1
      sam.end.idx <- length(vcf.col.nam)
    } # Build vcf column name
  }
  ln=readLines(con,n=1)
  ln_num <- ln_num+1
}
close(con)
date()
