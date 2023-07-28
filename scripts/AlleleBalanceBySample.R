# AB filter
# 输入为vcf文件，按照format中的AD字段对应的位置，提取每个突变样本（0/1开始）的ad，计算ab = ad[2] / sum(ad)
# 如果ab满足条件，不进行改变，输出，如果ab不满足条件，将0/1变为./.，输出

library(stringr)
library(magrittr)
library(readr)
library(getopt)

# 参数描述
spec <- matrix(
  # 每行五个，第五个可选，也就是说第五列可以不写
  # byrow 按行填充矩阵的元素
  # ncol  每行填充五个元素
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
  # ... 这里你也可以自定义一些东放在里面
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
    vcf.col <- ln %>% str_split(., pattern = "\t") %>% unlist # 得到vcf文件的列
    format.col <- vcf.col[format.idx] %>% str_split(., pattern = ":") %>% unlist # 得到vcf文件中FORMAT字段内的列
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
        if((ab < minAB) | is.na(ab) | (ab > maxAB) ){ # 如果AB不满足0.25的条件，替换为./. ；如果0/1:.:0,0会报错，极特殊情况
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
    } # 构建vcf列名
  }
  ln=readLines(con,n=1)
  ln_num <- ln_num+1
}
close(con)
date()
