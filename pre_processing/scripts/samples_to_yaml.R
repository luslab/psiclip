library(data.table)
library(dplyr)

mdat = fread("../../metadata.tsv")
head(mdat)


for (row in 1:nrow(mdat)){
  r = mdat[row,]
  write(paste0(r$`Source Name`, ': ["data/fastq/', r$`Comment[ENA_RUN]`, '.fq.gz", ""]'), 
        file="TEMP_CONFIG.txt",
        append=TRUE)
  }


