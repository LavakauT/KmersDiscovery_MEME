# loop code for change Class annotation from 1,0 to pos,neg after pCRE-------
library(dplyr)
# library(pgirmess)
library(stringr)


dir <- 'path/to/output/downsizing/negative/gene/fasta/folder/'
middle <- 'neg'
ends <- '.txt.fa_pcre_df_p0.01.txt'

for ( i in 1:10){
  number <- i
  df <- read.delim(paste0(dir, paste(middle, number, sep = '_'), ends))
  df$Class <- str_replace_all(df$Class, '1', 'pos') # class 1 replace with postive
  df$Class <- str_replace_all(df$Class, '0', 'neg')# class 0 replace with negative
  
  # output result
  write.table(df,
              paste0(dir, paste(middle, number, sep = '_'), ends)),
  row.names = F,
  quote = F,
  sep = '\t'}