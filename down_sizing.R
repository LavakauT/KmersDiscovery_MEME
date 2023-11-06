# get randomly negative data set-------
library(dplyr)
# library(pgirmess)
library(stringr)


dir <- 'path/to/output/downsizing/negative/gene/list/folder/'
middle <- 'neg'
neg <- read.delim('path/to/full/negative/gene/list')
n.neg <- nrow(neg)
pos <- read.delim('path/to/full/positive/gene/list')
n.pos <- nrow(pos)

for ( i in 1:10){
  number <- i
  set.seed(number)
  random.neg <- sample(seq_len(n.neg), size = n.pos)
  df <- neg[random.neg,]
  
  # output result
  write.table(df, 
              paste0(dir, paste('neg', number, sep = '_'), '.txt'),
              row.names = F,
              quote = F,
              sep = '\t')
}
