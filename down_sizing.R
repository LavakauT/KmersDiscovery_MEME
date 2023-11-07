# get randomly negative data set-------
library(dplyr)
# library(pgirmess)
library(stringr)


dir <- '/RAID1/working/R425/lavakau/'
sub.dir <- 'sample_data/'
middle <- 'neg'
neg <- read.delim('/RAID1/working/R425/lavakau/pCRE/sample_data/NN.txt')
n.neg <- nrow(neg)
group <- list.files(paste0(dir, '/', 'sample_data'), pattern = '*.txt', full.names = TRUE)


pos <- read.delim('/RAID1/working/R425/lavakau/pCRE/sample_data/UU.txt')
n.pos <- nrow(pos)

for ( i in 1:10){
  number <- i
  set.seed(number)
  random.neg <- sample(seq_len(n.neg), size = n.pos)
  df <- neg[random.neg,]
  
  # output result
  write.table(df, 
              paste0(dir, sub.dir, paste('neg', number, sep = '_'), '.txt'),
              row.names = F,
              quote = F,
              sep = '\t')
}
