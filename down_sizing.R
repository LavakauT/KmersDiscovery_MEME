#### get randomly negative data set ####
#### load packages ####
library(dplyr)
library(stringr)


dir <- '/RAID1/working/R425/lavakau/pCRE/'
sub.dir <- 'sample_data/'
sub.dir2 <- 'sample_data'
middle <- 'neg'
neg <- read.delim(paste0(dir, sub.dir, 'other/', 'NN.txt'))
n.neg <- nrow(neg)
groups <- list.files(paste0(dir, sub.dir2), pattern = '*.txt', full.names = FALSE)


for (j in 1:length(groups)) {
  group <- groups[j]
  group_name <- str_replace(group, '.txt', '')
  pos <- read.delim(paste0(dir, sub.dir, group_name, '.txt'))
  n.pos <- nrow(pos)
  
  for ( i in 1:10){
    number <- i
    set.seed(number)
    random.neg <- sample(seq_len(n.neg), size = n.pos)
    df <- neg[random.neg,]
    
    # output result
    write.table(df, 
                paste0(dir, sub.dir, group_name, '/',
                       paste('neg', number, sep = '_'), '.txt'),
                row.names = F,
                quote = F,
                sep = '\t')
  }
}
