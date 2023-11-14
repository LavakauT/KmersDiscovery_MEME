# filter: one kmer match perfectly to the other one or pcc > 0.9-----
## ∆ CONSENSUS motif comparison------------
library(dplyr)
library(tibble)
library(universalmotif)
library(stringr)
library(tidyverse)
# library(ComplexHeatmap)
library(circlize)
library(gridtext)
library(pgirmess)


dir <- '/full/path/to/all/pCRE/folder'
dir2 <- '/full/path/to/all/ML_master/folder'
middle <- 'neg'
folder <- list.dirs(dir,
                    full.names = FALSE,
                    recursive = FALSE)

for (x in 1:length(folder)) {
  folder_name <- folder[x]
  filenames <- list.files(paste(dir, folder[x], sep = '/'),
                          pattern="*_FETresults.txt",
                          full.names=FALSE)
  
  filenames2 <- list.files(paste(dir, folder[x], sep = '/'),
                           pattern="*_df_p0.01.txt",
                           full.names=FALSE)
  
  # head(filenames)
  # head(filenames2)
  
  for (i in 1:10) {
    df2 <- read.delim(paste(dir, folder_name, filenames[i], sep = '/')) %>%
      select(1,4)
    names(df2) <- c('motif', 'pvalue')
    df2$times <- paste('neg', i, sep = '_')
    
    list.motif <- list()
    for (z in 1:nrow(df2)) {
      m <- create_motif(df2$motif[z], name = df2$motif[z], family = df2$times[z])
      list.motif <- c(list.motif, m)
    }
    
    comparisons <- compare_motifs(list.motif,
                                  method = "PCC",
                                  min.mean.ic = 0,
                                  score.strat = "a.mean")
    
    kmers <- colnames(comparisons)
    columns = c('motif', 'ppc', 'pvalue', 'group') 
    sub.com = data.frame(matrix(nrow = 0, ncol = length(columns))) 
    colnames(sub.com) = columns
    for (q in 1:nrow(df2)) {
      number <- q
      kmer <- kmers[q]
      sub.com2 <- data.frame(comparisons) %>%
        select(all_of(kmer)) %>% 
        rownames_to_column(var = 'kmer')
      
      colnames(sub.com2) <- c(kmer, 'ppc')
      sub.com2 <- sub.com2 %>% 
        filter(ppc > 0.9)
      colnames(sub.com2) <- c('motif', 'ppc')
      sub.com2 <- inner_join(sub.com2, df2, by = 'motif')
      sub.com2$group <- kmer
      sub.com <- rbind(sub.com, sub.com2)
    }
    
    sub.com3 <- sub.com %>% 
      group_by(group) %>% 
      summarise_all(min)
    
    sub.com3 <- sub.com3[,c('motif', 'pvalue')] %>% 
      group_by(motif) %>% 
      summarise(pvalue = min(pvalue))
    
    write.table(sub.com3,
                paste(dir2, folder_name,
                      paste0(middle, '_', folder_name, '_', i,'_distinct_pcc_enriched_kmer', '.txt'),
                      sep = '/'),
                row.names = F,
                quote = F,
                sep = '\t')
    
    df3 <- read.delim(paste(dir, folder_name, filenames2[i], sep = '/'))
    df3 <- cbind(df3[,1:2], df3[,sub.com3$motif])
    
    write.table(df3,
                paste(dir2, folder_name, filenames2[i], sep = '/'),
                row.names = F,
                quote = F,
                sep = '\t')
    
  }
  
}
