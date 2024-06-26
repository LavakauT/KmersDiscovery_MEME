# Motif Discovery
## repository source: Shiulab (https://github.com/ShiuLab/MotifDiscovery; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7671360/)
The latest pipeline uses SciPy to determine enrichment, Numpy and Pandas for dataframe management, and SciKit Learn for RandomForest (via ML-pipeline Github repository).

Original pipeline used python for processing dataframes and R for running Random Forest. For those instructions see bottom of the document (Python & R Pipeline).

## modification concept: Dr. Liu Ming-Jung and Dr. Wu Ting-Ying

If two k-mers were both significantly enriched and one k-mer sequence exactly matched the other one, only the one with the lower adjusted P value was retained (https://academic.oup.com/plcell/article/30/7/1445/6100070).

The similarity between two k-mers will be determined using the Pearson correlation coefficient (PCC), and k-mers with a PCC lower than 0.9 will be filtered out. Among these, only the k-mer with the lower adjusted P value will be retained.



# Environment Requirements
- biopython 1.78
- matplotlib 3.5.3
- numpy 1.21.5
- pandas 1.3.5
- python 3.7.0
- scikit-learn 1.0.2
- scipy 1.7.3
```
$ wget http://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh
$ bash ~/miniconda.sh -b -p $HOME/miniconda
$ export PATH="$HOME/miniconda/bin:$PATH"
# source conda.sh
$ conda create -n ml python==3.7.0 
$ conda install biopython
$ conda install matplotlib
$ conda install pandas
$ conda install scikit-learn

```

# Kmer finding and ML predictions
## Get positive and negative gene lists from log FC data set
  1. Get up-reg gene list (logFC >= 1, adj. p-val < 0.05), down-reg gene list (logFC <= -1, adj. p-val < 0.05), and negative gene list (-0.8 < logFC < 0.8, adj. p-val > 0.05).

You can adjust your negative data set with different ranges of logFC and adj. p-val. Select the range which is acceptable to ML results but please keep applying the same range in different datasets.

> Typically, the numbers of positive and negative genes in datasets are different, leading to an imbalance that can bias machine learning models toward the larger class (usually negative genes). To address this issue, we randomly downsize the dataset 10 times to achieve a balanced status before running machine learning algorithms. The downsizing process involves sampling with replacement. Please run [down_sizing.R](https://github.com/LavakauT/KmersDiscovery_MEME/blob/main/down_sizing.R) to get the subsets of negative gene lists.
![ml](https://github.com/LavakauT/KmersDiscovery_MEME/assets/132649549/a5d43f8b-a660-467c-a349-80005da40dfd)


```
$ conda info --envs
$ conda create -n R4.3.2
$ conda activate R4.3.2
$ conda install -c conda-forge r-base==4.3.2
$ conda install -c conda-forge r-biocmanager


$ conda install  -c conda-forge r-dplyr
$ conda install -c conda-forge r-stringr
$ conda install -c conda-forge r-tidyverse
$ conda install -c conda-forge r-circlize


# install constant package in conda env
$ R
>library(dplyr)
>library(stringr)
>library(tidyverse)
>library(tibble)
>library(BiocManager)

>BiocManager::install('universalmotif')
>BiocManager::install("ComplexHeatmap")
>BiocManager::install("gridtext")
```

  2. Download Genome fasta and annotation file

Please download the GFF file and genome FASTA corresponding to the version of your gene list.
> For Marchantia polymorpha, we can download MpTak_v6.1r1.gff and MpTak_v6.1r1.genome.fasta in MarpolBase(https://marchantia.info/download/MpTak_v6.1/).


  3. Get promoter coordinates based on TSS and 5' UTR or only TSS
     

![pro_cor](https://github.com/LavakauT/KmersDiscovery_MEME/assets/132649549/33cffc1a-255a-4dac-8369-12fe3b401b77)

The choice of promoter sequences depends on your preference. Please ensure that the same kind is used for both positive and negative datasets. In our dataset, we use promoter coordinates based on Transcription Start Site (TSS) and 5'UTR from Differentially Expressed Genes (DEGs) in RNA-seq. Additionally, please note that I have modified the original FastaManager.py into FastaManager_modified.py, which includes custom length settings.

get promoter coordinates based on transcription start site (TSS) and 5'UTR
```
python full/path/to/FastaManager_modified.py -f gff_prom_to_coord_5utr -gff [full/path/to/gff file]
# output: gff_file.coord
```

get promoter coordinates based on transcription start site (TSS)
```
python full/path/to/FastaManager_modified.py -f gff_prom_to_coord2 -gff [full/path/to/gff file]
# output: gff_file_prom-5utr.coord
```


  4. Get fasta sequence using coords file and genome fasta
```
# TSS and 5'UTR, output: gff_file_prom-5utr.coord.fa
python full/path/to/FastaManager_modified.py -f get_stretch4 -coords [full/path/to/prom-5utr.coord]  -fasta [full/path/to/genome.fa]

# TSS, output: gff_file.coord.fa
python full/path/to/FastaManager_modified.py -f get_stretch4 -coords [full/path/to/prom.coord]  -fasta [full/path/to/genome.fa]
```

To prevent the promoter sequences from mapping to incorrect gene annotation, which may include redundant words or other irrelevant information, it's important to clean up the gene annotation. Otherwise, it may lead to inaccurate mapping or interpretation of the data. Please enter R again and follow the script below to overwrite prom.coord.fa/prom-5utr.coord.fa:
```
# modify fasta file-------
library(Biostrings)
library(stringr)

pro <- readDNAStringSet('/full/path/to/MpTak_v6.1r1.gff_prom-5utr.coord.fa')
name <- names(pro)# extract gene name
name <- str_replace_all(name, "\\.\\w[:blank:].*", "") # process gene name
seq <- paste(pro) # extract promoter sequences
dfa <- data.frame(name, seq) # re-combine seq and name

# export fast function definition:
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

writeFasta(dfa, '/full/path/to/MpTak_v6.1r1.gff_prom-5utr.coord.fa')
```

5. Get fasta sequence to each gene list
```
# check before run
$ for inp in folder_to_all_gene_list/*.txt
> do 
> echo "python FastaManager_modified.py -f getseq2 -fasta full/path/to/gff_file_prom-5utr.coord.fa -name $inp"
> done

# then, run it
$ for inp in folder_to_all_gene_list/*.txt
> do 
> python FastaManager_modified.py -f getseq2 -fasta full/path/to/gff_file_prom-5utr.coord.fa -name $inp
> done
```

Till this steps, outputs with coordinated promtoer sequences should exist in where you put the gene list.

# Get enriched kmers

  5. Get enriched kmer dataframe using Fisher's Exact Test
you can put all negative fasta files in one folder to run loop
```
python full/path/to/pCRE_Finding_FET.py -pos <pos fasta> -neg <neg fasta> -k 6mer.txt -save <name of output files>

# loop
$ for inp in folder/to/all/*.txt.fa
> do
> neg=$inp
> python pCRE_Finding_FET.py -pos full/path/to/pos.txt.fa -neg $neg -k wgcna/6mer.txt FDR Y -save $inp.pcre
> done

for inp in tom/e_1kb/*.fa.out; do neg=$inp; python pCRE_Finding_FET.py -pos tom/e_1kb/e_1kb.fa.out -neg $neg -k tom/6mer.txt FDR Y -save $neg.pcre ; done
```
there're some parameters applied in this line, notice that FDR correction was run after enrichment test. Therefore, you should type -FDR Y:
```
 -pos_str  String for what codes for the positive example (Default = 1)
 -neg_str  String for what codes for the negative example (Default = 0)
 -k        List of kmers to start with (/mnt/home/azodichr/ML_Python/6mers.txt or 5mers.txt)
 -pval     P-value cut off for Fisher's exact test (Default = 0.01)
 -FDR      Default: N. Designate (Y/N) if you want to run FDR correction during enrichment test
```

the output:
```
 - save/file.txt.fa_df_p0.01.txt       Dataframe that goes into ML-pipeline
 - save/file.txt.fa_FETresults.txt     Dataframe that k-mers counts and p-value which can help us filtering exactly match ones
```

# Optional: 1,0 to pos,neg convertor
Please follow [pos_neg_conversion.R](https://github.com/LavakauT/KmersDiscovery_MEME/blob/main/pos_neg_conversion.R)


# Optional: PCC filtering
Please run [pcc_filtering.R](https://github.com/LavakauT/KmersDiscovery_MEME/blob/main/pcc_filtering.R) to rewrite all save/file.txt.fa_df_p0.01.txt files, then go to machine learning part.


# Machine learning
Use ML Pipeline (most recent version) here to get class predictions (see ML-pipeline https://github.com/bmmoore43/ML-Pipeline and https://github.com/LavakauT/ML-pipeline)
