#Looking for polymorphisms de novo without alignment ####
#currently using DiscoSNP++, though may also want to try ebwt2Indel (which requires BCR_LRP_GSA or similar to do a Burrows weaver transform on the data first)

#fuller file elsewhere, but general approach:

#First get the relevant samples locally (this is the kind of thing I used at a mac command line) ####
rsync -vzPR mcdssck2@rds-ssh.itservices.manchester.ac.uk:/mnt/cifs-KnightLab/Entire_Share/Common/SequenceData/SoilMetagenomics/Trimmed/Sample_75-185 ./Sample_75-185

#I've generally run discoSNP from R because it was simple, but it isn't necessary:

system("~/Desktop/DiscoSNP++-v2.5.4-Source/run_discoSnp++.sh -r DEallfof.txt -T -p 'DEall'")

#That took ~2 days to run on my (old but high memory) desktop for all 33 German samples
#That's more complex than what you, Kunal, need to do (where you only need to consider a single sample), but for Xuan, that input file (DEallfof.txt) and the files it points to were created using the following code: #####
setwd("~/Dropbox (The University of Manchester)/knightmareDocs/soil/OcÃ©ane/sequence/metagenomics/")
library(tidyverse)
library(fs)

key <- read_csv(file = "KeyToMetagenomicSamplesAug.csv")
dir_create("DEall")

keyDE <- key %>% 
  filter(Country == "DE")

findFiles <- function(FolderName, Label, ...){
  files <- list.files(path = paste("/Volumes/K4/Dropbox (The University of Manchester)/knightmareDocs/soil/metagenomicData/Trimmed/", folderName, sep = "" ), pattern = ".fastq.gz", recursive = TRUE, full.names = TRUE )
  files <- paste("\"", files, "\"", sep = "")
  fof <- paste("./DEall/", Label, ".txt", sep = "")
  cat(files, file = fof, sep = "\n")
  return(paste( Label, ".txt", sep = ""))}

keyDE <- keyDE %>% 
  mutate(fof = pmap_chr(., findFiles))

cat(paste("\"", "/Volumes/K4/Dropbox (The University of Manchester)/knightmareDocs/soil/OcÃ©ane/sequence/metagenomics/DEall/", keyDE$fof, "\"", sep = ""), file = "./DEall/DEallfof.txt", sep = "\n")





#This is then code to parse the output:
library(seqinr)

dfa <- read.fasta("DEall_k_31_c_3_D_100_P_3_b_0_coherent.fa", as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)

#first split up the fasta file to put each SNP into a single row
df <- tibble(mut = rep(1:(length(dfa)/2),each = 2), path = rep(c("higher", "lower"), length.out = length(dfa))) %>% 
  mutate(annotation = attr(dfa, "names"),
         sequence = unlist(dfa)) %>% 
  pivot_wider(names_from = path, values_from = c(annotation,sequence))

#then start parsing the annotation (NB this uses 33 because it's for 33 samples, for fewer samples it will need modifying/simplifying)
df <- df %>% 
  separate(annotation_higher, into = c("mut_higher", "polymorphism_higher", "complexity_higher", "npol_higher", "Lunitig_higher", "Runitig_higher", "Lcontig_higher", "Rcontig_higher",  paste("coverage", 1:33, "higher", sep = "_"),  paste("Q", 1:33, "higher", sep = "_"), paste("genotype", 1:33, "higher", sep = "_") , "rank"), sep = "\\|", remove = FALSE) %>% 
  separate(annotation_lower, into = c("mut_lower", "polymorphism_lower", "complexity_lower", "npol_lower", "Lunitig_lower", "Runitig_lower", "Lcontig_lower", "Rcontig_lower",  paste("coverage", 1:33, "lower", sep = "_"),  paste("Q", 1:33, "lower", sep = "_"), paste("genotype", 1:33, "lower", sep = "_") , "rank"), sep = "\\|", remove = FALSE) %>% 
  mutate(across(starts_with("coverage_"), function(x) str_split_fixed(x, "_", n =2)[,2] %>% as.numeric()))

#and save:
saveRDS(df, "DEallDf.rds")


#NB this can get too big to work with, e.g. to count up polymorphisms in different samples using these functions:####
WhichPoly <- function(x) {
  tb <- ifelse(x$higher > 0 & x$lower > 0, "poly", "mono")
  whichPoly <- paste(as.numeric(tb == "poly"), collapse ="")
  return(whichPoly)
}
Npoly3 <- function(x) {
  sum(as.numeric(unlist(strsplit(x, split = ""))))
}

# One way of cutting things down is to iterate through in chunks (this also fits a model to look for significant variation - might want to remove that bit)

N <- 100
for(i in 1:N){
  df <- readRDS("DEallDf.rds") %>% 
    slice(ceiling(((i-1)*nrow(.)/N)+1):ceiling(((i)*nrow(.)/N))) %>% 
     pivot_longer(starts_with("coverage"), names_to = c(NA, "file", ".value"), names_sep = "_") %>% 
    select(-starts_with("genotype"), -starts_with("Q_")) %>% 
    #left_join(keyDE) %>% 
    group_by(across(c(1:20))) %>% #this differs among files - should be the columns that define the SNP, not the sample
    nest() %>%   
    mutate(whichPoly = map_chr(data, WhichPoly),
           Npoly = map_dbl(whichPoly, Npoly3)) %>% 
   # mutate(m = map(data, ~glm(cbind(higher,lower)~file, data = ., family = binomial)),
           #p = map_dbl(m, ~anova(., test = "LRT")$'Pr(>Chi)'[2])) %>% 
    saveRDS(file = paste("DEalldf", i, ".rds", sep = ""))
  rm(df)
  gc()
}
