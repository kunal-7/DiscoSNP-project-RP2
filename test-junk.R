ssh kunal@137.205.69.144

run_discoSnp++.sh -r "path.txt" -p "ECOR_01"
(base) kunal@spectre:~/PROJECTS/trial$ run_discoSnp++.sh -r "/home/kunal/PROJECTS/trial/fof_path_for_files.txt"

/home/kunal/PROJECTS/trial

/media/storage/DATA/ECOR

017217_ECOR_01_1_trimmed.fastq.gz
017217_ECOR_01_2_trimmed.fastq.gz
017217_ECOR_01_U1_trimmed.fastq.gz
017217_ECOR_01_U2_trimmed.fastq.gz

ECOR_34 ECOR_48 ECOR_51 ECOR_62 not there

sink("mylist.txt")
print(mylist)
sink()

mkdir learning_c
cd learning_c
touch bspl{0001..0003}.c

rsync -av kunal@137.205.69.144:~/PROJECTS/filenames.txt ./CLIMB/"Filenames.txt"


DiscoSNP++    
  
  run_discoSnp++.sh -r "/home/kunal/PROJECTS/reads/path.txt"
system("~/DiscoSNP++-v2.5.4-bin-Linux/run_discoSnp++.sh -r fof.txt -T")

for(i in 1:50){
  cmd <- paste("~/DiscoSNP++-v2.5.4-bin-Linux/run_discoSnp++.sh -r", df[i,"filename"], "-T", sep = " ")
  system(cmd)
}


system("~/DiscoSNP++-v2.5.4-bin-Linux/run_discoSnp++.sh -r fof.txt -T")
for(i in 1:50){
  cmd <- paste("~/DiscoSNP++-v2.5.4-bin-Linux/run_discoSnp++.sh -r", df[i,"filename"], "-T", sep = " ")
  system(cmd)
}

cdsDiscoSNP++               
  
  
  
loopseq<- list.files("~/PROJECTS/reads",full.names=TRUE,pattern=".txt",include.dirs=TRUE)
write.csv2(seq,file = "loopseq.txt", quote =F, row.names = F)


getwd()
con <- readLines("loopseq.txt")
head(con)
for (i in 1:75){
  + command<-paste("run_discoSnp++.sh -r ",con[i]," -p sample_",as.character(i)," -T",sep = "")
  + print(command)
  
  
  datasets <- c("3mW.csv", "5mW.csv")
  
  d <- map_dfr(datasets, read_csv, .id = "experiment")
}

#######--------
dfa <- read.fasta("sample_69_k_31_c_3_D_100_P_3_b_0_coherent.fa", as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)

#first split up the fasta file to put each SNP into a single row
df <- tibble(mut = rep(1:(length(dfa)/2),each = 2), path = rep(c("higher", "lower"), length.out = length(dfa))) %>% 
  mutate(annotation = attr(dfa, "names"),
         sequence = unlist(dfa)) %>% 
  pivot_wider(names_from = path, values_from = c(annotation,sequence))

#then start parsing the annotation (NB this uses 33 because it's for 33 samples, for fewer samples it will need modifying/simplifying)
df <- df %>% 
  separate(annotation_higher, into = c("mut_higher", "polymorphism_higher", "complexity_higher", "npol_higher", "Lunitig_higher", "Runitig_higher", "Lcontig_higher", "Rcontig_higher",  paste("coverage",  "higher", sep = "_"),  paste("Q",  "higher", sep = "_"), paste("genotype",  "higher", sep = "_") , "rank"), sep = "\\|", remove = FALSE) %>% 
  separate(annotation_lower, into = c("mut_lower", "polymorphism_lower", "complexity_lower", "npol_lower", "Lunitig_lower", "Runitig_lower", "Lcontig_lower", "Rcontig_lower",  paste("coverage",  "lower", sep = "_"),  paste("Q",  "lower", sep = "_"), paste("genotype",  "lower", sep = "_") , "rank"), sep = "\\|", remove = FALSE) %>% 
  mutate(across(starts_with("coverage_"), function(x) str_split_fixed(x, "_", n =2)[,2] %>% as.numeric()))

#and save:
saveRDS(df, "DEallDf.rds")


