library(seqinr)
library(tidyverse)
library(fs)

datasets<-list.files(pattern="coherent.fa")
paste(datasets)

parsing<-function(x){
  dfa <- read.fasta(x, as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)
  
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
  }

df<-map_dfr(datasets, parsing, .id = "experiment")

#and save:
saveRDS(df, "allDf.rds")

df %>%
  ggplot(aes(coverage_higher+coverage_lower))+
  geom_histogram(bins = 100)

hist(df$coverage_higher+df$coverage_lower) 

df %>%
  mutate(total_coverage = coverage_higher + coverage_lower) %>%
  ggplot(aes(y = total_coverage)) +
  geom_histogram(bins = 50)


ggsave("totalCoverageHistogram.pdf") 

getwd()

data<-read_rds("allDF.rds")

data %>%
  ggplot(aes(coverage_higher+coverage_lower))+
  geom_histogram(bins = 100)+xlim(0,1250)


#total coverage
data %>%
  mutate(total_coverage = coverage_higher + coverage_lower) %>%
  filter(data$complexity_higher=="high" & data$complexity_lower=="high") %>%
  ggplot(aes(x = total_coverage)) +
  geom_histogram(bins = 100)+ xlim(0,1250)
  #geom_histogram(bins = 100)+ xlim(100,400) 


#minor allele frequency
data <- data %>%
  mutate(total_coverage = coverage_higher + coverage_lower, 
        higher_frequency=coverage_higher/total_coverage, 
        MAF=ifelse(higher_frequency>0.5,1-higher_frequency,higher_frequency))  

  data %>%
    filter(complexity_higher=="high" & complexity_lower=="high" & total_coverage>=100 & total_coverage<=400) %>%
    ggplot(aes(x = MAF)) +
    geom_histogram(bins = 100)

  
  
filtered_data<- data %>%
  filter(complexity_higher=="high" & complexity_lower=="high" & total_coverage>=100 & total_coverage<=400) %>%
  group_by(experiment)%>%
  #filter(complexity_higher=="high" & complexity_lower=="high" & total_coverage>=100 & total_coverage<=400 & MAF>0.015 & MAF<0.03) %>%  
  arrange(desc(MAF)) %>% 
  mutate(num_above=row_number())

filtered_data %>%
  #filter(experiment%in%c(as.character(1:75))& MAF>0.01 &MAF<0.035)%>%
  ggplot(aes(x=MAF,y=num_above,colour=experiment))+
  geom_line()+ scale_x_log10()+scale_y_log10()
  #facet_wrap(~experiment)





