# start some overall analysis using the DE samples - have all of them done together and a randomised 100k from that.

library(tidyverse)
setwd("~/Dropbox (The University of Manchester)/knightmareDocs/soil/OcÃ©ane/sequence/metagenomics/DiscoSNP/DEall/")

df <- readRDS("DE100kstk1.rds")


df %>% 
  ggplot(aes(x = Npoly)) + 
  # geom_histogram(bins =100, aes(y=after_stat(density))) +
  geom_histogram(bins =100) +
  # geom_freqpoly(bins = 100) +
  scale_x_continuous("Number of polymorphic samples") +
  scale_y_continuous("SNP count") +
  theme_classic()

ggsave("./plots/Npoly1e5.pdf", scale = 0.5)


# 2483 indels to 97517 SNPs i.e. ~2.5% indels
di <- df %>% 
  separate(mut_higher, sep ="_",into = c("MutationType", NA, NA, NA),remove = FALSE ) %>% 
  ggplot(aes(x = MutationType))+
  geom_bar()

# plan
# c) remove loci with alleles perfectly separated among soils (never compete)
# b) remove all low complexity SNP
# a) assume everything with both alleles present in >1 soil is inter-specific variation #Expectation that those present in >1 soil will be better conserved genes - could check by BLAST (e.g. number of different spp with 100% coverage)
# d2) check for significant variation among sample form single soil where polymorphic
# remove everything that isn't directly polymorphic in at least one single sample (could be conserved genes rare)
# d) check for significant variation within soils - remove those that don't vary
# e) focus on initial disturbance those where at least one treatment varies significantly from the control (?)
# f) calculate fitness of strain (both alleles combined) and relative fitness of alleles (sign doesn't matter) expect v weak U-shaped association - higher relative fitness of alleles (the causal ones) in strains that 

rs <- df$data[[1]]$Replicate_site
NpolySoils <- function(x) {
  polys <- as.numeric(unlist(strsplit(x, split = "")))
  tb <- data.frame(rs, polys) %>% 
    group_by(rs) %>% 
    summarise(Nps = sum(polys))
  
  tb %>% 
    summarise(Npoly = sum(Nps > 0)) %>% 
    unlist()
}

df <- df %>%
  mutate(NpolySoils = map_dbl(whichPoly, NpolySoils))

#see how many of the soils are polymorphic
df %>% 
  ggplot(aes(x = NpolySoils)) + 
  geom_histogram(bins =10) +
  scale_x_continuous("Number of polymorphic soils") +
  scale_y_continuous("SNP count") +
  theme_classic()

ggsave("./plots/NpolySoils1e5.pdf", scale = 0.5)


#for the SNPs that are polymorphic within only one soil, how many samples are they polymorphic in? (maximum of 10)
df %>% 
  filter(NpolySoils == 1) %>% 
  ggplot(aes(x = Npoly)) + 
  geom_histogram(bins =30) +
  scale_x_continuous("Number of polymorphic samples", breaks = c(1,5,9)) +
  scale_y_continuous("SNP count") +
  theme_classic()

ggsave("./plots/NpolyInSoil1e5.pdf", scale = 0.5)


#now focus down on those that are polymorphic within the initial sample
ini <- df$data[[1]]$Sampling_time
init <- function(x){
  polys <- as.numeric(unlist(strsplit(x, split = "")))
  ifelse(sum(polys > 0 & ini == "Initial")>0, "iniPoly", "iniMono")
}
whichInit <- function(x){
  polys <- as.numeric(unlist(strsplit(x, split = "")))
  paste("DE", which(polys > 0 & ini == "Initial")-30,sep = "")
}

#21k polymophic within one sample, of which ~6k are polymorphic in the initial sample
dfI <- df %>% 
  filter(NpolySoils == 1) %>%
  mutate(InitPoly = map_chr(whichPoly, init)) %>% 
  filter(InitPoly == "iniPoly") %>% 
  mutate(whichInitPoly = map_chr(whichPoly, whichInit))

#interestingly, different soils have different amounts of standing variation:
dfI %>% 
  ggplot(aes(x = whichInitPoly)) + 
  geom_bar() +
  scale_x_discrete("Initial Soil") +
  scale_y_continuous("SNP count") +
  theme_classic()

ggsave("./plots/NInitSNP1e5.pdf", scale = 0.5)

#Again, most of that standing variation only appears in that initial sample, dropping down to v few in 10 and none in all 11 samples.
dfI %>% 
  ggplot(aes(x = Npoly)) + 
  geom_histogram(bins =30) +
  scale_x_continuous("Number of polymorphic samples", breaks = c(1,5,9)) +
  scale_y_continuous("SNP count") +
  theme_classic()

ggsave("./plots/NpolyInit1e5.pdf", scale = 0.5)


#look at the fitness and fates of those SNPs

#can now look at r (selection rate coefficient, see http://myxo.css.msu.edu/ecoli/srvsrf.html) for these SNPs
#NB has units, so if comparing S1 to initial and S4 to S1, need to divide latter by 4 because it's over 4x the period (4wk vs 1wk)
#may want to add 1 to everything to stop negative infinite coefficients for everything that goes extinct

#slow way of going about things...!
rCtf <- function(dat, whichSoil){
  xI <- dat %>% 
    filter(Replicate_site == whichSoil) %>% 
    mutate(minor = if(higher[11] > lower[11]){lower}else{higher},
           major = if(higher[11] > lower[11]){higher}else{lower}) %>% 
    filter(Sampling_time == "Initial" | Treatment == "Ct")
  rCtS1 <- log2(xI$minor[1]/xI$minor[3]) - log2(xI$major[1]/xI$major[3])
  rCtS4 <- (log2(xI$minor[2]/xI$minor[1]) - log2(xI$major[2]/xI$major[1]))/4
  r1CtS1 <- log2((xI$minor[1]+1)/(xI$minor[3]+1)) - log2((xI$major[1]+1)/(xI$major[3]+1))
  r1CtS4 <- (log2((xI$minor[2]+1)/(xI$minor[1]+1)) - log2((xI$major[2]+1)/(xI$major[1]+1)))/4
  return(tibble(rCtS1, rCtS4, r1CtS1, r1CtS4))
}


dfI <- dfI %>% 
  mutate(rCt = map2(.x = data, .y = whichInitPoly, .f =rCtf)) %>% 
  unnest(rCt)

#problem is that, when both alleles go extinct, adding 1 will always make the minor allele fitter
#when only one goes extinct, still becomes infinite without adding one - if it's the major allele it becomes negative, if minor, positive
#nonethesless, gives some answers - most SNPs survive S1 in control, but not all the way to S4

table(!is.na(dfI$rCtS1))
table(!is.na(dfI$rCtS4))

#but of those that do, most are beneficial in S1, but become neutral after that

dfI %>% 
  ggplot(aes(x = rCtS1))+
  geom_density() +
  # geom_histogram(bins = 70) +
  scale_x_continuous("Selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)

ggsave("./plots/ControlSRCS11e5.pdf", scale = 0.5)


dfI %>% 
  ggplot(aes(x = rCtS4))+
  geom_density() +
  scale_x_continuous("Selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)

ggsave("./plots/CtSRCS41e5.pdf", scale = 0.5)


#but look at everything, using the +1:

dfI %>% 
  # filter(is.na(rCtS1)) %>%
  # filter(is.infinite(rCtS1)) %>%
  filter(is.infinite(rCtS1)|is.na(rCtS1)) %>%
  ggplot(aes(x = r1CtS1))+
  geom_density() +
  # geom_histogram(bins = 70) +
  scale_x_continuous("Selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)

ggsave("./plots/CtSRC1S11e5.pdf", scale = 0.5)


dfI %>% 
  ggplot(aes(x = r1CtS4))+
  # geom_histogram(bins = 70) +
  geom_density() +
  scale_x_continuous("Selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)

ggsave("./plots/CtSRCS41e5.pdf", scale = 0.5)


#may be misleading...look at fates in either S1 or S4 (doesn't say if it's the same ones at S4)
Fate <- function(dat, whichSoil){
  xI <- dat %>% 
    filter(Replicate_site == whichSoil) %>% 
    mutate(minor = if(higher[11] > lower[11]){lower}else{higher},
           major = if(higher[11] > lower[11]){higher}else{lower})
  xI1 <- xI %>% 
    filter(Sampling_time == "S1")
  ftS1 <- if(isTRUE(all.equal(xI1$major, rep(0,5))) & isTRUE(all.equal(xI1$minor, rep(0,5)))){"Both extinct"}else{
    if(isTRUE(all.equal(xI1$minor, rep(0,5)))){"Minor extinct"}else{
      if(isTRUE(all.equal(xI1$major, rep(0,5)))){"Major extinct"}else{
        if(sum(as.logical(xI1$major) & as.logical(xI1$minor))>0){"Polymorphism persists"}else{
          "Segregation"
        }
      }
    }
  }
  xI4 <- xI %>% 
    filter(Sampling_time == "S4")
  ftS4 <- if(isTRUE(all.equal(xI4$major, rep(0,5))) & isTRUE(all.equal(xI4$minor, rep(0,5)))){"Both extinct"}else{
    if(isTRUE(all.equal(xI4$minor, rep(0,5)))){"Minor extinct"}else{
      if(isTRUE(all.equal(xI4$major, rep(0,5)))){"Major extinct"}else{
        if(sum(as.logical(xI4$major) & as.logical(xI4$minor))>0){"Polymorphism persists"}else{
          "Segregation"
        }
      }
    }
  }
  return(tibble(ftS1, ftS4))
}

dfI <- dfI %>% 
  mutate(fate = map2(data, whichInitPoly, Fate)) %>% 
  unnest(fate)

dfI %>% 
  pivot_longer(c(ftS1, ftS4), names_to = "Sampling_time", values_to = "fate") %>% 
  mutate(Sampling_time = recode(Sampling_time, ftS1 = "S1", ftS4 = "S4")) %>% 
  ggplot(aes(Sampling_time, fill = fate))+
  geom_bar() +
scale_y_continuous("Number of SNPs") +
  scale_x_discrete("Sampling time") +
  theme_classic()

ggsave("./plots/Fates1e5.pdf", scale = 0.5)


#also try for 4 treatments:
rDf <- function(dat, whichSoil){
  xI <- dat %>% 
    filter(Replicate_site == whichSoil) %>% 
    mutate(minor = if(higher[11] > lower[11]){lower}else{higher},
           major = if(higher[11] > lower[11]){higher}else{lower}) %>% 
    filter(Sampling_time == "Initial" | Treatment == "D")
  rDS1 <- log2(xI$minor[1]/xI$minor[3]) - log2(xI$major[1]/xI$major[3])
  rDS4 <- (log2(xI$minor[2]/xI$minor[1]) - log2(xI$major[2]/xI$major[1]))/4
  r1DS1 <- log2((xI$minor[1]+1)/(xI$minor[3]+1)) - log2((xI$major[1]+1)/(xI$major[3]+1))
  r1DS4 <- (log2((xI$minor[2]+1)/(xI$minor[1]+1)) - log2((xI$major[2]+1)/(xI$major[1]+1)))/4
  return(tibble(rDS1, rDS4, r1DS1, r1DS4))
}

rFdf <- function(dat, whichSoil){
  xI <- dat %>% 
    filter(Replicate_site == whichSoil) %>% 
    mutate(minor = if(higher[11] > lower[11]){lower}else{higher},
           major = if(higher[11] > lower[11]){higher}else{lower}) %>% 
    filter(Sampling_time == "Initial" | Treatment == "Fd")
  rFdS1 <- log2(xI$minor[1]/xI$minor[3]) - log2(xI$major[1]/xI$major[3])
  rFdS4 <- (log2(xI$minor[2]/xI$minor[1]) - log2(xI$major[2]/xI$major[1]))/4
  r1FdS1 <- log2((xI$minor[1]+1)/(xI$minor[3]+1)) - log2((xI$major[1]+1)/(xI$major[3]+1))
  r1FdS4 <- (log2((xI$minor[2]+1)/(xI$minor[1]+1)) - log2((xI$major[2]+1)/(xI$major[1]+1)))/4
  return(tibble(rFdS1, rFdS4, r1FdS1, r1FdS4))
}

rFzf <- function(dat, whichSoil){
  xI <- dat %>% 
    filter(Replicate_site == whichSoil) %>% 
    mutate(minor = if(higher[11] > lower[11]){lower}else{higher},
           major = if(higher[11] > lower[11]){higher}else{lower}) %>% 
    filter(Sampling_time == "Initial" | Treatment == "Fz")
  rFzS1 <- log2(xI$minor[1]/xI$minor[3]) - log2(xI$major[1]/xI$major[3])
  rFzS4 <- (log2(xI$minor[2]/xI$minor[1]) - log2(xI$major[2]/xI$major[1]))/4
  r1FzS1 <- log2((xI$minor[1]+1)/(xI$minor[3]+1)) - log2((xI$major[1]+1)/(xI$major[3]+1))
  r1FzS4 <- (log2((xI$minor[2]+1)/(xI$minor[1]+1)) - log2((xI$major[2]+1)/(xI$major[1]+1)))/4
  return(tibble(rFzS1, rFzS4, r1FzS1, r1FzS4))
}

rHf <- function(dat, whichSoil){
  xI <- dat %>% 
    filter(Replicate_site == whichSoil) %>% 
    mutate(minor = if(higher[11] > lower[11]){lower}else{higher},
           major = if(higher[11] > lower[11]){higher}else{lower}) %>% 
    filter(Sampling_time == "Initial" | Treatment == "H")
  rHS1 <- log2(xI$minor[1]/xI$minor[3]) - log2(xI$major[1]/xI$major[3])
  rHS4 <- (log2(xI$minor[2]/xI$minor[1]) - log2(xI$major[2]/xI$major[1]))/4
  r1HS1 <- log2((xI$minor[1]+1)/(xI$minor[3]+1)) - log2((xI$major[1]+1)/(xI$major[3]+1))
  r1HS4 <- (log2((xI$minor[2]+1)/(xI$minor[1]+1)) - log2((xI$major[2]+1)/(xI$major[1]+1)))/4
  return(tibble(rHS1, rHS4, r1HS1, r1HS4))
}

#put together in one go from df
dfI <- df %>% 
  filter(NpolySoils == 1) %>%
  mutate(InitPoly = map_chr(whichPoly, init)) %>% 
  filter(InitPoly == "iniPoly") %>% 
  mutate(whichInitPoly = map_chr(whichPoly, whichInit)) %>% 
  # slice(1:10) %>% 
##########ä¸é¢å¯ä»¥æ¿æ¢æ####
  dfI <- dfIni %>%
  mutate(rCt = map2(.x = data, .y = whichInitPoly, .f =rCtf),
         rD = map2(.x = data, .y = whichInitPoly, .f =rDf),
         rFd = map2(.x = data, .y = whichInitPoly, .f =rFdf),
         rFz = map2(.x = data, .y = whichInitPoly, .f =rFzf),
         rH = map2(.x = data, .y = whichInitPoly, .f =rHf)) %>%  
  unnest(c(rCt,rD, rFd, rFz, rH)) %>% 
  mutate(rdDS1 = rDS1 - rCtS1,
         rdFdS1 = rFdS1 - rCtS1,
         rdFzS1 = rFzS1 - rCtS1,
         rdHS1 = rHS1 - rCtS1,
         rdDS4 = rDS4 - rCtS4,
         rdFdS4 = rFdS4 - rCtS4,
         rdFzS4 = rFzS4 - rCtS4,
         rdHS4 = rHS4 - rCtS4) 

saveRDS(dfI, file = "dfI1e5.rds")

dfI <- readRDS(file = "dfI1e5.rds")

dfI %>% 
  pivot_longer(c(r1CtS1, r1DS1, r1FdS1, r1FzS1, r1HS1), names_to = "treatment", values_to = "r1S1") %>% 
  ggplot(aes(x = r1S1, fill = treatment))+
geom_density(alpha = 0.2) +
scale_x_continuous("Selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)


ggsave("./plots/SRC1S11e5.pdf", scale = 0.5)

dfI %>% 
  pivot_longer(c(rCtS1, rDS1, rFdS1, rFzS1, rHS1), names_to = "treatment", values_to = "rS1") %>% 
  ggplot(aes(x = rS1, fill = treatment))+
  geom_density(alpha = 0.2) +
  scale_x_continuous("Selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)


ggsave("./plots/SRCS11e5.pdf", scale = 0.5)


#look at all the differences from control:
dfI %>% 
  # mutate(rdDS1 = rDS1 - rCtS1,
  #        rdFdS1 = rFdS1 - rCtS1,
  #        rdFzS1 = rFzS1 - rCtS1,
  #        rdHS1 = rHS1 - rCtS1) %>% 
  pivot_longer(c(rdDS1, rdFdS1, rdFzS1, rdHS1), names_to = "treatment", values_to = "rdS1") %>% 
  ggplot(aes(x = rdS1, fill = treatment))+
  geom_density(alpha = 0.2) +
  scale_x_continuous("Selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)


ggsave("./plots/dSRCS11e5.pdf", scale = 0.5)


#look at the distributions - makes it clear that there are very different numbers shared between pairs of treatments - heat with freezing and flooding
library(GGally)
dfI %>% 
  select(starts_with("rd",ignore.case = FALSE) & ends_with("S1")) %>% 
  # ggpairs(columnLabels = c("Drought", "Flooding", "Freezing","Heat"), upper = list(continuous = wrap(ggally_cor, use = "pairwise.complete.obs" , method = "kendall", stars = FALSE)))
  ggpairs(columnLabels = c("Drought", "Flooding", "Freezing","Heat"), upper = "blank", lower = list(continuous = wrap(ggally_points, alpha = 0.2 )))
  
ggsave("./plots/S1pairs1e5.pdf")

dfI %>% 
  select(starts_with("rd",ignore.case = FALSE) & (ends_with("S1")|ends_with("S4"))) %>% 
  ggpairs(columnLabels = c("Drought S1", "Flooding S1", "Freezing S1","Heat S1", "Drought S4", "Flooding S4", "Freezing S4","Heat S4"), upper = "blank", lower = list(continuous = wrap(ggally_points, alpha = 0.2 )))

ggsave("./plots/S1S4pairs1e5.pdf")

#looks like antagonistic pleiotropy between S1 and S4 phases, particularly in freezing treatment, where 56 samples with good values
dfI %>% 
  filter(!is.na(rdFzS1), !is.na(rdFzS4), is.finite(rdFzS1), is.finite(rdFzS4)) %>% 
  ggplot(aes(x = rdFzS1, y = rdFzS4, colour = whichInitPoly))+
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
scale_x_continuous("Selection during freezing") +
  scale_y_continuous("Selection during recovery") +
  theme_classic()

ggsave("./plots/freezingS1S4.pdf",  scale = 0.5)



#might also want to summarise the number of initial SNPs polymorphic in each treatment

# di <- dfI %>% 
  dfI %>%
  # pivot_longer(c(rdDS1, rdFdS1, rdFzS1, rdHS1),  names_to = "treatment", values_to = "rdS1") %>%
  pivot_longer(c(starts_with("rd",ignore.case = FALSE)),  names_to = c("treatment", ".value"), names_sep = "S") %>%
  group_by(treatment) %>% 
  summarise(nS1 = sum(!is.na(`1`) & !is.infinite(`1`)), nS4 = sum(!is.na(`4`) & !is.infinite(`4`))) %>% 
  # summarise(nS1 = sum(!is.na(rdS1) & !is.infinite(rdS1)))
  pivot_longer(c(nS1, nS4), names_to = "Sampling time", values_to = "N SNPs") %>% 
  mutate(treatment = recode(treatment, rdD = "Drought", rdFd = "Flooding", rdFz = "Freezing", rdH = "Heat")) %>% 
  ggplot(aes(x = `Sampling time`, y = `N SNPs`, fill = treatment)) +
  geom_col()+
  scale_y_continuous("Number of SNPs with fitness estimate") +
  scale_x_discrete("Sampling time") +
  theme_classic() +
  scale_fill_brewer(palette = "Set1")

ggsave("./plots/NSNPr1e5.pdf")

dfI %>% 
  mutate(r1dDS1 = r1DS1 - r1CtS1,
         r1dFdS1 = r1FdS1 - r1CtS1,
         r1dFzS1 = r1FzS1 - r1CtS1,
         r1dHS1 = r1HS1 - r1CtS1) %>% 
  pivot_longer(c(r1dDS1, r1dFdS1, r1dFzS1, r1dHS1), names_to = "treatment", values_to = "r1dS1") %>% 
  ggplot(aes(x = r1dS1, fill = treatment))+
  geom_density(alpha = 0.2) +
  scale_x_continuous("Selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)

#Want to blast sequences to ID. Options:
# BlastSequences::annotate
# rBlast (github)
# gschofl/blastr
#blast::BoSSA (not obvious it's in there)

#can't persuade bioconductor to use http rather than https repositories (except for BioManager), so failing to install on linux machine
# looks like it didn't work, using this:
BiocManager::repositories()
#but actually, this seems to add a usable option to the RStudio gui installer:
options(repos = c(BioCsoft = "http://bioconductor.org/packages/3.10/bioc"))
#but need to switch back:
options(repos = c(CRAN = "http://cloud.r-project.org/"))
#still has problems with certificates
library(annotate)
blastSequences("GGCCTTCATTTACCCAAAATG")

#could resort to website (blastX): 
#look at extreme sequences that are positively selected in freezing and negatively in recovery, or vice versa.
library(seqinr)
difz<- dfI %>% 
  filter(!is.na(rdFzS1), !is.na(rdFzS4), is.finite(rdFzS1), is.finite(rdFzS4),
         rdFzS1>2, rdFzS4 < -0.5) %>% 
  select(sequence_higher, polymorphism_higher, mut)

  write.fasta(sequences = as.list(unlist(difz$sequence_higher)), names = paste(unlist(di$mut), unlist(difz$polymorphism_higher), sep ="_"), as.string = TRUE, file.out = "FreezeNotRecov.fa" )

    direc<- dfI %>% 
    filter(!is.na(rdFzS1), !is.na(rdFzS4), is.finite(rdFzS1), is.finite(rdFzS4),
           rdFzS1 < -2, rdFzS4 > 0.5) %>% 
    select(sequence_higher, polymorphism_higher, mut)
    
    write.fasta(sequences = as.list(unlist(direc$sequence_higher)), names = paste(unlist(direc$mut), unlist(direc$polymorphism_higher), sep ="_"), as.string = TRUE, file.out = "RecovNotFreeze.fa" )
    
#also want to look at the organisms
Nseq <- dfI$data[[1]]$Nseq

#read counts do vary somewhat among samples
dfI$data[[1]] %>% 
  ggplot(aes(x = Nseq)) + 
  geom_histogram(bins =50) +
  scale_x_continuous("Read count") +
  scale_y_continuous("Number of samples") +
  theme_classic()
ggsave("./plots/readCountsDE.pdf",  scale = 0.5)


rOf <- function(dat, whichSoil){
  xI <- dat %>% 
    filter(Replicate_site == whichSoil) %>% 
    mutate(total = higher+lower,
           Nother = Nseq - total) 
  xICt <- xI %>% filter(Treatment == "Ct"| Sampling_time == "Initial")
  xID <- xI %>% filter(Treatment == "D"| Sampling_time == "Initial")
  xIFd <- xI %>% filter(Treatment == "Fd"| Sampling_time == "Initial")
  xIFz <- xI %>% filter(Treatment == "Fz"| Sampling_time == "Initial")
  xIH <- xI %>% filter(Treatment == "H"| Sampling_time == "Initial")
  rOCtS1 <- log2(xICt$total[1]/xICt$total[3]) - log2(xICt$Nother[1]/xICt$Nother[3])
  rOCtS4 <- (log2(xICt$total[2]/xICt$total[1]) - log2(xICt$Nother[2]/xICt$Nother[1]))/4
  rODS1 <- log2(xID$total[1]/xID$total[3]) - log2(xID$Nother[1]/xID$Nother[3])
  rODS4 <- (log2(xID$total[2]/xID$total[1]) - log2(xID$Nother[2]/xID$Nother[1]))/4
  rOFdS1 <- log2(xIFd$total[1]/xIFd$total[3]) - log2(xIFd$Nother[1]/xIFd$Nother[3])
  rOFdS4 <- (log2(xIFd$total[2]/xIFd$total[1]) - log2(xIFd$Nother[2]/xIFd$Nother[1]))/4
  rOFzS1 <- log2(xIFz$total[1]/xIFz$total[3]) - log2(xIFz$Nother[1]/xIFz$Nother[3])
  rOFzS4 <- (log2(xIFz$total[2]/xIFz$total[1]) - log2(xIFz$Nother[2]/xIFz$Nother[1]))/4
  rOHS1 <- log2(xIH$total[1]/xIH$total[3]) - log2(xIH$Nother[1]/xIH$Nother[3])
  rOHS4 <- (log2(xIH$total[2]/xIH$total[1]) - log2(xIH$Nother[2]/xIH$Nother[1]))/4
  return(tibble(rOCtS1, rOCtS4, rODS1, rODS4, rOFdS1, rOFdS4, rOFzS1, rOFzS4, rOHS1, rOHS4))
}

dfI <- dfI %>% 
  mutate(rO = map2(.x = data, .y = whichInitPoly, .f =rOf)) %>% 
  unnest(rO) %>% 
    mutate(
      rOdDS1 = rODS1 - rOCtS1,
      rOdFdS1 = rOFdS1 - rOCtS1,
      rOdFzS1 = rOFzS1 - rOCtS1,
      rOdHS1 = rOHS1 - rOCtS1,
      rOdDS4 = rODS4 - rOCtS4,
      rOdFdS4 = rOFdS4 - rOCtS4,
      rOdFzS4 = rOFzS4 - rOCtS4,
      rOdHS4 = rOHS4 - rOCtS4)
 
saveRDS(dfI, file = "dfI1e5.rds")
setwd("~/Dropbox (The University of Manchester)/knightmareDocs/soil/OcÃ©ane/sequence/metagenomics/DiscoSNP/DEall/")

dfI <- readRDS(file = "dfI1e5.rds")

#look at organism frequency distributions
dfI %>% 
  # pivot_longer(c(starts_with("rOd",ignore.case = FALSE)),  names_to = c("treatment", ".value"), names_sep = "S") %>%
  pivot_longer(c(starts_with("rOd",ignore.case = FALSE)),  names_to = "time_trt", values_to = "rOd") %>%
  mutate(samplingTime = str_trunc(time_trt, 2, ellipsis = "", side = "left"),
         treatment = str_remove_all(time_trt, "[rOdS14]")) %>% 
  ggplot(aes(x = `rOd`, fill = treatment))+
  geom_density(alpha = 0.2) +
  facet_wrap(~samplingTime) +
  scale_x_continuous("Organismal selection rate coefficient (doublings per week)") +
  scale_y_continuous("density") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_fill_brewer(palette = "Set1")
ggsave("./plots/OrganismSelection.pdf",  scale = 0.5)

#want to plot organism against snp fitness should be able to split names using pivot-longer, but don't have consistent length or consistent letter beyond 'd' which appears in Fd and messes stuff up
# di <- dfI %>% 
dfI %>%   
pivot_longer(c(starts_with("rOd",ignore.case = FALSE),starts_with("rd",ignore.case = FALSE)), names_to = "thing", values_to = "selectionRate")%>%
  mutate(samplingTime = str_trunc(thing, 2, ellipsis = "", side = "left"),
         treatment = str_remove_all(thing, "[rOdS14]"),
         what = ifelse(str_detect(thing, "O"), "Organism", "SNP")) %>% 
  pivot_wider(id_cols = c(mut, treatment, samplingTime), names_from = what, values_from = selectionRate) %>% 
  ggplot(aes(x = Organism, y = SNP, colour = treatment))+
  geom_point(alpha = 0.2) +
  facet_grid(treatment~samplingTime) +
  scale_x_continuous("Organismal selection rate coefficient (doublings per week)") +
  scale_y_continuous("SNP selection rate coefficient (doublings per week)") +
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_colour_brewer(palette = "Set1")
ggsave("./plots/OrganismVsSNP.pdf",  scale = 0.8)

 
#bring in models:
setwd("~/Dropbox (The University of Manchester)/knightmareDocs/soil/OcÃ©ane/sequence/metagenomics/DiscoSNP/DEall/")

dfmods <- readRDS("DEallstk2.rds")
dfI <- readRDS(file = "dfI1e5.rds") %>% 
  left_join(dfmods) %>% 
  mutate(fdr = p.adjust(p, method = "fdr"))

dfI %>% 
  ggplot(aes(x = p)) +
  geom_histogram()+
  scale_x_log10("P value")+
  geom_vline(xintercept = c(0.05, 0.01, 0.001), linetype = c(4, 2, 3))+
  scale_y_continuous("SNP count") +
  theme_classic()
  
ggsave("./plots/Pval.pdf",  scale = 0.5)


dfI %>% 
  ggplot(aes(x = fdr)) +
  geom_histogram()+
  scale_x_log10("P value")+
  geom_vline(xintercept = c(0.05, 0.01, 0.001), linetype = c(4, 2, 3))+
  scale_y_continuous("SNP count") +
  theme_classic()

ggsave("./plots/PvalFDR.pdf",  scale = 0.5)

dfI %>% 
  ggplot(aes(x = dfnull, y = nulldev)) +
  geom_point(alpha = 0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_x_continuous("null degrees of freedom") +
  scale_y_continuous("null deviance") +
  theme_classic()
ggsave("./plots/nullDeviance.pdf",  scale = 0.5)


saveRDS(dfI, file = "dfI1e5.rds")
