library(tidyverse)

###Clean up ensembl gene data

#Obtain data from ensembl BioMart (Ensembl Genes v108)
bed <- read.delim("Ensembl gene data.txt", header = T)
colnames(bed) <- c("GeneID", "Chromosome", "Start", "End", "GeneType", "GeneName", "Strand")

#Only keep chromosomes 1-19, X, Y, and mitochondria
#Also, clean up chromosome names
chrom <- c(1:19, "X", "Y", "MT")
bed <- filter(bed, Chromosome %in% chrom) %>% 
  mutate(Strand = ifelse(Strand == "1", "+", "-"))

#Official bed format has the start 0-indexed, meaning that the first base of the chromosome is numbered 0
bed <- mutate(bed, Start = Start - 1)

#Rearrange columns to bed format and collapse info into one column 
bed <- mutate(bed, info = paste(GeneID, GeneName, GeneType, Strand, sep = ";"))
bed <- select(bed, Chromosome, Start, End, info)

#Now load in tRNA data
tRNA <- read.table("mm10-tRNAs.bed")
tRNA <- mutate(tRNA,
               V1 = sub("chr", "", V1),
               V13 = paste0(V4, ";", V6)) %>% 
  filter(V1 %in% chrom) %>%
  select(V1, V2, V3, V13) 
colnames(tRNA) <- colnames(bed)
bed <- rbind(bed, tRNA)

bed <- mutate(bed, Chromosome = factor(Chromosome, levels = chrom)) %>%
  arrange(Chromosome, Start)

fbed <- bed[grep("\\;\\+", bed$info),]
rbed <- bed[grep("\\;\\-", bed$info),]

write.table(bed, "genes.bed", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(fbed, "fgenes.bed", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(rbed, "rgenes.bed", sep = "\t", row.names = F, col.names = F, quote = F)

