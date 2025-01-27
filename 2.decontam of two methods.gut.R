rm(list = ls())
####get path#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#path <- getwd()
#path
library(phyloseq)
library(ggplot2)
library(decontam)
library(openxlsx)
library(readxl)
library(dplyr, warn.conflicts = F)
####Built-in data that can be used to imitate data types####
#ps1 <- readRDS(system.file("extdata", "MUClite.rds", package="decontam"))
#sample_data(ps1)
#otu_table(ps1, taxa_are_rows = F)
####DNA.concen.method####
asv.sample <- read.delim('DNA.concen.method.sample.txt', row.names = 1)
DNA.concen <- read.delim('DNA.concen.group.txt', row.names = 1)
tax <- read.delim('DNA.concen.method.sample.txt', row.names = 1)
tax <- tax[, c('Taxonomy','Phylum')]
tax$Phylum <- NULL
asv.sample <- asv.sample[, 1:(which(names(asv.sample) == "Taxonomy") - 1), drop = FALSE]
#####phyloseq####
otu_table <- as.data.frame(t(asv.sample))
otu_table <- otu_table(otu_table, taxa_are_rows = F)
sample_data <- sample_data(DNA.concen)
ps <- phyloseq(otu_table,sample_data)
#head(sample_data(ps))
#head(otu_table(ps))
#####threshold = 0.1####
contamdf.freq.0.1 <- isContaminant(ps, 
                               method = "frequency", 
                               conc = "quant_reading", 
                               threshold = 0.1)
contamdf.freq.0.1$freq0.1 <- ifelse(contamdf.freq.0.1$contaminant == FALSE, 0, 1)
#write.table(contamdf.freq.0.1, "contamdf_freq0.1.txt", row.names = F, sep = "\t", quote = FALSE)
#head(contamdf.freq.0.1)
#head(contamdf.freq.0.1$freq)
#head(which(contamdf.freq.0.1$contaminant))
which(contamdf.freq.0.1$contaminant)
table(contamdf.freq.0.1$contaminant)
contamdf.freq.0.1$ASV0.1 <- rownames(contamdf.freq.0.1)
ASV0.1 <- 'ASV0.1'
contamdf.freq.0.1 <- contamdf.freq.0.1[,c(ASV0.1, setdiff(names(contamdf.freq.0.1), ASV0.1))]
rownames(contamdf.freq.0.1) <- NULL
contamdf.freq.0.1$ASV <- rownames(contamdf.freq.0.1)
ASV.name <- 'ASV'
contamdf.freq.0.1 <- contamdf.freq.0.1[,c(ASV.name, setdiff(names(contamdf.freq.0.1), ASV.name))]
ASV <- contamdf.freq.0.1[contamdf.freq.0.1$freq0.1 == 1, "ASV"]
freq.0.1 <- as.data.frame(ASV)
freq.0.1$ASV0.1 <- contamdf.freq.0.1[contamdf.freq.0.1$freq0.1 == 1, "ASV0.1"]
head(freq.0.1,10)
######plot of some contam####
taxa_to_plot1 <- freq.0.1$ASV0.1[1:24]
taxa_to_plot1
p <- plot_frequency(ps, taxa_to_plot1, 
               conc = "quant_reading") +
  xlab("DNA Concentration") 
p
ggsave("frequence.DNA.conc.0.1.pdf", width = 12, height = 10)
ggsave('frequence.DNA.conc.0.1.png', p, width = 12, height = 10)
#p
######plot of some un.contam####
taxa_to_plot1 <- taxa_names(ps)[c(1,2,3,4,5,6,7,8,9)]
p.p <- plot_frequency(ps, taxa_to_plot1, 
                    conc = "quant_reading") +
  xlab("DNA Concentration") 
p.p
ggsave("frequence.DNA.conc.0.1.uncon.pdf", width = 8, height = 6)
ggsave('frequence.DNA.conc.0.1.uncon.png', p.p, width = 8, height = 6)
#####threshold = 0.5####
contamdf.freq.0.5 <- isContaminant(ps, 
                                   method = "frequency", 
                                   conc = "quant_reading", 
                                   threshold = 0.5)
contamdf.freq.0.5$freq0.5 <- ifelse(contamdf.freq.0.5$contaminant == FALSE, 0, 1)
#write.table(contamdf.freq.0.5, "contamdf_freq0.5.txt", row.names = T, sep = "\t", quote = FALSE)
contamdf.freq.0.1$freq0.5 <- contamdf.freq.0.5$freq0.5
#write.table(contamdf.freq.0.1, "contamdf_freq0.1.txt", row.names = T, sep = "\t", quote = FALSE)
contamdf.freq.0.5$ASV0.5 <- rownames(contamdf.freq.0.5)
ASV0.5 <- 'ASV0.5'
contamdf.freq.0.5 <- contamdf.freq.0.5[,c(ASV0.5, setdiff(names(contamdf.freq.0.5), ASV0.5))]
rownames(contamdf.freq.0.5) <- NULL
contamdf.freq.0.5$ASV5 <- rownames(contamdf.freq.0.5)
ASV.name.0.5 <- 'ASV5'
contamdf.freq.0.5 <- contamdf.freq.0.5[,c(ASV.name.0.5, setdiff(names(contamdf.freq.0.5), ASV.name.0.5))]
ASV5 <- contamdf.freq.0.5[contamdf.freq.0.5$freq0.5 == 1, "ASV5"]
freq.0.5 <- as.data.frame(ASV5)
freq.0.5$ASV0.5 <- contamdf.freq.0.5[contamdf.freq.0.5$freq0.5 == 1, "ASV0.5"]
#freq.0.5
#head(contamdf.freq.0.5)
#head(contamdf.freq.0.5$freq)
#head(which(contamdf.freq.0.5$contaminant))
#which(contamdf.freq.0.5$contaminant)
table(contamdf.freq.0.5$contaminant)
table(contamdf.freq.0.1$contaminant)
######plot of some contam####
taxa_to_plot2 <- freq.0.5$ASV0.5[1:16]
taxa_to_plot2
#taxa_to_plot2 <- taxa_names(ps)[c(13,21,27,29,30,31,35,38,43)]
p1.0.5 <- plot_frequency(ps, taxa_to_plot2, 
                     conc = "quant_reading") +
                    xlab("DNA Concentration") 
p1.0.5
ggsave("frequence.DNA.conc.0.5.pdf", width = 10, height = 8)
ggsave('frequence.DNA.conc.0.5.png', p1.0.5, width = 10, height = 8)
######plot of some un.contam####
taxa_to_plot1 <- taxa_names(ps)[c(1,5,10,11,12,14,15,16,17)]
pp <- plot_frequency(ps, taxa_to_plot1, 
                      conc = "quant_reading") +
  xlab("DNA Concentration") 
pp
ggsave("frequence.DNA.conc.0.5.uncon.pdf", width = 8, height = 6)
ggsave('frequence.DNA.conc.0.5.uncon.png', pp, width = 8, height = 6)
####NC.method####
asv_p <- read.delim('NC.method.sample.txt', row.names = 1)
asv_p <- asv_p[, 1:(which(names(asv_p) == "Taxonomy") - 1), drop = FALSE]
tax_p <- read.delim('NC.method.sample.txt', row.names = 1)
tax_p <- tax_p[, c('Taxonomy','Phylum')]
tax$Phylum <- NULL
map <- read.delim("NC.method.sample.group.txt", row.names = 1)
#relative
conlum_sum <- colSums(asv_p[,1:NCOL(asv_p)])
#conlum_sum
asv_p.re <- asv_p[, 1:NCOL(asv_p)] / rep(conlum_sum, each = nrow(asv_p))
ASV.re <- otu_table(asv_p.re, taxa_are_rows = TRUE)
tax_p.re <- cbind(asv_p.re, Taxonomy = tax_p$Taxonomy)
tax_p.re <- as.matrix(tax_p.re)
tax_p <- as.matrix(tax_p)
TAX.re <- tax_table(tax_p.re)
#####phyloseq####
ASV <- otu_table(asv_p, taxa_are_rows = TRUE)
TAX <- tax_table(tax_p)
map <- sample_data(map)
#head(map)
ps2 <- phyloseq(ASV, TAX, map)
ps2
sample_data(ps2)$is.neg <- sample_data(ps2)$Group == "control"
contamdf.prev0.1 <- isContaminant(ps2, 
                               method = "prevalence", 
                               neg = "is.neg",
                               threshold = 0.1)
contamdf.prev0.5 <- isContaminant(ps2, 
                               method = "prevalence", 
                               neg = "is.neg",
                               threshold = 0.5)
physeq.contam0.1 <- prune_taxa(contamdf.prev0.1$contaminant, ps2)
physeq.contam0.5 <- prune_taxa(contamdf.prev0.5$contaminant, ps2)
table(contamdf.prev0.1$contaminant) #show results
table(contamdf.prev0.5$contaminant)
contamdf.prev0.1$abprev0.1 <- ifelse(contamdf.prev0.1$contaminant == FALSE, 0, 1)
contamdf.prev0.5$abprev0.5 <- ifelse(contamdf.prev0.5$contaminant == FALSE, 0, 1)
contamdf.freq.0.1$abprev0.1 <- contamdf.prev0.1$abprev0.1
contamdf.freq.0.1$abprev0.5 <- contamdf.prev0.5$abprev0.5
####relative abundance####
ps3 <- phyloseq(ASV.re, TAX.re, map)
ps3
sample_data(ps3)$is.neg <- sample_data(ps3)$Group == "control"
contamdf.prev0.1.re <- isContaminant(ps3, 
                                  method = "prevalence", 
                                  neg = "is.neg",
                                  threshold = 0.1)
contamdf.prev0.5.re <- isContaminant(ps3, 
                                  method = "prevalence", 
                                  neg = "is.neg",
                                  threshold = 0.5)
physeq.contam0.1.re <- prune_taxa(contamdf.prev0.1.re$contaminant, ps3)
physeq.contam0.5.re <- prune_taxa(contamdf.prev0.5.re$contaminant, ps3)
table(contamdf.prev0.1.re$contaminant)
table(contamdf.prev0.5.re$contaminant)
contamdf.prev0.1.re$reprev0.1 <- ifelse(contamdf.prev0.1.re$contaminant == FALSE, 0, 1)
contamdf.prev0.5.re$reprev0.5 <- ifelse(contamdf.prev0.5.re$contaminant == FALSE, 0, 1)
contamdf.freq.0.1$reprev0.1 <- contamdf.prev0.1.re$reprev0.1
contamdf.freq.0.1$reprev0.5 <- contamdf.prev0.5.re$reprev0.5
Sum <- rowSums(contamdf.freq.0.1[, (ncol(contamdf.freq.0.1) - 5):ncol(contamdf.freq.0.1)])
contamdf.freq.0.1$Sum <- Sum
contam.result <- contamdf.freq.0.1[, c(2, (ncol(contamdf.freq.0.1) - 6):ncol(contamdf.freq.0.1))]
colnames(contam.result)[1] <- "ASV_num"
contam.result.filter <- contamdf.freq.0.1[, c(1,2,9:15)]
write.table(contam.result, "contam.result.txt", row.names = F, sep = "\t", quote = FALSE)
write.table(contam.result.filter, "contamdf.result.No.txt", row.names = F, sep = "\t", quote = FALSE)
write.table(contamdf.freq.0.5, "contamdf.freq.0.5.txt", row.names = F, sep = "\t", quote = FALSE)
write.table(contamdf.prev0.1, "contamdf.abprev0.1.txt", row.names = F, sep = "\t", quote = FALSE)
write.table(contamdf.prev0.5, "contamdf.abprev0.5.txt", row.names = F, sep = "\t", quote = FALSE)
write.table(contamdf.prev0.1.re, "contamdf.reprev0.1.txt", row.names = F, sep = "\t", quote = FALSE)
write.table(contamdf.prev0.5.re, "contamdf.reprev0.5.txt", row.names = F, sep = "\t", quote = FALSE)
#####plot of threshold=0.1####
physeq_p.pa <- transform_sample_counts(ps2,function(abund) 100*(abund > 0))
physeq_p.pa.neg <- prune_samples(sample_data(physeq_p.pa)$Group == "control", physeq_p.pa)
physeq_p.pa.pos <- prune_samples(sample_data(physeq_p.pa)$Group == "real", physeq_p.pa)
df.pa0.1 <- data.frame(pa.pos = taxa_sums(physeq_p.pa.pos), pa.neg = taxa_sums(physeq_p.pa.neg), 
                    contaminant = contamdf.prev0.1$contaminant)
p2 <- ggplot(data = df.pa0.1, aes(x = pa.neg, y = pa.pos, color = contaminant)) + 
       geom_point() + 
       xlab("Prevalence (Negative Controls)  0.1") + 
       ylab("Prevalence (True Samples)")
p2
ggsave("Prevalence0.1.pdf", width = 6, height = 5)
ggsave('Prevalence0.1.png', p2, width = 6, height = 5)
#####plot of threshold=0.5####
df.pa0.5 <- data.frame(pa.pos = taxa_sums(physeq_p.pa.pos), pa.neg = taxa_sums(physeq_p.pa.neg), 
                       contaminant = contamdf.prev0.5$contaminant)
p3 <- ggplot(data = df.pa0.5, aes(x = pa.neg, y = pa.pos, color = contaminant)) + 
  geom_point() + 
  xlab("Prevalence (Negative Controls) 0.5") + 
  ylab("Prevalence (True Samples)")
p3
ggsave("Prevalence0.5.pdf", width = 6, height = 5)
ggsave('Prevalence0.5.png', p3, width = 6, height = 5)
####out put####
rownames(otu_table(physeq.contam0.1))  #show contaminant ASV
rownames(otu_table(physeq.contam0.5))  #show contaminant ASV
physeq.noncontam0.1 <- prune_taxa(!contamdf.prev0.1$contaminant, ps2)
physeq.noncontam0.5 <- prune_taxa(!contamdf.prev0.5$contaminant, ps2)
otu_final0.1 = otu_table(physeq.noncontam0.1)
otu_final0.5 = otu_table(physeq.noncontam0.5)
otu_final0.1 = as.data.frame(otu_final0.1)
otu_final0.5 = as.data.frame(otu_final0.5)
otu_colname0.1 = colnames(otu_final0.1)
otu_colname0.5 = colnames(otu_final0.5)
otu_final0.1 = cbind(rownames(otu_final0.1),otu_final0.1)
otu_final0.5 = cbind(rownames(otu_final0.5),otu_final0.5)
colnames(otu_final0.1) = c("ASV_num", otu_colname0.1)
colnames(otu_final0.5) = c("ASV_num", otu_colname0.5)
tax_final0.1 = tax_table(physeq.noncontam0.1)
tax_final0.5 = tax_table(physeq.noncontam0.5)
#ASV taxa table OTU_ID
tax_final0.1 = as.data.frame(tax_final0.1)
tax_final0.5 = as.data.frame(tax_final0.5)
tax_colname0.1 = colnames(tax_final0.1)
tax_colname0.5 = colnames(tax_final0.5)
tax_final0.1 = cbind(rownames(tax_final0.1),tax_final0.1)
tax_final0.5 = cbind(rownames(tax_final0.5),tax_final0.5)
#colnames(tax_final) = c("ASV_num", tax_colname)
write.table(otu_final0.1, "decon_asv0.1.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(tax_final0.1, "decon_asv0.1_tax.txt", sep = "\t", quote = F, col.names = T, row.names = F) #decontam.asv.table
write.table(otu_final0.5, "decon_asv0.5.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(tax_final0.5, "decon_asv0.5_tax.txt", sep = "\t", quote = F, col.names = T, row.names = F) #decontam.asv.table
#Notice####
##Decontam results just provide a reference for the ASV table.
##The author should use these results carefully.
##I recommend checking the results for each ASV feature, and I do the same.