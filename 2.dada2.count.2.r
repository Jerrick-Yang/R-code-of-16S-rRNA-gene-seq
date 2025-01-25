#rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#path <- getwd()
#path
library(dada2)
library(Biostrings)
library(openxlsx)
library(dplyr)
#读取特征序列
# 设置输入和输出文件路径
# 定义前缀
prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
silva_ref <- readDNAStringSet("silva_nr99_v138.2_toSpecies_trainset.fa")
#count.2####
asv_fasta <- "feature.sequence.2.fasta"  # 特征序列文件路径
seqs.2 <- readDNAStringSet("feature.sequence.2.fasta")
#seqs.2
#silva_ref.s <- readDNAStringSet("silva_v138.2_assignSpecies.fa")
#silva_ref
taxa.count.2 <- assignTaxonomy(seqs.2, 
                       refFasta = "silva_nr99_v138.2_toSpecies_trainset.fa", 
                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                       minBoot = 50,
                       tryRC = TRUE,
                       outputBootstraps = TRUE, 
                       multithread = TRUE)
taxa.2 <- as.data.frame(taxa.count.2$tax)
taxa.3 <- taxa.2
taxa.3$Species <- NULL
#taxa.2 <- as.data.frame(taxa.2)
taxa.2$species.2 <- ifelse(is.na(taxa.2$Species), NA, 
                             paste0(taxa.2$Genus, "_", taxa.2$Species))
taxa.2$Species <- taxa.2$species.2
taxa.2$species.2 <- NULL
# 遍历每一列，添加对应的前缀
for (i in 1:ncol(taxa.2)) {
  taxa.2[[i]] <- ifelse(is.na(taxa.2[[i]]), NA, 
                          paste0(prefixes[i], taxa.2[[i]]))
}
taxa.count.2.boot <- as.data.frame(taxa.count.2$boot)
#####addspecies#####
taxa.2.species <- addSpecies(taxa.3,
                            refFasta = 'silva_v138.2_assignSpecies_modified.fa',
                            allowMultiple = TRUE,
                            tryRC = TRUE,
                            n = 2000,
                            verbose = TRUE)
taxa.2.species <- as.data.frame(taxa.2.species)
# 遍历每一列，添加对应的前缀
for (i in 1:ncol(taxa.2.species)) {
  taxa.2.species[[i]] <- ifelse(is.na(taxa.2.species[[i]]), NA, 
                        paste0(prefixes[i], taxa.2.species[[i]]))
}
# 提取特征序列的名称
seq_names.2 <- names(seqs.2)
# 将序列名称添加到 taxa.2 数据框中
taxa.2 <- cbind(ASV_num = seq_names.2, taxa.2)
taxa.2.species <-  cbind(ASV_num = seq_names.2, taxa.2.species)
taxa.count.2.boot <-  cbind(ASV_num = seq_names.2, taxa.count.2.boot)
#####Taxonomy####
taxa.2 <- taxa.2 %>%
  mutate(Taxonomy = case_when(
    is.na(Species) & !is.na(Genus) ~ paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = ";"),
    is.na(Genus) & !is.na(Family) ~ paste(Kingdom, Phylum, Class, Order, Family, sep = ";"),
    is.na(Family) & !is.na(Order) ~ paste(Kingdom, Phylum, Class, Order, sep = ";"),
    is.na(Order) & !is.na(Class) ~ paste(Kingdom, Phylum, Class, sep = ";"),
    is.na(Class) & !is.na(Phylum) ~ paste(Kingdom, Phylum, sep = ";"),
    is.na(Phylum) & !is.na(Kingdom) ~ Kingdom,
    TRUE ~ paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";")
  )) %>%
  relocate(Taxonomy, .before = Kingdom)
taxa.2.asv <- read.delim('ASV_table_2.txt')
taxa.2.asv <- taxa.2.asv %>% rename(Taxonomy.2 = Taxonomy)
taxa.2.merged <- taxa.2.asv %>% left_join(taxa.2, by = "ASV_num")
#taxa.2.species
taxa.2.species <- taxa.2.species %>%
  mutate(Taxonomy = case_when(
    is.na(Species) & !is.na(Genus) ~ paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = ";"),
    is.na(Genus) & !is.na(Family) ~ paste(Kingdom, Phylum, Class, Order, Family, sep = ";"),
    is.na(Family) & !is.na(Order) ~ paste(Kingdom, Phylum, Class, Order, sep = ";"),
    is.na(Order) & !is.na(Class) ~ paste(Kingdom, Phylum, Class, sep = ";"),
    is.na(Class) & !is.na(Phylum) ~ paste(Kingdom, Phylum, sep = ";"),
    is.na(Phylum) & !is.na(Kingdom) ~ Kingdom,
    TRUE ~ paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";"))) %>%
  relocate(Taxonomy, .before = Kingdom)
taxa.2.species.merged <- taxa.2.asv %>% left_join(taxa.2.species, by = "ASV_num")
# 将 Taxonomy 赋值给 Taxonomy.3，并将 Taxonomy.3 插入到 Taxonomy.2 后面
taxa.2.merged <- taxa.2.merged %>%
                 mutate(Taxonomy.muti = taxa.2.species.merged$Taxonomy,
                        Species.muti = taxa.2.species.merged$Species,
                        Diff = case_when(
                          is.na(Species.muti) & !is.na(Species) ~ 1,
                          !is.na(Species.muti) & is.na(Species) ~ 1,
                          Species.muti != Species ~ 1,
                          is.na(Species.muti) & is.na(Species) ~ 0,
                          TRUE ~ 0)) %>%
                 relocate(Taxonomy.muti, .after = Taxonomy.2)
#####excel####
# 添加每个数据框到工作簿中，每个数据框一个工作表
wb <- createWorkbook()
#sheet
addWorksheet(wb, "taxa.2")
writeData(wb, sheet = "taxa.2", taxa.2.merged, rowNames = FALSE)
#sheet
addWorksheet(wb, "taxa.2.species")
writeData(wb, sheet = "taxa.2.species", taxa.2.species.merged, rowNames = FALSE)
#sheet
addWorksheet(wb, "taxa.count.2.boot")
writeData(wb, sheet = "taxa.count.2.boot", taxa.count.2.boot, rowNames = FALSE)
#xlsx
saveWorkbook(wb, "ASV.count.2.tax.xlsx", overwrite = TRUE)
####count.1####
asv_fasta.1 <- "feature.sequence.1.fasta"  # 特征序列文件路径
seqs.1 <- readDNAStringSet("feature.sequence.1.fasta")
taxa.count.1 <- assignTaxonomy(seqs.1, 
                               refFasta = "silva_nr99_v138.2_toSpecies_trainset.fa", 
                               taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                               minBoot = 50,
                               tryRC = TRUE,
                               outputBootstraps = TRUE, 
                               multithread = TRUE)
taxa.1 <- as.data.frame(taxa.count.1$tax)
taxa.0 <- taxa.1
taxa.0$Species <- NULL
#taxa.2 <- as.data.frame(taxa.2)
taxa.1$species.2 <- ifelse(is.na(taxa.1$Species), NA, 
                           paste0(taxa.1$Genus, "_", taxa.1$Species))
taxa.1$Species <- taxa.1$species.2
taxa.1$species.2 <- NULL
# 遍历每一列，添加对应的前缀
for (i in 1:ncol(taxa.1)) {
  taxa.1[[i]] <- ifelse(is.na(taxa.1[[i]]), NA, 
                        paste0(prefixes[i], taxa.1[[i]]))
}
taxa.count.1.boot <- as.data.frame(taxa.count.1$boot)
#####addspecies#####
taxa.1.species <- addSpecies(taxa.0,
                             refFasta = 'silva_v138.2_assignSpecies_modified.fa',
                             allowMultiple = TRUE,
                             tryRC = TRUE,
                             n = 2000,
                             verbose = TRUE)
taxa.1.species <- as.data.frame(taxa.1.species)
# 遍历每一列，添加对应的前缀
for (i in 1:ncol(taxa.1.species)) {
  taxa.1.species[[i]] <- ifelse(is.na(taxa.1.species[[i]]), NA, 
                                paste0(prefixes[i], taxa.1.species[[i]]))
}
# 提取特征序列的名称
seq_names.1 <- names(seqs.1)
# 将序列名称添加到 taxa.1数据框中
taxa.1 <- cbind(ASV_num = seq_names.1, taxa.1)
taxa.1.species <-  cbind(ASV_num = seq_names.1, taxa.1.species)
taxa.count.1.boot <-  cbind(ASV_num = seq_names.1, taxa.count.1.boot)
#####Taxonomy####
taxa.1 <- taxa.1 %>%
  mutate(Taxonomy = case_when(
    is.na(Species) & !is.na(Genus) ~ paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = ";"),
    is.na(Genus) & !is.na(Family) ~ paste(Kingdom, Phylum, Class, Order, Family, sep = ";"),
    is.na(Family) & !is.na(Order) ~ paste(Kingdom, Phylum, Class, Order, sep = ";"),
    is.na(Order) & !is.na(Class) ~ paste(Kingdom, Phylum, Class, sep = ";"),
    is.na(Class) & !is.na(Phylum) ~ paste(Kingdom, Phylum, sep = ";"),
    is.na(Phylum) & !is.na(Kingdom) ~ Kingdom,
    TRUE ~ paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";")
  )) %>%
  relocate(Taxonomy, .before = Kingdom)
taxa.1.asv <- read.delim('ASV_table_1.txt')
taxa.1.asv <- taxa.1.asv %>% rename(Taxonomy.2 = Taxonomy)
taxa.1.merged <- taxa.1.asv %>% left_join(taxa.1, by = "ASV_num")
#taxa.2.species
taxa.1.species <- taxa.1.species %>%
  mutate(Taxonomy = case_when(
    is.na(Species) & !is.na(Genus) ~ paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = ";"),
    is.na(Genus) & !is.na(Family) ~ paste(Kingdom, Phylum, Class, Order, Family, sep = ";"),
    is.na(Family) & !is.na(Order) ~ paste(Kingdom, Phylum, Class, Order, sep = ";"),
    is.na(Order) & !is.na(Class) ~ paste(Kingdom, Phylum, Class, sep = ";"),
    is.na(Class) & !is.na(Phylum) ~ paste(Kingdom, Phylum, sep = ";"),
    is.na(Phylum) & !is.na(Kingdom) ~ Kingdom,
    TRUE ~ paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";"))) %>%
  relocate(Taxonomy, .before = Kingdom)
taxa.1.species.merged <- taxa.1.asv %>% left_join(taxa.1.species, by = "ASV_num")
# 将 Taxonomy 赋值给 Taxonomy.3，并将 Taxonomy.3 插入到 Taxonomy.2 后面
taxa.1.merged <- taxa.1.merged %>%
  mutate(Taxonomy.muti = taxa.1.species.merged$Taxonomy,
         Species.muti = taxa.1.species.merged$Species,
         Diff = case_when(
           is.na(Species.muti) & !is.na(Species) ~ 1,
           !is.na(Species.muti) & is.na(Species) ~ 1,
           Species.muti != Species ~ 1,
           is.na(Species.muti) & is.na(Species) ~ 0,
           TRUE ~ 0)) %>%
  relocate(Taxonomy.muti, .after = Taxonomy.2)
#####excel####
# 添加每个数据框到工作簿中，每个数据框一个工作表
wb.1 <- createWorkbook()
#sheet
addWorksheet(wb.1, "taxa.1")
writeData(wb.1, sheet = "taxa.1", taxa.1.merged, rowNames = FALSE)
#sheet
addWorksheet(wb.1, "taxa.1.species")
writeData(wb.1, sheet = "taxa.1.species", taxa.1.species.merged, rowNames = FALSE)
#sheet
addWorksheet(wb.1, "taxa.count.1.boot")
writeData(wb.1, sheet = "taxa.count.1.boot", taxa.count.1.boot, rowNames = FALSE)
#xlsx
saveWorkbook(wb.1, "ASV.count.1.tax.xlsx", overwrite = TRUE)
#finall####
sessionInfo()
