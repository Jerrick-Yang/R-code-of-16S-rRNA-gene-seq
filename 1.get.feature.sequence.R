#rm(list = ls())
####get path#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#path <- getwd()
#path
library(dplyr, warn.conflicts = F)
####data####
ASV.all <- read.delim('ASV_table.txt', header = TRUE, sep = '\t')
####feature.sequence####
#feature sequences
ASV <- ASV.all
ASV$ASVs <- paste(">", ASV$ASV_num, sep = "") 
feature.sequence <- read.delim('feature.sequence.2.fasta', header = F,  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
colnames(feature.sequence)[colnames(feature.sequence) == 'V1'] <- 'Sequence'
feature.sequence$Sequence <- sub("^>(.*?) .*", ">\\1", feature.sequence$Sequence)
feature.sequence$ASVs <- feature.sequence$Sequence
#feature.sequence$ID <- rownames(feature.sequence)
feature.sequence$ID <- as.numeric(rownames(feature.sequence)) 
#Copies the value of each odd row to the even row that follows it
feature.sequence$ASVs[seq(2, nrow(feature.sequence), by = 2)] <- feature.sequence$Sequence[seq(1, nrow(feature.sequence), by = 2)]
#clean.feature.sequence
result <- inner_join(feature.sequence, ASV, by = "ASVs")
clean.result <- result[, c("Sequence")]
write.table(clean.result, "feature.sequence.fasta", row.names = F, col.names = F,  sep = "\t", quote = FALSE)
#session
sessionInfo()
