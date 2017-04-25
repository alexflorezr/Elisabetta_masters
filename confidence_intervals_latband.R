######## CREATION OF CONFIDENCE INTERVALS FOR LATITUDINAL BANDS #######
library(readxl)

cytb_matlab <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/cytb_inside_range/matlab_latband_result.csv", sep = "\t", header = F)
cytb_matlab$V8 = paste(cytb_matlab$V2, cytb_matlab$V3, sep = " - ")
colnames(cytb_matlab) <- c("Sp_name", "Latband1", "Latband2", "Num_seqs", "Nuc_div", "Num_mut_bp", "Tot_bp", "Latband")

cytb_sorted <- cytb_matlab[order(cytb_matlab$Latband),]
latband <- unique(cytb_sorted$Latband)

i <- 0
for(l in latband){
  subset_db <- subset(cytb_sorted, cytb_sorted$Latband == l)
  count <- nrow(subset_db)
  avg <- mean(sample(subset_db$Nuc_div, count, replace = T))
  i[l] <- avg
}
