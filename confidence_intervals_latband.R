######## CREATION OF CONFIDENCE INTERVALS FOR LATITUDINAL BANDS #######

avg_stdev_f <- function(lband, db){
  errorbar_db <- as.data.frame(matrix(ncol = 3, nrow = length(unique(db$Latband))))
  colnames(errorbar_db) <- c("Latband", "mean", "stdev")
  avgs <- c()
  stddevs <- c()
  for(l in lband){
    temp <- c()
    subset_db <- subset(db, db$Latband == l)
    count <- nrow(subset_db)
    for(iter in 1:100){
      temp[iter] <- mean(sample(subset_db$Nuc_div, count, replace = T))
      # temp has 100 random avgs
    }
    avgs[l] <- mean(temp)
    stddevs[l] <- sd(temp)
  }
  errorbar_db[,1] <- lband
  errorbar_db[,2] <- avgs
  errorbar_db[,3] <- stddevs
  return(errorbar_db)
}

######## CYTB ########
library(readxl)

cytb_matlab <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/cytb_inside_range/matlab_latband_result.csv", sep = "\t", header = F)
cytb_matlab$V8 = paste(cytb_matlab$V2, cytb_matlab$V3, sep = " - ")
colnames(cytb_matlab) <- c("Sp_name", "Latband1", "Latband2", "Num_seqs", "Nuc_div", "Num_mut_bp", "Tot_bp", "Latband")

cytb_sorted <- cytb_matlab[order(cytb_matlab$Latband),]
cytb_sorted <- cytb_sorted[-which(is.nan(cytb_sorted$Nuc_div)),]
latband <- unique(cytb_sorted$Latband)

cytb_errorbar <- avg_stdev_f(latband, cytb_sorted)

######## CO1 ##########
library(readxl)

co1_matlab <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/co1_inside_range/matlab_latband_result.csv", sep = "\t", header = F)
co1_matlab[,8] <- paste(co1_matlab$V2, co1_matlab$V3, sep = " - ")
colnames(co1_matlab) <- c("Sp_name", "Latband1", "Latband2", "Num_seqs", "Nuc_div", "Num_mut_bp", "Tot_bp", "Latband")

co1_sorted <- co1_matlab[order(co1_matlab$Latband),]
co1_sorted <- co1_sorted[-which(is.nan(co1_sorted$Nuc_div)),]
latb <- unique(co1_sorted$Latband)

co1_errorbar <- avg_stdev_f(latb, co1_sorted)



