####### CORRELATION BETWEEN DATABASE FOR VALIDATION OF LATITUDINAL PATTERN #######
## CYTB
# load the file
library(readxl)
cytb_all <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/Birds_cytb_Latband.xlsx")
# filter teh database
cytb_all <- cytb_all[,c(1,2,3,9,12,13)]
colnames(cytb_all) <- c("latband", "all_num_sps", "all_num_seqs", "all_pair_muts_bp", "all_tot_bp", "all_GD")
#load second file
cytb_ir <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/cytb_latband_insiderange.xlsx")
colnames(cytb_ir) <- c("latband", "ir_num_sps", "ir_num_seqs", "ir_pair_muts_bp", "ir_tot_bp", "ir_GD")
#merge
cytb <- merge(cytb_all, cytb_ir, by = "latband")
cytb2 <- cbind(cytb_ir, cytb_all)
cytb2 <- cytb2[, -7]
#correlation test and plot
cor_obj <- cor.test(cytb2$ir_GD, cytb2$all_GD)
plot(cytb2$ir_GD, cytb2$all_GD, xlab = "GD breeding range", ylab = "GD all sequences")

## CO1
co1_all <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_latband_allseqs.xlsx")
colnames(co1_all) <- c("latband", "all_num_sps", "all_num_seqs", "all_pair_muts_bp", "all_tot_bp", "all_GD")
co1_ir <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_latband_insiderange.xlsx")
colnames(co1_ir) <- c("latband", "ir_num_sps", "ir_num_seqs", "ir_pair_muts_bp", "ir_tot_bp", "ir_GD")
co1 <- cbind(co1_ir,co1_all)
co1 <- co1[, -7]
cor_obj_co1 <- cor.test(co1$ir_GD, co1$all_GD)
plot(co1$ir_GD, co1$all_GD, xlab = "GD breeding range", ylab = "GD all sequences")
