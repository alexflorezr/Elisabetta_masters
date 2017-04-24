####### CORRELATION BETWEEN DATASETS FOR VALIDATION OF LATITUDINAL PATTERN #######
####### CYTB ###

## ---- CYTB_CORRELATION ----
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
cytb2$latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 -60", "60 - 70", "70 - 80")

#correlation test and plot
cor_obj <- cor.test(cytb2$ir_GD, cytb2$all_GD)
pv <- cor_obj$p.value
r_square <- cor_obj$estimate^2
par(mar=c(5,5,4,2))
plot(cytb2$ir_GD, cytb2$all_GD, xlab = "GD breeding range", ylab = "GD all sequences", xlim=c(0, 0.01), ylim=c(0, 0.01), main = "Cytochrome-b")
lm_obj <- lm(cytb2$all_GD ~ cytb2$ir_GD)
abline(lm_obj, col = "red")
text(0.00712, 0.0085, labels = bquote(paste(R^2~"=", ~.(r_square))), cex = 1)
text(0.0075, 0.009, labels = bquote(paste(p-value~"=", ~.(pv))), cex = 1)

###### CO1 ######
## ---- CO1_CORRELATION ----
# load the files
library(readxl)
co1_all <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_latband_allseqs.xlsx")
colnames(co1_all) <- c("latband", "all_num_sps", "all_num_seqs", "all_pair_muts_bp", "all_tot_bp", "all_GD")
co1_ir <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_latband_insiderange.xlsx")
colnames(co1_ir) <- c("latband", "ir_num_sps", "ir_num_seqs", "ir_pair_muts_bp", "ir_tot_bp", "ir_GD")
co1 <- cbind(co1_ir,co1_all)
co1 <- co1[, -7]

# correlation test and plot
cor_obj_co1 <- cor.test(co1$ir_GD, co1$all_GD)
pv2 <- cor_obj_co1$p.value
r_square_co1 <- cor_obj_co1$estimate^2
lm_obj_co1 <- lm(co1$all_GD ~ co1$ir_GD)

plot(co1$ir_GD, co1$all_GD, xlab = "GD breeding range", ylab = "GD all sequences", main = "CO1")
abline(lm_obj_co1, col = "red")
text(0.006, 0.003, labels = bquote(paste(R^2~"=", ~.(r_square_co1))), cex = 1)
text(0.006, 0.0025, labels = bquote(paste(p-value~"=", ~.(pv2))), cex = 1)
