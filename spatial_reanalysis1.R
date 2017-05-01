######## SPATIAL RE-ANALYSIS: EXCLUSION OF POORLY SAMPLED SPECIES #####
######## CYTB ######
library(readxl)

cytb_all <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/cytb_grid_insiderange.xlsx")
cytb_more_than_5 <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/cytb_grid_more5seqs.csv")

cytb_merged <- merge(cytb_all, cytb_more_than_5, by = 'Grid_ID')
cytb_merged <- cytb_merged[, c(1,6,11)]


par(mar = c(5,5,4,2))
plot(cytb_merged$GD.x, cytb_merged$GD.y, xlab = "GD all sequences", ylab = "GD more than 5 sequences", main = "Cytochrome-b")
cytb_cor <- cor.test(x = cytb_merged$GD.x, y = cytb_merged$GD.y)
r_square <- cytb_cor$estimate^2
l <- lm(cytb_merged$GD.y ~ cytb_merged$GD.x)
abline(l, col = "red")

######### CO1 ########
library(readxl)

co1_all <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_grid_insiderange.xlsx")
co1_more_than_5 <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_grid_more5seqs.csv")
colnames(co1_all) <- c("Grid_ID", "Num_sps", "Num_seqs", "Nuc_div", "Tot_bp", "GD")

co1_merged <- merge(co1_all, co1_more_than_5, by = "Grid_ID")
co1_merged <- co1_merged[, c(1,6,11)]

par(mar = c(5,5,4,2))
plot(co1_merged$GD.x, co1_merged$GD.y, xlab = "GD all sequences", ylab = "GD more than 5 sequences", main = "CO1")
co1_cor <- cor.test(x = co1_merged$GD.x, y = co1_merged$GD.y)
r_square2 <- co1_cor$estimate^2
l2 <- lm(co1_merged$GD.y ~ co1_merged$GD.x)
abline(l2, col = "red")
