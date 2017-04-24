###### RESULT PLOTS CO1 ####### LATITUDINAL BANDS ######
## ---- CO1_IR_LATBAND ----
# load the file --> inside the breeding range
library(readxl)

ir_co1 <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_latband_insiderange.xlsx", col_names=TRUE)
ir_co1$Latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")

# create the plot -> barplot with genetic diversity and number of sequences
layout(rbind(1,2), heights=c(7,1))
par(mar=c(5, 6, 3, 4))

barplot(ir_co1$`GD value`, horiz = T, cex.axis = 1, xlim = c(0, 0.008), cex.lab = 1, col="#FF149390", names.arg = ir_co1$Latband, las = 1)

par(new=T)

barplot(ir_co1$`num seqs`, horiz = T, col="#69696990", las = 1, axes = F)
axis(3, cex.axis = 1, at = c(0, 300, 600, 900, 1200))
title(sub = "Inside breeding range only (CO1)", xlab="Genetic diversity", font.sub=2)
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("GD", "Number of sequences"), fill = c("#FF149390", "#69696990"), cex= 0.8, bty = "n", ncol = 2)

##################################
## ---- CO1_ALL_LATBAND ----
# load the file --> all sequences
library(readxl)
all_co1 <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_latband_allseqs.xlsx")
all_co1$Latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")

# create the plot -> barplot with genetic diversity and number of sequences
layout(rbind(1,2), heights=c(7,1))
par(mar=c(5, 6, 3, 4))

barplot(all_co1$`GD value`, horiz = T, cex.axis = 1, xlim = c(0, 0.008), cex.lab = 1, col="#FF149390", names.arg = all_co1$Latband, las = 1)

par(new=T)

barplot(all_co1$`num seqs`, horiz = T, col="#69696990", las = 1, axes = F)
axis(3, cex.axis = 1, at = c(0, 500, 1000, 1500, 2000))
title(sub = "All sequences (CO1)", xlab="Genetic diversity", font.sub=2)
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("GD", "Number of sequences"), fill = c("#FF149390", "#69696990"), cex= 0.8, bty = "n", ncol = 2)

#################################
# barplot with genetic diversity of the two datasets on top of each other
## ---- CO1_INSOUT_LATBAND ----
library(readxl)
all_co1 <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_latband_allseqs.xlsx")
all_co1$Latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")
ir_co1 <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_latband_insiderange.xlsx", col_names=TRUE)
ir_co1$Latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")

layout(rbind(1,2), heights=c(7,1))
par(mar=c(5, 6, 3, 4) + 0.5)
barplot(names.arg=ir_co1$Latband, ir_co1$`GD value`, col = "#00BFFF90", las = 1, horiz = T, xlim = c(0,0.008))
barplot(names.arg=ir_co1$Latband, all_co1$`GD value`, add = T, col = "#F0808090", las = 1, horiz = T, xlim = c(0,0.008))
title("Genetic diversity over latitudinal bands (CO1)", xlab = "Genetic Diversity")
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("inside breeding range", "all sequences"), fill = c("#00BFFF90", "#F0808090"), cex= 0.8, bty = "n", ncol = 2)
