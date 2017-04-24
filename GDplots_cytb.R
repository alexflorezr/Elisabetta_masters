### RESULT PLOTS CYTOCHROME-B ####### LATITUDINAL BANDS ########

## ---- LATBAND_IR ----
# load the file --> inside the breeding range
library(readxl)

gd_cytb <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/cytb_latband_insiderange.xlsx", col_names=TRUE)
gd_cytb$Latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")

# barplot with genetic diversity and number of sequences
layout(rbind(1,2), heights=c(7,1))
par(mar=c(5, 6, 3, 4))

barplot(gd_cytb$`GD value`, horiz = T, cex.axis = 1, xlim = c(0, 0.01), cex.lab = 1, col="#76EE0090", names.arg = gd_cytb$Latband, las = 1)

par(new=T)

barplot(gd_cytb$`num seqs`, horiz = T, col="#69696990", las = 1, axes = F)
axis(3, cex.axis = 1, at = c(0, 300, 600, 900, 1200))
title(sub = "Inside breeding range only (Cytochrome-b)", xlab="Genetic diversity", font.sub=2)
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("GD", "Number of sequences"), fill = c("#76EE0090", "#69696990"), cex= 0.8, bty = "n", ncol = 2)

##################################
## ---- LATBAND_ALL ----
# load the file --> all sequences
library(readxl)
all_cytb <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/Birds_cytb_Latband.xlsx")
all_cytb$latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")

# barplot with genetic diversity and number of sequences
layout(rbind(1,2), heights=c(7,1))
par(mar=c(5, 6, 3, 4))

barplot(all_cytb$GD_divided_by_sp, horiz = T, cex.axis = 1, xlim = c(0, 0.008), cex.lab = 1, col="#76EE0090", names.arg = all_cytb$Latband, las = 1)

par(new=T)

barplot(all_cytb$num_seqs, horiz = T, col="#69696990", las = 1, axes = F)
axis(3, cex.axis = 1, at = c(0, 300, 600, 900, 1200))
title(sub = "All sequences (Cytochrome-b)", xlab="Genetic diversity", font.sub=2)
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("GD", "Number of sequences"), fill = c("#76EE0090", "#69696990"), cex= 0.8, bty = "n", ncol = 2)

##################################
# barplot with genetic diversity of the two datasets on top of each other
## ---- LATBAND_INSOUT ----
library(readxl)
all_cytb <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/Birds_cytb_Latband.xlsx")
all_cytb$latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")
gd_cytb <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/cytb_latband_insiderange.xlsx", col_names=TRUE)
gd_cytb$Latband <- c("-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")

layout(rbind(1,2), heights=c(7,1))
par(mar=c(5, 6, 3, 4) + 0.5)
barplot(names.arg=gd_cytb$Latband, gd_cytb$`GD value`, col = "#00BFFF90", las = 1, horiz = T, xlim = c(0,0.01))
barplot(names.arg=gd_cytb$Latband, all_cytb$GD_divided_by_sp, add = T, col = "#F0808090", las = 1, horiz = T, xlim = c(0,0.01))
title("Genetic diversity over latitudinal bands (Cytochrome-b)", xlab = "Genetic Diversity")
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("inside breeding range", "all sequences"), fill = c("#00BFFF90", "#F0808090"), cex= 0.8, bty = "n", ncol = 2)

