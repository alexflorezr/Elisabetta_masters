###### NUCLEOTIDE DIVERSITY AMONG IUCN CATEGORIES ####
## ---- cytb_iucn ----
# load the datasets
cytb_gd <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/cytb_all_species/matlab_output_table.csv", header = T)
cytb_db <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/cytb_all_species/Cytb_database.csv", sep = "\t", header = F)
# subset only the names as they are in gd table and iucn categories
cytb_db <- cytb_db[,c(6,16)]
colnames(cytb_db) <- c("Species", "IUCN")
cytb_db[,1] <- as.character(cytb_db[,1])
cytb_db[,2] <- as.character(cytb_db[,2])

# create a database with only the species name and the respective IUCN category
species <- unique(cytb_db$Species)
iucn_db <- as.data.frame(matrix(ncol = 2, nrow = length(species)), stringsAsFactors = F)
colnames(iucn_db) <- c("Species", "IUCN")
i <- 1
for(s in seq_along(species)){
  subdb <- subset(cytb_db, cytb_db$Species == species[s])
  iucn_db[i,] <- c(species[s], subdb[1,2])
  i <- i+1
}

# merge the two dataframes 
cytb_merged <- merge(cytb_gd, iucn_db, by = "Species")
# remove NaNs
cytb_merged <- cytb_merged[-which(is.nan(cytb_merged$Nuc_div)),]
cytb_merged <- cytb_merged[-which(is.na(cytb_merged$IUCN)),]
cytb_merged_more5seqs <- subset(cytb_merged, cytb_merged$Num_seqs > 5)
cytb_merged_more10seqs <- subset(cytb_merged, cytb_merged$Num_seqs > 10)
cytb_merged_corrected <- cytb_merged[-which(cytb_merged$Nuc_div > 0.1),]

# merging some of the categories and excluding NR (= not recognized)
cytb_less_categories <- cytb_merged_more5seqs
#unique(cytb_less_categories$IUCN) # look which categories are present in order to know what to merge

cytb_less_categories[which(cytb_less_categories$IUCN == "CR (PEW)"),5] <- "CR"
cytb_less_categories <- cytb_less_categories[-which(cytb_less_categories$IUCN == "NR"),]
cytb_less_categories$IUCN <- factor(cytb_less_categories$IUCN, levels = c("LC", "NT", "VU", "EN", "CR", "EW"))

## ---- boxplot1 ----
# create boxplot
par(mar = c(5,5,4,2))
b <- boxplot(Nuc_div ~ IUCN, data = cytb_merged_more5seqs, xlab = "IUCN cateogries", ylab = "Nucleotide diversity", cex.axis = 0.8, las =2, main = "Species with more than 5 sequences")

## ---- boxplot2 ----
# create boxplot with less categories
b2 <- boxplot(Nuc_div ~ IUCN, data = cytb_less_categories, 
              xlab = "IUCN cateogries", ylab = "Nucleotide diversity", 
              col = c("#228B22","#C0FF3E","#FFD700","#EE7600","#CD2626","#551A8B"),
              main = "Cytochrome-b")

## ---- Kruskal-Wallis test and post hoc Dunn test ----
library(PMCMR)
kruskal.test(Nuc_div ~ IUCN, data = cytb_less_categories)
posthoc.kruskal.dunn.test(Nuc_div ~ factor(IUCN), data = cytb_less_categories, p.adjust.method = 'holm')
#posthoc.kruskal.dunn.test(Nuc_div ~ factor(IUCN), data = cytb_less_categories, p.adjust.method = 'bonferroni')

## ---- categories ----
LC <- cytb_merged_more5seqs[which(cytb_merged_more5seqs$IUCN == 'LC'),2]
NT <- cytb_merged_more5seqs[which(cytb_merged_more5seqs$IUCN == 'NT'),2]
VU <- cytb_merged_more5seqs[which(cytb_merged_more5seqs$IUCN == 'VU'),2]
EN <- cytb_merged_more5seqs[which(cytb_merged_more5seqs$IUCN == 'EN'),2]
CR <- cytb_merged_more5seqs[which(cytb_merged_more5seqs$IUCN == 'CR' | cytb_merged_more5seqs$IUCN == 'CR (PE)' | cytb_merged_more5seqs$IUCN == 'CR (PEW)'),2]
EW <- cytb_merged_more5seqs[which(cytb_merged_more5seqs$IUCN == 'EW'),2]

grouped <- list(g1=c(LC, NT), g2=VU, g3=EN, g4=CR, g5=EW)
kruskal.test(grouped)
posthoc.kruskal.dunn.test(grouped, p.adjust.method = 'holm')




