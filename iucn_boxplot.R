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

# create boxplot
par(mar = c(5,5,4,2))
boxplot(Nuc_div ~ IUCN, data = cytb_merged, xlab = "IUCN cateogries", ylab = "Nucleotide diversity")

# merging some of the categories and excluding NR (= not recognized)
cytb_merged[which(cytb_merged$IUCN == "CR (PE)"),4] <- "CR"
cytb_merged[which(cytb_merged$IUCN == "CR (PEW)"),4] <- "CR"
cytb_merged <- cytb_merged[-which(cytb_merged$IUCN == "NR"),]
cytb_merged$IUCN <- factor(cytb_merged$IUCN, levels = c("DD", "LC", "NT", "VU", "EN", "CR", "EW"))

boxplot(Nuc_div ~ IUCN, data = cytb_merged, xlab = "IUCN cateogries", ylab = "Nucleotide diversity", col = c("#DEDEDE","#228B22","#C0FF3E","#FFD700","#EE7600","#CD2626","#551A8B"))




