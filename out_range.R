##################### CALCULATING NUMBER OF SEQUENCES OUTSIDE THE BREEDING RANGE ######################

# load raster package and create empty raster
rm(list=ls())
library(raster)
empty_raster <- raster()

# open breeding range database
b_ranges <- read.table("/Users/Elisabetta/Documents/UCPH/Thesis/Data/WorldBreedingBirdsRanges_CMEC_2016-09-12.txt", header = F, stringsAsFactors = F, sep = ",") 

# open and clean cytochrome-b database
library(readxl)
cytb_db <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Cytb_database.xlsx") 
cytb_db <- cytb_db[-1,]

# subset original database and remove lines with string "NA" in coordinates columns
cytb_subset <- cytb_db[,c(2,4,5,12,11)]
colnames(cytb_subset) <- c("Accession_number", "IOC_name", "CMEC_name", "longitude", "latitude")
cytb_subset <- cytb_subset[-which(cytb_subset$latitude == "NA"),]
cytb_subset[,4] <- as.numeric(cytb_subset$longitude)
cytb_subset[,5] <- as.numeric(cytb_subset$latitude)

# loop over ranges to find cells and sequences of species that fall outside the range.
# I use my database with the species, I go through the names, take the correspondent binomial that is supposed to be 
# in the range database and use that name to compare with the species range.
# The final results includes in the first column the names of the species in my database
outside_range_f <- function(mydb, rangedb){
  sp_names <- unique(mydb$IOC_name)
  out_of_range <- as.data.frame(matrix(nrow=length(sp_names), ncol=6))
  colnames(out_of_range) <- c("IOC_name", "In_ranges_file", "Cells_out","Total_range_cells", "Seqs_out", "Total_seqs")
  for(s in seq_along(sp_names)){
    temp_db <- subset(mydb, mydb$IOC_name == sp_names[s])
    to_check_name <- unique(temp_db$CMEC_name) 
    if(is.element(to_check_name, rangedb$V1)){
      temp_range_sp <- subset(rangedb, rangedb$V1 == to_check_name)
      cells_range <- cellFromXY(empty_raster, cbind(temp_range_sp$V4, temp_range_sp$V5))
      cells_samples <- cellFromXY(empty_raster, cbind(as.numeric(temp_db$longitude), as.numeric(temp_db$latitude)))
      how_many_seqs <- sum(is.na(match(cells_samples, unique(cells_range))))
      how_many_cells <- sum(is.na(match(unique(cells_samples), unique(cells_range))))
      out_of_range[s,] <- c(sp_names[s], "YES", how_many_cells, length(unique(cells_range)), how_many_seqs, dim(temp_db)[1])
    } else {
      out_of_range[s,] <- c(sp_names[s], "NO",NA, NA, NA,dim(temp_db)[1])
    }
    print(s)
  }
  return(out_of_range)
}

cytb_result <- outside_range_f(cytb_subset, b_ranges)

# backup
cytb_backup <- cytb_result

# subset the database to see which species are present in the range database and which are not 
cytb_missing_name <- subset(cytb_result, cytb_result$In_ranges_file == "NO")
cytb_present_name <- subset(cytb_result, cytb_result$In_ranges_file == "YES")

# calculate percentage of sequences out of range over total sequences
cytb_present_name$pct_out <- NA
pct_f <- function(db){
  seq_out <- db$Seqs_out
  tot_seq <- db$Total_seqs
  r <- (100*as.numeric(seq_out))/as.numeric(tot_seq)
  return(r)
}

pct <- pct_f(cytb_present_name)
cytb_present_name$pct_out <- as.integer(pct)

# subset the database to have included only the species with sequences outside the range
more_than_one <- subset(cytb_present_name, cytb_present_name$Cells_out >= 1)
sort(as.numeric(more_than_one$Seqs_out))
sum(as.numeric(more_than_one$Seqs_out))
# 3253 total sequences out of range

# panel plot
pdf(file="cells_outside_range_cytb.pdf")
#par(mfrow=c(2,1))
# plot the barplot for the number of cells
oor_table <- table(as.numeric(cytb_present_name$Cells_out))
cytb_bp <- barplot(oor_table, ylim=c(0, 1700))
mtext("Number of cells outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of species =", dim(cytb_present_name)[1], sep=" "), side=3, line=2, cex=2)
mtext(paste("Total number of cells out of range =", sum(as.numeric(cytb_present_name$Cells_out)), sep = " "), side = 3, cex = 1)
text(cytb_bp, oor_table+50, labels = oor_table, cex = 0.8)
dev.off()

# plot the barplot for the number of sequences
pdf("seq_outside_range_cytb.pdf")
seq_table <- table(as.numeric(cytb_present_name$Seqs_out))
seq_table_20_plus <- c(seq_table[1:20], sum(seq_table[21:length(seq_table)]))
names(seq_table_20_plus)[21] <- "20+"
cytb_bp_20 <- barplot(seq_table_20_plus, ylim=c(0, 1700))
mtext("Number of sequences outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of sequences out of range =", sum(as.numeric(cytb_present_name$Seqs_out)), sep = " "), side = 3, cex = 1.5)
text(cytb_bp_20, seq_table_20_plus + 50, labels = seq_table_20_plus, cex = 0.6)
dev.off()


#################### CO1 calculation #######################

# open and clean the dataset
co1_db <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/CO1_database.xlsx")
co1_db <- co1_db[-1,]

co1_subset <- co1_db[,c(2,4,5,12,11)]
colnames(co1_subset) <- c("Accession_number", "IOC_name", "CMEC_name", "longitude", "latitude")
co1_subset <- co1_subset[-which(co1_subset$latitude == "NA"),]
co1_subset[,4] <- as.numeric(co1_subset$longitude)
co1_subset[,5] <- as.numeric(co1_subset$latitude)

# run the function
co1_result <- outside_range_f(co1_subset, b_ranges)

# subset the resulting database
co1_missing_name <- subset(co1_result, co1_result$In_ranges_file == "NO")
co1_present_name <- subset(co1_result, co1_result$In_ranges_file == "YES")
co1_cells_out <- subset(co1_present_name, co1_present_name$Cells_out >= 1)

# calculate how many sequences are outside the range
sort(as.numeric(co1_cells_out$Seqs_out))
sum(as.numeric(co1_cells_out$Seqs_out))
# 4219 sequences out of range

# calculate percentage of sequences out of range over total seqs per species
co1_present_name$pct_out <- NA
co1_pct <- pct_f(co1_present_name)
co1_present_name$pct_out <- as.integer(co1_pct)

# panel plot
pdf(file="cells_outside_the_range_co1.pdf")
#par(mfrow=c(2,1))
# plot the barplot for the number of cells
co1_table <- table(as.numeric(co1_present_name$Cells_out))
co1_bp <- barplot(co1_table, ylim=c(0, 2000))
mtext("Number of cells outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of species =", dim(co1_present_name)[1], sep=" "), side=3, line=2, cex=2)
mtext(paste("Total number of cells out of range =", sum(as.numeric(co1_present_name$Cells_out)), sep = " "), side = 3, cex = 1)
text(co1_bp, co1_table + 50, labels = co1_table, cex = 0.8)
dev.off()

# plot the barplot for the number of sequences
pdf("seq_outside_the_range.pdf")
co1_seq_table <- table(as.numeric(co1_present_name$Seqs_out))
co1_table_20_plus <- c(co1_seq_table[1:20], sum(co1_seq_table[21:length(co1_seq_table)]))
names(co1_table_20_plus)[21] <- "20+"
co1_bp_20 <- barplot(co1_table_20_plus, ylim=c(0, 1700))
mtext("Number of sequences outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of sequences out of range =", sum(as.numeric(co1_present_name$Seqs_out)), sep = " "), side = 3, cex = 1.5)
text(co1_bp_20, co1_table_20_plus+50, labels = co1_table_20_plus, cex = 0.6)
dev.off()

