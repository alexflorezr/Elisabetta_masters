##################### CYTB OUT OF RANGE ######################

rm(list=ls())
library(raster)
empty_raster <- raster()

# open and clean database
library(readxl)
cytb_db <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Cytb_database.xlsx") 
cytb_db <- cytb_db[-1,]

# subset original database and remove lines with string "NA" in coordinates columns
little_db <- cytb_db[,c(4,5,11,12)]
colnames(little_db) <- c("IOC_name", "CMEC_name", "latitude", "longitude")
without_NA <- little_db[-which(little_db$latitude == "NA"),]

# open breeding range database
b_ranges <- read.table("/Users/Elisabetta/Documents/UCPH/Thesis/Data/WorldBreedingBirdsRanges_CMEC_2016-09-12.txt", header = F, stringsAsFactors = F, sep = ",") 

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

cytb_result <- outside_range_f(without_NA, b_ranges)

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
barplot(oor_table, ylim=c(0, 1700))
mtext("Number of cells outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of species =", dim(cytb_present_name)[1], sep=" "), side=3, line=2, cex=2)
mtext(paste("Total number of cells out of range =", sum(as.numeric(cytb_present_name$Cells_out)), sep = " "), side = 3, cex = 1)
text(0.8, 1120, labels = "1078", cex = 0.8)
text(2, 750, labels = "704", cex = 0.8)
text(3.1, 150, labels = "114", cex = 0.8)
text(4.3, 70, labels = "34", cex = 0.8)
text(5.5, 55, labels = "15", cex = 0.8)
text(6.7, 40, labels = "6", cex = 0.8)
text(7.9, 40, labels = "5", cex = 0.8)
text(9.1, 35, labels = "3", cex = 0.8)
text(10.3, 35, labels = "3", cex = 0.8)
text(11.5, 35, labels = "3", cex = 0.8)
text(12.7, 35, labels = "1", cex = 0.8)
dev.off()

# plot the barplot for the number of sequences
pdf("seq_outside_range_cytb.pdf")
seq_table <- table(as.numeric(cytb_present_name$Seqs_out))
seq_table_20_plus <- c(seq_table[1:20], sum(seq_table[21:length(seq_table)]))
names(seq_table_20_plus)[21] <- "20+"
barplot(seq_table_20_plus, ylim=c(0, 1700))
mtext("Number of sequences outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of sequences out of range =", sum(as.numeric(cytb_present_name$Seqs_out)), sep = " "), side = 3, cex = 1.5)
text(0.7, 1115, labels = "1078", cex = 0.6)
text(1.9, 590, labels = "560", cex = 0.6)
text(3.1, 170, labels = "143", cex = 0.6)
text(4.3, 82, labels = "55", cex = 0.6)
text(5.5, 53, labels = "26", cex = 0.6)
text(6.7, 40, labels = "18", cex = 0.6)
text(7.9, 45, labels = "21", cex = 0.6)
text(9.1, 32, labels = "11", cex = 0.6)
text(10.3, 30, labels = "8", cex = 0.6)
text(11.5, 27, labels = "3", cex = 0.6)
text(12.7, 27, labels = "4", cex = 0.6)
text(13.9, 27, labels = "2", cex = 0.6)
text(15.1, 27, labels = "3", cex = 0.6)
text(16.3, 27, labels = "1", cex = 0.6)
text(17.5, 27, labels = "3", cex = 0.6)
text(18.7, 27, labels = "3", cex = 0.6)
text(19.9, 27, labels = "1", cex = 0.6)
text(21.1, 27, labels = "3", cex = 0.6)
text(22.3, 27, labels = "3", cex = 0.6)
text(23.5, 27, labels = "1", cex = 0.6)
text(24.7, 41, labels = "17", cex = 0.6)
dev.off()


#################### CO1 calculation #######################

# open and clean the dataset
co1_db <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/CO1_database.xlsx")
co1_db <- co1_db[-1,]

ldb <- co1_db[,c(4,5,11,12)]
colnames(ldb) <- c("IOC_name", "CMEC_name", "latitude", "longitude")
ldb <- ldb[-which(ldb$latitude == "NA"),]

# run the function
co1_result <- outside_range_f(ldb, b_ranges)

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
barplot(co1_table, ylim=c(0, 2000))
mtext("Number of cells outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of species =", dim(co1_present_name)[1], sep=" "), side=3, line=2, cex=2)
mtext(paste("Total number of cells out of range =", sum(as.numeric(co1_present_name$Cells_out)), sep = " "), side = 3, cex = 1)
text(0.7, 1570, labels = "1516", cex = 0.8)
text(1.9, 930, labels = "886", cex = 0.8)
text(3.1, 320, labels = "273", cex = 0.8)
text(4.3, 140, labels = "104", cex = 0.8)
text(5.5, 100, labels = "55", cex = 0.8)
text(6.7, 80, labels = "37", cex = 0.8)
text(7.9, 65, labels = "16", cex = 0.8)
text(9.1, 51, labels = "8", cex = 0.8)
text(10.3, 38, labels = "2", cex = 0.8)
text(11.5, 38, labels = "2", cex = 0.8)
text(12.7, 38, labels = "1", cex = 0.8)
text(13.9, 38, labels = "1", cex = 0.8)
dev.off()

# plot the barplot for the number of sequences
pdf("seq_outside_the_range.pdf")
co1_seq_table <- table(as.numeric(co1_present_name$Seqs_out))
co1_table_20_plus <- c(co1_seq_table[1:20], sum(co1_seq_table[21:length(co1_seq_table)]))
names(co1_table_20_plus)[21] <- "20+"
barplot(co1_table_20_plus, ylim=c(0, 1700))
mtext("Number of sequences outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of sequences out of range =", sum(as.numeric(co1_present_name$Seqs_out)), sep = " "), side = 3, cex = 1.5)
text(0.7, 1550, labels = "1516", cex = 0.6)
text(1.9, 570, labels = "530", cex = 0.6)
text(3.1, 370, labels = "330", cex = 0.6)
text(4.3, 235, labels = "203", cex = 0.6)
text(5.5, 150, labels = "115", cex = 0.6)
text(6.7, 100, labels = "72", cex = 0.6)
text(7.9, 65, labels = "36", cex = 0.6)
text(9.1, 51, labels = "24", cex = 0.6)
text(10.3, 47, labels = "19", cex = 0.6)
text(11.5, 40, labels = "12", cex = 0.6)
text(12.7, 35, labels = "8", cex = 0.6)
text(13.9, 35, labels = "5", cex = 0.6)
text(15.1, 38, labels = "9", cex = 0.6)
text(16.3, 27, labels = "1", cex = 0.6)
text(17.5, 27, labels = "2", cex = 0.6)
text(18.7, 27, labels = "3", cex = 0.6)
text(19.9, 27, labels = "1", cex = 0.6)
text(21.1, 27, labels = "2", cex = 0.6)
text(22.3, 27, labels = "2", cex = 0.6)
text(23.5, 27, labels = "2", cex = 0.6)
text(24.7, 38, labels = "9", cex = 0.6)
dev.off()


