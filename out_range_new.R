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
oor_backup <- out_of_range

# subset the database to see which species are present in the range database and which are not 
oor_missing_name <- subset(out_of_range, out_of_range$In_ranges_file == "NO")
oor_present_name <- subset(out_of_range, out_of_range$In_ranges_file == "YES")

# subset the database to have included only the species with sequences outside the range
more_than_one <- subset(oor_present_name, oor_present_name$Cells_out >= 1)
sort(as.numeric(more_than_one$Seqs_out))
sum(as.numeric(more_than_one$Seqs_out))

# panel plot
pdf(file="outside_the_range_cytb.pdf")
par(mfrow=c(2,1))
# plot the barplot for the number of cells
oor_table <- table(as.numeric(oor_present_name$Cells_out))
barplot(oor_table, ylim=c(0, 1700))
mtext("Number of cells outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of species =", dim(oor_present_name)[1], sep=" "), side=3, line=2, cex=2)
text(0.8, 1120, labels = "1078", cex = 0.8)
text(2, 600, labels = "560", cex = 0.8)
text(3.1, 180, labels = "143", cex = 0.8)
text(4.3, 90, labels = "55", cex = 0.8)
text(5.5, 50, labels = "26", cex = 0.8)
text(6.7, 45, labels = "18", cex = 0.8)
text(7.9, 45, labels = "21", cex = 0.8)
text(9.1, 39, labels = "11", cex = 0.8)
text(10.3, 32, labels = "8", cex = 0.8)
text(11.5, 29, labels = "5", cex = 0.8)

# plot the barplot for the number of sequences
oor_table <- table(as.numeric(oor_present_name$Seqs_out))
oor_table_20_plus <- c(oor_table[1:20], sum(oor_table[21:length(oor_table)]))
names(oor_table_20_plus)[21] <- "20+"
barplot(oor_table_20_plus, ylim=c(0, 1700))
mtext("Number of sequences outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
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

# panel plot
pdf(file="outside_the_range_co1.pdf")
par(mfrow=c(2,1))
# plot the barplot for the number of cells
co1_table <- table(as.numeric(co1_present_name$Cells_out))
barplot(co1_table, ylim=c(0, 2000))
mtext("Number of cells outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of species =", dim(co1_present_name)[1], sep=" "), side=3, line=2, cex=2)
text(0.7, 1550, labels = "1516", cex = 0.8)
text(1.9, 560, labels = "530", cex = 0.8)
text(3.1, 360, labels = "330", cex = 0.8)
text(4.3, 230, labels = "203", cex = 0.8)
text(5.5, 145, labels = "115", cex = 0.8)
text(6.7, 100, labels = "72", cex = 0.8)
text(7.9, 65, labels = "36", cex = 0.8)
text(9.1, 52, labels = "24", cex = 0.8)
text(10.3, 45, labels = "19", cex = 0.8)
text(11.5, 38, labels = "12", cex = 0.8)

# plot the barplot for the number of sequences
co1_table <- table(as.numeric(co1_present_name$Seqs_out))
co1_table_20_plus <- c(co1_table[1:20], sum(co1_table[21:length(co1_table)]))
names(co1_table_20_plus)[21] <- "20+"
barplot(co1_table_20_plus, ylim=c(0, 1700))
mtext("Number of sequences outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
dev.off()


