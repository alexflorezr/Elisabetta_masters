#clean the environment
rm(list=ls())

# set the working directory
setwd("/Users/afr/Desktop/Postdoc/Elisabetta/Data/Birds_data/")

# empty raster file 
library(raster)
empty_raster <- raster()

# upload, clean and separate coordinates to diferent columns
sample_points <- read.table("Master_Database.txt", header = T, stringsAsFactors = F, sep = "\t")
geo_sample_points <- sample_points[!is.na(sample_points$Geonames.coordinates),]
        # If you have the coordinates separated already you can skip this part
        geo_sample_points <- geo_sample_points[!geo_sample_points$Geonames.coordinates == ", ",]
        geo_sample_points <- geo_sample_points[!geo_sample_points$Geonames.coordinates == "insufficient_info",]
        XY <- as.data.frame(matrix(nrow = dim(geo_sample_points)[1], ncol = 2))
        colnames(XY) <- c("x_coor", "y_coor")
        temp_unlist <- unlist(strsplit(geo_sample_points$Geonames.coordinates, split = ","))
        XY$y_coor <- temp_unlist[seq(1,length(temp_unlist), 2)]
        XY$x_coor <- temp_unlist[seq(2,length(temp_unlist), 2)]
        geo_sample_points_XY <- cbind(geo_sample_points, XY)       

# upload ranges
###############
ranges <- read.table("WorldBreedingBirdsRanges_CMEC_2016-09-12.txt", header = F, stringsAsFactors = F, sep = ",")
        
# map the points for georefereced samples
#########################################
library(rworldmap)
map <- getMap("world")
plot(map)
        # this is an example for "Phylloscopus forresti"
        # all range points
        temp_ranges <- subset(ranges, ranges$V1 == "Phylloscopus forresti")
        points(cbind(temp_ranges$V4, temp_ranges$V5))
        # sampling points
        temp_points <- subset(geo_sample_points_XY, geo_sample_points_XY$IOC.Sp..Name == "Phylloscopus forresti")
        points(cbind(temp_points$x_coor, temp_points$y_coor), pch=21, bg="#EE2C2C80", col="#47474790")


# create a new variable using only the genera and the sp name (bi)
##################################################################
geo_sample_points_XY$sp_name_bi <- NA
for (sp in seq_along(geo_sample_points_XY$sp_name_bi)){
        temp_sp <- strsplit(geo_sample_points_XY$IOC.Sp..Name[sp], split = " ")[[1]]
        if(length(temp_sp) > 2){
                temp_sp <- temp_sp[c(1,2)] 
        }
        geo_sample_points_XY$sp_name_bi[sp] <- paste(temp_sp[1], temp_sp[2], sep=" ")
}

# Loop over ranges using the bi names in sample points to estimate
# the number of cells and seqs outside the range
##################################################################
species_bi <- unique(geo_sample_points_XY$sp_name_bi)
out_of_range <- as.data.frame(matrix(nrow=length(species_bi), ncol=6))
colnames(out_of_range) <- c("Species", "In_ranges_file", "Cells_out","Total_range_cells", "Seqs_out", "Total_seqs")
for(s in seq_along(species_bi)){
        if(is.element(species_bi[s],ranges$V1)){
                temp_range_sp <- subset(ranges, ranges$V1 == species_bi[s])
                cells_range <- cellFromXY(empty_raster, cbind(temp_range_sp$V4, temp_range_sp$V5))
                temp_sample_sp <- subset(geo_sample_points_XY,  geo_sample_points_XY$sp_name_bi == species_bi[s])
                cells_samples <- cellFromXY(empty_raster, cbind(as.numeric(temp_sample_sp$x_coor), as.numeric(temp_sample_sp$y_coor)))
                how_many_seqs <- sum(is.na(match(cells_samples, unique(cells_range))))
                how_many_cells <- sum(is.na(match(unique(cells_samples), unique(cells_range))))
                out_of_range[s,] <- c(species_bi[s], "YES", how_many_cells, length(unique(cells_range)), how_many_seqs, dim(temp_sample_sp)[1])
        }else{
                out_of_range[s,] <- c(species_bi[s], "NO",NA, NA, NA,dim(temp_sample_sp)[1])
        }
        print(s)
}
# this is a backup
out_of_range_bck <- out_of_range
# out of range missing name
OOR_missing_name <- subset(out_of_range, out_of_range$In_ranges_file == "NO")
# out of range present name
OOR_present_name <- subset(out_of_range, out_of_range$In_ranges_file == "YES")

more_than_one <- subset(OOR_present_name, OOR_present_name$Cells_out >= 1)
sort(as.numeric(more_than_one$Seqs_out))
sum(as.numeric(more_than_one$Seqs_out))
# panel plot
pdf(file="outside_the_range.pdf")
par(mfrow=c(2,1))
# plot the barplot for the number of cells
OOR_table <- table(as.numeric(OOR_present_name$Cells_out))
barplot(OOR_table, ylim=c(0, 1700))
mtext("Number of cells outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
mtext(paste("Total number of species =", dim(OOR_present_name)[1], sep=" "), side=3, line=2, cex=2)


# plot the barplot for the number of sequences
OOR_table <- table(as.numeric(OOR_present_name$Seqs_out))
OOR_table_20_plus <- c(OOR_table[1:20], sum(OOR_table[21:length(OOR_table)]))
names(OOR_table_20_plus)[21] <- "20+"
barplot(OOR_table_20_plus, ylim=c(0, 1700))
mtext("Number of sequences outside the range", side=1, line=2.5, cex=1.5)
mtext("Number of species", side=2, line=2.5, cex=1.5)
dev.off()

##### this is for the next step, for now DO NOT pay attention
# estimate the distance between points outside and the range
library(geosphere)
distm()
