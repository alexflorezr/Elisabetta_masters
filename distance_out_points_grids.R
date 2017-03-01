################## CHECK IF GENBANK COORDINATES FALL INSIDE THE SPECIES RANGE ##############

cytb_db #starting database

db <- cytb_db[,c(2,4,5,7,8,9,10)]
colnames(db) <- c("Acc_num", "IOC_name", "CMEC_name", "Geonames_lat", "Geonames_long", "Genbank_lat", "Genbank_long")

# I remove lines with no coordinates in the Genbank columns (clean the database)
final_db <- db[which(db$Genbank_lat != "NA"), ]
final_db <- final_db[which(final_db$Genbank_lat != "xxx"), ]

sp_names <- unique(final_db$IOC_name)
out_of_range <- as.data.frame(matrix(nrow=length(sp_names), ncol=6))
colnames(out_of_range) <- c("IOC_name", "In_ranges_file", "Cells_out","Total_range_cells", "Seqs_out", "Total_seqs")
for(s in seq_along(sp_names)){
    temp_db <- subset(final_db, final_db$IOC_name == sp_names[s])
    to_check_name <- unique(temp_db$CMEC_name) 
    if(is.element(to_check_name, b_ranges$V1)){
        temp_range_sp <- subset(b_ranges, b_ranges$V1 == to_check_name)
        cells_range <- cellFromXY(empty_raster, cbind(temp_range_sp$V4, temp_range_sp$V5))
        cells_samples <- cellFromXY(empty_raster, cbind(as.numeric(temp_db$Genbank_long), as.numeric(temp_db$Genbank_lat)))
        how_many_seqs <- sum(is.na(match(cells_samples, unique(cells_range))))
        how_many_cells <- sum(is.na(match(unique(cells_samples), unique(cells_range))))
        out_of_range[s,] <- c(sp_names[s], "YES", how_many_cells, length(unique(cells_range)), how_many_seqs, dim(temp_db)[1])
    }else{
        out_of_range[s,] <- c(sp_names[s], "NO",NA, NA, NA,dim(temp_db)[1])
    }
    print(s)
}

oor_present <- subset(out_of_range, out_of_range$In_ranges_file == "YES")
nrow(oor_present) # 186 species
seq_out <- subset(oor_present, oor_present$Seqs_out >= 1)
nrow(seq_out) # 70 species with coordinates out of range  
par(mfrow=c(4,1))
pdf("Genbank_seq_out_cytb.pdf")
hist(as.numeric(seq_out$Seqs_out), xlim = c(0, 250), xlab = "Number of sequences out of range", ylab = "Number of species", main = "Number of Genbank sequences out of range" )
less_50 <- seq_out[which(seq_out$Seqs_out <= 20), ]
hist(as.numeric(less_50$Seqs_out), xlab = "Number of sequences out of range", ylab = "Number of species", main = "Number of Genbank sequences out of range", ylim = c(0, 50))
mtext("Interval between 0 and 20")
seq_out$pct_out <- NA
pct <- pct_f(seq_out)
seq_out$pct_out <- as.integer(pct)
hist(seq_out$pct_out, main = "Percentage of sequences out of range per species", xlab = "Percentage of sequences")
no_all_seq <- subset(seq_out, seq_out$Seqs_out != seq_out$Total_seqs)
hist(no_all_seq$pct_out, main = "Percentage of sequences out of range per species", xlab = "Percentage of sequences")
mtext("Number of total sequences does not equal the number of sequences out of range", side = 3)
dev.off()

##### CALCULATING DISTANCE BETWEEN GENBANK POINTS AND GEONAMES ESTIMATES ##############
final_db #starting database

# I create a function to calculate euclidean distance that uses both lat and long
euc.dist <- function(x,y,a,b){
  sqrt((x-a)^2 + (y-b)^2)
}

# check that is correct
a <- euc.dist(as.numeric(final_db$Geonames_long[1]), as.numeric(final_db$Geonames_lat[1]), as.numeric(final_db$Genbank_long[1]), as.numeric(final_db$Genbank_lat[1]))

# I get the distance for each row of my database
b <- c()
for(r in 1:nrow(final_db)){
  c <- euc.dist(as.numeric(final_db$Geonames_long[r]), as.numeric(final_db$Geonames_lat[r]), as.numeric(final_db$Genbank_long[r]), as.numeric(final_db$Genbank_lat[r]))
  print(r)
  print(c)
  b <- append(b, c)
}

# I calculate the mean distance
mean(b, na.rm = TRUE)

# Alternative way: use pointDistance() function of raster package
p1 <- cbind(as.numeric(final_db$Geonames_long[1]), as.numeric(final_db$Geonames_lat[1]))
p2 <- cbind(as.numeric(final_db$Genbank_long[1]), as.numeric(final_db$Genbank_lat[1]))

d <- pointDistance(p1, p2, lonlat=F)

d <- c()
for(l in 1:nrow(final_db)){
  p1 <- cbind(as.numeric(final_db$Geonames_long[l]), as.numeric(final_db$Geonames_lat[l]))
  p2 <- cbind(as.numeric(final_db$Genbank_long[l]), as.numeric(final_db$Genbank_lat[l]))
  r <- pointDistance(p1, p2, lonlat = F)
  print(l)
  print(r)
  d <- append(d, r)
}

# b and d are the same
par(mfrow=c(2,1))
pdf("distance_between_points_cytb.pdf")
hist(d, ylim = c(0, 1000), main = "Histogram of distance between points", xlab = "Distance", xlim = c(0, 250))
m <- mean(d, na.rm=T)
abline(v = m, col = "red")
mtext("mean value = 6.956998 (red line)", cex = 0.8, side = 3)
subset <- subset(d, d <= 25)
hist(subset, main = "Histogram of distance between points", xlab = "Distance")
mtext("Values interval between 0 and 25, mean = 6.956998", cex = 0.8, side = 3)
dev.off()

