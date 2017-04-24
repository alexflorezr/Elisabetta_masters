# CREATING DATABASE WITH ONLY SEQUENCES INSIDE RANGE #
#database containing species of which all sequences fall within the range
noout_cytb <- cytb_present_name[which(cytb_present_name$Seqs_out == 0),]
ir_cytb <- cytb_subset[which(cytb_subset$IOC_name %in% noout_cytb$IOC_name),]
ir_cytb <- ir_cytb[,-3]
colnames(ir_cytb) <- c("Accession_number", "Name", "Long", "Lat")

#database containing species which have part of the total sequences falling inside the range
partout_cytb <- cytb_present_name[which(cytb_present_name$Total_seqs != cytb_present_name$Seqs_out),]
partout_cytb <- partout_cytb[which(partout_cytb$Seqs_out != 0),]

inrange_db_f <- function(partout_db, mydb, rangedb){
  inrange_db <- as.data.frame(matrix(ncol = 4))
  colnames(inrange_db) <- c("Accession_number", "Name", "Long", "Lat")
  i <- 1
  for(s in partout_db$IOC_name){
    temp_db <- subset(mydb, mydb$IOC_name == s)
    check_name <- unique(temp_db$CMEC_name)
    temp_range_sp <- subset(rangedb, rangedb$V1 == check_name)
    sample_coor <- matrix(as.numeric(cbind(temp_db$longitude, temp_db$latitude)), ncol = 2)
    range_coor <- cbind(temp_range_sp$V4, temp_range_sp$V5)
    cells_range <- cellFromXY(empty_raster, range_coor)
    cells_samples <- cellFromXY(empty_raster, sample_coor)
    oor <- sample_coor[-which(is.na(match(cells_samples, unique(cells_range)))),]
    if(is.vector(oor) == TRUE){
      ir_temp_db <- temp_db[which(temp_db$longitude == oor[1]), ]
      inrange_db[i,] <- c(ir_temp_db[1], s, oor[1], oor[2])
      i <- i+1
    }else{
      for(c in 1:nrow(oor)){
        match_var <- which(temp_db$longitude == oor[c,1])
        ir_temp_db <- temp_db[match_var,]
        if(length(match_var) >= 2){
          for(r in 1:nrow(ir_temp_db)){
            inrange_db[i,] <- c(ir_temp_db[r,1], s, ir_temp_db[r,4], ir_temp_db[r,5])
            i <- i+1
          }
        }else{  
        ir_temp_db <- temp_db[which(temp_db$longitude == oor[c,1]), ]
        inrange_db[i,] <- c(ir_temp_db[c,1], s, oor[c,1], oor[c,2])
        i <- i+1
        print(c)
        }
      }
    }
  }
  return(inrange_db)
}

inrange_cytb <- inrange_db_f(partout_cytb, cytb_subset, b_ranges)
final_inrange_cytb <- rbind(ir_cytb, inrange_cytb)
final_inrange_cytb <- final_inrange_cytb[order(final_inrange_cytb$Name),]
write.csv(final_inrange_cytb, file = "/Users/Elisabetta/Documents/UCPH/Thesis/Data/cytb_inside_range.csv")

######### CO1
#database containing species of which all sequences fall within the range
noout_co1 <- co1_present_name[which(co1_present_name$Seqs_out == 0),]
ir_co1 <- co1_subset[which(co1_subset$IOC_name %in% noout_co1$IOC_name),]
ir_co1 <- ir_co1[,-3]
colnames(ir_co1) <- c("Accession_number", "Name", "Long", "Lat")

#database containing species which have part of the total sequences falling inside the range
partout_co1 <- co1_present_name[which(co1_present_name$Total_seqs != co1_present_name$Seqs_out),]
partout_co1 <- partout_co1[which(partout_co1$Seqs_out != 0),]

inrange_co1 <- inrange_db_f(partout_co1, co1_subset, b_ranges)
final_inrange_co1 <- rbind(ir_co1, inrange_co1)
final_inrange_co1 <- final_inrange_co1[order(final_inrange_co1$Name),]
write.csv(final_inrange_co1, file = "/Users/Elisabetta/Documents/UCPH/Thesis/Data/co1_inside_range.csv")


##### DATABASE WITH ONLY SEQUENCES OUTSIDE THE RANGE ########
#database containing species of which all sequences fall within the range
library(raster)
empty_raster <- raster()

out_cytb <- cytb_present_name[which(cytb_present_name$Seqs_out >= 1),]
onlyout_cytb <- out_cytb[which(out_cytb$Total_seqs == out_cytb$Seqs_out),]
or_cytb <- cytb_subset[which(cytb_subset$IOC_name %in% onlyout_cytb$IOC_name),]
or_cytb <- or_cytb[,-3]
colnames(or_cytb) <- c("Accession_number", "Name", "Long", "Lat")

#database containing species which have part of the total sequences falling inside the range
partout_cytb <- out_cytb[which(out_cytb$Total_seqs != out_cytb$Seqs_out),]

outrange_db_f <- function(partout_db, mydb, rangedb){
  outrange_db <- as.data.frame(matrix(ncol = 4))
  colnames(outrange_db) <- c("Accession_number", "Name", "Long", "Lat")
  i <- 1
  for(s in partout_db$IOC_name){
    temp_db <- subset(mydb, mydb$IOC_name == s)
    check_name <- unique(temp_db$CMEC_name)
    temp_range_sp <- subset(rangedb, rangedb$V1 == check_name)
    sample_coor <- matrix(as.numeric(cbind(temp_db$longitude, temp_db$latitude)), ncol = 2)
    range_coor <- cbind(temp_range_sp$V4, temp_range_sp$V5)
    cells_range <- cellFromXY(empty_raster, range_coor)
    cells_samples <- cellFromXY(empty_raster, sample_coor)
    oor <- sample_coor[which(is.na(match(cells_samples, unique(cells_range)))),]
    if(is.vector(oor) == TRUE){
      or_temp_db <- temp_db[which(temp_db$latitude == oor[2]), ]
      outrange_db[i,] <- c(or_temp_db[1], s, oor[1], oor[2])
      i <- i+1
    }else{
      #print(oor)
      for(c in 1:nrow(oor)){
        match_var <- which(temp_db$latitude == oor[c,2])
        or_temp_db <- temp_db[match_var,]
        if(length(match_var) >= 2){
          for(r in 1:nrow(or_temp_db)){
            outrange_db[i,] <- c(or_temp_db[r,1], s, or_temp_db[r,4], or_temp_db[r,5])
            i <- i+1
          }
        }else{
        outrange_db[i,] <- c(or_temp_db[c,1], s, oor[c,1], oor[c,2])
        i <- i+1
        print(c)
        }
      }
    }
  }
  return(outrange_db)
}

outrange_cytb <- outrange_db_f(partout_cytb, cytb_subset, b_ranges)
final_outrange_cytb <- rbind(or_cytb, outrange_cytb)
final_outrange_cytb <- final_outrange_cytb[order(final_outrange_cytb$Name),]
write.csv(final_outrange_cytb, file = "/Users/Elisabetta/Documents/UCPH/Thesis/Data/cytb_outside_range.csv")

final_with_latband <- final_outrange_cytb[order(final_outrange_cytb$Lat),]
final_with_latband[,5] <- NA
colnames(final_with_latband) <- c("Accession_number", "Name", "Long", "Lat","Latband")
final_with_latband[which(final_with_latband$Lat <= -70),5] <- "-80 - -70" 
final_with_latband[which(final_with_latband$Lat > -70 & final_with_latband$Lat <= -60), 5] <- "-70 - -60" 
# do the same for the other latbands

# create database to plot frequency in each latband
freq <- as.data.frame(matrix(ncol = 2, nrow = 16))
colnames(freq) <- c("Num_seqs", "Latband")
freq[,2] <- c("-80 - -70", "-70 - - 60", "-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80")
# I had to insert latitudinal bamds manually because one is missing values, so it's not reported in final_with_latband$Latband
freq[1,1] <- nrow(final_with_latband[which(final_with_latband$Latband == "-80 - -70"),])
# do the same for all latbands
par(mar = c(6,6,4,4))
b <- barplot(freq$Num_seqs, names.arg = freq$Latband, horiz = T, las = 1, xlim = c(0,1000), xlab = "Number of sequences")
title("Number of sequences outside the breeding range per latitudinal band", sub = "(Cytochrome-b)")
text(freq$Num_seqs + 20, b, labels = freq$Num_seqs)

###### CO1 #######

out_co1 <- co1_present_name[which(co1_present_name$Seqs_out >= 1),]
onlyout_co1 <- out_co1[which(out_co1$Total_seqs == out_co1$Seqs_out),]
or_co1 <- co1_subset[which(co1_subset$IOC_name %in% onlyout_co1$IOC_name),]
or_co1 <- or_co1[,-3]
colnames(or_co1) <- c("Accession_number", "Name", "Long", "Lat")

partout_co1 <- out_co1[which(out_co1$Total_seqs != out_co1$Seqs_out),]
outrange_co1 <- outrange_db_f(partout_co1, co1_subset, b_ranges)
# during the compilation it gave me an error because JN801425 doesn't have any geographic coordinates so I removed the row
final_outrange_co1 <- rbind(or_co1, outrange_co1)
final_outrange_co1 <- final_outrange_co1[order(final_outrange_co1$Name),]
write.csv(final_outrange_co1, file = "/Users/Elisabetta/Documents/UCPH/Thesis/Data/co1_outside_range.csv")

# add latband column
final_outrange_co1$Lat <- as.numeric(final_outrange_co1$Lat)
final_outrange_co1$Long <- as.numeric(final_outrange_co1$Long)
with_latband <- final_outrange_co1[order(final_outrange_co1$Lat),]
with_latband[,5] <- NA
colnames(with_latband) <- c("Accession_number", "Name", "Long", "Lat","Latband")
with_latband[which(with_latband$Lat <= -60), 5] <- "-70 - -60"
with_latband[which(with_latband$Lat > -60 & with_latband$Lat <= -50), 5] <- "-60 - -50" 
# frequency database
freq_co1 <- as.data.frame(matrix(ncol=2, nrow = 17))
colnames(freq_co1) <- c("Num_seqs", "Latband")
freq_co1[,2] <- c("-80 - -70", "-70 - - 60", "-60 - -50", "-50 - -40", "-40 - -30", "-30 - -20", "-20 - -10", "-10 - 0", "0 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "60 - 70", "70 - 80", "80 - 90")
freq_co1[1,1] <- nrow(with_latband[which(with_latband$Latband == "-80 - -70"),])
par(mar = c(6,6,4,4))
b2 <- barplot(freq_co1$Num_seqs, names.arg = freq_co1$Latband, horiz = T, las = 1, xlim = c(0,2500), xlab = "Number of sequences")
title("Number of sequences outside the breeding range per latitudinal band", sub = "(CO1)")
text(freq_co1$Num_seqs + 70, b2, labels = freq_co1$Num_seqs)
