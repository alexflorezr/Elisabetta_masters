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
        ir_temp_db <- temp_db[which(temp_db$longitude == oor[c,1]), ]
        inrange_db[i,] <- c(ir_temp_db[1], s, oor[c,1], oor[c,2])
        i <- i+1
        print(c)
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
