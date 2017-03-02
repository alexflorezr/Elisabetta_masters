############### USING AVERAGE ERROR DISTANCE ON SEQUENCES OUT OF RANGE ##############

more_than_one, without_NA, b_ranges # starting databases

# I use a test database that has only coords for Phylloscopus humei
test <- subset(without_NA, without_NA$IOC_name == "Phylloscopus humei")
test$latitude <- as.numeric(test$latitude)
test$longitude <- as.numeric(test$longitude)

for(s in seq_along(more_than_one$IOC_name)){
  temp_db <- subset(without_NA, without_NA$IOC_name == s)
  check_name <- unique(temp_db$CMEC_name)
  temp_range_sp <- subset(b_ranges, b_ranges$V1 == check_name)
  cells_range <- cellFromXY(empty_raster, cbind(temp_range_sp$V4, temp_range_sp$V5))
  cells_samples <- cellFromXY(empty_raster, cbind(temp_db$longitude, temp_db$latitude))
  coor_new <- temp_db[which(is.na(match(cells_samples, unique(cells_range)))), c(3,4)] + m
  cell_new <- cellFromXY(empty_raster, cbind(coor_new[2], coor_new[1]))
  r_new <- match(cell_new, unique(cells_range))
  if(r_new != NA){
    print("The sequence has new coordinates!")
  }else{
    coor_new <- temp_db[which(is.na(match(cells_samples, unique(cells_range)))), c(3,4)] - m
    cell_new <- cellFromXY(empty_raster, cbind(coor_new[2], coor_new[1]))
    mtc <- match(cell_new, unique(cells_range))
    print(mtc)
  }
}
