############### CALCULATING AVERAGE MINIMUM DISTANCE OF SEQUENCES OUT OF RANGE PER SPECIES ##############
# cytb
more_than_one, without_NA, b_ranges # starting databases

without_NA[ ,3] <- as.numeric(without_NA$latitude)
without_NA[ ,4] <- as.numeric(without_NA$longitude)
cytb_dist_db <- as.data.frame(matrix(nrow = 3253, ncol = 4))
colnames(cytb_dist_db) <- c("Name", "Long", "Lat", "min_dist")
i <- 1
for(s in more_than_one$IOC_name){
  temp_db <- subset(without_NA, without_NA$IOC_name == s)
  check_name <- unique(temp_db$CMEC_name)
  temp_range_sp <- subset(b_ranges, b_ranges$V1 == check_name)
  sample_coor <- cbind(temp_db$longitude, temp_db$latitude)
  range_coor <- cbind(temp_range_sp$V4, temp_range_sp$V5)
  cells_range <- cellFromXY(empty_raster, range_coor)
  cells_samples <- cellFromXY(empty_raster, sample_coor)
  oor <- sample_coor[which(is.na(match(cells_samples, unique(cells_range)))),]
  dst <- c()
  if(is.vector(oor) == TRUE){
    d <- min(pointDistance(oor, range_coor, lonlat = T))/1000
    #dist_db[s,] <- c(s, oor[1], oor[2], d)
    cytb_dist_db[i,] <- c(s, oor[1], oor[2], d)
    i <- i+1
  }else{
    for(c in 1:nrow(oor)){
      d <- min(pointDistance(oor[c,], range_coor, lonlat = T))/1000
      #dist_db[s,] <- c(s, oor[c,1], oor[c,2], dst[c])
      cytb_dist_db[i,] <- c(s, oor[c,1], oor[c,2], d)
      i <- i+1
      dst <- append(dst, d)
      print(c)
    }
  }
}

cytb_mean_db <- as.data.frame(matrix(nrow = length(unique(cytb_dist_db$Name)), ncol = 2))
colnames(cytb_mean_db) <- c("Name", "Mean_distance")
cytb_mean_db$Name <- unique(cytb_dist_db$Name)
i <- 1
for(s in unique(cytb_dist_db$Name)){
  t_db <- subset(cytb_dist_db, cytb_dist_db$Name == s)
  av <- mean(as.numeric(t_db$min_dist))
  cytb_mean_db[i,2] <- av
  i <- i+1
}

nrow(cytb_mean_db)
hist(cytb_mean_db$Mean_distance)
hist(cytb_mean_db[which(cytb_mean_db$Mean_distance <= 1000),2])

# co1
co1_cells_out, ldb, b_ranges # starting databases

co1_dist_db <- as.data.frame(matrix(nrow = sum(as.numeric(co1_cells_out$Seqs_out)), ncol = 4))
colnames(co1_dist_db) <- c("Name", "Long", "Lat", "min_dist")
i <- 1
for(s in co1_cells_out$IOC_name){
  temp_db <- subset(ldb, ldb$IOC_name == s)
  check_name <- unique(temp_db$CMEC_name)
  temp_range_sp <- subset(b_ranges, b_ranges$V1 == check_name)
  sample_coor <- matrix(as.numeric(cbind(temp_db$longitude, temp_db$latitude)), ncol = 2)
  range_coor <- cbind(temp_range_sp$V4, temp_range_sp$V5)
  cells_range <- cellFromXY(empty_raster, range_coor)
  cells_samples <- cellFromXY(empty_raster, sample_coor)
  oor <- sample_coor[which(is.na(match(cells_samples, unique(cells_range)))),]
  dst <- c()
  if(is.vector(oor) == TRUE){
    d <- min(pointDistance(oor, range_coor, lonlat = T))/1000
    #dist_db[s,] <- c(s, oor[1], oor[2], d)
    co1_dist_db[i,] <- c(s, oor[1], oor[2], d)
    i <- i+1
  }else{
    for(c in 1:nrow(oor)){
      d <- min(pointDistance(oor[c,], range_coor, lonlat = T))/1000
      #dist_db[s,] <- c(s, oor[c,1], oor[c,2], dst[c])
      co1_dist_db[i,] <- c(s, oor[c,1], oor[c,2], d)
      i <- i+1
      dst <- append(dst, d)
      print(c)
    }
  }
}

nrow(co1_dist_db[which(co1_dist_db$min_dist <= 400),])

co1_mean_db <- as.data.frame(matrix(nrow = length(unique(co1_dist_db$Name)), ncol = 2))
colnames(co1_mean_db) <- c("Name", "Mean_distance")
co1_mean_db$Name <- unique(co1_dist_db$Name)
i <- 1
for(s in unique(co1_dist_db$Name)){
  t_db <- subset(co1_dist_db, co1_dist_db$Name == s)
  av <- mean(as.numeric(t_db$min_dist))
  co1_mean_db[i,2] <- av
  i <- i+1
}

nrow(co1_mean_db)
hist(co1_mean_db$Mean_distance)
nrow(co1_mean_db[which(co1_mean_db$Mean_distance <= 1000),])
hist(co1_mean_db[which(co1_mean_db$Mean_distance <= 1000),2])
nrow(co1_mean_db[which(co1_mean_db$Mean_distance <= 400),])
hist(co1_mean_db[which(co1_mean_db$Mean_distance <= 400),2])

## END HERE ##
########################################