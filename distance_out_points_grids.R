############ CALCULATE DISTANCE FROM POINTS TO GRID CELL ##########

# this is a test with just one species of my cytb database

library(raster)
empty_raster <- raster()
# subset the cytb database containing species name and coordinates that are not NA
test <- without_NA[which(without_NA$IOC_name == "Phylloscopus humei"), ]
subtest <- subset(b_ranges, b_ranges$V1 == "Phylloscopus humei")

library(rworldmap)
world <- getMap("world")
# plot the points to see where it falls out the range
plot(world)
points(test$longitude, test$latitude, col = "red") 
points(subtest$V4, subtest$V5, col = "blue")

# calculate the cell numbers for the sequences in my database
cells_range <- cellFromXY(empty_raster, cbind(subtest$V4, subtest$V5))
cells_samples <- cellFromXY(empty_raster, cbind(as.numeric(test$longitude), as.numeric(test$latitude)))
seqs_out <- match(cells_samples, unique(cells_range))

# fill a new empty raster with values
# bird range points have value 1
# the cell number of the sequence that falls outside the range gets value NA
r <- raster()
r[cells_range] <- 1
r[cells_samples[which(is.na(match(cells_samples, unique(cells_range))))]] <- NA
plot(r)

# calculate distance with distance() function
dist <- distance(r)
plot(dist/1000)
# I'm getting a raster file not a measure of distance
# And also all the cells that don't have value 1 have value NA so distance() calculates the distance also from these
# and not only from the cell out of range in my database

# calculate distance with pointDistance() function
NA_coor <- cbind(test$longitude[3], test$latitude[3])
range_coor <- xyFromCell(empty_raster, cells_range)
pointDistance(NA_coor, range_coor, lonlat = T)
min(pointDistance(NA_coor, range_coor, lonlat = T)) # 803770.7
#in km
min(pointDistance(NA_coor, range_coor, lonlat = T))/1000 # 803.7707

# adding cells adjacent to range ("buffer")
ad <- adjacent(empty_raster, cells_range, pairs = FALSE)
crange_new <- c(cells_range, ad)
r2 <- raster()
r2[crange_new] <- 1
plot(r2)
plot(r)
match(cells_samples, unique(crange_new))



