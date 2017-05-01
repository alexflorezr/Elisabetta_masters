################ RAREFACTION ANALYSIS ##############

# load the file
# the file already contains values for only species that have >5 seqs per grid cell
# if a species in a grid cell has <5 seqs, the values for that grid cell is NaN

cytb <- read.csv("/Users/Elisabetta/Desktop/rarefaction_cytb.csv", header = F)
colnames(cytb) <- c("species", "grid_id", "seq1", "seq2", "length_seq1", "length_seq2", "overlap", "commons", "num_mut_bp")
# remove the NaN, which removes also sequences that overlap <0.5
cytb <- cytb[-which(is.nan(cytb$num_mut_bp)),]

sp_names <- unique(cytb$species)
grid_cell <- unique(cytb$grid_id)

N <- 100
temp_db <- as.data.frame(matrix(ncol = 3))
colnames(temp_db) <- c("species", "grid", "avg5")
inter_db <- as.data.frame(matrix(ncol = 100, nrow = length(grid_cell)))
rownames(inter_db) <- grid_cell
for(n in 1:N){
  temp_db <- as.data.frame(matrix(ncol = 3))
  colnames(temp_db) <- c("species", "grid", "avg5")
  i <- 1
  for(s in sp_names){
    #subset the database for each species
    subdb.sp <- subset(cytb, cytb$species == s)
    temp_grid <- unique(subdb.sp$grid_id)
    for(g in seq_along(temp_grid)){
      # subset again for each  grid cell
      subdb.sp.cell <- subset(subdb.sp, subdb.sp$grid_id == temp_grid[g])
      # sample 5 random values of number of mutations and do the average --> pi-hat
      avg.5 <- mean(sample(subdb.sp.cell$num_mut_bp, 5))
      # the final database contains for each species and each cell of that species a pi-hat value
      temp_db[i,] <- c(s,temp_grid[g],avg.5)
      i <- i+1
    }
  }
  print(n)
  for(c in seq_along(unique(temp_db$grid))){
    # subset the final database for each grid cell
    subdb.cell <- subset(temp_db, temp_db$grid == unique(temp_db$grid)[c])
    # calculate GD for the cell by summing pi-hats of that cell and dividing by their frequence (number of species)
    gd.cell <- mean(as.numeric(subdb.cell$avg5))
    # the new database contains for each grid cell 100 values of GD
    inter_db[which(rownames(inter_db) == unique(temp_db$grid)[c]),n] <- gd.cell
  }
}

# calculate for each grid cell the mean and standard deviation using the 100 values
final_db <- as.data.frame(matrix(ncol = 3))
colnames(final_db) <- c("grid_ID", "avg", "stdev")
i <- 1
for(cell in seq_along(grid_cell)){
  temp_avg <- mean(as.numeric(inter_db[which(rownames(inter_db) == grid_cell[cell]),]))
  temp_stdev <- sd(as.numeric(inter_db[which(rownames(inter_db) == grid_cell[cell]),]))
  final_db[i,] <- c(grid_cell[cell], temp_avg, temp_stdev)
  i <- i+1 
}

