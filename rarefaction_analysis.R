################ RAREFACTION ANALYSIS ##############

rarefaction_f <- function(db, species, gridcells){
  N <- 100
  temp_db <- as.data.frame(matrix(ncol = 3))
  colnames(temp_db) <- c("species", "grid", "avg5")
  inter_db <- as.data.frame(matrix(ncol = 100, nrow = length(gridcells)))
  rownames(inter_db) <- gridcells
  for(n in 1:N){
    temp_db <- as.data.frame(matrix(ncol = 3))
    colnames(temp_db) <- c("species", "grid", "avg5")
    i <- 1
    for(s in species){
      #subset the database for each species
      subdb.sp <- subset(db, db$species == s)
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
  colnames(final_db) <- c("Grid_ID", "avg", "stdev")
  i <- 1
  for(cell in seq_along(gridcells)){
    temp_avg <- mean(as.numeric(inter_db[which(rownames(inter_db) == gridcells[cell]),]))
    temp_stdev <- sd(as.numeric(inter_db[which(rownames(inter_db) == gridcells[cell]),]))
    final_db[i,] <- c(gridcells[cell], temp_avg, temp_stdev)
    i <- i+1 
  }
  return(final_db)
}


#####################################################################################
## CYTB
# load the file
# the file already contains values for only species that have >5 seqs per grid cell
# if a species in a grid cell has <5 seqs, the values for that grid cell is NaN
library(readxl)

cytb <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Sensitivity analysis/CYTB/matlab_rarefaction_cytb.csv", header = F)
colnames(cytb) <- c("species", "grid_id", "seq1", "seq2", "length_seq1", "length_seq2", "overlap", "commons", "num_mut_bp")
# remove the NaN, which removes also sequences that overlap <0.5
cytb <- cytb[-which(is.nan(cytb$num_mut_bp)),]

sp_names <- unique(cytb$species)
grid_cell <- unique(cytb$grid_id)
cytb_rarefaction <- rarefaction_f(cytb, sp_names, grid_cell)
write.csv(cytb_rarefaction, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/Sensitivity analysis/CYTB/rarefaction_result_cytb.csv",
          row.names = F, quote = F)

cytb_gd <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Sensitivity analysis/CYTB/cytb_grid_more5seqs.csv")
cytb_merged <- merge(cytb_gd, cytb_rarefaction, by = "Grid_ID")

# correlation with average
plot(cytb_merged$avg, cytb_merged$GD, xlab = "Average GD from rarefaction analysis", ylab = "GD", main = "Cytochrome-b")
cytb_cor.avg <- cor.test(cytb_merged$avg, cytb_merged$GD)
r_square.avg <- cytb_cor.avg$estimate^2
cytb_lm.avg <- lm(cytb_merged$GD ~ cytb_merged$avg)
abline(cytb_lm.avg, col = "red")

# correlation with standard deviation
plot(cytb_merged$stdev, cytb_merged$GD, xlab = "Standard Deviation GD from rarefaction analysis", ylab = "GD", main = "Cytochrome-b")
cytb_cor.std <- cor.test(cytb_merged$stdev, cytb_merged$GD)
r_square.std <- cytb_cor.std$estimate^2
cytb_lm.std <- lm(cytb_merged$GD ~ cytb_merged$stdev)
abline(cytb_lm.std, col = "red")

####################################################################################
## CO1
library(readxl)

co1 <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Sensitivity analysis/CO1/matlab_rarefaction_co1.csv", header = F)
colnames(co1) <- c("species", "grid_id", "seq1", "seq2", "length_seq1", "length_seq2", "overlap", "commons", "num_mut_bp")
co1 <- co1[-which(is.nan(co1$num_mut_bp)),]

spname <- unique(co1$species)
grids <- unique(co1$grid_id)

co1_rarefaction <- rarefaction_f(co1, spname, grids)
write.csv(co1_rarefaction, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/Sensitivity analysis/CO1/rarefaction_result_co1.csv",
          row.names = F, quote = F)

co1_gd <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_grid_insiderange.xlsx")
colnames(co1_gd) <- c("Grid_ID","Num_sps", "Num_seqs", "Sum_nuc_div", "Tot_mut_bp", "GD")
co1_merged <- merge(co1_gd, co1_rarefaction, by = "Grid_ID")

# correlation with average
plot(co1_merged$avg, co1_merged$GD, xlab = "Average GD from rarefaction analysis", ylab = "GD", main = "CO1")
co1_cor.avg <- cor.test(co1_merged$avg, co1_merged$GD)
r_square.avg <- co1_cor.avg$estimate^2
co1_lm.avg <- lm(co1_merged$GD ~ co1_merged$avg)
abline(co1_lm.avg, col = "red")

# correlation with standard deviation
plot(co1_merged$stdev, co1_merged$GD, xlab = "Standard Deviation of GD from rarefaction analysis", ylab = "GD", main = "CO1")
co1_cor.std <- cor.test(co1_merged$stdev, co1_merged$GD)
r_square.std <- co1_cor.std$estimate^2
co1_lm.std <- lm(co1_merged$GD ~ co1_merged$stdev)
abline(co1_lm.std, col = "red")
