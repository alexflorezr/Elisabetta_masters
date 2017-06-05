############## CREATION OF DATABASE FOR TAXONOMIC COVERAGE ################
## THE FINAL DATABASE WILL BE IMPORTED IN ARCMAP TO PRODUCE THE FIGURES ##
## CYTB ##
library(readxl)

range_db <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/Map of Ignorance/ranges_gridded.txt", header = T, sep = ";")
range_db <- range_db[,c(3,6)] # keep only column with species name and grid ID
agg <- aggregate(. ~ FID_1, data = range_db, FUN = function(x){length(unique(x))})

colnames(agg) <- c("Grid_ID", "Num_species_range")

cytb_grid <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/cytb_grid_insiderange.xlsx")
cytb_grid <- cytb_grid[,1:2] # keep only column with grid ID and number of species

grid_new <- merge(agg, cytb_grid, by = "Grid_ID")
colnames(grid_new) <- c("Grid_ID", "Num_species_range", "Num_species_cytb")
final_grid <- grid_new
final_grid$Pct <- with(final_grid, 100/(Num_species_range/Num_species_cytb))

write.csv(final_grid, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/Map of Ignorance/taxonomic_coverage_cytb.csv",
          quote = F, row.names = F)

## CO1 ##
co1_grid <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_grid_insiderange.xlsx")
co1_grid <- co1_grid[,1:2] # keep only column with grid ID and number of species
colnames(co1_grid) <- c("Grid_ID", "Num_species")

co1_grid_new <- merge(agg, co1_grid, by = "Grid_ID")
colnames(co1_grid_new) <- c("Grid_ID", "Num_species_range", "Num_species_co1")

co1_grid_new$Pct <- with(co1_grid_new, 100/(Num_species_range/Num_species_co1))
write.csv(co1_grid_new, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/Map of Ignorance/taxonomic_coverage_co12.csv", quote = F, row.names = F)

################## CREATION OF DATABASE FOR TOTAL NUMBER OF BASEPAIRS PER GRID CELL ###################
## CYTB ##

cytb <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/cytb_inside_range/matlab_grid_result.csv", sep = "\t", header = F)
cytb <- cytb[,c(2,7)] # keep only columns with grid ID and number of basepairs
colnames(cytb) <- c("Grid_ID", "Num_bp")
cytb_bp <- aggregate(. ~ Grid_ID, data = cytb, FUN = sum)

write.csv(cytb_bp, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/Map of Ignorance/tot_bp_grid_cytb.csv", row.names = F, quote = F)

## CO1 ##

co1 <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/co1_inside_range/matlab_grid_result.csv", sep = "\t", header = F)
co1 <- co1[,c(2,7)] # keep only columns with grid ID and number of basepairs
colnames(co1) <- c("Grid_ID", "Num_bp")
co1_bp <- aggregate(. ~ Grid_ID, data = co1, FUN = sum)

write.csv(co1_bp, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/Map of Ignorance/tot_bp_grid_co1.csv", row.names = F, quote = F)



