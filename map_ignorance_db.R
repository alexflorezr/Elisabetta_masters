######### CREATION OF DATABASE FOR MAP OF IGNORANCE ######
library(readxl)
range_grid <- read.csv("/Users/Elisabetta/Downloads/ignorance map/num_species_grid.csv", header = F, sep = "\t")
cytb_grid <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CYTB/cytb_grid_insiderange.xlsx")
cytb_grid <- cytb_grid[,1:2]
colnames(range_grid) <- c("Grid_ID", "Num_species")

grid_new <- merge(range_grid, cytb_grid, by = "Grid_ID")
colnames(grid_new) <- c("Grid_ID", "Num_species_range", "Num_species_cytb")

final_grid <- grid_new
final_grid$Pct <- with(final_grid, 100/(Num_species_range/Num_species_cytb))
write.csv(final_grid, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/Map of Ignorance/num_species_grid_cytb.csv", quote = F, row.names = F)

##################################################
co1_grid <- read_excel("/Users/Elisabetta/Documents/UCPH/Thesis/Data/Results/CO1/co1_grid_insiderange.xlsx")
co1_grid <- co1_grid[,1:2]
colnames(co1_grid) <- c("Grid_ID", "Num_species")

grid_new2 <- merge(range_grid, co1_grid, by = "Grid_ID")
colnames(grid_new2) <- c("Grid_ID", "Num_species_range", "Num_species_co1")

grid_new2$Pct <- with(grid_new2, 100/(Num_species_range/Num_species_co1))
write.csv(grid_new2, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/Map of Ignorance/num_species_grid_co1.csv", quote = F, row.names = F)


