##### CREATE DATASETS WITH ONLY SPECIES WITH LESS THAN 5 SEQUENCES PER GRID CELL ####

### CYTB ###
cytb_matlab_grid <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/cytb_inside_range/matlab_grid_result.csv", 
                             sep = "\t", header = F)
colnames(cytb_matlab_grid) <- c("Sp_name", "Grid_ID", "grid", "Num_seqs", "Nuc_div", "Num_mut_bp", "Tot_bp")

less_than_5_cytb <- cytb_matlab_grid[which(cytb_matlab_grid$Num_seqs <5),]
write.csv(less_than_5_cytb, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/cytb_inside_range/cytb_zoogeo.csv", 
          col.names = F, row.names = F, sep = "\t", quote = F)

### CO1 ###
co1_matlab_grid <- read.csv("/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/co1_inside_range/matlab_grid_result.csv", 
                             sep = "\t", header = F)
colnames(co1_matlab_grid) <- c("Sp_name", "Grid_ID", "grid", "Num_seqs", "Nuc_div", "Num_mut_bp", "Tot_bp")

less_than_5_co1 <- co1_matlab_grid[which(co1_matlab_grid$Num_seqs <5),]
write.csv(less_than_5_co1, "/Users/Elisabetta/Documents/UCPH/Thesis/Data/GD calculation/co1_inside_range/cytb_zoogeo.csv", 
          col.names = F, row.names = F, sep = "\t", quote = F)
