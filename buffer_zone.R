# CALCULATING NUMBER OF SEQUENCES THAT FALL OUTSIDE THE RANGE WITH BUFFER ZONE OF 1 DEGREE #
library(readxl)
library(raster)

sp_names <- unique(co1_subset$IOC_name)
co1_out_of_new_range <- as.data.frame(matrix(nrow=length(sp_names), ncol=6))
colnames(co1_out_of_new_range) <- c("IOC_name", "In_ranges_file", "Cells_out","Total_range_cells", "Seqs_out", "Total_seqs")
for(s in seq_along(sp_names)){
  temp_db <- subset(co1_subset, co1_subset$IOC_name == sp_names[s])
  to_check_name <- unique(temp_db$CMEC_name) 
  if(is.element(to_check_name, b_ranges$V1)){
    temp_range_sp <- subset(b_ranges, b_ranges$V1 == to_check_name)
    cells_range <- cellFromXY(empty_raster, cbind(temp_range_sp$V4, temp_range_sp$V5))
    cells_samples <- cellFromXY(empty_raster, cbind(as.numeric(temp_db$longitude), as.numeric(temp_db$latitude)))
    cnear <- adjacent(empty_raster, cells_range, directions = 4, pairs = F)
    new_cells_range <- c(cells_range, cnear)
    how_many_seqs <- sum(is.na(match(cells_samples, unique(new_cells_range))))
    how_many_cells <- sum(is.na(match(unique(cells_samples), unique(new_cells_range))))
    co1_out_of_new_range[s,] <- c(sp_names[s], "YES", how_many_cells, length(unique(new_cells_range)), how_many_seqs, dim(temp_db)[1])
  } else {
    co1_out_of_new_range[s,] <- c(sp_names[s], "NO",NA, NA, NA,dim(temp_db)[1])
  }
  print(s)
}

co1_oonr_present_name <- subset(co1_out_of_new_range, co1_out_of_new_range$In_ranges_file == "YES")
oonr_seqs_co1 <- sum(as.numeric(co1_oonr_present_name$Seqs_out))
oonr_missing_name <- subset(out_of_new_range, out_of_new_range$In_ranges_file == "NO")
sum(as.numeric(oonr_missing_name$Total_seqs))