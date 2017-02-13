################# CALCULATION OF COVERAGE(%) OF SPECIES TOTAL RANGES ##################

#db_test <- head(result_db, n = 10)

range_coverage <- function(db1, db2){
  db_names <- db1$species
  mydb_names <- db2$species
  db_percentage <- as.data.frame(matrix(nrow = length(mydb_names), ncol = 2))
  colnames(db_percentage) <- c("species", "coverage_pct")
  db_percentage$species <- mydb_names
  for(sp in mydb_names){
    m <- match(sp, db_names)
    f <- db2[which(db2$species == sp), 2]
    temp_p <- f*100/db1[m,2]
    db_percentage[which(db_percentage$species == sp), 2] <- temp_p
  }
  return(db_percentage)
}

db_pct <- range_coverage(kk, result_db)

#db_pct[which(db_pct$coverage_pct == max(db_pct$coverage_pct)), 2]
# Some species have a coverage percentage that exceeds the 100%. This is probably due to extensive sampling that also
# may come from a single place. Since the function to find the area simply adds the cells relative to latitude, some
# cell areas may be summed multiple times, increasing the total number.

hist(db_pct$coverage_pct)
hist(db_pct$coverage_pct, breaks = seq(0, 19000, by = 100))
hist(db_pct$coverage_pct, breaks = seq(0, 19000, by = 10))



