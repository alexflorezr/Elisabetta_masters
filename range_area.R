findarea <- function(db){
  larus <- subset(db, db$species == "Larus dominicanus")
  for(l in larus$lat){
    n <- abs(l) + 0.5
    ar <- c(careas[n])
  }
  return(colSums(ar))
}

test <- function(larus){
  minidb <- head(larus)
  for(l in minidb$lat){
    n <- abs(l) + 0.5
    ar <- c(careas[n])
  }
  return(sum(ar))
}

test(larus)
