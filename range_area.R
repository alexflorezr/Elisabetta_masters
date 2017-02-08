db <- read.delim(file.choose(), header=F, sep = ",")
colnames(db) <- c("species", "C1", "C2", "long", "lat")
# create an empty data frame 
all_species <- unique(db$species)
db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
colnames(db_area) <- c("species", "range_area")
db_area$species <- all_species
#this function only works if the latitude values that you have in your data are in the form of "0.5"
findarea <- function(db, sp_names){
        all_species <- unique(sp_names)
        db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
        colnames(db_area) <- c("species", "range_area")
        db_area$species <- all_species
        for (sp in sp_names){
                temp_sp <- subset(db, db$species == sp)
                ar <- 0
                for(l in temp_sp$lat){
                        n <- abs(l) + 0.5
                        ar <- ar + as.numeric(careas[n])
                }
                db_area[which(db_area$species == sp), 2] <- ar
        }
        return(db_area)
}
kk <- findarea(db, unique(db$species))


kk$species[which(kk$range_area == min(kk$range_area))]

test(larus)
str(table(larus6$lat))
x <- abs(as.numeric(unlist(dimnames(table(larus6$lat)))))
y <- as.numeric(unlist(dimnames(cellareas)))
sum(table(larus6$lat)*cellareas[match(x,y)])

sp <- "Larus dominicanus"

kk3 <- kk[order(kk$species),]
plot(kk3$range_area, kk2$Freq)
tt <- lm(kk2$Freq~kk3$range_area)
abline(tt, col="red")
dev.off()


############# CYTOCHROME-B SPECIES RANGE AREA #####################
library(readxl)
cytb_db <- read_excel(file.choose())
cytb_db <- subset(cytb_db[-1,5:8]) # I remove the first row which is part of the header and doesn't contain data
colnames(cytb_db) <- c("spname", "location", "lat", "long")
table(is.na(small_db$spname))
which(cytb_db$lat == "NA")
sdb <- cytb_db[ ! small_db$spname %in% "NA", ] # I remove the rows that have the string "NA" in the names column 
which(is.na(sdb))
nrow(sdb)
# it gives me NA in positions that exceed the number of rows in my database
sdb2 <- sdb[-c(134530,134531,143539,143555),] # I remove anyway the rows that it tells me they have missing values
sdb_with_coord <- sdb2[ ! sdb2$lat %in% "NA", c(1,3,4)]
# I remove the rows that have the string "NA" in the lat column
which(is.na(sdb_with_coord)) # I still have NA
nrow(sdb_with_coord) 
# The number of rows in my dataframe is 7611 
# however it tells me that there are NAs in the position 18716, which shouldn't exist

# my data don't have latitude values in decimals every 0.5 so I need to change the function
findarea2 <- function(db, sp_names){
  all_species <- unique(sp_names)
  db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
  colnames(db_area) <- c("species", "range_area")
  db_area$species <- all_species
  for (sp in sp_names){
    temp_sp <- subset(db, db$spname == sp)
    ar <- 0
    for(l in temp_sp$lat){
      n <- floor(abs(as.numeric(l))) + 0.5
      ar <- ar + as.numeric(careas[n])
    }
    db_area[which(db_area$species == sp), 2] <- ar
  }
  return(db_area)
}
a1 <- findarea2(sdb_with_coord, unique(sdb_with_coord$spname))
# run the function to find the area of the species that have geographic coordinates
# it gives me error:
#     Error in `[<-.data.frame`(`*tmp*`, which(db_area$species == sp), 2, value = numeric(0)) : 
#     replacement has length zero
# I think it's because of NAs in my database which I don't know how to remove.. it's weird that it gives me NAs in
# positions that don't exist (or shouldn't) in sdb_with_coord

############## END HERE #################
############## THE PART THAT FOLLOWS IS BASED ON WRONG CALCULATIONS ###########################

which(is.na(a1$range_area))
table(is.na(a1$range_area))
NAs <- a1[which(is.na(a1$range_area)),]
which(NAs$species == small_db$spname)

# therefore I need to remove the speices with NA as range size value
a2 <- a1[complete.cases(a1),]
a2$species[which(a2$range_area == max(a2$range_area))]
a2$species[which(a2$range_area == min(a2$range_area))]

pdf("Species range areas cytb.pdf")
hist(a2$range_area, xlab = "Species range area", main = "Cytochrome-b species area size")
dev.off()

# create a dataset with the species in first column and n of sequences in 2nd column (frequency)
nseq <- as.data.frame(table(sdb$spname))
nseq <- (nseq[-1,]) # I delete the first row because of a wrong species name

sort_a2 <- a2[order(a2$species),]

plot(a2$range_area, nseq$Freq)



str(table(sdb2$lat))
x <- floor(abs(as.numeric(unlist(dimnames(table(sdb$lat)))))) + 0.5
y <- as.numeric(unlist(dimnames(cellareas)))
sum(table(sdb2$lat)*cellareas[match(x,y)])




