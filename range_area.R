<<<<<<< HEAD
# read the file with the species ranges (in this case: WorldBreedingBirdsRanges_CMEC_2016-09-12.txt)
db_ranges <- read.delim(file.choose(), header=F, sep = ",", stringsAsFactors = F)
        colnames(db_ranges) <- c("species", "C1", "C2", "long", "lat")
# read the sequences database (in this case: Cytb_database.txt )
db_seqs <- read.delim(file.choose(), header=T, sep = "\t", stringsAsFactors = F)
        # check if there are NAs in the species names and latitude
        ifelse(sum(is.na(db_seqs$CMEC.Database.sp..Name)) >= 1, paste("There are", sum(is.na(db_seqs$IOC.sp.NAme)), "NAs in the names", sep = " "), "There are NOT Nas in the names")
        ifelse(sum(is.na(db_seqs$Geonames_Latitude)) >= 1, paste("There are", sum(is.na(db_seqs$Geonames_Latitude)), "NAs in the names", sep = " "), "There are NOT Nas in the names")
        # remove sp == NA and latitude == NA
        sdb_no_na_spname <- db_seqs[!is.na(db_seqs$CMEC.Database.sp..Name),]
        sdb_no_na_lat <- sdb_no_na_spname[!is.na(sdb_no_na_spname$Geonames_Latitude),]

=======
############################ RANGE SIZE CMEC DATABASE #############################

db <- read.delim(file.choose(), header=F, sep = ",")
colnames(db) <- c("species", "C1", "C2", "long", "lat")
# create an empty data frame 
all_species <- unique(db$species)
db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
colnames(db_area) <- c("species", "range_area")
db_area$species <- all_species
>>>>>>> 646e7a8316d61df8adb9138c0056f93600d166e5
#this function only works if the latitude values that you have in your data are in the form of "0.5"
findarea <- function(db, sp_names){
        all_species <- unique(sp_names)
        db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
        colnames(db_area) <- c("species", "range_area")
        db_area$species <- all_species
        for (sp in all_species){
                temp_sp <- subset(db, db$species == sp)
                db_area[db_area$species == sp, 2] <- sum(careas[abs(floor(temp_sp$lat)) + 0.5])
        }
        return(db_area)
}
<<<<<<< HEAD
species_area <- findarea(db_ranges, unique(sdb_no_na_lat$CMEC.Database.sp..Name))
# check that there is no NAs in the range area column
sum(is.na(species_area$range_area))
# species with estimates of area == zero
# I think is poecile montanus, not oecile montanus
# Be careful with the names you are including, neither "UNPUBLISHED NAME" nor "EXTINCT" should be in species names
species_area$species[which(species_area$range_area == 0)]
=======
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
cytb_db <- subset(cytb_db[-1, c(5,7,8)]) # I remove the first row which is part of the header and doesn't contain data
colnames(cytb_db) <- c("spname", "lat", "long")
#table(is.na(cytb_db$spname))
#which(cytb_db$spname == "NA")
#which(cytb_db$lat == "NA")
sdb <- cytb_db[-which(cytb_db$spname == "NA"), ] # I remove the rows that have the string "NA" in the "names" column
sdb2 <- sdb[-which(sdb$lat == "NA"), ] # I remove the rows that have the string "NA" in the "lat" column 
#which(is.na(sdb2))
#nrow(sdb2)
which(sdb2$spname == "EXTINCT")
which(sdb2$spname == "UNPUBLISHED NAME")
sdb4 <- sdb2[-c(530,531,2393,2394,2450,2451,4227), ] 
# I remove rows that have under the column "names" the strings
# "EXTINCT" and "UNPUBLISHED NAME"
nrow(sdb4)

# it gives me NA in positions that exceed the number of rows in my database
sdb3 <- sdb2[-c(18714,18715,20330,20343),] # I remove anyway the rows that it tells me they have missing values
which(is.na(sdb3)) # I still have NA
nrow(sdb3) 
# The number of rows in my dataframe is 7610 
# however it tells me that there are NAs in the position 18714, which shouldn't exist

# since I still don't know why it does like this, I decided to create a new table with set number of rows
test <- as.data.frame(matrix(nrow = 7602, ncol = 2))
colnames(test) <- c("sp_name", "lat")
test$sp_name <- sdb4$spname
test$lat <- sdb4$lat
which(is.na(test))

# my data don't have latitude values in decimals every 0.5 so I need to change the function
#findarea2 <- function(db, sp_names){
#  all_species <- unique(sp_names)
#  db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
#  colnames(db_area) <- c("species", "range_area")
#  db_area$species <- all_species
#  for (sp in sp_names){
#    temp_sp <- subset(db, db$sp_name == sp)
#    ar <- 0
#    for(l in temp_sp$lat){
#      n <- floor(abs(as.numeric(l))) + 0.5
#      ar <- ar + as.numeric(careas[n])
#    }
#    db_area[which(db_area$species == sp), 2] <- ar
#  }
#  return(db_area)
#}
#a1 <- findarea2(test, unique(test$sp_name))
# run the function to find the area of the species that have geographic coordinates
# it gives me error:
#     Error in `[<-.data.frame`(`*tmp*`, which(db_area$species == sp), 2, value = numeric(0)) : 
#     replacement has length zero
# I think it's because of NAs in my database which I don't know how to remove.. it's weird that it gives me NAs in
# positions that don't exist (or shouldn't) in sdb_with_coord


# try with another function
#findarea3 <- function(db, sp_names){
#  all_species <- unique(sp_names)
#  db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
#  colnames(db_area) <- c("species", "range_area")
#  db_area$species <- all_species
#  for (sp in sp_names){
#    temp_sp <- subset(db, db$sp_name == sp)
#    ar <- 0
#    for(l in temp_sp$lat){
#      x <- floor(abs(as.numeric(unlist(dimnames(table(temp_sp$lat)))))) + 0.5
#      y <- as.numeric(unlist(dimnames(cellareas)))
#      ar <- sum(table(temp_sp$lat)*cellareas[match(x,y)])
#    }
#    db_area[which(db_area$species == sp), 2] <- ar
#  }
#  return(db_area)
#}

#result_db <- findarea3(test, unique(test$sp_name))
# IT WORKS!!!!!

# improved function for unique cells (final one)
findarea4 <- function(db, sp_names){
  all_species <- unique(sp_names)
  db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
  colnames(db_area) <- c("species", "range_area")
  db_area$species <- all_species
  for (sp in sp_names){
    temp_sp <- subset(db, db$sp_name == sp)
    temp_lat <- unique(floor(as.numeric(temp_sp$lat)))
    ar <- 0
    for(l in temp_lat){
      x <- abs(temp_lat) + 0.5
      y <- as.numeric(unlist(dimnames(cellareas)))
      ar <- sum(table(temp_lat)*cellareas[match(x,y)])
    }
    db_area[which(db_area$species == sp), 2] <- ar
  }
  return(db_area)
}

result_db <- findarea4(test, unique(test$sp_name))


table(is.na(result_db)) 
result_db[which(result_db$range_area == max(result_db$range_area)), ]
result_db[which(result_db$range_area == min(result_db$range_area)), ]

pdf("cytb_sp_range_size.pdf")
hist(result_db$range_area, xlab = "Species range area", main = "Cytochrome-b species area size")
dev.off()
pdf("cytb_sp_range_size_high_def.pdf")
hist(result_db$range_area, breaks = seq(0, 245000, by = 1000), main = "", xlab = "Species range area")
mtext("Cytb species range size", line = 2, cex = 1.2, font = 2)
mtext("breaks by 1000", line = 0.5)
dev.off()

table(test$sp_name)
nseq <- as.data.frame(table(test$sp_name))
result_sorted <- result_db[order(result_db$species),]
cor(result_sorted$range_area, nseq$Freq)
tt <- lm(nseq$Freq ~ result_sorted$range_area)
pdf("Range_size_nseq.pdf")
plot(result_sorted$range_area, nseq$Freq, xlab = "Species range area", ylab = "Number of sequences", cex.lab = 0.9)
mtext("Range size over sequence availability", line = 1.7, cex = 1, font = 2)
abline(tt, col="red")
dev.off()

pdf("Logarithmic_scale_range_nseq.pdf")
ltt <- lm(log(nseq$Freq) ~ log(result_sorted$range_area))
plot(log(result_sorted$range_area), log(nseq$Freq), xlab = "log(range area)", ylab = "log(n° of sequences)", cex.lab = 0.9)
abline(ltt, col="red")
dev.off()
############## END HERE #################
>>>>>>> 646e7a8316d61df8adb9138c0056f93600d166e5
