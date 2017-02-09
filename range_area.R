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
species_area <- findarea(db_ranges, unique(sdb_no_na_lat$CMEC.Database.sp..Name))
# check that there is no NAs in the range area column
sum(is.na(species_area$range_area))
# species with estimates of area == zero
# I think is poecile montanus, not oecile montanus
# Be careful with the names you are including, neither "UNPUBLISHED NAME" nor "EXTINCT" should be in species names
species_area$species[which(species_area$range_area == 0)]