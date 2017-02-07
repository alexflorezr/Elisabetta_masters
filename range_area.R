db <- read.delim(file.choose(), header=F, sep = ",")
colnames(db) <- c("species", "C1", "C2", "long", "lat")
# create an empty data frame 
all_species <- unique(db$species)
db_area <- data.frame(matrix(nrow=length(all_species), ncol=2))
colnames(db_area) <- c("species", "range_area")
db_area$species <- all_species
findarea <- function(db, sp_names){
        all_species <- unique(db$species)
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


############# CYTOCHROME-B SPECIES RANGE AREA #####################
library(readxl)
cytb_db <- read_excel(file.choose())
small_db <- subset(cytb_db[-1,5:8])
colnames(small_db) <- c("spname", "location", "lat", "long")
sdb <- small_db[ ! small_db$spname %in% "NA", ]
# I create another database that has no "NA" as values in the column spname

a1 <- findarea(db, unique(sdb$spname))
# some of the species have missing results (NA) which is strange since I'm sure the names in my database match the ones of the ranges
# therefore I need to remove the speices with NA as range size value
# CONTINUE HERE


a1$species[which(a1$range_area == max(a1$range_area))]
a1$species[which(a1$range_area == min(a1$range_area))]
