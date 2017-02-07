# Uses the latitude of the midpoint of a cell to estimate its area
cellsurfarea <- function(latitude, cellsize.lat = 1, cellsize.long = cellsize.lat) {
        # Transform degrees into radians
        # You need this transformation to estimate lamda and phi
        degtorad <- function(deg){
                return(deg*pi/180)
        }
        # Estimates the area surface of a grid cell with size equals to cellsize.lat, in this case 1
        # The values for lambda1, lambda2,phi1 and phi2 will be estimated using the variable latitude, the function degtorad, the cellsize.lat and the cellsize.long
        surfarea <- function(R,lambda1, lambda2,phi1, phi2){
                return(R^2 * (lambda2-lambda1) * (sin(phi2) - sin(phi1)))
        }
        # The radius of the earth
        R <- 6731.007178
        # This is the output of the cellsurfarea function
        # It uses the function surfarea
        # lambda1 == degtorad(1 * cellsize.long)
        # lamba2 == degtorad(0)
        # phi1 == degtorad(latitude + 0.5 * cellsize.lat)
        # phi2 == degtorad(latitude - 0.5 * cellsize.lat)
        
        # As you can see the values for lamba1 and lamda2 are constant, which means that we are only estimating the area of the cell in the longitudinal band 0
        return(surfarea(R, degtorad(1 * cellsize.long), degtorad(0), degtorad(latitude + 0.5 * cellsize.lat), degtorad(latitude - 0.5 * cellsize.lat))) 
}

# define a very small latitudinal band, think in the concept of integration (summing up infinitesimal parts to get the area under a curve)
cellsize.lat <- 0.01
# Create a vector the cells for which we want to estimate the area, can you see that it only goes from 0 to 90 (northern hemisphere)
midpoints <- seq(0, 90, by = cellsize.lat)
# Add 0.5 * cellsize.lat, this is basically the midpoint of the cell. You are adding the half of the size (0.5 * cellsize.lat) to each value to get the midpoint
midpoints <- midpoints[- 9001] + 0.5 * cellsize.lat
# Now you will apply to each one of the midpoints the surfareacell function (see above), that will give you the area for each 0.01 cell
areas <- sapply(midpoints, cellsurfarea, cellsize.lat = cellsize.lat, cellsize.long = 1)
# To obtain the total area of the 1x1 degree you just have to sum the area of every 100 cells with 0.01 size
cellareas <- tapply(areas, rep(0:89 + 0.5, each = 100), sum)
careas <- unname(cellareas)
str(careas)
