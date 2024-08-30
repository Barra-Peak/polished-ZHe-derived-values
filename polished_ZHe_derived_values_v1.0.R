#### Derived Intermediate Zircon (U-Th)/He values for polished grains
## B. Peak June 27, 2024

# PURPOSE: This R script calculates volume, surface area, and related values for Ft corrections and eU calculation for 
  # polished partial zircon grains using standard 2D grain measurements.
# Volume, surface area, and Ft equations from Ketcham et al. (2011), modified to account for partial crystals
# Ft equations from Cooperdock et al. (2019)
# Volume, Rft, and Ft uncertainties applied following recommendations in Zeigler et al. (2024) where applicable

# INPUTS: standard CU TraIL data sheet as a .csv file with the addition of 4 columns labeled as follows:
  # Run: index for corresponding in-situ U-Pb data collection if applicable, if NA, fill with numeric 1
  # Grind.Depth: the depth removed from the grains by polishing, can be estimated with glass beads or other
    # object of known size polished alongside the grains
  # Geom: numeric grain geometry index:
            # 1 = ellipsoid
            # 2 = cylindrical
            # 3 = tetragonal/orthorhombic
            # 4 = hexagonal (unsupported)
  # Grind.Orientation.to.C.axis: orientation of the grain polished surface relative to the crystal c-axis:
            # 1 = polished perpendicular to c-axis
            # 2 = polished parallel to c-axis
# Note that importing the data file as something other than a .csv file (e.g., .xlsx) changes the data frame
  # column names which are referenced throughout the script

# OUTPUTS: default outputs are two .xlsx files:
    # 1. Full data sheet with all derived values and original measurements, analytical measurements 
      # and uncertainties
    # 2. Sheet formatted for use with the HeCalc python program for calculating (U-Th)/He dates and
      # uncertainty (Martin et al., 2022, 2023)

# Output file names can be changed at the end of the script.

# -------------------------------------------------------------------------------------------------------------

#### Load necessary packages ----
library(openxlsx) # read in and export microsoft excel files
library(dplyr) # easier data manipulation
library(tidyverse) # for manipulating datatables

#### Import data and set as numeric values ----
setwd("~/Desktop") # set working directory. This is where output files will be saved
he_data <- read.csv(
  "/Users/barrapeak/Dropbox/Thermo-Detrital/Data/Jacobsville_all_He_06_10_24.csv") # load He data
# he_data <- read.csv(file.choose()) # popup file picker use in published version
he_data$length.1..µm...c. <- as.numeric(he_data$length.1..µm...c.) 
he_data$width.1..µm...d. <- as.numeric(he_data$width.1..µm...d.) 
he_data$length.2..µm...c. <- as.numeric(he_data$length.2..µm...c.) 
he_data$width.2..µm...d. <- as.numeric(he_data$width.2..µm...d.) 

#### Calculate volume, surface area, Rsv and FT values for all geometries ----

# set up dataframe columns to store values
he_data$Run <- as.factor(he_data$Run)
he_data$Vol.Ellipse <- c(rep(0, nrow(he_data))) # ellipsoid volume
he_data$Vol.Ortho <- c(rep(0, nrow(he_data))) # tetragonal volume
he_data$Vol.Cyl <- c(rep(0, nrow(he_data))) # cylindrical volume
he_data$SA.Ellipse <- c(rep(0, nrow(he_data))) # ellipsoid surface area
he_data$SA.Ortho <- c(rep(0,nrow(he_data))) # tetragonal surface area
he_data$SA.Cyl <- c(rep(0, nrow(he_data))) # cylindrical surface area
he_data$Rsv.Ellipse <- c(rep(0, nrow(he_data))) # ellipsoid equivalent spherical radius
he_data$Rsv.Ortho <- c(rep(0, nrow(he_data))) # tetragonal equivalent spherical radius
he_data$Rsv.Cyl <- c(rep(0, nrow(he_data))) # cylindrical equivalent spherical radius
he_data$Ft_238U.e <- c(rep(0, nrow(he_data))) # 238U Ft value
he_data$Ft_235U.e <- c(rep(0, nrow(he_data))) # 235U Ft value
he_data$Ft_232Th.e <- c(rep(0, nrow(he_data))) # 232Th Ft value
he_data$Ft_147Sm.e <- c(rep(0, nrow(he_data))) # 147Sm Ft value
he_data$Ft_238U.c <- c(rep(0, nrow(he_data))) # 238U Ft value
he_data$Ft_235U.c <- c(rep(0, nrow(he_data))) # 235U Ft value
he_data$Ft_232Th.c <- c(rep(0, nrow(he_data))) # 232Th Ft value
he_data$Ft_147Sm.c <- c(rep(0, nrow(he_data))) # 147Sm Ft value
he_data$Ft_238U.o <- c(rep(0, nrow(he_data))) # 238U Ft value
he_data$Ft_235U.o <- c(rep(0, nrow(he_data))) # 235U Ft value
he_data$Ft_232Th.o <- c(rep(0, nrow(he_data))) # 232Th Ft value
he_data$Ft_147Sm.o <- c(rep(0, nrow(he_data))) # 147Sm Ft value

p <- 1.6075 # ellipsoid Ft coefficient from Ketcham et al. 2011
SD.U235 <- 18.05 # stopping distance for U235 from Cooperdock et al. 2019
SD.U238 <- 15.55 # stopping distance for U238 from Cooperdock et al. 2019
SD.Th232 <- 18.43 # stopping distance for Th232 from Cooperdock et al. 2019
SD.Sm147 <- 4.76 # stopping distance for S147 from Cooperdock et al. 2019

# Calculate values
for(i in 1:as.numeric(nrow(he_data))){
  g <- he_data$Grind.Depth[i] # depth of amount removed through polishing
  Np <- he_data$Np..f.[i] # number of tetragonal terminations
  if(he_data$Grind.Orientation.to.C.axis[i] == 1){ # polished perpendicular to c-axis
    # Ellipsoid 
    a.e <- (he_data$length.1..µm...c.[i] + he_data$width.2..µm...d.[i])/4
    b.e <- he_data$width.1..µm...d.[i]/2
    c.e <- he_data$length.2..µm...c.[i]
    
    he_data$Vol.Ellipse[i] <- (2/3)*pi*a.e*b.e*c.e
    he_data$SA.Ellipse[i] <- 2*pi*((a.e^p*b.e^p + b.e^p*c.e^p + c.e^p*a.e^p)/3)^(1/p)
    he_data$Rsv.Ellipse[i] <- 3*he_data$Vol.Ellipse[i]/he_data$SA.Ellipse[i] # ellipsoid equivalent spherical radius
    he_data$Ft_238U.e[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.U238/he_data$Rsv.Ellipse[i])^3 # 238U Ft value
    he_data$Ft_235U.e[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.U235/he_data$Rsv.Ellipse[i])^3 # 235U Ft value
    he_data$Ft_232Th.e[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.Th232/he_data$Rsv.Ellipse[i])^3 # 232Th Ft value
    he_data$Ft_147Sm.e[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.Sm147/he_data$Rsv.Ellipse[i])^3 # 147Sm Ft value
    
    # Cylinder
    h <- (he_data$length.2..µm...c.[i] + he_data$width.1..µm...d.[i])/2 # original height is average of 2 lengths + g
    r <- (he_data$width.1..µm...d.[i] + he_data$width.2..µm...d.[i])/4 # radius is average of 2 widths divided by 2
    
    he_data$Vol.Cyl[i] <- pi*r^2*h
    he_data$SA.Cyl[i] <- (2*pi*r*h) + (pi*r^2) # one end is polished so don't count
    he_data$Rsv.Cyl[i] <- 3*he_data$Vol.Cyl[i]/he_data$SA.Cyl[i]
    he_data$Ft_238U.c[i] <- 1 - (1/2)*((r+h)*SD.U238/(r*h)) + 0.2122*(SD.U238^2)/(r*h) + 0.0153*(SD.U238^3/r^3) # 238U Ft value
    he_data$Ft_235U.c[i] <- 1 - (1/2)*((r+h)*SD.U235/(r*h)) + 0.2122*(SD.U235^2)/(r*h) + 0.0153*(SD.U235^3/r^3) # 235U Ft value
    he_data$Ft_232Th.c[i] <- 1 - (1/2)*((r+h)*SD.Th232/(r*h)) + 0.2122*(SD.Th232^2)/(r*h) + 0.0153*(SD.Th232^3/r^3) # 232Th Ft value
    he_data$Ft_147Sm.c[i] <- 1 - (1/2)*((r+h)*SD.Sm147/(r*h)) + 0.2122*(SD.Sm147^2)/(r*h) + 0.0153*(SD.Sm147^3/r^3) # 147Sm Ft value
    
    # Tetragonal
    a.o <- min(he_data$width.1..µm...d.[i], he_data$width.2..µm...d.[i])
    b.o <- max(he_data$width.1..µm...d.[i], he_data$width.2..µm...d.[i])
    c.o <- (he_data$length.1..µm...c.[i] + he_data$length.2..µm...c.[i])/2
    
    he_data$Vol.Ortho[i] <- a.o*b.o*c.o - Np*(a.o/4)*(b.o^2 + (a.o^2)/3)
    he_data$SA.Ortho[i] <- 2*(a.o*b.o+b.o*c.o+a.o*c.o) - Np*((a.o^2-b.o^2)/2 + (2-sqrt(2))*a.o*b.o) - a.o*b.o # a*b is the area of the polished side
    he_data$Rsv.Ortho[i] <- 3*he_data$Vol.Ortho[i]/he_data$SA.Ortho[i]
    he_data$Ft_238U.o[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ortho[i]) + 
      (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))*
         (a.o + b.o)*Np)*(SD.U238^2)/he_data$Vol.Ortho[i] # 238U Ft value
    he_data$Ft_235U.o[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ortho[i]) + 
      (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))*
         (a.o + b.o)*Np)*(SD.U235^2)/he_data$Vol.Ortho[i] # 235U Ft value
    he_data$Ft_232Th.o[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ortho[i]) + 
      (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))*
         (a.o + b.o)*Np)*(SD.Th232^2)/he_data$Vol.Ortho[i] # 232Th Ft value
    he_data$Ft_147Sm.o[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ortho[i]) + 
      (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))*
         (a.o + b.o)*Np)*(SD.Sm147^2)/he_data$Vol.Ortho[i] # 147Sm Ft value
    
  } else { # polished parallel to c-axis
    # Ellipsoid, grind more than halfway 
    a.e <- he_data$width.2..µm...d.[i]
    b.e <- he_data$width.1..µm...d.[i]/2
    c.e <- (he_data$length.1..µm...c.[i] + he_data$length.2..µm...c.[i])/4
    he_data$Vol.Ellipse[i] <- (2/3)*pi*a.e*b.e*c.e
    he_data$SA.Ellipse[i] <-  2*pi*((a.e^p*b.e^p + a.e^p*c.e^p + b.e^p*c.e^p)/3)^(1/p)
    he_data$Rsv.Ellipse[i] <- 3*he_data$Vol.Ellipse[i]/he_data$SA.Ellipse[i] # ellipsoid equivalent spherical radius
    he_data$Ft_238U.e[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.U238/he_data$Rsv.Ellipse[i])^3 # 238U Ft value
    he_data$Ft_235U.e[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.U235/he_data$Rsv.Ellipse[i])^3 # 235U Ft value
    he_data$Ft_232Th.e[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.Th232/he_data$Rsv.Ellipse[i])^3 # 232Th Ft value
    he_data$Ft_147Sm.e[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.Sm147/he_data$Rsv.Ellipse[i])^3 # 147Sm Ft value
    
    # Cylinder
    r_ideal <- ((min(he_data$width.1..µm...d.[i],he_data$width.2..µm...d.[i]) + g) + max(he_data$width.1..µm...d.[i],he_data$width.2..µm...d.[i]))/4
    h <- (he_data$length.1..µm...c.[i] + he_data$length.2..µm...c.[i])/2
    if(r_ideal < g){ # ground more than halfway
      r <- min(he_data$width.1..µm...d.[i], he_data$width.2..µm...d.[i])
      b <- sqrt(r_ideal^2-(r_ideal-r)^2)
      he_data$Vol.Cyl[i] <- h*(1/2)*pi*r*b
      he_data$SA.Cyl[i] <- pi*r*b + h*pi/2*(3*(r + b) - sqrt((3*r+b)*(r+3*b)))
      he_data$Rsv.Cyl[i] <- 3*he_data$Vol.Cyl[i]/he_data$SA.Cyl[i]
      he_data$Ft_238U.c[i] <- 1 - (1/2)*((r+h)*SD.U238/(r*h)) + 0.2122*(SD.U238^2)/(r*h) + 0.0153*(SD.U238^3/r^3) # 238U Ft value
      he_data$Ft_235U.c[i] <- 1 - (1/2)*((r+h)*SD.U235/(r*h)) + 0.2122*(SD.U235^2)/(r*h) + 0.0153*(SD.U235^3/r^3) # 235U Ft value
      he_data$Ft_232Th.c[i] <- 1 - (1/2)*((r+h)*SD.Th232/(r*h)) + 0.2122*(SD.Th232^2)/(r*h) + 0.0153*(SD.Th232^3/r^3) # 232Th Ft value
      he_data$Ft_147Sm.c[i] <- 1 - (1/2)*((r+h)*SD.Sm147/(r*h)) + 0.2122*(SD.Sm147^2)/(r*h) + 0.0153*(SD.Sm147^3/r^3) # 147Sm Ft value
    } else if(r_ideal > g){ # ground less than halfway
      b <- sqrt(r_ideal^2 - (r_ideal - g)^2)
      r <- (2*r_ideal - g)/2
      he_data$Vol.Cyl[i] <- h*pi*(r_ideal^2 - (1/2)*g*b)
      he_data$SA.Cyl[i] <-  2*pi*(r_ideal^2 - (1/2)*g*b) + h*(2*pi*r - (pi/2)*(3*(g + b) - sqrt((3*g+b)*(g+3*b))))
      he_data$Rsv.Cyl[i] <- 3*he_data$Vol.Cyl[i]/he_data$SA.Cyl[i]
      he_data$Ft_238U.c[i] <- 1 - (1/2)*((r+h)*SD.U238/(r*h)) + 0.2122*(SD.U238^2)/(r*h) + 0.0153*(SD.U238^3/r^3) # 238U Ft value
      he_data$Ft_235U.c[i] <- 1 - (1/2)*((r+h)*SD.U235/(r*h)) + 0.2122*(SD.U235^2)/(r*h) + 0.0153*(SD.U235^3/r^3) # 235U Ft value
      he_data$Ft_232Th.c[i] <- 1 - (1/2)*((r+h)*SD.Th232/(r*h)) + 0.2122*(SD.Th232^2)/(r*h) + 0.0153*(SD.Th232^3/r^3) # 232Th Ft value
      he_data$Ft_147Sm.c[i] <- 1 - (1/2)*((r+h)*SD.Sm147/(r*h)) + 0.2122*(SD.Sm147^2)/(r*h) + 0.0153*(SD.Sm147^3/r^3) # 147Sm Ft value
    }
    
    # Tetragonal
    a.o <- min(he_data$width.1..µm...d.[i], he_data$width.2..µm...d.[i])
    b.o <- max(he_data$width.1..µm...d.[i], he_data$width.2..µm...d.[i])
    c.o <- (he_data$length.1..µm...c.[i] + he_data$length.2..µm...c.[i])/2
    a_1 <- a.o + g # original width pre-polishing
    if(a_1/2 < g){ # ground more than halfway
      he_data$Vol.Ortho[i] <- (2*a.o*b.o*c.o - Np*(2*a.o/4)*sqrt(b.o^2+((2*a.o)^2)/3))/2
      he_data$SA.Ortho[i] <- (2*(2*a.o*b.o+b.o*c.o+2*a.o*c.o) - Np*(((2*a.o)^2-b.o^2)/2 + (2-sqrt(2))*2*a.o*b.o))/2
      he_data$Rsv.Ortho[i] <- 3*he_data$Vol.Ortho[i]/he_data$SA.Ortho[i]
      he_data$Ft_238U.o[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ortho[i]) + 
        (0.2095* (2*a.o + b.o +c.o) - (0.096 - 0.013*((4*a.o^2 + b.o^2)/c.o^2))*
           (2*a.o + b.o)*Np)*(SD.U238^2)/he_data$Vol.Ortho[i] # 238U Ft value
      he_data$Ft_235U.o[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ortho[i]) + 
        (0.2095* (2*a.o + b.o +c.o) - (0.096 - 0.013*((4*a.o^2 + b.o^2)/c.o^2))*
           (2*a.o + b.o)*Np)*(SD.U235^2)/he_data$Vol.Ortho[i] # 235U Ft value
      he_data$Ft_232Th.o[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ortho[i]) + 
        (0.2095* (2*a.o + b.o +c.o) - (0.096 - 0.013*((4*a.o^2 + b.o^2)/c.o^2))*
           (2*a.o + b.o)*Np)*(SD.Th232^2)/he_data$Vol.Ortho[i] # 232Th Ft value
      he_data$Ft_147Sm.o[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ortho[i]) + 
        (0.2095* (2*a.o + b.o +c.o) - (0.096 - 0.013*((4*a.o^2 + b.o^2)/c.o^2))*
           (2*a.o + b.o)*Np)*(SD.Sm147^2)/he_data$Vol.Ortho[i] # 147Sm Ft value
    } else if (a_1/2 > g){ # ground less than halfway
      gc <- c.o - 2*sqrt((a.o/2)^2 + (b.o/2)^2)
      he_data$Vol.Ortho[i] <- (a_1*b.o*c.o - Np*(a_1/4)*sqrt(b.o^2+((a_1)^2)/3)) - (1/2)*(2*g*b.o*gc - Np*(g/2)*sqrt(b.o^2 + ((2*g)^2)/3))
      he_data$SA.Ortho[i] <- 2*(a_1*b.o+b.o*c.o+a_1*c.o) - Np*((a_1^2-b.o^2)/2 + (2-sqrt(2))*a_1*b.o) - (1/2)*(2*(2*g*b.o + b.o*gc + 2*g*gc) - Np*((4*g^2 + b.o^2)/2 + (2 - sqrt(2))*2*g*b.o))
      he_data$Rsv.Ortho[i] <- 3*he_data$Vol.Ortho[i]/he_data$SA.Ortho[i]
      he_data$Ft_238U.o[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ortho[i]) + 
        (0.2095* (a_1 + b.o +c.o) - (0.096 - 0.013*((a_1^2 + b.o^2)/c.o^2))*
           (a_1 + b.o)*Np)*(SD.U238^2)/he_data$Vol.Ortho[i] - (1 - (3/4)*(SD.U238/he_data$Rsv.Ortho[i]) + 
                                                                 (0.2095* (2*g + b.o + gc) - (0.096 - 0.013*((4*g^2 + b.o^2)/gc^2))*
                                                                    (2*g + b.o)*Np)*(SD.U238^2)/he_data$Vol.Ortho[i]) # 238U Ft value
      he_data$Ft_235U.o[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ortho[i]) + 
        (0.2095* (a_1 + b.o +c.o) - (0.096 - 0.013*((a_1^2 + b.o^2)/c.o^2))*
           (a_1 + b.o)*Np)*(SD.U235^2)/he_data$Vol.Ortho[i] - (1 - (3/4)*(SD.U235/he_data$Rsv.Ortho[i]) + 
                                                                 (0.2095* (2*g + b.o + gc) - (0.096 - 0.013*((4*g^2 + b.o^2)/gc^2))*
                                                                    (2*g + b.o)*Np)*(SD.U235^2)/he_data$Vol.Ortho[i]) # 235U Ft value
      he_data$Ft_232Th.o[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ortho[i]) + 
        (0.2095* (a_1 + b.o +c.o) - (0.096 - 0.013*((a_1^2 + b.o^2)/c.o^2))*
           (a_1 + b.o)*Np)*(SD.Th232^2)/he_data$Vol.Ortho[i] - (1 - (3/4)*(SD.Th232/he_data$Rsv.Ortho[i]) + 
                                                                  (0.2095* (2*g + b.o + gc) - (0.096 - 0.013*((4*g^2 + b.o^2)/gc^2))*
                                                                     (2*g + b.o)*Np)*(SD.Th232^2)/he_data$Vol.Ortho[i]) # 232Th Ft value
      he_data$Ft_147Sm.o[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ortho[i]) + 
        (0.2095* (a_1 + b.o +c.o) - (0.096 - 0.013*((a_1^2 + b.o^2)/c.o^2))*
           (a_1 + b.o)*Np)*(SD.Sm147^2)/he_data$Vol.Ortho[i] - (1 - (3/4)*(SD.Sm147/he_data$Rsv.Ortho[i]) + 
                                                                  (0.2095* (2*g + b.o + gc) - (0.096 - 0.013*((4*g^2 + b.o^2)/gc^2))*
                                                                     (2*g + b.o)*Np)*(SD.Sm147^2)/he_data$Vol.Ortho[i]) # 147Sm Ft value
      
    } 
  }
}

#### Calculate uncertainties ----
vu.e <- 0.21 # volume uncertainty percentage for ellipsoidal grains from Zeigler et al. 2024
vu.o <- 0.13 # volume uncertainty percentage for tetragonal grains from Zeigler et al. 2024
vu.c <- 0.21 # applied volume uncertainty percentage for cylindrical grains - the larger value from Zeigler et al. 2024. Cylindrical volume uncertainty currently unquantified.

he_data$Vol.Ellipse.Unc <- vu.e*he_data$Vol.Ellipse # ellipsoid volume uncertainty
he_data$Vol.Ortho.Unc <- vu.o*he_data$Vol.Ortho # tetragonal volume uncertainty
he_data$Vol.Cyl.Unc <- vu.c*he_data$Vol.Cyl # cylindrical volume uncertainty

# clean table of recalculated derived values
he_data_rc <- he_data %>%
  select(Run, Sample.Name.and.Aliquot, length.1..µm...c., width.1..µm...d., length.2..µm...c., width.2..µm...d., Grind.Depth, Geom, Np..f.,
         X4He..fmol...g., X...h., U..ng...i., X...h..1, Th..ng...j., X...h..2, X147Sm..ng...k., X...h..3, Alpha.Ejection.Applied.)

# set up value columns
he_data_rc$Volume <- c(rep(0, nrow(he_data))) # volume to use
he_data_rc$Volume.unc <- c(rep(0, nrow(he_data))) # volume uncertainty to use
he_data_rc$Rsv <- c(rep(0, nrow(he_data))) # Rsv to use
he_data_rc$Rsv.unc <- c(rep(0, nrow(he_data))) # Rsv uncertainty to use
he_data_rc$Ft.U238 <- c(rep(0, nrow(he_data))) # U238 Ft value
he_data_rc$Ft.U238.unc <- c(rep(0, nrow(he_data))) # U238 Ft value uncertainty
he_data_rc$Ft.U235 <- c(rep(0, nrow(he_data))) # U235 Ft value
he_data_rc$Ft.U235.unc <- c(rep(0, nrow(he_data))) # U235 Ft value uncertainty
he_data_rc$Ft.Th232 <- c(rep(0, nrow(he_data))) # Th232 Ft value
he_data_rc$Ft.Th232.unc <- c(rep(0, nrow(he_data))) # Th232 Ft value uncertainty
he_data_rc$Ft.Sm147 <- c(rep(0, nrow(he_data))) # Sm147 Ft value
he_data_rc$Ft.Sm147.unc <- c(rep(0, nrow(he_data))) # Sm147 Ft value uncertainty

# assign values based on geometry
for(i in 1:(as.numeric(nrow(he_data_rc)))){
  if(he_data_rc$Geom[i] == 1){ # ellipsoid
    he_data_rc$Volume[i] <- he_data$Vol.Ellipse[i]
    he_data_rc$Volume.unc[i] <- he_data$Vol.Ellipse.Unc[i]
    he_data_rc$Rsv[i] <- he_data$Rsv.Ellipse[i]
    he_data_rc$Rsv.unc[i] <- 3*he_data$Vol.Ellipse.Unc[i]/he_data$SA.Ellipse[i]
    he_data_rc$Ft.U238[i] <- he_data$Ft_238U.e[i]
    he_data_rc$Ft.U238.unc[i] <- he_data_rc$Ft.U238[i]*0.03
    he_data_rc$Ft.U235[i] <- he_data$Ft_235U.e[i]
    he_data_rc$Ft.U235.unc[i] <- he_data_rc$Ft.U235[i]*0.04
    he_data_rc$Ft.Th232[i] <- he_data$Ft_232Th.e[i]
    he_data_rc$Ft.Th232.unc[i] <- he_data_rc$Ft.Th232[i]*0.04
    he_data_rc$Ft.Sm147[i] <- he_data$Ft_147Sm.e[i]
    he_data_rc$Ft.Sm147.unc[i] <- he_data_rc$Ft.Sm147[i]*0.01
  } else if(he_data_rc$Geom[i] == 3){ # tetragonal
    he_data_rc$Volume[i] <- he_data$Vol.Ortho[i]
    he_data_rc$Volume.unc[i] <- he_data$Vol.Ortho.Unc[i]
    he_data_rc$Rsv[i] <- he_data$Rsv.Ortho[i]
    he_data_rc$Rsv.unc[i] <- 3*he_data$Vol.Ortho.Unc[i]/he_data$SA.Ortho[i]
    he_data_rc$Ft.U238[i] <- he_data$Ft_238U.o[i]
    he_data_rc$Ft.U238.unc[i] <- he_data_rc$Ft.U238[i]*0.03
    he_data_rc$Ft.U235[i] <- he_data$Ft_235U.o[i]
    he_data_rc$Ft.U235.unc[i] <- he_data_rc$Ft.U235[i]*0.04
    he_data_rc$Ft.Th232[i] <- he_data$Ft_232Th.o[i]
    he_data_rc$Ft.Th232.unc[i] <- he_data_rc$Ft.Th232[i]*0.05
    he_data_rc$Ft.Sm147[i] <- he_data$Ft_147Sm.o[i]
    he_data_rc$Ft.Sm147.unc[i] <- he_data_rc$Ft.Sm147[i]*0.01
  } else if(he_data_rc$Geom[i] == 2){ # cylindrical
    he_data_rc$Volume[i] <- he_data$Vol.Cyl[i]
    he_data_rc$Volume.unc[i] <- he_data$Vol.Cyl.Unc[i]
    he_data_rc$Rsv[i] <- he_data$Rsv.Cyl[i]
    he_data_rc$Rsv.unc[i] <- 3*he_data_rc$Volume.unc[i]/he_data$SA.Cyl[i]
    he_data_rc$Ft.U238[i] <- he_data$Ft_238U.c[i]
    he_data_rc$Ft.U238.unc[i] <- he_data_rc$Ft.U238[i]*0.03
    he_data_rc$Ft.U235[i] <- he_data$Ft_235U.c[i]
    he_data_rc$Ft.U235.unc[i] <- he_data_rc$Ft.U235[i]*0.04
    he_data_rc$Ft.Th232[i] <- he_data$Ft_232Th.c[i]
    he_data_rc$Ft.Th232.unc[i] <- he_data_rc$Ft.Th232[i]*0.05
    he_data_rc$Ft.Sm147[i] <- he_data$Ft_147Sm.c[i]
    he_data_rc$Ft.Sm147.unc[i] <- he_data_rc$Ft.Sm147[i]*0.01
  }
}

#### Calculate masses and parent isotope concentrations in ppm ----
d <- 4.65*10^(-12) # density of zircon in gram/micrometer^3

he_data_rc$Mass <- d*he_data_rc$Volume # calculate mass in g
he_data_rc$Mass.unc <- d*he_data_rc$Volume.unc # calculate mass uncertainty
he_data_rc$He.nmol.g <- (he_data$X4He..fmol...g./he_data_rc$Mass)*10^(-6) # calculate He in nmol/g # calculate 4He concentration
he_data_rc$He.nmol.g.unc <- he_data_rc$He.nmol.g*sqrt((he_data_rc$X...h./he_data_rc$X4He..fmol...g.)^2 +
                                                        (he_data_rc$Mass.unc/he_data_rc$Mass)^2) # calculate He uncertainty
he_data_rc$U.ppm <- (he_data_rc$U..ng...i./1000)/he_data_rc$Mass # U concentration in ppm
he_data_rc$U.ppm.unc <- he_data_rc$U.ppm*sqrt((he_data_rc$X...h..1/he_data_rc$U..ng...i.)^2 + 
                                                (he_data_rc$Mass.unc/he_data_rc$Mass)^2)
he_data_rc$Th.ppm <- (he_data_rc$Th..ng...j./1000)/he_data_rc$Mass # Th concentration in ppm
he_data_rc$Th.ppm.unc <- he_data_rc$Th.ppm*sqrt((he_data_rc$X...h..2/he_data_rc$Th..ng...j.)^2 + 
                                                  (he_data_rc$Mass.unc/he_data_rc$Mass)^2)
he_data_rc$Sm.ppm <- case_when(he_data_rc$X147Sm..ng...k. == "n.m." ~ 0, # Sm concentration in ppm (n.m. converted to zeros for ease of data manipulation)
                               he_data_rc$X147Sm..ng...k. == "0" ~ 0,
                               .default = (as.numeric(he_data_rc$X147Sm..ng...k.)/1000)/he_data_rc$Mass)
he_data_rc$Sm.ppm.unc <- case_when(he_data_rc$X...h..3 == "n.m." ~ 0,
                                   he_data_rc$X...h..3 == "0" ~ 0,
                                   .default = as.numeric(he_data_rc$Sm.ppm)*
                                     sqrt((as.numeric(he_data_rc$X...h..3)/as.numeric(he_data_rc$X147Sm..ng...k.))^2 + 
                                            (he_data_rc$Mass.unc/he_data_rc$Mass)^2))

#### Calculate eU ----
he_data_rc$eU.ppm <- he_data_rc$U.ppm + 0.238*he_data_rc$Th.ppm + 0.0083*he_data_rc$Sm.ppm
he_data_rc$eU.ppm.unc <- sqrt(he_data_rc$U.ppm.unc^2 + (0.238)*he_data_rc$Th.ppm.unc^2 + (0.0083)*he_data_rc$Sm.ppm.unc^2)

#### Calculate combined Ft and uncertainty ----
Ft.c.calc <- data.frame(Run = he_data_rc$Run, 
                        Sample = he_data_rc$Sample.Name.and.Aliquot, 
                        U = he_data_rc$U..ng...i.,
                        sUf = he_data_rc$X...h..1/he_data_rc$U..ng...i., # fractional uncertainty
                        Th = he_data_rc$Th..ng...j.,
                        sThf = he_data_rc$X...h..2/he_data_rc$Th..ng...j., # fractional uncertainty
                        Th.U = he_data_rc$Th..ng...j./he_data_rc$U..ng...i.,
                        Ft.238 = he_data_rc$Ft.U238,
                        sFt238f = he_data_rc$Ft.U238.unc/he_data_rc$Ft.U238,
                        Ft.235 = he_data_rc$Ft.U235,
                        sFt235f = he_data_rc$Ft.U235.unc/he_data_rc$Ft.U235,
                        Ft.232 = he_data_rc$Ft.Th232,
                        sFt232f = he_data_rc$Ft.Th232.unc/he_data_rc$Ft.Th232)
Ft.c.calc <- Ft.c.calc %>%
  mutate(sTh.U.f = sqrt(sThf^2 + sUf^2)) %>%
  mutate(A238 = (1.04 + 0.247*Th.U)^(-1)) %>%
  mutate(sA238f = 0.247*sTh.U.f) %>%
  mutate(A232 = (1 + 4.21/Th.U)^(-1)) %>%
  mutate(sA232f = (1/4.21)*sTh.U.f) %>%
  mutate(A235 = 1 - A238 - A232) %>%
  mutate(sA235f = sqrt((sA238f*A238)^2 + (sA232f*A232)^2)/A235) %>%
  mutate(t1 = A238*Ft.238) %>% # first term in combined sum
  mutate(st1a = t1*sqrt(sA238f^2 + sFt238f^2)) %>% # first term absolute uncertainty
  mutate(t2 = A232 *Ft.232) %>% # second term in combined sum
  mutate(st2a = t2*sqrt(sA232f^2 + sFt232f^2)) %>% # second term absolute uncertainty
  mutate(t3 = A235 *Ft.235) %>% # third term in combined sum
  mutate(st3a = t3*sqrt(sA235f^2 + sFt235f^2)) %>% # third term absolute uncertainty
  mutate(Ft.c = t1 + t2 + t3) %>%
  mutate(Ft.c.unc.1s = sqrt(st1a^2 + st2a^2 + st3a^2))

# add columns to combined data frame
he_data_rc$Ft.c <- Ft.c.calc$Ft.c
he_data_rc$Ft.c.unc.1s <- Ft.c.calc$Ft.c.unc.1s

#### Calculate RFT and uncertainty ----
Rft.calc <- data.frame(Run = he_data_rc$Run, 
                        Sample = he_data_rc$Sample.Name.and.Aliquot, 
                        U = he_data_rc$U..ng...i.,
                        sUf = he_data_rc$X...h..1/he_data_rc$U..ng...i., # fractional uncertainty
                        Th = he_data_rc$Th..ng...j.,
                        sThf = he_data_rc$X...h..2/he_data_rc$Th..ng...j., # fractional uncertainty
                        Th.U = he_data_rc$Th..ng...j./he_data_rc$U..ng...i.,
                        Ft.c = he_data_rc$Ft.c,
                        Ft.c.unc = he_data_rc$Ft.c.unc.2s/2)
Rft.calc <- Rft.calc %>%
  mutate(sTh.U.f = sqrt(sThf^2 + sUf^2)) %>%
  mutate(A238 = (1.04 + 0.247*Th.U)^(-1)) %>%
  mutate(sA238f = 0.247*sTh.U.f) %>%
  mutate(A232 = (1 + 4.21/Th.U)^(-1)) %>%
  mutate(sA232f = (1/4.21)*sTh.U.f) %>%
  mutate(A235 = 1 - A238 - A232) %>%
  mutate(sA235f = sqrt((sA238f*A238)^2 + (sA232f*A232)^2)/A235) %>%
  mutate(t1 = A238*SD.U238) %>% # first term in combined sum
  mutate(st1a = sA238f*SD.U238) %>% # first term absolute uncertainty
  mutate(t2 = A232 *SD.Th232) %>% # second term in combined sum
  mutate(st2a = sA232f*SD.Th232) %>% # second term absolute uncertainty
  mutate(t3 = A235 *SD.U235) %>% # third term in combined sum
  mutate(st3a = sA235f*SD.U235) %>% # third term absolute uncertainty
  mutate(S.c = t1 + t2 + t3) %>%
  mutate(S.R = 1.681 - 2.428*Ft.c + 1.153*Ft.c^2 - 0.406*Ft.c^3) %>%
  mutate(R.ft = S.c/S.R) %>%
  mutate(R.ft.unc = 0.08*R.ft)

# add columns to combined data frame
he_data_rc$R.ft <- Rft.calc$R.ft
he_data_rc$R.ft.unc <- Rft.calc$R.ft.unc


#### hecalc input table ----
hecalc_table <- he_data_rc %>%
  mutate(Sample = paste(as.character(Run), Sample.Name.and.Aliquot, sep = "_")) %>%
  mutate(X147Sm..ng...k. = case_when(X147Sm..ng...k. == "n.m." ~ 0, # Sm concentration in ppm (n.m. converted to zeros for ease of data manipulation)
                                     he_data_rc$X147Sm..ng...k. == "0" ~ 0,
                                     .default = as.numeric(X147Sm..ng...k.))) %>%
  mutate(X...h..3 = case_when(X...h..3 == "n.m." ~ 0, # Sm concentration in ppm (n.m. converted to zeros for ease of data manipulation)
                              X...h..3 == "0" ~ 0,
                              .default = as.numeric(X...h..3))) %>%
  select(Sample, X4He..fmol...g., X...h., U..ng...i., X...h..1, Th..ng...j., X...h..2, X147Sm..ng...k.,
         X...h..3, Ft.U238, Ft.U238.unc, Ft.U235, Ft.U235.unc, Ft.Th232, Ft.Th232.unc, Ft.Sm147, Ft.Sm147.unc) %>%
  mutate(He.mol = 10^(-15)*X4He..fmol...g., .after = X...h.) %>%
  mutate(He.mol.unc = 10^(-15)*X...h., .after = He.mol) %>%
  mutate(mol.238 = (U..ng...i.*10^(-9)/238), .before = Th..ng...j.) %>%
  mutate(mol.238.unc = (X...h..1*10^(-9)/238), .before = Th..ng...j.) %>%
  mutate(mol.232 = (Th..ng...j.*10^(-9)/232), .before = X147Sm..ng...k.) %>%
  mutate(mol.232.unc = (X...h..2*10^(-9)/232), .before = X147Sm..ng...k.) %>%
  mutate(mol.147 = as.numeric(X147Sm..ng...k.)*10^(-9)/147, .before = Ft.U238) %>%
  mutate(mol.147.unx = as.numeric(X...h..3)*10^(-9)/147, .before = Ft.U238)
hecalc_header <- c("Sample", "fmol 4He", "±", "mol 4He", "±",
                   "ng 238U", "±", "mol 238U", "±","ng 232Th", "±", "mol 232Th", "±",
                   "ng 147Sm", "±", "mol 147Sm", "±",
                   "238Ft" ,"±", "235Ft", "±", "232Ft", "±", "147Ft", "±")
colnames(hecalc_table) <- hecalc_header

### Save data as microsoft excel files ----
write.xlsx(he_data_rc, "all_recalculated_values.xlsx")
write.xlsx(hecalc_table, "helcalc_input_file.xlsx")
