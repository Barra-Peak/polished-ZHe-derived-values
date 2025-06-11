#### Derived Intermediate Zircon (U-Th)/He values for polished grains
## B.A. Peak May 02, 2025

# PURPOSE: This R script calculates volume, alpha-ejection surface area, and related values for Ft corrections and eU calculation using 2D measurements of polished zircon grains
# Equations are modified from Ketcham et al. (2011) to account for partial crystals and range of starting grain geometries
# Volume, Ft and RFt uncertainties applied following recommendations in Zeigler et al. (2024) where applicable based on geometry classification

# INPUTS: .xlsx file with required headers, see example. Input file may have additional columns.
# Note that required column names are referenced throughout the script and should not be changed

# OUTPUTS: default outputs are two .xlsx files:
    # 1. Full data sheet with all derived values and original measurements and analytical measurements and uncertainties
    # 2. Sheet formatted for use with the HeCalc python program for calculating (U-Th)/He dates and uncertainty (Martin et al., 2022, 2023)

# Output file names are set by the user at the beginning of the script.

# -------------------------------------------------------------------------------------------------------------

#### Load necessary packages ----
# these packages must be installed prior to running the script
library(openxlsx) # read in and export microsoft excel files
library(dplyr) # easier data manipulation
library(tidyverse) # for manipulating data tables

#### Import data  and settings for output files ----
setwd("~/Desktop") # set working directory. This is where output files will be saved
filename_root <- "2025_05_14_real_zircon_g0_" # root of filename to be saved. e.g. "Sample_1"

he_data_import <- read.xlsx(file.choose(), startRow = 1, colNames = TRUE, check.names = TRUE, sep.names = "_") # interactive popup file picker
# startRow tells R what row data starts on, if there is no table title, set to 2
# he_data_import <- read.xlsx("/Filepath/example_input_file.xlsx",
#                   startRow = 2, colNames = TRUE, check.names = TRUE, sep.names = "_") #load data using full file path

#### Clean data ----
# the next two lines clean the data table to remove the unit header and footnotes
cut_index <- which(is.na(he_data_import$Grain), arr.ind = TRUE)
he_data <- he_data_import[-cut_index, ]

# the next line converts all values that should be numeric to numeric
# suppressWarnings() supresses an NA message (NAs are from analyses with no Sm)
suppressWarnings(he_data <- he_data %>%
  mutate(across(-c("Grain", "Crystal_Fragment."), as.numeric)) %>%
    mutate_at(c("X147Sm", "X._1σ.3"), ~replace(., is.na(.), 0)))

#### Calculate volume, surface area, Rsv and FT values for all geometries ----

# set up data frame columns to store values
he_data <- he_data %>%
  mutate(Analysis_Session = as.factor(Analysis_Session)) %>% # convert to factor
  mutate(Vol.Ellipse = 0, # ellipsoid volume
         Vol.Ortho = 0, # tetragonal volume
         Vol.Cyl = 0, # cylindrical volume
         SA.Ellipse = 0, # ellipsoid surface area
         SA.Ortho = 0, # tetragonal surface area
         SA.Cyl = 0, # cylindrical surface area
         Rsv.Ellipse = 0, # ellipsoid equivalent spherical radius
         Rsv.Ortho = 0, # tetragonal equivalent spherical radius
         Rsv.Cyl = 0, # cylindrical equivalent spherical radius
         Ft_238U.e = 0, # isotope-specific FT values for each geometry
         Ft_235U.e = 0,
         Ft_232Th.e = 0,
         Ft_174Sm.e = 0,
         Ft_238U.o = 0,
         Ft_235U.o = 0,
         Ft_232Th.o = 0,
         Ft_174Sm.o = 0,
         Ft_238U.c = 0,
         Ft_235U.c = 0,
         Ft_232Th.c = 0,
         Ft_174Sm.c = 0)

### Define constants ----
p <- 1.6075 # ellipsoid coefficient used in Knudson's equation for surface area of an ellipsoid from Ketcham et al. 2011
SD.U235 <- 18.05 # stopping distance for U235 from Cooperdock et al. 2019
SD.U238 <- 15.55 # stopping distance for U238 from Cooperdock et al. 2019
SD.Th232 <- 18.43 # stopping distance for Th232 from Cooperdock et al. 2019
SD.Sm147 <- 4.76 # stopping distance for S147 from Cooperdock et al. 2019

# Calculate values
for(i in 1:as.numeric(nrow(he_data))){
  g <- he_data$Grind_Depth[i] # depth of amount removed through grinding and polishing
  Np <- he_data$Np[i] # number of tetragonal terminations
  if(g == 0) { # whole grain, no polishing
    # Ellipsoid
    a.e <- he_data$Width_1[i]/2
    b.e <- he_data$Width_2[i]/2
    c.e <- (he_data$Length_1[i] + he_data$Length_2[i])/4
    
    he_data$Vol.Ellipse[i] <- (4/3)*pi*a.e*b.e*c.e
    
    he_data$SA.Ellipse[i] <- 4*pi*((((a.e^p*b.e^p) + (b.e^p*c.e^p) + (c.e^p*a.e^p))/3)^(1/p))
    
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
    h <- (he_data$Length_1[i] + he_data$Length_2[i])/2 # original height is average of 2 lengths + g
    r <- (he_data$Width_1[i] + he_data$Width_2[i])/4 # radius is average of 2 widths divided by 2
    a = he_data$Width_1[i]/2 # ellipsoidal termination values
    b = he_data$Width_2[i]/2 # ellipsoidal termination values
    k = ((a-b)^2)/((a+b)^2) # ellipsoid intermediate value
    
    he_data$Vol.Cyl[i] <- pi*a*b*h
    
    he_data$SA.Cyl[i] <- (2*pi*a*b) + (h*pi*(a+b)*(1 + (3*k/(10 + sqrt(4-(3*k))))))
   
    he_data$Rsv.Cyl[i] <- 3*he_data$Vol.Cyl[i]/he_data$SA.Cyl[i]
    
    he_data$Ft_238U.c[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Cyl[i]) + # 238U Ft value
      (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.U238/he_data$Rsv.Cyl[i])^2 +
      ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.U238/he_data$Rsv.Cyl[i])^3
    he_data$Ft_235U.c[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Cyl[i]) + # 238U Ft value
      (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.U235/he_data$Rsv.Cyl[i])^2 +
      ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.U235/he_data$Rsv.Cyl[i])^3
    he_data$Ft_232Th.c[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Cyl[i]) + # 238U Ft value
      (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.Th232/he_data$Rsv.Cyl[i])^2 +
      ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.Th232/he_data$Rsv.Cyl[i])^3
    he_data$Ft_147Sm.c[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Cyl[i]) + # 238U Ft value
      (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.Sm147/he_data$Rsv.Cyl[i])^2 +
      ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.Sm147/he_data$Rsv.Cyl[i])^3
    
    # Tetragonal
    a.o <- min(he_data$Width_1[i], he_data$Width_2[i])
    b.o <- max(he_data$Width_1[i], he_data$Width_2[i])
    c.o <- (he_data$Length_1[i] + he_data$Length_2[i])/2
    
    he_data$Vol.Ortho[i] <- a.o*b.o*c.o - Np*(a.o/4)*(b.o^2 + ((a.o^2)/3))
    
    he_data$SA.Ortho[i] <- 2*(a.o*b.o+b.o*c.o+a.o*c.o) - Np*((a.o^2-b.o^2)/2 + (2-sqrt(2))*a.o*b.o)
    
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
 
  
  } else if(g != 0 & he_data$Orientation[i] == 1){ # polished perpendicular to c-axis
    if(Np == 2){
      Np.p <- 1
    }else{Np.p <- Np}
    # Ellipsoid 
    a.e <- he_data$Width_1[i]/2
    b.e <- (he_data$Length_1[i] + he_data$Width_2[i])/4
    c.e <- he_data$Length_2[i]
    
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
    h <- (he_data$Length_1[i] + he_data$Length_2[i])/2 # original height is average of 2 lengths + g
    r <- (he_data$Width_1[i] + he_data$Width_2[i])/4 # radius is average of 2 widths divided by 2
    a = he_data$Width_1[i]/2 # ellipsoidal termination values
    b = he_data$Width_2[i]/2 # ellipsoidal termination values
    k = ((a-b)^2)/((a+b)^2)
    
    he_data$Vol.Cyl[i] <- pi*a*b*h
    
    he_data$SA.Cyl[i] <- (pi*a*b) + (h*pi*(a+b)*(1 + (3*k/(10 + sqrt(4-(3*k)))))) # one end is polished so don't count
    
    he_data$Rsv.Cyl[i] <- 3*he_data$Vol.Cyl[i]/he_data$SA.Cyl[i]
    
    he_data$Ft_238U.c[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Cyl[i]) + # 238U Ft value
      (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.U238/he_data$Rsv.Cyl[i])^2 +
      ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.U238/he_data$Rsv.Cyl[i])^3
    he_data$Ft_235U.c[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Cyl[i]) + # 238U Ft value
      (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.U235/he_data$Rsv.Cyl[i])^2 +
      ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.U235/he_data$Rsv.Cyl[i])^3
    he_data$Ft_232Th.c[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Cyl[i]) + # 238U Ft value
      (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.Th232/he_data$Rsv.Cyl[i])^2 +
      ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.Th232/he_data$Rsv.Cyl[i])^3
    he_data$Ft_147Sm.c[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Cyl[i]) + # 238U Ft value
      (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.Sm147/he_data$Rsv.Cyl[i])^2 +
      ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.Sm147/he_data$Rsv.Cyl[i])^3
    
    # Tetragonal
    a.o <- min(he_data$Width_1[i], he_data$Width_2[i])
    b.o <- max(he_data$Width_1[i], he_data$Width_2[i])
    c.o <- (he_data$Length_1[i] + he_data$Length_2[i])/2
    
    he_data$Vol.Ortho[i] <- a.o*b.o*c.o - Np.p*(a.o/4)*(b.o^2 + ((a.o^2)/3))
    
    he_data$SA.Ortho[i] <- 2*(a.o*b.o+b.o*c.o+a.o*c.o) - Np.p*((a.o^2-b.o^2)/2 + (2-sqrt(2))*a.o*b.o) - (a.o*b.o) # subtract polished face
    
    he_data$Rsv.Ortho[i] <- 3*he_data$Vol.Ortho[i]/he_data$SA.Ortho[i]
    
    he_data$Ft_238U.o[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ortho[i]) + 
      (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))*
         (a.o + b.o)*Np.p)*(SD.U238^2)/he_data$Vol.Ortho[i] # 238U Ft value
    he_data$Ft_235U.o[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ortho[i]) + 
      (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))*
         (a.o + b.o)*Np.p)*(SD.U235^2)/he_data$Vol.Ortho[i] # 235U Ft value
    he_data$Ft_232Th.o[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ortho[i]) + 
      (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))*
         (a.o + b.o)*Np.p)*(SD.Th232^2)/he_data$Vol.Ortho[i] # 232Th Ft value
    he_data$Ft_147Sm.o[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ortho[i]) + 
      (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))*
         (a.o + b.o)*Np.p)*(SD.Sm147^2)/he_data$Vol.Ortho[i] # 147Sm Ft value
    
  } else if(g != 0 & he_data$Orientation[i] == 2) { # polished parallel to c-axis
    # Ellipsoid
    if(g >= (he_data$Width_2[i] + g)/2){ # ground halfway or more
    a.e <- he_data$Width_1[i]/2
    b.e <- he_data$Width_2[i]
    c.e <- (he_data$Length_1[i] + he_data$Length_2[i])/4
    
    he_data$Vol.Ellipse[i] <- (2/3)*pi*a.e*b.e*c.e
    
    he_data$SA.Ellipse[i] <-  2*pi*(((a.e^p*b.e^p + a.e^p*c.e^p + b.e^p*c.e^p)/3)^(1/p))
    
    he_data$Rsv.Ellipse[i] <- 3*he_data$Vol.Ellipse[i]/he_data$SA.Ellipse[i] # ellipsoid equivalent spherical radius
    
    he_data$Ft_238U.e[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.U238/he_data$Rsv.Ellipse[i])^3 # 238U Ft value
    he_data$Ft_235U.e[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.U235/he_data$Rsv.Ellipse[i])^3 # 235U Ft value
    he_data$Ft_232Th.e[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.Th232/he_data$Rsv.Ellipse[i])^3 # 232Th Ft value
    he_data$Ft_147Sm.e[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ellipse[i]) + 
      ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.Sm147/he_data$Rsv.Ellipse[i])^3 # 147Sm Ft value
    
    } else if(g < (he_data$Width_2[i] + g)/2){ # ground less than halfway
      a.e <- he_data$Width_1[i]/2
      b.e <- (he_data$Width_2[i] + g)/2
      c.e <- (he_data$Length_1[i] + he_data$Length_2[i])/4
      a.p <- he_data$Width_P[i]/2
      b.p <- g
      c.p <- he_data$Length_P[i]/2
      
      he_data$Vol.Ellipse[i] <- (4/3)*pi*a.e*b.e*c.e - (2/3)*pi*a.p*b.p*c.p
      
      he_data$SA.Ellipse[i] <-  4*pi*((a.e^p*b.e^p + a.e^p*c.e^p + b.e^p*c.e^p)/3)^(1/p) -
        2*pi*((a.p^p*b.p^p + a.p^p*c.p^p + b.p^p*c.p^p)/3)^(1/p)
      
      he_data$Rsv.Ellipse[i] <- 3*he_data$Vol.Ellipse[i]/he_data$SA.Ellipse[i] # ellipsoid equivalent spherical radius
      
      he_data$Ft_238U.e[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ellipse[i]) + 
        ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.U238/he_data$Rsv.Ellipse[i])^3 # 238U Ft value
      he_data$Ft_235U.e[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ellipse[i]) + 
        ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.U235/he_data$Rsv.Ellipse[i])^3 # 235U Ft value
      he_data$Ft_232Th.e[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ellipse[i]) + 
        ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.Th232/he_data$Rsv.Ellipse[i])^3 # 232Th Ft value
      he_data$Ft_147Sm.e[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ellipse[i]) + 
        ((1/16) + 0.1686*(1 - a.e/he_data$Rsv.Ellipse[i])^2)*(SD.Sm147/he_data$Rsv.Ellipse[i])^3 # 147Sm Ft value
    
      } 
    
    # Cylinder
    r_ideal <- (g + he_data$Width_2[i])/2
    h <- (he_data$Length_1[i] + he_data$Length_2[i])/2
    if(he_data$Width_2[i] <= g){ # ground halfway or more
      a = he_data$Width_1[i]/2
      b = he_data$Width_2[i]
      k = ((a-b)^2)/((a+b)^2)
      
      he_data$Vol.Cyl[i] <- (pi/2)*a*b*h
      
      he_data$SA.Cyl[i] <- (pi*a*b) + (((pi*(a + b))/2)*h*(1 + ((3*k)/(10 + sqrt(4-3*k)))))
      
      he_data$Rsv.Cyl[i] <- 3*he_data$Vol.Cyl[i]/he_data$SA.Cyl[i]
      
      he_data$Ft_238U.c[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Cyl[i]) + # 238U Ft value
        (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.U238/he_data$Rsv.Cyl[i])^2 +
        ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.U238/he_data$Rsv.Cyl[i])^3
      he_data$Ft_235U.c[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Cyl[i]) + # 238U Ft value
        (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.U235/he_data$Rsv.Cyl[i])^2 +
        ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.U235/he_data$Rsv.Cyl[i])^3
      he_data$Ft_232Th.c[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Cyl[i]) + # 238U Ft value
        (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.Th232/he_data$Rsv.Cyl[i])^2 +
        ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.Th232/he_data$Rsv.Cyl[i])^3
      he_data$Ft_147Sm.c[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Cyl[i]) + # 238U Ft value
        (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.Sm147/he_data$Rsv.Cyl[i])^2 +
        ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.Sm147/he_data$Rsv.Cyl[i])^3
      
    } else if(he_data$Width_2[i] > g){ # ground less than halfway
      a <- he_data$Width_1[i]/2
      b <- (he_data$Width_2[i] + g)/2
      k <- ((a-b)^2)/((a+b)^2)
      
      a.p <- he_data$Width_P[i]/2
      b.p <- g
      k.p <- ((a.p-b.p)^2)/((a.p+b.p)^2)
      
      he_data$Vol.Cyl[i] <- pi*h*(a*b - 0.5*a.p*b.p)
      
      he_data$SA.Cyl[i] <-  2*pi*(a*b - 0.5*a.p*b.p) +
        pi*h*(((a + b)*(1 + (3*k/(10 + sqrt(4-(3*k)))))) - 
                 (0.5*(a.p + b.p)*(1 + (3*k.p/(10 + sqrt(4 - (3*k.p)))))))
      
      he_data$Rsv.Cyl[i] <- 3*he_data$Vol.Cyl[i]/he_data$SA.Cyl[i]
      
      he_data$Ft_238U.c[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Cyl[i]) + # 238U Ft value
        (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.U238/he_data$Rsv.Cyl[i])^2 +
        ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.U238/he_data$Rsv.Cyl[i])^3
      he_data$Ft_235U.c[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Cyl[i]) + # 238U Ft value
        (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.U235/he_data$Rsv.Cyl[i])^2 +
        ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.U235/he_data$Rsv.Cyl[i])^3
      he_data$Ft_232Th.c[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Cyl[i]) + # 238U Ft value
        (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.Th232/he_data$Rsv.Cyl[i])^2 +
        ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.Th232/he_data$Rsv.Cyl[i])^3
      he_data$Ft_147Sm.c[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Cyl[i]) + # 238U Ft value
        (((3/2)*0.2112*he_data$Rsv.Cyl[i])/(h + (((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))))*(SD.Sm147/he_data$Rsv.Cyl[i])^2 +
        ((0.0153*(he_data$Rsv.Cyl[i])^3)/((((2/3)*he_data$Rsv.Cyl[i]*h)/(h - ((2/3)*he_data$Rsv.Cyl[i])))^3))*(SD.Sm147/he_data$Rsv.Cyl[i])^3
      } 
    
    # Tetragonal
    c.o <- (he_data$Length_1[i] + he_data$Length_2[i])/2
    b_1 <- he_data$Width_2[i] + g # original width pre-polishing
    if(b_1/2 <= g){ # ground halfway or more
      a.o <- min(he_data$Width_1[i], 2*he_data$Width_2[i])
      b.o <- max(he_data$Width_1[i], 2*he_data$Width_2[i])
      
      he_data$Vol.Ortho[i] <- (a.o*b.o*c.o - Np*(a.o/4)*(b.o^2+((a.o^2)/3)))/2
      
      he_data$SA.Ortho[i] <- (2*(a.o*b.o+b.o*c.o+a.o*c.o) - Np*(((a.o^2-b.o^2)/2) + (2-sqrt(2))*a.o*b.o))/2
      
      he_data$Rsv.Ortho[i] <- 3*he_data$Vol.Ortho[i]/he_data$SA.Ortho[i]
      he_data$Ft_238U.o[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ortho[i]) + 
        ((0.2095*(a.o + b.o +c.o)) - (Np*(a.o + b.o)*(0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))))*
           ((SD.U238^2)/(2*he_data$Vol.Ortho[i])) # 238U Ft value
      he_data$Ft_235U.o[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ortho[i]) + 
        ((0.2095*(a.o + b.o +c.o)) - (Np*(a.o + b.o)*(0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))))*
        ((SD.U235^2)/(2*he_data$Vol.Ortho[i])) # 235U Ft value
      he_data$Ft_232Th.o[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ortho[i]) + 
        ((0.2095*(a.o + b.o +c.o)) - (Np*(a.o + b.o)*(0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))))*
        ((SD.Th232^2)/(2*he_data$Vol.Ortho[i])) # 232Th Ft value
      he_data$Ft_147Sm.o[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ortho[i]) + 
        ((0.2095*(a.o + b.o +c.o)) - (Np*(a.o + b.o)*(0.096 - 0.013*((a.o^2 + b.o^2)/c.o^2))))*
        ((SD.Sm147^2)/(2*he_data$Vol.Ortho[i])) # 147Sm Ft value
      
    } else if (b_1/2 > g){ # ground less than halfway
      a.o <- min(he_data$Width_1[i], he_data$Width_2[i] + g)
      b.o <- max(he_data$Width_1[i], he_data$Width_2[i] + g)
      a.g <- min(he_data$Width_1[i], 2*g)
      b.g <- max(he_data$Width_1[i], 2*g)
      c.g <- c.o - Np*((he_data$Width_2[i] - g)/2)
      if(a.o == he_data$Width_1[i]){
        a_g <- a.o
        b_g <- he_data$Width_2[i]
      } else if(a.o == he_data$Width_2[i] + g){
        a_g <- he_data$Width_2[i]
        b_g <- b.o
      }
      
      he_data$Vol.Ortho[i] <- a.o*b.o*c.o - Np*(a.o/4)*(b.o^2+((a.o^2)/3)) - 
        ((1/2)*(a.g*b.g*c.g - Np*(a.g/4)*(b.g^2+((a.g^2)/3))))
      
      he_data$SA.Ortho[i] <- 2*(a.o*b.o+b.o*c.o+a.o*c.o) - Np*(((a.o^2-b.o^2)/2) + (2-sqrt(2))*a.o*b.o) - 
        (1/2)*(2*(a.g*b.g+b.g*c.g+a.g*c.g) - Np*(((a.g^2-b.g^2)/2) + (2-sqrt(2))*a.g*b.g))
      
      he_data$Rsv.Ortho[i] <- 3*he_data$Vol.Ortho[i]/he_data$SA.Ortho[i]
      
      he_data$Ft_238U.o[i] <- 1 - (3/4)*(SD.U238/he_data$Rsv.Ortho[i]) + 
        (0.2095* (a.o + b.o + c.o) - (0.096 - 0.013*((a_g^2 + b_g^2)/c.o^2))*
           (a_g + b_g)*Np)*(SD.U238^2)/he_data$Vol.Ortho[i]
      he_data$Ft_235U.o[i] <- 1 - (3/4)*(SD.U235/he_data$Rsv.Ortho[i]) + 
        (0.2095* (a.o + b.o +c.o) - (0.096 - 0.013*((a_g^2 + b_g^2)/c.o^2))*
           (a_g + b_g)*Np)*(SD.U235^2)/he_data$Vol.Ortho[i]
      he_data$Ft_232Th.o[i] <- 1 - (3/4)*(SD.Th232/he_data$Rsv.Ortho[i]) + 
        (0.2095* (a.o + b.o + c.o) - (0.096 - 0.013*((a_g^2 + b_g^2)/c.o^2))*
           (a_g + b_g)*Np)*(SD.Th232^2)/he_data$Vol.Ortho[i]
      he_data$Ft_147Sm.o[i] <- 1 - (3/4)*(SD.Sm147/he_data$Rsv.Ortho[i]) + 
        (0.2095* (a.o + b.o + c.o) - (0.096 - 0.013*((a_g^2 + b_g^2)/c.o^2))*
           (a_g + b_g)*Np)*(SD.Sm147^2)/he_data$Vol.Ortho[i]
      
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
  select(Analysis_Session, Grain, Length_1, Width_1, Length_2, Width_2, Geometry, Np,
         Length_P, Width_P, Grind_Depth, Orientation, Crystal_Fragment.,
         X4He, X._1σ, U, X._1σ.1, Th, X._1σ.2, X147Sm, X._1σ.3) %>%
  mutate(Volume = 0, # volume to use
         Volume.unc = 0, # volume uncertainty to use
         Rsv = 0, # Rsv to use
         Rsv.unc = 0, # Rsv uncertainty to use
         Ft.U238 = 0, # isotope-specific FT values and uncertanties
         Ft.U238.unc = 0,
         Ft.U235 = 0,
         Ft.U235.unc = 0,
         Ft.Th232 = 0,
         Ft.Th232.unc = 0,
         Ft.Sm147 = 0,
         Ft.Sm147.unc = 0)

# assign values based on geometry
for(i in 1:(as.numeric(nrow(he_data_rc)))){
  if(he_data_rc$Geometry[i] == 1){ # ellipsoid
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
  } else if(he_data_rc$Geometry[i] == 3){ # tetragonal
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
  } else if(he_data_rc$Geometry[i] == 2){ # cylindrical
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
he_data_rc$He.nmol.g <- (he_data$X4He/he_data_rc$Mass)*10^(-6) # calculate He in nmol/g # calculate 4He concentration
he_data_rc$He.nmol.g.unc <- he_data_rc$He.nmol.g*sqrt((he_data_rc$X._1σ/he_data_rc$X4He)^2 +
                                                        (he_data_rc$Mass.unc/he_data_rc$Mass)^2) # calculate He uncertainty
he_data_rc$U.ppm <- (he_data_rc$U/1000)/he_data_rc$Mass # U concentration in ppm
he_data_rc$U.ppm.unc <- he_data_rc$U.ppm*sqrt((he_data_rc$X._1σ.1/he_data_rc$U)^2 + 
                                                (he_data_rc$Mass.unc/he_data_rc$Mass)^2)
he_data_rc$Th.ppm <- (he_data_rc$Th/1000)/he_data_rc$Mass # Th concentration in ppm
he_data_rc$Th.ppm.unc <- he_data_rc$Th.ppm*sqrt((he_data_rc$X._1σ.2/he_data_rc$Th)^2 + 
                                                  (he_data_rc$Mass.unc/he_data_rc$Mass)^2)
he_data_rc$Sm.ppm <- (he_data_rc$X147Sm/1000)/he_data_rc$Mass # Sm concentration in ppm
he_data_rc$Sm.ppm.unc <- case_when(he_data_rc$X147Sm != 0 ~
  he_data_rc$X147Sm*sqrt((he_data_rc$X._1σ.3/he_data_rc$X147Sm)^2 + 
                                                  (he_data_rc$Mass.unc/he_data_rc$Mass)^2),
                         .default = 0)

#### Calculate eU ----
he_data_rc$eU.ppm <- he_data_rc$U.ppm + 0.238*he_data_rc$Th.ppm + 0.0083*he_data_rc$Sm.ppm
he_data_rc$eU.ppm.unc <- sqrt(he_data_rc$U.ppm.unc^2 + (0.238)*he_data_rc$Th.ppm.unc^2 + (0.0083)*he_data_rc$Sm.ppm.unc^2)

#### Calculate combined Ft and uncertainty ----
Ft.c.calc <- data.frame(Run = he_data_rc$Analysis_Session, 
                        Sample = he_data_rc$Grain, 
                        U = he_data_rc$U,
                        sUf = he_data_rc$X._1σ.1/he_data_rc$U, # fractional uncertainty
                        Th = he_data_rc$Th,
                        sThf = he_data_rc$X._1σ.2/he_data_rc$Th, # fractional uncertainty
                        Th.U = he_data_rc$Th/he_data_rc$U,
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
Rft.calc <- data.frame(Run = he_data_rc$Analysis_Session, 
                        Sample = he_data_rc$Grain, 
                        U = he_data_rc$U,
                        sUf = he_data_rc$X._1σ.1/he_data_rc$U, # fractional uncertainty
                        Th = he_data_rc$Th,
                        sThf = he_data_rc$X._1σ.2/he_data_rc$Th, # fractional uncertainty
                        Th.U = he_data_rc$Th/he_data_rc$U,
                        Ft.c = he_data_rc$Ft.c,
                        Ft.c.unc = he_data_rc$Ft.c.unc.1s)
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

#### Format full data table for export ----
all_data_header <- c("Analysis Session", "Grain", "Length 1", "Width 1", "Length 2", "Width 2", "Geometry", "Np", "Length P", "Width P", "Grind Depth", "Orientation",
"Crystal Fragment?", "4He", "± 1σ", "U", "± 1σ", "Th", "± 1σ", "147Sm", "± 1σ",
"Volume", "± 1σ", "Mass", "± 1σ","4He", "± 1σ", "U", "± 1σ", "Th", "± 1σ", "147Sm", "± 1σ", "eU", "± 1σ", "RSV", "± 1σ", 
"FT238", "± 1σ", "FT235", "± 1σ", "FT232", "± 1σ", "FT147", "± 1σ",  "Combined FT", "± 1σ", "RFT", "± 1σ")

unit_row <- data.frame("Analysis_Session" = "", "Grain" = "", "Length_1" = "(µm)", "Width_1" = "(µm)", "Length_2" = "(µm)",
                       "Width_2" = "(µm)", "Geometry" = "", "Np" = "", "Length_P" = "(µm)", "Width_P" = "(µm)", "Grind_Depth" = "(µm)", 
                       "Orientation" = "", "Crystal_Fragment." = "", "X4He" = "(fmol)",
                       "X._1σ" = "(fmol)", "U" = "(ng)", "X._1σ.1" = "(ng)", "Th" = "(ng)", "X._1σ.2" = "(ng)", "X147Sm" = "(ng)", "X._1σ.3" = "(ng)",
                       "Volume" = "(µm3)", "Volume.unc" = "(µm3)", "Rsv" = "(µm)", "Rsv.unc" = "(µm)", "Ft.U238" = "", "Ft.U238.unc" = "", 
                       "Ft.U235" = "", "Ft.U235.unc" = "", "Ft.Th232" = "", "Ft.Th232.unc" = "", "Ft.Sm147" = "", "Ft.Sm147.unc" = "",
                       "Mass" = "(g)", "Mass.unc" = "(g)", "He.nmol.g" = "(nmol/g)", "He.nmol.g.unc" = "(nmol/g)", "U.ppm" = "(ppm)",
                       "U.ppm.unc" = "(ppm)", "Th.ppm" = "(ppm)", "Th.ppm.unc" = "(ppm)", "Sm.ppm" = "(ppm)", "Sm.ppm.unc" = "(ppm)",
                       "eU.ppm" = "(ppm)", "eU.ppm.unc" = "(ppm", "Ft.c" = "", "Ft.c.unc.1s" = "", "R.ft" = "(µm)", "R.ft.unc" = "(µm)")

all_vals <- rbind(unit_row, he_data_rc) %>%
  relocate(Mass, .after = Volume.unc) %>%
  relocate(Mass.unc, .after = Mass) %>%
  relocate(He.nmol.g, .after = Mass.unc) %>%
  relocate(He.nmol.g.unc, .after = He.nmol.g) %>%
  relocate(U.ppm, .after = He.nmol.g.unc) %>%
  relocate(U.ppm.unc, .after = U.ppm) %>%
  relocate(Th.ppm, .after = U.ppm.unc) %>%
  relocate(Th.ppm.unc, .after = Th.ppm) %>%
  relocate(Sm.ppm, .after = Th.ppm.unc) %>%
  relocate(Sm.ppm.unc, .after = Sm.ppm) %>%
  relocate(eU.ppm, .after = Sm.ppm.unc) %>%
  relocate(eU.ppm.unc, .after = eU.ppm)

colnames(all_vals) <- all_data_header

#### HeCalc input table ----
hecalc_table <- he_data_rc %>%
  mutate(Sample = paste(as.character(Analysis_Session), Grain, sep = "_")) %>%
  select(Sample, X4He, X._1σ, U, X._1σ.1, Th, X._1σ.2, X147Sm, X._1σ.3, 
         Ft.U238, Ft.U238.unc, Ft.U235, Ft.U235.unc, Ft.Th232, Ft.Th232.unc, 
         Ft.Sm147, Ft.Sm147.unc) %>%
  mutate(He.mol = 10^(-15)*X4He, .after = X._1σ) %>%
  mutate(He.mol.unc = 10^(-15)*X._1σ, .after = He.mol) %>%
  mutate(mol.238 = (U*10^(-9)/238), .before = Th) %>%
  mutate(mol.238.unc = (X._1σ.1*10^(-9)/238), .before = Th) %>%
  mutate(mol.232 = (Th*10^(-9)/232), .before = X147Sm) %>%
  mutate(mol.232.unc = (X._1σ.2*10^(-9)/232), .before = X147Sm) %>%
  mutate(mol.147 = (X147Sm*10^(-9)/147), .before = Ft.U238) %>%
  mutate(mol.147.unc = (X._1σ.3*10^(-9)/147), .before = Ft.U238)
hecalc_header <- c("Sample", "fmol 4He", "±", "mol 4He", "±",
                   "ng 238U", "±", "mol 238U", "±","ng 232Th", "±", "mol 232Th", "±",
                   "ng 147Sm", "±", "mol 147Sm", "±",
                   "238Ft" ,"±", "235Ft", "±", "232Ft", "±", "147Ft", "±")
colnames(hecalc_table) <- hecalc_header

### Save data as microsoft excel files ----
write.xlsx(all_vals, paste(filename_root, "all_values.xlsx", sep = ""))
write.xlsx(hecalc_table, paste(filename_root, "helcalc_input_file.xlsx", sep = ""))

