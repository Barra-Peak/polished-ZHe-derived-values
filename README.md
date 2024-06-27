# polished-ZHe-derived-values
This R script calculates volume, surface area, and related values for Ft corrections and eU calculation for polished partial zircon grains using standard 2D grain measurements.

Volume, surface area, and Ft equations from Ketcham et al. (2011), modified to account for partial crystals
Ft equations from Cooperdock et al. (2019)
Volume, Rft, and Ft uncertainties applied following recommendations in Zeigler et al. (2024) where applicable

INPUTS: standard CU TraIL data sheet as a .csv file with the addition of 4 columns labeled as follows:
  Run: index for corresponding in-situ U-Pb data collection if applicable, if NA, fill with numeric 1
  Grind.Depth: the depth removed from the grains by polishing, can be estimated with glass beads or other
  object of known size polished alongside the grains
  Geom: numeric grain geometry index:
            # 1 = ellipsoid
            # 2 = cylindrical
            # 3 = tetragonal/orthorhombic
            # 4 = hexagonal (unsupported)
  Grind.Orientation.to.C.axis: orientation of the grain polished surface relative to the crystal c-axis:
            # 1 = polished perpendicular to c-axis
            # 2 = polished parallel to c-axis

Note that importing the data file as something other than a .csv file (e.g., .xlsx) changes the data frame
column names which are referenced throughout the script

OUTPUTS: default outputs are two .xlsx files:
    # 1. Full data sheet with all derived values and original measurements, analytical measurements 
    and uncertainties
    # 2. Sheet formatted for use with the HeCalc python program for calculating (U-Th)/He dates and
    uncertainty (Martin et al., 2022, 2023)

Output file names can be changed at the end of the script.
