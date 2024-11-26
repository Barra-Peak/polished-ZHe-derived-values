# polished-ZHe-derived-values

This R script calculates volume, surface area, and related values for Ft corrections and eU calculation for polished partial zircon grains using standard 2D grain measurements.

Volume, surface area, and Ft equations from Ketcham et al. (2011), modified to account for partial crystals.
Ft equations from Cooperdock et al. (2019).
Volume, Rft, and Ft uncertainties applied following recommendations in Zeigler et al. (2024) where applicable.

INPUTS: .xlsx file with required headers, see example_input_file.xlsx in the repository. Input file may have additional columns beyond those required. Note that required column names are referenced throughout the script and should not be changed.

OUTPUTS: default outputs are two .xlsx files:

1. Full data sheet with all derived values and original measurements, analytical measurements and uncertainties
    
2. Sheet formatted for use with the HeCalc python program for calculating (U-Th)/He dates and uncertainty (Martin et al., 2022, 2023)

Output file names can be changed at the beginning of the script.

<a href="https://doi.org/10.5281/zenodo.14219282"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14219282.svg" alt="DOI"></a>
