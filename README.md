# polished-ZHe-derived-values

This R script calculates volume, surface area, and related values for Ft corrections and eU calculation for polished partial zircon grains using standard 2D grain measurements.

Volume, surface area, and Ft equations from Ketcham et al. (2011) and Cooperdock et al. (2019), modified to account for partial crystals.
Volume, Rft, and Ft uncertainties applied following recommendations in Zeigler et al. (2024) where applicable.

INPUTS: .xlsx file with required headers, see example_input_file.xlsx in the repository. Input file may have additional columns beyond those required. Note that required column names are referenced throughout the script and should not be changed.

OUTPUTS: default outputs are two .xlsx files:

1. Full data sheet with all derived values and original measurements, analytical measurements and uncertainties
    
2. Sheet formatted for use with the HeCalc python program for calculating (U-Th)/He dates and uncertainty (Martin et al., 2022, 2023)

Output file names can be changed at the beginning of the script.

<a href="https://doi.org/10.5281/zenodo.15642289"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.15642289.svg" alt="DOI"></a>

References:
Cooperdock, E. H. G., Ketcham, R. A., and Stockli, D. F.: Resolving the effects of 2-D versus 3-D grain measurements on apatite (U-Th)ĝ/ĝHe age data and reproducibility, Geochronology, 1, 17–41, https://doi.org/10.5194/gchron-1-17-2019, 2019.![image](https://github.com/user-attachments/assets/d079d31d-b884-484f-9322-4f0e764272d1)

Ketcham, R. A., Gautheron, C., and Tassan-Got, L.: Accounting for long alpha-particle stopping distances in (U–Th–Sm)/He geochronology: Refinement of the baseline case, Geochim. Cosmochim. Acta, 75, 7779–7791, https://doi.org/10.1016/j.gca.2011.10.011, 2011.![image](https://github.com/user-attachments/assets/f0669ea1-4006-4806-af3d-0a3d9f1c61ea)

Martin, P.: HeCalc (1.0.1), Zenodo, https://doi.org/10.5281/zenodo.7453426, 2022. 

Martin, P. E., Metcalf, J. R., and Flowers, R. M.: Calculation of uncertainty in the (U–Th) ∕ He system, Geochronology, 5, 91–107, https://doi.org/10.5194/gchron-5-91-2023, 2023.




