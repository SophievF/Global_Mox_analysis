# Global_Mox_analysis
Global analysis of SOC abundance and oxalate-extractable metals (Mox, Alox and Feox) across climate regions and soil depth.

Finishing date: May 2024

Author: Sophie von Fromm

This repository contains all the code to reproduce the analysis and all figures in the publication von Fromm et al. (2024), to be submitted to Biogeochemistry.

The folder 'Code' contains all the R code, the folder 'Data' contains all the data needed to run the R Scripts and the folder 'Output' contains all the figures and tables that are produced with the R code.

Only the data file 'Database_all_merged_2024-05-27.csv' is needed to reproduce the analysis and figures in the manuscript. All other files are can be generated with the corresponding R scripts. 


Folder Data:
- Database_all_merged_2024-05-27.csv: Compiled dataset (can also be accessed via zenodo, doi: TO BE ADDED)
- Database_HLZ_2024-05-27.csv: Entire database with extracted climate data (MAP, MAT, HLZ). Can be generated with Mox_SOC_GlobalData.R
- Database_HLZ_grp_2024-05-27.csv: Final database with grouped HLZ. Can be generated with Mox_SOC_HLZ_Grouping.R
  
Folder Code:
- Mox_SOC_GlobalData.R: Extract climate data (MAP, MAT, HLZ). Needed to generate Database_HLZ_2024-05-27.csv
- Mox_SOC_HLZ_Grouping.R Group HLZ. Needed to generate Database_HLZ_grp_2024-05-27.csv
- Mox_SOC_DataDistribution.R: Data distribution analysis. Needed to generate Figures 1, 2, A1-A7, Table A1
- Mox_SOC_LinearMixedEffectsModels_allHLZ.R: Linear mixed-effects models with all HLZ. Needed to generate Figure 5, and Table A8 and A9
- Mox_SOC_LinearMixedEffectsModels_moistHLZ.R: Linear mixed-effects models with HLZ grouped by moisture. Needed to generate Figures 3, A8, and A9; Tables 1, and A2-A4.
- Mox_SOC_LinearMixedEffectsModels_temperatureHLZ.R: Linear mixed-effects models with HLZ grouped by temperature. Needed to generate Figures 4, A12, and A13; Tables 2, and A5-A7.
