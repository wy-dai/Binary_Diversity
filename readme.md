# Binary_Diversity

Additional supporting information. MATLAB scripts in Dai, W.-Y. et al. Diverse configurations of binary asteroids explained by multi-generation satellites (Nature Communications).

## Project Structure

### `/binary_asteroid_data/`
Contains observational data for binary asteroids:
- **johnstons archive binary [0304].xlsx** - Binary asteroid catalog data from Johnston's Archive

### `/DEM_simulation/`
Contains MATLAB scripts used for the input data generation and post-processing of DEMBody simulations:

#### `/input/`
Input generation scripts:
- **diskGen.m** - Disk generation script for creating initial configurations
- **satPlusDisk.m** - Script for assembling the disk and the satellite

#### `/post-processing/`
Scripts for cluster searching and analysis:
- **clusterAnalysis_no_sat.m** - Cluster properties analysis for simulations without pre-existing satellite
- **clusterAnalysis_sat.m** - Cluster properties analysis for simulations with pre-existing satellite
- **clusterInfoFinder_no_sat.m** - Searching clusters for simulations without pre-existing satellite
- **clusterInfoFinder_pile_sat.m** - Cluster information finder for imulations with pre-existing satellite

## Citation
If you use this repository or its data, please cite:

> Dai, W.-Y. *et al.* (2026) Diverse configurations of binary asteroids explained by multi-generation satellites. *Nature Communications*.

Data and code: **Binary_Diversity** repository.

For questions regarding DEMBody datasets, please contact the corresponding author listed in the publication.