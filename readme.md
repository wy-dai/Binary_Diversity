# Binary_Diversity

Additional supporting information. Source data and MATLAB scripts for producing figures in Dai, W.-Y. et al. Diverse configurations of binary asteroids explained by multi-generation satellites (Nature Communications).

The raw simulation snapshot data (approx. 50 GB) is not included due to size constraints. For questions regarding the simulation datasets, please contact the corresponding author listed in the publication.

## Project Structure

### `/binary_asteroid_data/`
Contains observational data for binary asteroids:
- **johnstons archive binary [0304].xlsx** - Binary asteroid catalog data from Johnston's Archive (date: 4 Mar 2025). Only binaries with a D \< 20 km primary are included, and the TNO and Trojan binaries are excluded. Applied in Figure 1, 9 and S2.

### `/source_data/`
Contains raw .fig MATLAB figure files of figures in the paper, in which all the underlying data are accessible:
- **diagram.fig** - Raw .fig file for Figure 2.
- **simulations.fig** - Raw .fig file for Figure 3.
- **diagram_simu.fig** - Raw .fig file for Figure 6.
- **binaryPrimarySpin.fig** - Raw .fig file for Figure S2.
- **simu_isolation.fig** - Raw .fig file for Figure S3.
- **diskMass.fig** - Raw .fig file for Figure S4.

### `/DEM_simulation/`
Contains MATLAB scripts used for the input data generation and post-processing of DEMBody simulations:

#### `/input/`
Input generation scripts:
- **diskGen.m** - Disk generation script for creating initial configurations.
- **satPlusDisk.m** - Script for assembling the disk and the satellite.

#### `/post-processing/`
Scripts for cluster searching and analysis:
- **clusterInfoFinder_pile_sat.m** - Searching cluster information for DEM simulations with pre-existing satellite.
- **clusterAnalysis_sat.m** - Cluster properties analysis for DEM simulations with pre-existing satellite. Applied for producing Figure 3 and S3.
- **clusterInfoFinder_no_sat.m** - Searching clusters for DEM simulations without pre-existing satellite.
- **diskAnalysis.m** - Cluster information finder for DEM simulations with pre-existing satellite. Applied for producing Figure S4.

## Citation
If you use this repository or its data/scripts, please cite:

> Dai, W.-Y. *et al.* (2026) Diverse configurations of binary asteroids explained by multi-generation satellites. *Nature Communications*.

Data and code: **Binary_Diversity** repository.
