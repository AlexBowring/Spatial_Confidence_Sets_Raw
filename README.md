# Spatial Confidence Sets for Raw Effect Size Images
Supporting code to perform the simulations and analyses of the manuscript with the same title available at [https://doi.org/10.1101/631473](https://doi.org/10.1101/631473).

## Table of contents
   * [How to cite](#how-to-cite)
   * [Dependencies](#dependencies)
   * [Simulations](#simulations)
      * [Sim_01](#sim_01)
      * [Sim_02](#sim_02)
      * [Sim_03](#sim_03)
      * [Sim_04](#sim_04)
      * [Sim_05](#sim_05)
      * [Sim_06](#sim_06)
      * [Sim_07](#sim_07)
      * [Sim_08](#sim_08)
      * [Sim_09](#sim_09)
      * [Sim_10](#sim_10)
      * [Sim_11](#sim_11)
      * [Sim_12](#sim_12)
      * [Sim_13](#sim_13)
   * [HCP Analyses](#hcp-analyses)
   * [Analyse Your Own Data](#analyse-your-own-data)

## How to cite

To cite this repository, please cite the corresponding manuscript: 

"Spatial Confidence Sets for Raw Effect Size Images" Alexander Bowring, Fabian Telschow\*, Armin Schwartzman\*, Thomas Nichols\*. bioRxiv 631473; doi: [https://doi.org/10.1101/631473](https://doi.org/10.1101/631473) 

## Dependencies
* `MATLAB2016b` (or later) for the simulations.
* `SPM8` for the simulations.
* `SPM12` to run the permutation test for the traditional inference results displayed in the supplementary figures of the manuscript.

## Simulations
The `Sim_??_scripts` in this repository contain all the code used to obtain the 2D & 3D simulation results presented in Section 4.1, 4.2 & 4.3 of the manuscript. 

In particular, the `Sim_[01,..,04]_scripts` directories contain the code used for the 2D simulations in Section 4.2 of the manuscript. The `Sim_[05,..,11]_scripts` directories contain the code used for the 3D simulations in Section 4.3 of the manuscript. The `Sim_[12,13]_scripts` directories contain the code used for the 2D simulations applying the Gaussian wild bootstrap method from SSS, for the methodological comparisons made in Section 4.1 of the manuscript.

In each `Sim_??_scripts` directory, the `Sim_??.m` file contains the code to generate the synthetic data, apply the confidence sets method, and then assess the empirical coverage. The `Sim_??_subjects.m` and `Sim_??_subjects.sh` files contain code to run multiple simulation jobs in parallel on a high performance computing cluster, saving each job as a `.mat` file. **Users will need to change the SPM8 directory path in the `Sim_??_subjects.m` file to the appropriate directory on their own machine.** Finally, the `concatenate_??_subject.m` files contain code to concatenate the `.mat` results files obtained from each simulation job into a single `.mat` file.

### Sim_01
The `Sim_01_scripts` directory contains the code for the 2D Linear Ramp signal simulations with homogeneous variance, corresponding with the red curve results in Figure 7.1 of the manuscript.

### Sim_02
The `Sim_02_scripts` directory contains the code for the 2D Linear Ramp signal simulations with heterogeneous variance, corresponding with the blue curve results in Figure 7.1 of the manuscript.

### Sim_03
The `Sim_03_scripts` directory contains the code for the 2D Circle signal simulations with homogeneous variance, corresponding with the red curve results in Figure 7.2 of the manuscript, as well as the multiplier t-bootstrap results in Figure 6.1 of the manuscript.

### Sim_04
The `Sim_04_scripts` directory contains the code for the 2D Circle signal simulations with heterogeneous variance, corresponding with the red curve results in Figure 7.2 of the manuscript.

### Sim_05
The `Sim_05_scripts` directory contains the code for the 3D Small Sphere signal simulations with homogeneous variance, corresponding with the red curve results in Figure 8.1 of the manuscript.

### Sim_06
The `Sim_06_scripts` directory contains the code for the 3D Small Sphere signal simulations with heterogeneous variance, corresponding with the blue curve results in Figure 8.1 of the manuscript.

### Sim_07
The `Sim_07_scripts` directory contains the code for the 3D Large Sphere signal simulations with homogeneous variance, corresponding with the red curve results in Figure 8.2 of the manuscript, as well as the multiplier t-bootstrap results in Figure 6.2 of the manuscript.

### Sim_08
The `Sim_08_scripts` directory contains the code for the 3D Large Sphere signal simulations with heterogeneous variance, corresponding with the blue curve results in Figure 8.2 of the manuscript.

### Sim_09
The `Sim_09_scripts` directory contains the code for the 3D Multi Sphere signal simulations with homogeneous variance, corresponding with the red curve results in Figure 8.3 of the manuscript.

### Sim_10
The `Sim_10_scripts` directory contains the code for the 3D Multi Sphere signal simulations with heterogeneous variance, corresponding with the blue curve results in Figure 8.3 of the manuscript.

### Sim_11
The `Sim_11_scripts` directory contains the code for the Biobank signal simulations with Biobank variance, corresponding with the red curve results in Figure 8.4 of the manuscript.

### Sim_12
The `Sim_12_scripts` directory contains the code for the 2D Circle signal simulations with the gaussian wild bootstrap, corresponding with the gaussian wild bootstrap results in Figure 6.1 of the manuscript.

### Sim_13
The `Sim_13_scripts` directory contains the code for the 2D Large Sphere signal simulations with the gaussian wild bootstrap, corresponding with the gaussian wild bootstrap results in Figure 6.2 of the manuscript.

## HCP Analyses
The `HCP_analysis_scripts` directory contains the code used to analyse the HCP data. The `HCP_contour_inf.m` file contains the code used to obtain the confidence sets (and yellow point estimate set) for the HCP data displayed in Figure 9 and Figure 10 of the manuscript. The `permutation_test.m` file is the `SPM12` script used to apply the traditional inference procedure corresponding to the results in the supplementary figures for the manuscript. 

## Analyse Your Own Data
The `Confidence_Sets.m` function contains code to obtain Confidence Sets on your own data. The function has 4 inputs: `4D_COPES`, a 4D volume containing the (subject-level) effect estimate images which you would like to obtain confidence sets for. Individual effect estimate images can be concatenated into a 4D image using the `fslmerge` command-line tool that comes with `FSL`:
`fslmerge -t 4D_FILENAME EFFECT_ESTIMATE_IMAGE_1 EFFECT_ESTIMATE_IMAGE_2 ... EFFECT_ESTIMATE_IMAGE_N`.
`GROUP_MASK_IMAGE` is the group mask. `THRESH` is the threshold c in raw change units. `OUT` is the output directory where the CS images are saved.  