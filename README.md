

# Auditory Model Bifurcation Analysis using DDE-Biftool

This repository contains MATLAB code for performing bifurcation computations in an auditory model using **DDE-Biftool**.

## Requirements

Ensure that the following folders are added to your MATLAB path:

```matlab
   base=[pwd(),filesep(),'ddebiftool_snapshot_23October2022',filesep()];
    addpath([base,'ddebiftool'],...
            [base,'ddebiftool_extra_psol'],...
            [base,'ddebiftool_utilities'],...
            [base,'ddebiftool_extra_rotsym'],...
            [base,'ddebiftool_extra_nmfm'],...
            [base,'ddebiftool_extra_symbolic'],...
            [base,'ddebiftool_coco'],...
            [base,'Supporting_functions']);
```

### Note: 
You **do not** need to run any scripts from the `ddebiftool_snapshot_23October2022`. They simply need to be in the MATLAB path in all the scripts in this repository.

## Script Execution

The scripts can be divided into **sets**, and as each script saves its results with the `save(...)` command, which are loaded by subsequent scripts using the `load(...)` command. they must be run in a sequence. 
 
**Run the scripts in following order:**

### Set 1: Bifurcation Analysis Based on Experimental Data
1. `branch_of_sympos_original_case_and_threshold_crossing.m`
2. `branching_towrds_asymmetric_sols_original_case.m`
3. `sym_breaking_original_case.m`
4. `tracking_threshold_crossing_asymmetric.m`
5. `identifying_Sols_in_different_regions.m`
6. `Plot_auditory_model_original_case_submissionversion.m`

### Set 2: One-Parameter computation for different values of the tone duration (TD) and delay (D) parameters 
1. `br_of_POs_different_values_TD_D_in_loop.m`
2. `asymmetric_branches_case12.m`
3. `plot_one_parameter_for_all_caes.m`

### Set 3: Computation in (rp, df)-parameter space (Large TD)
1. `touching_theta_case_TD.m`
2. `symmetry_breaking_for_case_TD.m`

### Set 4: Computation in (rp, df)-parameter space (D = TD = 0.05)
1. `touching_theta_case_TD_and_D.m`
2. `symmetry_breaking_for_case_TD_and_D.m`

### Set 5: Computation in (rp, df)-parameter space (Large D = 0.05)
1. `touching_theta_case_D.m`
2. `symmetry_breaking_for_case_D.m`

