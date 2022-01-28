README

Master_Script.m runs the full analysis for the entire paper, while various functions
(described in Master_Script) are related to the different figures/sections. 

You will need the following toolboxes on your path: 
+microstate
Fieldtrip
MVPA-Lite

Data to run the analysis can be downloaded from https://osf.io/db9u4/

File HCP230.mat was downloaded from 
https://github.com/lukewtait/evaluate_inverse_methods/tree/master/reduced_atlas
and contains the HCP230 volumetric atlas used for defining ROIs in both simulated
and source reconstructed data sets. 

File PNAS_Smith09_rsn10.nii.gz was downloaded from
https://www.fmrib.ox.ac.uk/datasets/brainmap+rsns/PNAS_Smith09_rsn10.nii.gz
and contains the 10 well-matched fMRI-derived RSNs from Smith et al. (2009)
PNAS 106:13040-13045, and was used for generating ground truth maps when simulating
microstates. 

File layout.mat was created using +microstate.functions.layout_creator, and 
additionally is available with tutorial 5 of the toolbox. 
