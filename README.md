# extended-cessations
Code related to Shinozuka et al. (2025), "Neuroelectrophysiological correlates of extended cessation of consciousness in advanced meditators: A multimodal EEG and MEG study" - https://www.biorxiv.org/content/10.1101/2025.09.19.677455v2

Note that the abbreviation "NS" is used throughout the code to refer to extended cessations of consciousness, since nirodha samapatti is a specific type of extended cessation.

This repository contains the following code:

Source reconstruction:
- common_SR.m runs source reconstruction on participants using the following fucntions:
- GetSourceAnalysis_combined_NS.m, which computes spatial filters for joint EEG+MEG data (n = 3)
- GetSourceAnalysis_eeg_NS.m, which computes spatial filters on just EEG data (n = 1)
- GetSourceAnalysis_meg_NS.m, which computes spatial filters on just MEG data (n = 1)

Power spectrum:
- PSD_NS.m computes the power spectrum, runs within-subject statistical tests both globally and regionally, and plots the power spectrum and spatial distribution of power

Complexity:
- LZc_NS.m computes Lempel-Ziv complexity (both with and without normalization by phase-shuffling) using the LZ76 algorithm, available here: https://github.com/pmediano/EntRate, runs
  within-subject statistical tests both globally and regionally, and plots LZc both globally and regionally
- PE_NS.m computes permutation entropy using this algorithm: https://www.mathworks.com/matlabcentral/fileexchange/44161-permutation-entropy-fast-algorithm, runs within-subject
  statistical tests (only globally), and plots PE globally
- compute_phi_star_simple.m computes integrated information using the method described by Barrett & Seth (2011), doi.org/10.1371/journal.pcbi.1001052
- compare_phi_star_simple.m runs the above function, performs within-subject statistical tests (only globally), and plots integrated information globally

Functional connectivity:
- compute_wPLI_epoch.m computes the weighted phase-lag index using the method described by Vinck et al. (2011), doi.org/10.1016/j.neuroimage.2011.01.055
- compare_wPLI_and_graph_measures.m runs the above function, performs within-subject statistical tests between pairs of Yeo networks, and plots the wPLI matrix. It also computes several
graph-theoretic measures using functions in the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/).
