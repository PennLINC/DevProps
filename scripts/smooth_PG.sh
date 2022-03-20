# converting 15 mm FWHM to sigma: https://brainder.org/2011/08/20/gaussian-kernels-convert-fwhm-to-sigma/
wb_command -metric-smoothing ~/dropbox/fsaverage5_lh.pial_avg.surf.gii hcp.gradients_L_10k.func.gii 6.37 hcpd_gradients_L_10k_15smooth.func.gii
wb_command -metric-smoothing ~/dropbox/fsaverage5_rh.pial_avg.surf.gii hcp.gradients_R_10k.func.gii 6.37 hcpd_gradients_R_10k_15smooth.func.gii
