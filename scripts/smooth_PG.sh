# converting 6 FWHM to sigma: https://brainder.org/2011/08/20/gaussian-kernels-convert-fwhm-to-sigma/
wb_command -metric-smoothing ~/dropbox/fsaverage5_lh.pial_avg.surf.gii hcp.gradients_L_10k.func.gii 2.55 hcpd_gradients_L_10k_6smooth.func.gii
wb_command -metric-smoothing ~/dropbox/fsaverage5_rh.pial_avg.surf.gii hcp.gradients_R_10k.func.gii 2.55 hcpd_gradients_R_10k_6smooth.func.gii
