# PWs

# 1. Pre-processing
[preProc_PW.m](https://github.com/PennLINC/PWs/blob/main/scripts/preProc_PW.m) will motion mask subject time series, and generate indices of TRs in 10+ length continuous windows uninterrupted by FD > 0.2 / outlier frames. These windows (.txt files) are saved to projects/hcpd, not projects/pinesParcels at the moment. This script then calls [downsample_TS.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_TS.sh) to downsample this timeseries from the original spatial resolution (fslr). Finally, this script will compute optical flow on the freesurfer sphere (with [OpFl_Sph_CompVer.m](https://github.com/PennLINC/PWs/blob/main/scripts/OpFl_Sph_CompVer.m) ) , calculate angular distance from the PGG at each face for each vector (with [PGG_AngDistCalc.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc.m) ), and then re-mask out the medial wall for the subject (with [mask_mw_faces](https://github.com/PennLINC/PWs/blob/main/scripts/mask_mw_faces.m) )

items to validate in preProc_PW:

#### [apply_motion_mask.m](https://github.com/PennLINC/PWs/blob/main/scripts/apply_motion_mask.m): txt files generated for continuous segments is accurate

#### [OpFl_Sph_fs4.m](https://github.com/PennLINC/PWs/blob/main/scripts/OpFl_Sph_fs4.m): angular transformations are accurate

#### [PGG_AngDistCalc4.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc4_CompVer.m): Angular distance calculated accurately

###### note 1: separate processing stream exists for carit data: should be exactly equivalent, but have "\_c\" after script names to differentiate
###### note 2: OpFl_Sph and PGG_AngDistCalc are run as compiled progams (more specifically, from [here](https://github.com/PennLINC/PWs/blob/main/scripts/run_OpFl_Sph_CompVer.sh) and [here](https://github.com/PennLINC/PWs/blob/main/scripts/run_PGG_AngDistCalc4_CompVer.sh), respectively.) They are compiled directly from the linked files, but note the scripts used for large SGE submissions are derived from those linked. Expect non-compiled versions to run 3-5 times slower.
 
# 2. Reference map feature extraction
The PGG needs to be derived from a downsampled PG, and 1000 null PGGs need to be derived from spun PGs. The PG is first downsampled to fsaverage4 ([downsample_gPG4.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_gPG4.sh)). Next, the downsampled PG is converted from a gifti to .mat ([PGfuncgii_2_mat.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGfuncgii_2_mat.m)), as subsequent compiled matlab code cannot utilize the gifti() command properly without xml error. The equivalent script is run for the curvature map in [Cfuncgii_2_mat.m](https://github.com/PennLINC/PWs/blob/main/scripts/Cfuncgii_2_mat.m)

items to validate in section 2:

#### No L/R mislabels

# 3. Reference map spinning + feature extraction: spatial null
To create null maps, the [spin test](https://github.com/spin-test/spin-test) is adapted: maps are spun and then gradients are calculated on spun maps to obtain spun vector fields. See [SpinPG.m](https://github.com/PennLINC/PWs/blob/main/scripts/spin_pg.m) for implementaiton of spins, and [AggSpins_TakeGrad.m](https://github.com/PennLINC/PWs/blob/main/scripts/AggSpins_TakeGrad.m) for conversion of spun PGs to spun PGGs.

# 4. Extract Directional info from OpFlow output
[Extract_BUTD_ResultantVecs.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs.m) will calculate the TRs for which signal is hierarchically ascending vs. descending, and saveout subj metrics. [BUTD_to_Rformat.m](https://github.com/PennLINC/PWs/blob/main/scripts/BUTD_to_Rformat.m) will convert these values to .csv files for R to leverage.

# 5. Group-level directionality

[Bin_And_Aggregate_PGGDistribution_180.m](https://github.com/PennLINC/PWs/blob/main/scripts/Bin_And_Aggregate_BuProp_180.m) segments subject-level angular distributions into 18 angular distance bins. This reduction is used to aggregate a group-level angular distance histogram. Further, this script saves out normative, binned directionality for each cortical face.

[Bin_PGGDistributions360.m](https://github.com/PennLINC/PWs/blob/main/scripts/Bin_PGGDistributions360.m) segments subject-level angular distributions into 36 angular distance bins. This reduction is used to aggregate a group-level angular distance histogram in [Aggregate_PGGDistributions360.m](https://github.com/PennLINC/PWs/blob/main/scripts/Aggregate_PGGDistributions360.m). The second script also saves out normative, binned directionality for each cortical face.

# 6. Task effects

facewise_ttest_([L](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_L.R)/[R](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_R.R)).R: this will conduct mass univariate tests on the extracted top-down proportion metrics. Follow it up with [FDR_facewise_t.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise_t.R) to correct for multiple comparisons

- taskEffects.md - from .rmd

# 6. Age effects

facewise_stats_([L](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_stats_L.R)/[R](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_stats_R.R)).R: this will conduct mass univariate tests on the extracted top-down proportion metrics. Follow it up with [FDR_facewise.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise.R) to correct for multiple comparisons
 
- ageEffects.md - from .rmd

###### note 1: separate processing stream for carit is facewise_ttest_(L/R) (and FDR), paired t-test for rest vs. task carit















