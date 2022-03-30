# PWs

# 1. Pre-processing
[preProc_PW.m](https://github.com/PennLINC/PWs/blob/main/scripts/preProc_PW.m) will motion mask subject time series, and generate indices of TRs in 10+ length continuous windows uninterrupted by FD > 0.2 / outlier frames. These windows (.txt files) are saved to projects/hcpd, not projects/pinesParcels at the moment. This script then calls [downsample_TS.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_TS.sh) to downsample this timeseries from the original spatial resolution (fslr). Finally, this script will compute optical flow on the freesurfer sphere (with [OpFl_Sph_CompVer.m](https://github.com/PennLINC/PWs/blob/main/scripts/OpFl_Sph_CompVer.m) ) , calculate angular distance from the PGG at each face for each vector (with [PGG_AngDistCalc.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc.m) ), and then re-mask out the medial wall for the subject (with [mask_mw_faces](https://github.com/PennLINC/PWs/blob/main/scripts/mask_mw_faces.m) )

items to validate in preProc_PW:

#### [apply_motion_mask.m](https://github.com/PennLINC/PWs/blob/main/scripts/apply_motion_mask.m): txt files generated for continuous segments is accurate

#### [OpFl_Sph_fs4.m](https://github.com/PennLINC/PWs/blob/main/scripts/OpFl_Sph_fs4.m): angular transformations are accurate

#### [PGG_AngDistCalc4.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc4_CompVer.m): Angular distance calculated accurately

###### note 1: separate processing stream exists for carit data: should be exactly equivalent, but have "\_c\" after script names to differentiate
###### note 2: OpFl_Sph and PGG_AngDistCalc are run as compiled progams (more specifically, from [here](https://github.com/PennLINC/PWs/blob/main/scripts/run_OpFl_Sph_CompVer.sh) and [here](https://github.com/PennLINC/PWs/blob/main/scripts/run_PGG_AngDistCalc4_CompVer.sh), respectively.) They are compiled directly from the linked files, but note the scripts used for large SGE submissions are derived from those linked. Expect non-compiled versions to run 3-5 times slower.
 
# 2.1 Reference map feature extraction
The PGG needs to be derived from a downsampled PG, and 1000 null PGGs need to be derived from spun PGs. The PG is first downsampled to fsaverage4 ([downsample_gPG4.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_gPG4.sh)). Next, the downsampled PG is converted from a gifti to .mat ([PGfuncgii_2_mat.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGfuncgii_2_mat.m)), as subsequent compiled matlab code cannot utilize the gifti() command properly without xml error. The equivalent script is run for the curvature map in [Cfuncgii_2_mat.m](https://github.com/PennLINC/PWs/blob/main/scripts/Cfuncgii_2_mat.m)

items to validate in 2.1:

#### No L/R mislabels

# 2.2 Reference map spinning + feature extraction: spatial null
To create null maps, the [spin test](https://github.com/spin-test/spin-test) is adapted: maps are spun and then gradients are calculated on spun maps to obtain spun vector fields. See [SpinPG.m](https://github.com/PennLINC/PWs/blob/main/scripts/SpinPG.m) for implementaiton of spins, and [AggSpins_TakeGrad.m](https://github.com/PennLINC/PWs/blob/main/scripts/AggSpins_TakeGrad.m) for conversion of spun PGs to spun PGGs.

items to validate in 2.2:

#### No medial wall mis-indexing

# 3.1 Extract Directional info from OpFlow output
[Extract_BUTD_ResultantVecs.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs.m) will calculate the TRs for which signal is hierarchically ascending vs. descending, and saveout subj metrics. Task data can be processed the same way with an equivalent script, [Extract_BUTD_ResultantVecs_c.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs_c.m). *Curvature* data can also be processed the same way, also with an equivalent script [Extract_BUTD_ResultantVecs_curv.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs_curv.m). 

For replicating the subject data depicted in the schematic figure, we can use [viz_vec_fields.m](https://github.com/PennLINC/PWs/blob/main/scripts/viz_vec_fields4.m).

items to validate in 3.1:

#### No L/R mislabels

# 3.2 Extract Directional info for dip test, compare to null models

Here we need to compute the dip statistic (Hartigan's Dip Statistic, see [This scipt](https://github.com/diazlab/chance/blob/master/hartigansdiptest.m), which is a matlab implementation of [this paper](https://www.jstor.org/stable/2347485?seq=1)). Because the angular distributions over each face over each TR over each subject and over spun reference angles is too much to store in one data frame, we will run the dip test in matlab iteratively as it runs over each spin. The "true" dip statistic is recorded last. 

[PGG_AngDistCalc_snull_CompVer.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc_snull_CompVer.m) runs on individual participants. You can call an altered form of [preProc_PW.m](https://github.com/PennLINC/PWs/blob/main/scripts/preProc_PW2.m) through an altered [JobSubmitta.py](https://github.com/PennLINC/PWs/blob/main/scripts/JobSubmitta.py) to qsub individual null-vs.-true dip 
Equivalent script for curvature: [CG_AngDistCalc_snull_CompVer.m](https://github.com/PennLINC/PWs/blob/main/scripts/CG_AngDistCalc_snull_CompVer.m)

Finally, this can be visualized with a few [gglines of code](https://github.com/PennLINC/PWs/blob/main/scripts/NullPlotting.md) in R

# 4. Group-level directionality

[aggregate_OpFlDistrModes.m](https://github.com/PennLINC/PWs/blob/main/scripts/aggregate_OpFlDistrModes.m) segments subject-level angular distributions into 18 angular distance bins. This reduction is needed to reduce subject-level angular distance data to a granularity feasible for aggregation into a single group-level distribution (across faces, TRs, and subjects). This script also captures the cross-subject mode of these out directionality bins for each cortical face, allowing for projection of these data onto the cortex. Note that for the modal plot, prominence of the top mode (height relative to non-neighboring modes) is also fed into [Vis_FaceVec_modes.m](https://github.com/PennLINC/PWs/blob/main/scripts/Vis_FaceVec_modes.m). Finally, the mean proportion of propagations that are top-down vs. bottom-up at each face are fed into [Vis_FaceVec_Bup.m)](https://github.com/PennLINC/PWs/blob/main/scripts/Vis_FaceVec_Bup.m)

The above script also supports setting "subject" as a variable (character string), and running these visualization scripts on individual-level data for modes and proportions.

# 5. Task effects

[BUTD_to_Rformat_c.m](https://github.com/PennLINC/PWs/blob/main/scripts/BUTD_to_Rformat_c.m)

facewise_ttest_([L](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_L.R)/[R](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_R.R)).R: this will conduct mass univariate tests on the extracted top-down proportion metrics. Follow it up with [FDR_facewise_t.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise_t.R) to correct for multiple comparisons

- taskEffects.md - from .rmd

# 6. Age effects

[Bin_And_Aggregate_BuProp_180.m](https://github.com/PennLINC/PWs/blob/e9dafb2118b3814b971d74d9cab4c35c51916aa2/scripts/Bin_And_Aggregate_BuProp_180.m) will derive individual-level whole-cortex metrics for proportion bottom-up (or top-down). 

[BUTD_to_Rformat.m](https://github.com/PennLINC/PWs/blob/main/scripts/BUTD_to_Rformat.m) will convert these values to .csv files for R to leverage

facewise_stats_([L](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_stats_L.R)/[R](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_stats_R.R)).R: this will conduct mass univariate tests on the extracted top-down proportion metrics. Follow it up with [FDR_facewise.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise.R) to correct for multiple comparisons
 
- ageEffects.md - from .rmd

###### note 1: separate processing stream for carit is facewise_ttest_(L/R) (and FDR), paired t-test for rest vs. task carit















