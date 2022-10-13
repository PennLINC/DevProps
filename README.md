# Development of Propagations: Replication Guide

# 1. Pre-processing
[preProc_PW.m](https://github.com/PennLINC/PWs/blob/main/scripts/preProc_PW.m) will motion mask subject time series, and generate indices of TRs in 10+ length continuous windows uninterrupted by FD > 0.2 / outlier frames. These windows (.txt files) are saved to projects/hcpd, not projects/pinesParcels at the moment. This script then calls [downsample_TS.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_TSfs4.sh) to downsample this timeseries from the original spatial resolution (fslr). Finally, this script will compute optical flow on the freesurfer sphere (with [OpFl_Sph_CompVer.m](https://github.com/PennLINC/PWs/blob/main/scripts/OpFl_Sph_CompVer.m) ) , calculate angular distance from the PGG at each face for each vector (with [PGG_AngDistCalc.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc.m) ), and then re-mask out the medial wall for the subject (with [mask_mw_faces](https://github.com/PennLINC/PWs/blob/main/scripts/mask_mw_faces_4.m) )

###### note 1: separate processing stream exists for carit data: should be exactly equivalent, but have "\_c\" after script names to differentiate
###### note 2: OpFl_Sph and PGG_AngDistCalc are run as compiled progams (more specifically, from [here](https://github.com/PennLINC/PWs/blob/main/scripts/run_OpFl_Sph_CompVer.sh) and [here](https://github.com/PennLINC/PWs/blob/main/scripts/run_PGG_AngDistCalc4_CompVer.sh), respectively.) They are compiled directly from the linked files, but note the scripts used for large SGE submissions are derived from those linked. Expect non-compiled versions to run 3-5 times slower.
 
# 2.1 Reference map feature extraction
The PGG needs to be derived from a downsampled PG, and 1000 null PGGs need to be derived from spun PGs. The PG is first downsampled to fsaverage4 ([downsample_gPG4.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_gPG4.sh)). Next, the downsampled PG is converted from a gifti to .mat ([PGfuncgii_2_mat.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGfuncgii_2_mat.m)), as subsequent compiled matlab code cannot utilize the gifti() command properly without xml error.

# 2.2 Reference map spinning + feature extraction, prep for spatial null
To create null maps, the [spin test](https://github.com/spin-test/spin-test) is adapted: maps are spun and then gradients are calculated on spun maps to obtain spun vector fields. See [SpinPG.m](https://github.com/PennLINC/PWs/blob/main/scripts/SpinPG.m) for implementation of spins, and [AggSpins_TakeGrad.m](https://github.com/PennLINC/PWs/blob/main/scripts/AggSpins_TakeGrad.m) for conversion of spun PGs to spun PGGs.

# 3.1 Spatial Nulls
Because the angular distributions over each face over each TR over each subject and over spun reference angles is too much to store in one data frame, we will run the dip test in matlab iteratively as it runs over each spin. The "true" dip statistic is recorded last. Spatial null models are evaluated at the single-subject level, but because the PGG is a group-level vector field, we can use the 1000 rotations generated above for all subjects. Specifically, we find the angular distribution of optical flow vectors versus spun PGGs and evaluate [the dip statistic](https://github.com/diazlab/chance/blob/master/hartigansdiptest.m) (from [this paper](https://www.jstor.org/stable/2347485?seq=1)) from those 1000 null distributions. Our spatial permutation test can then be evaluated as the dip statistic from the true distribution (angles relative to PGG) versus from the 1000 spatial permutations of the PGG.

This is a computationally intensive script. Consequently, the code is compiled. You can observe the pre-compilation code [here](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc_snull_CompVer.m). This script was compiled with this command (from the scripts directory):

> mcc -m PGG_AngDistCalc_snull_CompVer.m -R -singleCompThread -R -nodisplay -R -nojvm -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/lib/freesurfer/matlab/read_surf -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/ofd/util/grad -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/src/read_medial_wall_label -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/hermans-rasson/HartigansDipTest

Subsequently, individual subjects can be qsub'ed with the resulting run_PGG_AngDistCalc_snull_CompVer.sh command. Specifically, you will want to open a text file, print run_PGG_AngDistCalc_snull_CompVer.sh $MATLAB_DIR {subject ID} to a .sh file, and then qsub that .sh file with at least 13GB of virtual memory. See this [script](https://github.com/PennLINC/PWs/blob/main/scripts/preProc_PW2.m) fo ran example.

After the qsub has completed, results can be easily plotted in R. Note that the subject ID is not saved in [this script](https://github.com/PennLINC/PWs/blob/main/scripts/NullPlotting.Rmd), you have to enter it manually (to avoid uploading subj IDs). For individual-level plotting you can start at line 133.

# 3.2 Temporal Nulls
Unfortunately there is no shared temporal aspect between subjects. This means to generate the temporal null, each operation needs to be entirely executed on a subject-by-subject basis. It will take a day or two (maybe three) to run for a single subject, even when compiled/qsub'ed.

The idea here is we shuffle the cleaned time series for each subject. We can then see if resultant optical flow vectors adhere to hierarchical bimodality more than expected when we have a bunch of images of the same brain together, irrespective of when in time they occurred. 

So for each subject, we shuffle their cleaned time series in 100 different ways, run optical flow on each, and similarly compare the resultant dip statistics from those angular distributions versus the true dip statistic.

Because the saveout for 9k faces x ~1000 TR pairs x 100 temporal permutations is huge (hundreds of GB per subject), we run this whole thing in one shot rather than save intermediate files. So "combined" here refers to lumping optical flow calculation and angular-distance-from-pgg-calculation together into [one script](https://github.com/PennLINC/PWs/blob/main/scripts/tnull_comb_CompVer.m), in parts 1 and 2. If you've followed this guide chronologically, the only new component should be the shuffle vectors and shuffle matrix, from lines 71-105.

This one was compiled with:
>mcc -m tnull_comb_CompVer.m -R -singleCompThread -R -nodisplay -R -nojvm -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/lib/freesurfer/matlab/MRIread -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/lib/freesurfer/matlab/read_surf -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/ofd/of.m -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/ofd/util/linearsystem.m -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/src/read_medial_wall_label -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/hermans-rasson/HartigansDipTest.m -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/ofd/util/grad.m

This output can be qsubbed the same way as above, but requires 100 GB of virtual memory (unless subjects are at the lower end of TRs retained). It can also be visualized using the same r code as above, but note that we have 100 permutations and value 101 is the true dip stat, rather than 1000/1001 with the spatial permutations.

# 4.1 Extract Directional info from OpFlow output
[Extract_BUTD_ResultantVecs.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs.m) will calculate the TRs for which signal is hierarchically ascending vs. descending, and saveout subj metrics. Task data can be processed the same way with an equivalent script, [Extract_BUTD_ResultantVecs_c.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs_c.m). 

For replicating the subject data depicted in the schematic figure, we can use [Viz_VecFields.m](https://github.com/PennLINC/PWs/blob/main/scripts/Viz_VecFields.m).

# 4.2 Group-level aggregation and map visualization
[aggregate_OpFlDistrModes.m](https://github.com/PennLINC/PWs/blob/main/scripts/aggregate_OpFlDistrModes.m) segments subject-level angular distributions into 18 angular distance bins. This reduction is needed to reduce subject-level angular distance data to a granularity feasible for aggregation into a single group-level distribution (across faces, TRs, and subjects). This script also captures the cross-subject mode of these out directionality bins for each cortical face, allowing for projection of these data onto the cortex.

Finally, the mean proportion of propagations that are top-down vs. bottom-up at each face are fed into [Vis_FaceVec.m)](https://github.com/PennLINC/PWs/blob/main/scripts/Vis_FaceVec.m)

The above script also supports setting "subject" as a variable (character string), and running these visualization scripts on individual-level data for modes and proportions.

# 5. Task effects
First, we convert the output of Extract_BUTD_ResultantVecs.m to .csv's rather than .m files so R can use 'em: [BUTD_to_Rformat_c.m](https://github.com/PennLINC/PWs/blob/main/scripts/BUTD_to_Rformat_c.m)

Next, we run facewise t-tests to look at rest vs. task proportions. This is seperated into 2, qsubable scripts ([left](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_L.R) and [right](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_R.R) hemispheres).

Follow it up with [FDR_facewise_t.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise_t.R) to correct for multiple comparisons.

The resulting (masked) delta r^2 values can then be visualized with Vis_FaceVec.m.

Note the sample size here is different (smaller). Not all subjects had high-quality Carit data, see lines 121:132 of this R [script](https://github.com/PennLINC/PWs/blob/main/scripts/Group_level_analysis.rmd) for more sample construction information.

# 6.1 Age effects: Whole-cortex
[Bin_And_Aggregate_BuProp_180.m](https://github.com/PennLINC/PWs/blob/e9dafb2118b3814b971d74d9cab4c35c51916aa2/scripts/Bin_And_Aggregate_BuProp_180.m) will derive individual-level whole-cortex metrics for proportion bottom-up (or top-down). This output can read into R and plotted starting at line 206 of [this script](https://github.com/PennLINC/PWs/blob/main/scripts/Group_level_analysis.rmd), and within [this markdown](https://github.com/PennLINC/PWs/blob/main/scripts/Group_level_analysis.md). Code for reported stats is directly adjacent to plotting code.

# 6.2 Age effects: Mass univariate
As above, [BUTD_to_Rformat.m](https://github.com/PennLINC/PWs/blob/main/scripts/BUTD_to_Rformat.m) will convert extracted values to .csv files for R to leverage.

facewise_stats_([L](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_stats_L.R)/[R](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_stats_R.R)).R: this will conduct mass univariate tests on the extracted top-down proportion metrics. Follow it up with [FDR_facewise.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise.R) to correct for multiple comparisons

# 6.3 Age effects: angular distributions
Here, we create the permuted old - young PGG angular distance difference histogram (fig. 3d). These 1000 permutations are ran through [this script](https://github.com/PennLINC/PWs/blob/e2258b09c2ae78bc3a998faa59a283723d3cb085/scripts/Bin_And_Aggregate_PGGDistribution_180_permSubjs.m), which uses the old and young tertile splits (for the true histogram differences only) constructed in this R [script](https://github.com/PennLINC/PWs/blob/main/scripts/Group_level_analysis.rmd) starting at line 139. The resulting figure is plotted with the same .rmd, which can be run after permSubjs.m has ran.

# 7. Simulated posterior-anterior data
Here, we simulate a posterior-to-anterior traveling wave across the cortical sheet, and recover the directionality of this wave with optical flow. First, to generate the wave, we use [syn_wave_prop1](https://github.com/PennLINC/PWs/blob/main/scripts/syn_wave_prop1.m) (thank you Manish!) to generate the wave. This script pulls from On the Stability of BOLD fMRI Correlations by [Laumann et al., 2017](https://cpb-us-w2.wpmucdn.com/sites.wustl.edu/dist/7/1080/files/2018/06/simulate_BOLD_timecourse_func_v2-1eoyny0.m) to ensure simulated data resembles real BOLD signal. Second, we downsample this data to fsaverage4 resolution with [downsample_SIM_TSfs4.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_SIM_TSfs4.sh). Third, we run optical flow on this downsampled data with [OpFl_Sph_fs4_SIM](https://github.com/PennLINC/PWs/blob/main/scripts/OpFl_Sph_fs4_SIM.m) : note that the subject "name" is "sub-Sim", as set in the downsample script. Fourth, we calculate the posterior-anterior directionality by use of the 3d coordinate system innate to freesurfer. The coordinates are extracted and converted to a .mat file in [APfuncgii_2_mat.m](https://github.com/PennLINC/PWs/blob/main/scripts/APfuncgii_2_mat.m). Fifth, the surface gradient (i.e., nabla Posterior Anterior) is extracted and angular distance between said gradient and the synthetic data is calculated in [APG_AngDistCalc4_CompVer.m](https://github.com/PennLINC/PWs/blob/main/scripts/APG_AngDistCalc4_CompVer.m). Sixth, the face-wise medial wall mask is dropped from this data in [mask_mw_faces_4_ap](https://github.com/PennLINC/PWs/blob/main/scripts/mask_mw_faces_4_ap.m). Seventh (lastly!), we visualize the angular distribution output. More specifically, we observe the angular distance between nabla posterior anterior and the intended wave direction for each point on the cortex. This is done in R, see [lines 157-191](https://github.com/PennLINC/PWs/blob/main/scripts/opfl_vis.Rmd).

# 8.1 Midnight scan club replication: Preprocessing

Here we're going to run the midnight scan club data through an analagous pipeline. The only real difference is a TR of 2.2 seconds rather than 0.8 seconds, but we'll concatenate across resting-state scans so we'll still have plenty of TRs for each participant.

The parent script for this step (8.1) is [preProc_PWs_msc](https://github.com/PennLINC/PWs/blob/main/scripts/preProc_PW_msc.m). The first script called within is [apply_motion_mask_msc](https://github.com/PennLINC/PWs/blob/main/scripts/apply_motion_mask_msc.m), which applies and records the extended temporal censoring used to only analyze continuous segments of uninterrupted fMR images (> 10 TRs). The only distinct between this script and the motion masking used for HCPD is the FD threshold. Because longer TRs yield higher FDs, the FD threshold was relaxed to .4mm FD.

After motion masking and saving out the valid TR indices, the time series are downsampled with [downsample_TSfs4_msc](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_TSfs4_msc.sh).

The downsampled time series are then ready for optical flow. This [script](https://github.com/PennLINC/PWs/blob/main/scripts/OpFl_Sph_fs4_runSweep_msc.m) will run optical flow on individual scans (preproc_PWs leverages the C-compiled verison, see section 3 for more info on how matlab code was compiled). Concatenation of resting-state data is performed _after_ optical flow with this [script](https://github.com/PennLINC/PWs/blob/main/scripts/Concat_mscFl.m).

# 8.2 Midnight scan club replication: Post-processing

[This script](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc4_CompVer_msc.m) will calculate angular distances on MSC data (angular distances from nabla PG). However, the stats we get and report come from [comparing angular sitances among msc optical flow output to null models (snull)](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc_snull_CompVer_msc.m), which is all lumped together for computational efficiency. Note that as for other CompVer scripts, the above script is meant to be [run as its c-compiled version](https://github.com/PennLINC/PWs/blob/main/scripts/run_PGG_AngDistCalc_snull_CompVer_msc.sh). You will have to uncomment the addpath command up at the top of .m file for line-by-line use. Also note that more virtual memory is needed for qsubbing MSC snull jobs, as the concatenated time series is longer.

That's about it for MSC. Now we get to return to the warm comfort of r/rstudio to plot the output, labeled as SpunDips4_cnct.csv for each MSC participant. This only takes a few lines of code from [NullPlotting](https://github.com/PennLINC/PWs/blob/main/scripts/NullPlotting.Rmd). Specifically, lines 102-108.

# 9.1 Alternative Hierarchy: Dip statistics versus null

Here we use an independently-derived hierarchy (nabla myelin) to substantiate the claim that alignment of propgagtions with nabla PG reflect bottom-up and top-down directions. The original myelin map is downloaded from [BALSA](https://balsa.wustl.edu/study/show/mDBP0), which was uploaded as part of the efforts behind [this paper](https://www.sciencedirect.com/science/article/pii/S1053811922004797). Specifically, we used the transmit-bias corrected map from HCP-D. This processing stream is almost entirely parallel to that described above using the principal gradient: the exception is that whereas higher loadings are higher-order cortices for the principal gradient, lower values are higher-order for the myelin map. This only requires a trivial adjustment later on.

The first step to using this map is to downsample it to fsaverage4 space with [downsample_gMM.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_gMM.sh). Next, this downsampled version [is converted to a .mat file](https://github.com/PennLINC/PWs/blob/main/scripts/myGfuncgii_2_mat.m) for easy loading within the compiled scripts below. 

Because we are only changing the set of _reference_ directions for measurement of hierarchical directionality, optical flow does not need to be recalculated. The next steps then becomes to measure angular distance relative to the alternative hierarchy rather than the originally utilized one. This is done in [MyG (Myelin gradient) AngDistCalc4_CompVer.m](https://github.com/PennLINC/PWs/blob/main/scripts/MyG_AngDistCalc4_CompVer.m). As prior, this script will calculate the angular distances, but we are more immediately concerned with the angular distances relative to conservative spatial null models. That is conducted with the compiled version of [this script](https://github.com/PennLINC/PWs/blob/main/scripts/MyG_AngDistCalc_snull_CompVer.m), using [this script](https://github.com/PennLINC/PWs/blob/main/scripts/run_MyG_AngDistCalc_snull_CompVer.sh).

After running spatiall null testing on all participants, we can use r code to plot the observed values and null distributions at the subject and group-level. For a single subject lines 101-113 of [NullPlotting](https://github.com/PennLINC/PWs/blob/main/scripts/NullPlotting.Rmd) will do the trick, but to plot at the group-level we have to go through and load each subject's data, null distribution, and distance of the observed dip statistic from the null distribution. This can be accomplished with lines 7-88 of [NullPlotting](https://github.com/PennLINC/PWs/blob/main/scripts/NullPlotting.Rmd).

# 9.2 Alternative Hierarchy: Task effects

As prior, we need to calculate angular distances on optical flow from the carit task sessions. That is done with [the compiled version of this script](https://github.com/PennLINC/PWs/blob/main/scripts/MyG_AngDistCalc_c_CompVer.m) using [this script](https://github.com/PennLINC/PWs/blob/main/scripts/run_MyG_AngDistCalc_c_CompVer.sh) 

Second, we need to extract the proportion of propagations noted to be top-down at each point on the cortex over this task. We can do this with [Extract_BUTD_ResultantVecs_My_c.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs_My_c.m). Note that this script is NOT compiled; you will be subjected to license limitations when using it. Fortunately it is quick, so even though you can only run it on ~6 participants at a time, it's not a big rate-limiting step.

Third, we need to convert the output of Extract_BUTD_ResultantVecs_My_c.m to .csv's rather than .m files so R can use 'em. This is what [BUTD_to_Rformat_My_c.m](https://github.com/PennLINC/PWs/blob/main/scripts/BUTD_to_Rformat_My_c.m) is for.

Fourth, we run facewise t-tests to look at rest vs. task proportions. This is seperated into 2, qsubable scripts ([left](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_L_My.R) and [right](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_R_My.R) hemispheres).

Fifth, use [FDR_facewise_t_My.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise_t_My.R) to correct for multiple comparisons.

Sixth, visualize the output with the general-purpose [Vis_Facevec.m script](https://github.com/PennLINC/PWs/blob/main/scripts/Vis_FaceVec.m). Use readtable to load in the output of FDR-correction. Note there's no header to these files! Accordingly, 'ReadVariableNames` should be set to false. Finally, recall that where higher loadings are higher-order cortices for the principal gradient, lower values are higher-order for the myelin map. Correct for this flip of hierarchical directionality with a simple (*-1) to the corrected t-test vectors prior to inputting into the visualization script.

Within [Vis_Facevec.m](https://github.com/PennLINC/PWs/blob/main/scripts/Vis_FaceVec.m), mincol should be set to -10 and maxcol should be set to 10. See line 85 for the blue-orange colormapping specifically.


# 9.3 Alternative Hierarchy: Developmental effects















