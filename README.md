# PWs

# 1. Pre-processing
[preProc_PW.m](https://github.com/PennLINC/PWs/blob/main/scripts/preProc_PW.m) will motion mask subject time series, and generate indices of TRs in 10+ length continuous windows uninterrupted by FD > 0.2 / outlier frames. These windows (.txt files) are saved to projects/hcpd, not projects/pinesParcels at the moment. This script then calls [downsample_TS.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_TS.sh) to downsample this timeseries from the original spatial resolution (fslr). Finally, this script will compute optical flow on the freesurfer sphere (with [OpFl_Sph_CompVer.m](https://github.com/PennLINC/PWs/blob/main/scripts/OpFl_Sph_CompVer.m) ) , calculate angular distance from the PGG at each face for each vector (with [PGG_AngDistCalc.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc.m) ), and then re-mask out the medial wall for the subject (with [mask_mw_faces](https://github.com/PennLINC/PWs/blob/main/scripts/mask_mw_faces.m) )

###### note 1: separate processing stream exists for carit data: should be exactly equivalent, but have "\_c\" after script names to differentiate
###### note 2: OpFl_Sph and PGG_AngDistCalc are run as compiled progams (more specifically, from [here](https://github.com/PennLINC/PWs/blob/main/scripts/run_OpFl_Sph_CompVer.sh) and [here](https://github.com/PennLINC/PWs/blob/main/scripts/run_PGG_AngDistCalc4_CompVer.sh), respectively.) They are compiled directly from the linked files, but note the scripts used for large SGE submissions are derived from those linked. Expect non-compiled versions to run 3-5 times slower.
 
# 2.1 Reference map feature extraction
The PGG needs to be derived from a downsampled PG, and 1000 null PGGs need to be derived from spun PGs. The PG is first downsampled to fsaverage4 ([downsample_gPG4.sh](https://github.com/PennLINC/PWs/blob/main/scripts/downsample_gPG4.sh)). Next, the downsampled PG is converted from a gifti to .mat ([PGfuncgii_2_mat.m](https://github.com/PennLINC/PWs/blob/main/scripts/PGfuncgii_2_mat.m)), as subsequent compiled matlab code cannot utilize the gifti() command properly without xml error. The equivalent script is run for the curvature map in [Cfuncgii_2_mat.m](https://github.com/PennLINC/PWs/blob/main/scripts/Cfuncgii_2_mat.m)

# 2.2 Reference map spinning + feature extraction, prep for spatial null
To create null maps, the [spin test](https://github.com/spin-test/spin-test) is adapted: maps are spun and then gradients are calculated on spun maps to obtain spun vector fields. See [SpinPG.m](https://github.com/PennLINC/PWs/blob/main/scripts/SpinPG.m) for implementation of spins, and [AggSpins_TakeGrad.m](https://github.com/PennLINC/PWs/blob/main/scripts/AggSpins_TakeGrad.m) for conversion of spun PGs to spun PGGs.

# 3.1 Spatial Nulls
Because the angular distributions over each face over each TR over each subject and over spun reference angles is too much to store in one data frame, we will run the dip test in matlab iteratively as it runs over each spin. The "true" dip statistic is recorded last. Spatial null models are evaluated at the single-subject level, but because the PGG is a group-level vector field, we can use the 1000 rotations generated above for all subjects. Specifically, we find the angular distribution of optical flow vectors versus spun PGGs and evaluate [the dip statistic](https://github.com/diazlab/chance/blob/master/hartigansdiptest.m) (from [this paper](https://www.jstor.org/stable/2347485?seq=1)) from those 1000 null distributions. Our spatial permutation test can then be evaluated as the dip statistic from the true distribution (angles relative to PGG) versus from the 1000 spatial permutations of the PGG.

This is a computationally intensive script. Consequently, the code is compiled. You can observe the pre-compilation code [here](https://github.com/PennLINC/PWs/blob/main/scripts/PGG_AngDistCalc_snull_CompVer.m). This script was compiled with this command (from the scripts directory):

> mcc -m PGG_AngDistCalc_snull_CompVer.m -R -singleCompThread -R -nodisplay -R -nojvm -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/lib/freesurfer/matlab/read_surf -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/ofd/util/grad -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/src/read_medial_wall_label -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/hermans-rasson/HartigansDipTest

Subsequently, individual subjects can be qsub'ed with the resulting run_PGG_AngDistCalc_snull_CompVer.sh command. Specifically, you will want to open a text file, print run_PGG_AngDistCalc_snull_CompVer.sh $MATLAB_DIR {subject ID} to a .sh file, and then qsub that .sh file with at least 13GB of virtual memory. See this [script](https://github.com/PennLINC/PWs/blob/main/scripts/preProc_PW2.m) fo ran example.

After the qsub has completed, results can be easily plotted in R. Note that the subject ID is not saved in [this script](https://github.com/PennLINC/PWs/blob/main/scripts/NullPlotting.Rmd), you have to enter it manually (to avoid uploading subj IDs). For individual-level plotting you can start at line 133.

A parallel set of scripts exists for the curvature gradient. [This](https://github.com/PennLINC/PWs/blob/main/scripts/CG_AngDistCalc_snull_CompVer.m) is the key one, analagous to PGG_ANgDistCalc_snull_CompVer.m.

items to validate in 3.1:
#### No L/R mislabels

# 3.2 Temporal Nulls
Unfortunately there is no shared temporal aspect between subjects. This means to generate the temporal null, each operation needs to be entirely executed on a subject-by-subject basis. It will take a day or two (maybe three) to run for a single subject, even when compiled/qsub'ed.

The idea here is we shuffle the cleaned time series for each subject. We can then see if resultant optical flow vectors adhere to hierarchical bimodality more than expected when we have a bunch of images of the same brain together, irrespective of when in time they occurred. 

So for each subject, we shuffle their cleaned time series in 100 different ways, run optical flow on each, and similarly compare the resultant dip statistics from those angular distributions versus the true dip statistic.

Because the saveout for 9k faces x ~1000 TR pairs x 100 temporal permutations is huge (hundreds of GB per subject), we run this whole thing in one shot rather than save intermediate files. So "combined" here refers to lumping optical flow calculation and angular-distance-from-pgg-calculation together into [one script](https://github.com/PennLINC/PWs/blob/main/scripts/tnull_comb_CompVer.m), in parts 1 and 2. If you've followed this guide chronologically, the only new component should be the shuffle vectors and shuffle matrix, from lines 71-105.

This one was compiled with:
>mcc -m tnull_comb_CompVer.m -R -singleCompThread -R -nodisplay -R -nojvm -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/lib/freesurfer/matlab/MRIread -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/lib/freesurfer/matlab/read_surf -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/ofd/of.m -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/ofd/util/linearsystem.m -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/Code_mvNMF_l21_ard_v3_release/src/read_medial_wall_label -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/hermans-rasson/HartigansDipTest.m -a /cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/ofd/util/grad.m

This output can be qsubbed the same way as above, but requires 100 GB of virtual memory (unless subjects are at the lower end of TRs retained). It can also be visualized using the same r code as above, but note that we have 100 permutations and value 101 is the true dip stat, rather than 1000/1001 with the spatial permutations.

items to validate in 3.2:
#### shuffle behaving properly, no L/R flips

# 4.1 Extract Directional info from OpFlow output
[Extract_BUTD_ResultantVecs.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs.m) will calculate the TRs for which signal is hierarchically ascending vs. descending, and saveout subj metrics. Task data can be processed the same way with an equivalent script, [Extract_BUTD_ResultantVecs_c.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs_c.m). *Curvature* data can also be processed the same way, also with an equivalent script [Extract_BUTD_ResultantVecs_curv.m](https://github.com/PennLINC/PWs/blob/main/scripts/Extract_BUTD_ResultantVecs_curv.m). 

For replicating the subject data depicted in the schematic figure, we can use [Viz_VecFields.m](https://github.com/PennLINC/PWs/blob/main/scripts/Viz_VecFields.m).

# 4 Group-level directionality
[aggregate_OpFlDistrModes.m](https://github.com/PennLINC/PWs/blob/main/scripts/aggregate_OpFlDistrModes.m) segments subject-level angular distributions into 18 angular distance bins. This reduction is needed to reduce subject-level angular distance data to a granularity feasible for aggregation into a single group-level distribution (across faces, TRs, and subjects). This script also captures the cross-subject mode of these out directionality bins for each cortical face, allowing for projection of these data onto the cortex. Note that for the modal plot, prominence of the top mode (height relative to non-neighboring modes) is also fed into [Vis_FaceVec_modes.m](https://github.com/PennLINC/PWs/blob/main/scripts/Vis_FaceVec_modes.m). 

Finally, the mean proportion of propagations that are top-down vs. bottom-up at each face are fed into [Vis_FaceVec.m)](https://github.com/PennLINC/PWs/blob/main/scripts/Vis_FaceVec.m)

The above script also supports setting "subject" as a variable (character string), and running these visualization scripts on individual-level data for modes and proportions.

items to validate in 4:
#### Only top-down v. bottom-up proportion: extracted from calculation in Extract_BUTD_ResultantVecs.m

# 5. Task effects
First, we convert the output of Extract_BUTD_ResultantVecs.m to .csv's rather than .m files so R can use 'em: [BUTD_to_Rformat_c.m](https://github.com/PennLINC/PWs/blob/main/scripts/BUTD_to_Rformat_c.m)

Next, we run facewise t-tests to look at rest vs. task proportions. This is seperated into 2, qsubable scripts ([left](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_L.R) and [right](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_ttest_R.R) hemispheres).

Follow it up with [FDR_facewise_t.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise_t.R) to correct for multiple comparisons.

The resulting (masked) delta r^2 values can then be visualized with Vis_FaceVec.m.

Note the sample size here is different (smaller). Not all subjects had high-quality Carit data, see lines 114:156 of this R [script](https://github.com/PennLINC/PWs/blob/main/scripts/Group_level_analysis.rmd) for more sample construction information.

items to validate in 4:
#### T-test directionality: i.e., more top-down w/ task

# 6. Age effects

[Bin_And_Aggregate_BuProp_180.m](https://github.com/PennLINC/PWs/blob/e9dafb2118b3814b971d74d9cab4c35c51916aa2/scripts/Bin_And_Aggregate_BuProp_180.m) will derive individual-level whole-cortex metrics for proportion bottom-up (or top-down). Note this is in addition to mass univariate testing.

[BUTD_to_Rformat.m](https://github.com/PennLINC/PWs/blob/main/scripts/BUTD_to_Rformat.m) will convert these values to .csv files for R to leverage

facewise_stats_([L](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_stats_L.R)/[R](https://github.com/PennLINC/PWs/blob/main/scripts/facewise_stats_R.R)).R: this will conduct mass univariate tests on the extracted top-down proportion metrics. Follow it up with [FDR_facewise.R](https://github.com/PennLINC/PWs/blob/main/scripts/FDR_facewise.R) to correct for multiple comparisons

FINALLY, we need to create the permuted old - young PGG angular distance difference histogram (fig. 3d). These 1000 permutations are ran through [this script](https://github.com/PennLINC/PWs/blob/e2258b09c2ae78bc3a998faa59a283723d3cb085/scripts/Bin_And_Aggregate_PGGDistribution_180_permSubjs.m), which uses the old and young tertile splits (for the true histogram differences only) constructed in this R [script](https://github.com/PennLINC/PWs/blob/main/scripts/Group_level_analysis.rmd) from lines 101 to 110. The resulting figure is plotted with the same .rmd in the chunk beneath (lines 163:231), which can be run after permSubjs.m has ran.


And that's all I have to say about that














