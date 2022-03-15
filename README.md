# PWs

# 1. Pre-processing
preProc_PW.m will motion mask subject time series, and generate indices of TRs in 10+ length continuous windows uninterrupted by FD > 0.2 / outlier frames. These windows (.txt files) are saved to projects/hcpd, not projects/pinesParcels at the moment. Finally, this script will compute optical flow on the fs4 sphere, calculate angular distance from the PGG at each face for each vector, and then mask out the medial wall for the subject.

items to validate in preProc_PW:

a. apply_motion_mask: txt files generated for continuous segments is accurate
b. OpFl_Sph_fs4: angular transformations are accurate
c. PGG_AngDistCalc4: Angular distance calculated accurately

 -note sep. proc. stream exists for carit data: should be exactly equivalent, but have _c after filenames to differentiate

# 2. Extract Directional info from OpFlow output
Extract_BUTD_ResultantVecs.m will calculate the TRs for which signal is hierarchically ascending vs. descending, and saveout subj metrics. BUTD_to_Rformat.m will convert these values to .csv files for R to leverage.

# 3. Group-level directionality

Bin_And_Aggregate_PGGDistribution_180.m segments subject-level angular distributions into 18 angular distance bins. This reduction is used to aggregate a group-level angular distance histogram. Further, this script savesout normative, binned directionality for each cortical face.

Bin_And_Aggregate_PGGDistribution_360.m segments subject-level angular distributions into 36 angular distance bins. This reduction is used to aggregate a group-level angular distance histogram. Further, this script savesout normative, binned directionality for each cortical face.

# 4. Task effects

facewise_ttest_(L/R).R: this will conduct mass univariate tests on the extracted top-down proportion metrics. Follow it up with FDR_facewise_t.R to correct for multiple comparisons

- taskEffects.md - from .rmd

# 5. Age effects

facewise_stats_(L/R).R: this will conduct mass univariate tests on the extracted top-down proportion metrics. Follow it up with FDR_facewise.R to correct for multiple comparisons
 
- ageEffects.md - from .rmd

-note sep. proc stream for carit is facewise_ttest_(L/R) (and FDR), paired t-test for rest vs. task carit















