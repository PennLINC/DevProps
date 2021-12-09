# PWs
Pws in HCPD

1. preProc_PW.m will motion mask subject time series, and generate indices of TRs in long (long is user-defined) continuous windows uninterrupted by FD > 0.2 / outlier frames. Next, it will downsample the functional time series for feasible diffusion map embedding, which is then the subsequent step. As is, PGs are upsampled back to original resolution. The gradient (calculus) is then taken of the PG (dmap) to determine directionality of hierarchical ascent. Finally, this script uses OpFl.m to delineate phase directionality on downsampled flatmaps of the time series.

2. permut_PGs.m will create 1000 spun PGs from the subject's delineated PG. The gradient of spun PGs will then be derived for null models of hierarchical directionality. This is a lot of memory to save (10GB per subject), so the 1000 brain maps will be subsequently deleted. Spun PG gradients will be stored in a 1000 * 64 * 89 * 2 array. 1000 for spins, 64 and 89 for current flatmap res, 2 for x and y components. This script invokes the same upsampling and downsampling procedures.

3. Test_AngDist.m will evaluate angular distance of optical flow directionality versus null permutations. 
