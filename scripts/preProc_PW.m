function preProc_PWs(subj)
% Apply motion and outlier masks, bandpass continuous segments, make PG from continuous segments, run optical flow

% print subject being ran
subj

% add matlab path for used functions
addpath(genpath('/cbica/projects/hcpd/scripts/code_nmf_cifti/tool_folder'));

% make output folder
direcString=['/cbica/projects/pinesparcels/results/wave_output/' subj];
mkdirCommand=['mkdir ' direcString];
system(mkdirCommand)

% apply motion masks and extract TR indcies of continuous segments
apply_motion_mask_extractGS_genTRinds(subj)

% bandpass the global signal and time series to isolate freqs of interest
BandPass_ts(subj)

% don't mess about with this stuff if it's already been ran
% fp=['/cbica/projects/abcdfnets/results/wave_output/' subj];
% PGfp=[fp '/' subj '_PG_LR_32k.func.gii'];
% spin tests require 10k surfaces in parent fp, so need to re-run this block on all subjs
% if ~exist(PGfp)

% derive FC - grayOrd level and network level
deriveFCcommand=['python derive_fc.py ' subj];
system(deriveFCcommand)

% downsample aggregated TS 
dsCommand=['~/scripts/PWs/PWs/scripts/downsample_FC.sh ' subj];
system(dsCommand)

% derive personalized PG
derivePGcommand=['python derive_pg.py ' subj];
system(derivePGcommand)

%%% Run neuropattTB
OpFl(subj)

