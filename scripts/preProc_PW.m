function preProc_PWs(subj)
% Apply motion and outlier masks, bandpass continuous segments, make PG from continuous segments, run optical flow

% print subject being ran
subj

% add matlab path for used functions
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));

% make output folder
direcString=['/cbica/projects/pinesParcels/results/wave_output/' subj];
mkdirCommand=['mkdir ' direcString];
system(mkdirCommand)

% apply motion masks and extract TR indcies of continuous segments
apply_motion_mask_genTRinds(subj)

%%%%%%%%%%%%%%%
% don't mess about with this stuff if it's already been ran
% fp=['/cbica/projects/abcdfnets/results/wave_output/' subj];
% PGfp=[fp '/' subj '_PG_LR_32k.func.gii'];
% spin tests require 10k surfaces in parent fp, so need to re-run this block on all subjs
% if ~exist(PGfp)
%%%%%%%%%%%%%%%

% downsample aggregated TS 
dsCommand=['~/PWs/scripts/downsample_TS.sh ' subj];
system(dsCommand)

% make a "Proced" dir
direcString=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj];
mkdirCommand=['mkdir ' direcString];
system(mkdirCommand)

% derive personalized PG
derivePGcommand=['/cbica/projects/pinesParcels/miniconda3/envs/mv_preds/bin/python derive_pg.py ' subj];
system(derivePGcommand)

% upsample for equivalent grid resampling
dsCommand=['~/PWs/scripts/upsample_PG.sh ' subj];
system(dsCommand)

% note distinct interpolation for time series for opfl (within opflow) and time series to fsaverage5 (within downsample_ts)
% consider equivalent interpolation 

%%% Run neuropattTB
OpFl(subj)

