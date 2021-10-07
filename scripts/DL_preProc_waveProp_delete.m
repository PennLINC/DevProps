function DL_preProc_waveProp_delete(subj)
%%% This function will take a single subject's NDAR name, download their fMRI data and motion masks, concatenate the fMRI data, motion mask at 2mm FD, bandpass filter the global signal and the proc'ed data, derive global wave properties, and delete the input fMRI data.

% print subject being ran
subj

% add matlab path for used functions
addpath(genpath('/cbica/projects/abcdfnets/scripts/code_nmf_cifti/tool_folder'));

% tell it where AWS tools are for downloads
system('export PATH=/cbica/projects/abcdfnets/aws/dist/:$PATH')

% echo subj name into a .txt
subjTxtCommand=['echo ' subj ' >> /cbica/projects/abcdfnets/nda-abcd-s3-downloader/' subj '.txt'];
system(subjTxtCommand)

% make output folder
direcString=['/cbica/projects/abcdfnets/results/wave_output/' subj];
mkdirCommand=['mkdir ' direcString];
system(mkdirCommand)

% download that one subject's data
subjDlCommand=['python3 /cbica/projects/abcdfnets/nda-abcd-s3-downloader/download.py -i /cbica/projects/abcdfnets/nda-abcd-s3-downloader/datastructure_manifest_10_2.txt -o /scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/ -s /cbica/projects/abcdfnets/nda-abcd-s3-downloader/' subj '.txt -l /cbica/projects/abcdfnets/nda-abcd-s3-downloader/March_2021_DL/dl_logs -d /cbica/projects/abcdfnets/nda-abcd-s3-downloader/func_subsets.txt &']

% note: downloader tool does not seem to communicate when it is done to matlab
% added '&' and 'pause' so that matlab waits 5 minutes to proceed rather than getting caught up indefinitely
system(subjDlCommand)
pause(320)

% now the matlab portions. Apply the motion mask to the downloaded data and extract global signal from the non-proc'ed scans
apply_motion_mask_extractGS_genTRinds(subj)

% bandpass the global signal and time series to isolate freqs of interest
BandPass_ts(subj)

% don't mess about with this stuff if it's already been ran
fp=['/cbica/projects/abcdfnets/results/wave_output/' subj];
PGfp=[fp '/' subj '_PG_LR_32k.func.gii'];
if ~exist(PGfp)
	% derive FC - grayOrd level and network level
	deriveFCcommand=['python derive_fc.py ' subj];
	system(deriveFCcommand)

	% downsample aggregated TS 
	dsCommand=['~/scripts/PWs/PWs/scripts/downsample_FC.sh ' subj];
	system(dsCommand)

	% derive personalized PG
	derivePGcommand=['python derive_pg.py ' subj];
	system(derivePGcommand)

	% upsample derived principal gradient
	usCommand=['~/scripts/PWs/PWs/scripts/upsample_PG.sh ' subj];
	system(usCommand)
end

% spin the pg
spin_pg(subj);

% load in bigRot for spin delin
parentfp=['/scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/'];
spunFn=[parentfp subj '_spunions.mat']; 
bigRot=load(spunFn);

% read template PG for iteratively replacing with spun, subj-specific PG
PGL10k=gifti(['/cbica/projects/abcdfnets/data/hcp.gradients_L_10k.func.gii']);
PGR10k=gifti(['/cbica/projects/abcdfnets/data/hcp.gradients_R_10k.func.gii']);

% convert spins to dscalar and derive waveProps iteratively
for s=1:100
	% extract spin s
	rotL=bigRot.bigrotl(s,:);
	rotR=bigRot.bigrotr(s,:);
	PGL10k.cdata(:,1)=rotL';
	PGR10k.cdata(:,1)=rotR';
	write_cifti(PGL10k,[parentfp subj '_spunPGL.func.gii');
	write_cifti(PGR10k,[parentfp subj '_spunPGR.func.gii');
	% combinedRot=[rotL rotR];
	% upsample
	usSpunCommand=['~/scripts/PWs/PWs/scripts/upsample_PG_spun.sh ' subj];	
	% derive waveprops on spin
	SpinwavePropCommand=['python derive_WaveProps_Spun.py ' subj];
	system(SpinwavePropCommand)
	% record results of interest from this spin
end

% derive wave properties w/ python
wavePropCommand=['python derive_WaveProps.py ' subj];
system(wavePropCommand)

% deactivation propagation 
wavePropCommandS=['python derive_WaveProps_shadow.py ' subj];
system(wavePropCommandS)

% aggregate subject-level measures
MagSpeedPhaseTimeCommand=['Rscript derive_MagSpeedPhaseTimeinPW.R ' subj];
system(MagSpeedPhaseTimeCommand)

% shadow
MagSpeedPhaseTimeCommandS=['Rscript derive_MagSpeedPhaseTimeinPW_shadow.R ' subj];
system(MagSpeedPhaseTimeCommandS)

% group PG rather than individualized - not updated
% wavePropCommandG=['python derive_WaveProps_gPG.py ' subj];
% system(wavePropCommandG)

% basis time series wave prop - coarse
% wavePropCommandBTSC=['python derive_WaveProps_BasisTS.py ' subj];
% system(wavePropCommandBTSC)

% basis time series wave prop - fine
% wavePropCommandBTSF=['python derive_WaveProps_BasisTS_7.py ' subj];
% system(wavePropCommandBTSF)

% delete input data
Delete_input_data_w(subj)
