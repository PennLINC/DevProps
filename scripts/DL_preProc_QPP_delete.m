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
subjDlCommand=['python3 /cbica/projects/abcdfnets/nda-abcd-s3-downloader/download.py -i /cbica/projects/abcdfnets/nda-abcd-s3-downloader/datastructure_manifest_722.txt -o /scratch/abcdfnets/nda-abcd-s3-downloader/August_2021_DL/ -s /cbica/projects/abcdfnets/nda-abcd-s3-downloader/' subj '.txt -l /cbica/projects/abcdfnets/nda-abcd-s3-downloader/March_2021_DL/dl_logs -d /cbica/projects/abcdfnets/nda-abcd-s3-downloader/func_subsets.txt &']

% note: downloader tool does not seem to communicate when it is done to matlab
% added '&' and 'pause' so that matlab waits 5 minutes to proceed rather than getting caught up indefinitely
system(subjDlCommand)
pause(320)

% now the matlab portions. Apply the motion mask to the downloaded data and extract global signal from the non-proc'ed scans
apply_motion_mask_extractGS_genTRinds(subj)

% bandpass the global signal and time series to isolate freqs of interest
BandPass_ts(subj)

% derive FC - grayOrd level and network level
deriveFCcommand=['python derive_fc.py ' subj];
system(deriveFCcommand)

% downsample aggregated TS 
dsCommand=['~/scripts/PWs/PWs/scripts/downsample_FC.sh ' subj];
system(dsCommand)

% derive personalized PG
derivePGcommand=['python derive_pg.py ' subj];
system(derivePGcommand)

% derive qpp
deriveQPPcommand=['python derive_pw.py' subj];
system(deriveQPPcommand)

% upsample QPP frames
usQPPcommand=['~/scripts/PWs/PWs/scripts/upsample_PW.sh' subj];
system(usQPPcommand)

% concat frames
cQPPcommand=['python concat_pw.py' subj];
system(cQPPcommand)

% print out downsampled instances
vQPPcommand=['python ~/scripts/PWs/PWs/scripts/Viz_PWs_surf.py' subj];
system(vQPPcommand)

% upsample derived principal gradient
usCommand=['~/scripts/PWs/PWs/scripts/upsample_PG.sh ' subj];
system(usCommand)

% derive wave properties w/ python
wavePropCommand=['python derive_WaveProps.py ' subj];
system(wavePropCommand)

% group PG rather than individualized
wavePropCommandG=['python derive_WaveProps_gPG.py ' subj];
system(wavePropCommandG)

% basis time series wave prop - coarse
wavePropCommandBTSC=['python derive_WaveProps_BasisTS.py ' subj];
system(wavePropCommandBTSC)

% basis time series wave prop - fine
wavePropCommandBTSF=['python derive_WaveProps_BasisTS_7.py ' subj];
system(wavePropCommandBTSF)

% delete input data
Delete_input_data_w(subj)
