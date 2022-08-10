addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
% COMPLETE FILE PATH UPON OUTPUT
outdirgp='/cbica/projects/pinesParcels/results/PWs/FC/gro_CFC_';
% tertile list
tertiles={'young','old','mid'}

surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;

% for each tertile
for T=1:3;
	tertile=tertiles{T}
	subjs=fileread(['~/PWs/' tertile '_subs.txt']);
	% matlab is picky eater
	subjs=strsplit(subjs);
	subjs=subjs(1:(length(subjs)-1));
	% Make empty face-level CFC data matrix (# faces in mask,fsaverage4)
	CFCmatrix_l=zeros(4589,4589);
	CFCmatrix_r=zeros(4595,4595);
	for s=1:length(subjs);
		subj=subjs(s)
		% load circular FC matrix
		fw_CFCfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj{:} '/' subj{:} '_CircFC.mat'];
		fw_CFC=load(fw_CFCfp);
		% add em in
		CFCmatrix_l=CFCmatrix_l+fw_CFC.adjMats.L;
		CFCmatrix_R=CFCmatrix_r+fw_CFC.adjMats.R;	
	end
	% divide by subject number to average
	CFCmatrix_l=CFCmatrix_l./length(subjs);
	CFCmatrix_r=CFCmatrix_r./length(subjs);
	% write em out
	outdirgp_t_l=[outdirgp tertile '_l.csv'];
	outdirgp_t_r=[outdirgp tertile '_r.csv'];
	% save tertile FC average to subjdir
	csvwrite(outdirgp_t_l,CFCmatrix_l)
	csvwrite(outdirgp_t_r,CFCmatrix_r)
end
