addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
ProjectFolder = '/cbica/projects/pinesParcels/data/SingleParcellation';
% COMPLETE FILE PATH UPON OUTPUT
outdirgp='/cbica/projects/pinesParcels/results/PWs/FC/gro_FC_';
% tertile list
tertiles={'young','old','mid'}
% load in medial walls
surfML = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "isn't medial wall" vector for vertices
mw_L=ones(1,10242);
mw_L(mwIndVec_l)=0;
mw_R=ones(1,10242);
mw_R(mwIndVec_r)=0;
% initialize All-FC matrix
FCmatrix=zeros(20484,20484);
% initialize subj count
SubjNum=0;
% for each tertile
for T=1:3;
	tertile=tertiles{T}
	subjs=fileread(['~/PWs/' tertile '_subs.txt'])
	% matlab is picky eater
	subjs=strsplit(subjs);
	subjs=subjs(1:(length(subjs)-1));
	% Make empty vertex-level FC data matrix (# vert in mask,fsaverage5)
	FCmatrices=zeros(20484,20484);
	for s=1:length(subjs);
		SubjNum=SubjNum+1;
		subj=subjs(s)
	        vw_ts_l_p=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_L_10k.mgh'];
		vw_ts_r_p=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_R_10k.mgh'];
		% this is for sure some fussy matlab bs
		vw_ts_l_p=strcat(vw_ts_l_p{1},vw_ts_l_p{2},vw_ts_l_p{3},vw_ts_l_p{4},vw_ts_l_p{5});
		vw_ts_r_p=strcat(vw_ts_r_p{1},vw_ts_r_p{2},vw_ts_r_p{3},vw_ts_r_p{4},vw_ts_r_p{5});
		vw_ts_l=MRIread(vw_ts_l_p);
		vw_ts_r=MRIread(vw_ts_r_p);
		% series length
		vw_ts_l_SIZE=size(vw_ts_l.vol);
		vw_ts_l_LEN=vw_ts_l_SIZE(4)
		vw_ts_l=vw_ts_l.vol;
		vw_ts_r=vw_ts_r.vol;
		% stacking matrices so vertex number is doubled (not timepoints)
		vw_ts_both=[vw_ts_l vw_ts_r];
		vw_ts_bothrw=zeros(vw_ts_l_LEN,20484);
		% calculate FC
		for x=1:length(vw_ts_bothrw)
			vw_ts_bothrw(:,x)=vw_ts_both(1,x,1,:);
		end
		% bigass connectivity matrix, takes 5 seconds or so to calc
		ba_conmat=corrcoef(vw_ts_bothrw);
		FCmatrices=FCmatrices+ba_conmat;
		FCmatrix=FCmatrix+ba_conmat;	
		% save subject specific fc matrix
		outfn=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj{:} '/' subj{:} '_FC.csv'];
		csvwrite(outfn,ba_conmat);
	end
	% divide by subject number to average
	FCmatrices=FCmatrices./length(subjs);
	outdirgp_t=[outdirgp tertile '.csv'];
	% save tertile FC average to subjdir
	csvwrite(outdirgp_t,FCmatrices)
end
% average group FC matrix
FCmatrix=FCmatrix./SubjNum;
csvwrite('/cbica/projects/pinesParcels/results/PWs/FC/gro_FC.csv',FCmatrix);

