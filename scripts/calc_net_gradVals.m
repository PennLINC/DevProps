% add matlab path for used functions
addpath(genpath('/cbica/projects/abcdfnets/scripts/code_nmf_cifti/tool_folder'));
% load in gradients
ProjectFolder = '~/data';
pgl = gifti([ProjectFolder '/Gradients.lh.fsaverage5.func.gii']);
pgr = gifti([ProjectFolder '/Gradients.rh.fsaverage5.func.gii']);
% extract unimodal-transmodal gradient
grad_lh = pgl.cdata(:,1);
grad_rh = pgr.cdata(:,1);
% SNR Labels
lh_SNR_mask=read_label([],'~/data/lh.Mask_SNR.label');
rh_SNR_mask=read_label([],'~/data/rh.Mask_SNR.label');
% get group atlases from here
atlasdir='/cbica/projects/abcdfnets/data/cui2020_groupConsensus/';
Loading_Mat=load([atlasdir 'Group_AtlasLoading.mat']);
% convert to hard parcels
% left
V_l=Loading_Mat.sbj_AtlasLoading_lh';
% apply ZC code from step 6th - subject AtlasLabel will have parcels for viz
V_l_Max = max(V_l);
trimInd = V_l ./ max(repmat(V_l_Max, size(V_l, 1), 1), eps) < 5e-2;
V_l(trimInd) = 0;
sbj_AtlasLoading_NoMedialWall_l = V_l;
[~, sbj_AtlasLabel_NoMedialWall_l] = max(sbj_AtlasLoading_NoMedialWall_l, [], 2);
% right
V_r=Loading_Mat.sbj_AtlasLoading_rh';
% apply ZC code from step 6th - subject AtlasLabel will have parcels for viz
V_r_Max = max(V_r);
trimInd = V_r ./ max(repmat(V_r_Max, size(V_r, 1), 1), eps) < 5e-2;
V_r(trimInd) = 0;
sbj_AtlasLoading_NoMedialWall_r = V_r;
[~, sbj_AtlasLabel_NoMedialWall_r] = max(sbj_AtlasLoading_NoMedialWall_r, [], 2);

% remove PG and netlabs at snr verts
sbj_AtlasLabel_NoMedialWall_l(lh_SNR_mask(:,1)+1)=[];     
sbj_AtlasLabel_NoMedialWall_r(rh_SNR_mask(:,1)+1)=[];
grad_lh(lh_SNR_mask(:,1)+1)=[];
grad_rh(rh_SNR_mask(:,1)+1)=[];

% initialize hierarchy vector
HierarchyVec=zeros(17,1);

% for each network
for N=1:17
	leftN=find(sbj_AtlasLabel_NoMedialWall_l==N);
	PGValsL=grad_lh(leftN);	
	rightN=find(sbj_AtlasLabel_NoMedialWall_r==N);
	PGValsR=grad_rh(rightN);
	PGVals=[PGValsL' PGValsR'];
	HierarchyVec(N)=mean(PGVals);
end

% save out hierarchy vec
writetable(table(HierarchyVec),'~/results/hierarchyVec.csv')
