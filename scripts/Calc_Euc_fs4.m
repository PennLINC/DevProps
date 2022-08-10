% load in fs4 surface, saved out euclidean distance matrices for both hemispheres (pial-based)

%add needed paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% add surfaces (pial for equiv distances)
[vx_l, faces_l]=read_surf('/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4/surf/lh.pial');
[vx_r, faces_r]=read_surf('/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4/surf/rh.pial');
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;

% get center-of-face xyz coordinates
TR = TriRep(faces_l, vx_l);
cof_l = TR.incenters;
TR = TriRep(faces_r, vx_r);
cof_r = TR.incenters;

% use native freesurfer command for mw mask indices
surfML = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);

% apply to faces
faces_l=cof_l(g_noMW_combined_L,:);
faces_r=cof_r(g_noMW_combined_R,:);

% big diStance matrix
bdsml=zeros(length(faces_l),length(faces_l));
bdsmr=zeros(length(faces_r),length(faces_r));

% for all left faces
for F=1:length(bdsml);
		% vertex props (x,y,z coords)
		xVL=faces_l(F,1);
		yVL=faces_l(F,2);
		zVL=faces_l(F,3);
		% search through all of them for eucl. dist. calc.
		for i=1:length(bdsml);
			xL=faces_l(i,1);
			yL=faces_l(i,2);
			zL=faces_l(i,3);
			% calculate distance in 3d
			eucld_L=sqrt((xL-xVL)^2+(yL-yVL)^2+(zL-zVL)^2);
			bdsml(F,i)=eucld_L;
		end	
end

% for all right faces
for F=1:length(bdsmr);
                % vertex props (x,y,z coords)
                xVR=faces_r(F,1);
                yVR=faces_r(F,2);
                zVR=faces_r(F,3);
                % search through all of them for eucl. dist. calc.
                for i=1:length(bdsmr);
                        xR=faces_r(i,1);
                        yR=faces_r(i,2);
                        zR=faces_r(i,3);
                        % calculate distance in 3d
                        eucld_R=sqrt((xR-xVR)^2+(yR-yVR)^2+(zR-zVR)^2);
                        bdsmr(F,i)=eucld_R;
                end
end


save('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_left_fsaverage4_faces.mat','bdsml');
save('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_right_fsaverage4_faces.mat','bdsmr');

