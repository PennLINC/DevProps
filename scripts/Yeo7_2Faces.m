% add paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
% load in yeo7L
left=gifti('/cbica/projects/pinesParcels/data/y7_L_3k.label.gii');
%R
right=gifti('/cbica/projects/pinesParcels/data/y7_R_3k.label.gii');
%%% generate MW mask
% load surface data
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
surfL = [SubjectsFolder '/surf/lh.sphere'];
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
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
% initialize num faces x 7
membL=zeros(5120,7);
membR=zeros(5120,7);

% left hemisphere
for y=1:7
	% get numeric boolean of which vertices belong to this network
	bool_L_yInds=double(left.cdata==y);
	%L - for each face: interp
	% convert this network to faces, average value across three vertices defining face
	for f=g_noMW_combined_L;
		vert1 = F_L(f,1);
		vert2 = F_L(f,2);
		vert3 = F_L(f,3);
		BoolAv = mean([bool_L_yInds(vert1) bool_L_yInds(vert2) bool_L_yInds(vert3)]);
		% prop over 50?
		if BoolAv > .5
			membL(f,y)=1;
		else
		end
	end
end
			
%R
for y=1:7
        % get numeric boolean of which vertices belong to this network
        bool_R_yInds=double(right.cdata==y);
        %L - for each face: interp
        % convert this network to faces, average value across three vertices defining face
        for f=g_noMW_combined_R;
                vert1 = F_R(f,1);
                vert2 = F_R(f,2);
                vert3 = F_R(f,3);
                BoolAv = mean([bool_R_yInds(vert1) bool_R_yInds(vert2) bool_R_yInds(vert3)]);
                % prop over 50?
                if BoolAv > .5
                        membR(f,y)=1;
                else
                end
        end
end

% add last column to catch faces not assigned a network
membL(:,8)=.5;
membR(:,8)=.5;

% convert to facewise membership vector
[~, yeoLabel_L] = max(membL,[],2);
[~, yeoLabel_R] = max(membR,[],2);

% vis it
Vis_FaceVecY7(yeoLabel_L(g_noMW_combined_L),yeoLabel_R(g_noMW_combined_R),'yeo7.png')


% saveout as table
writetable(table(membL(g_noMW_combined_L,:)),'~/data/yeo7FaceBooleans_L.csv');
writetable(table(membR(g_noMW_combined_R,:)),'~/data/yeo7FaceBooleans_R.csv'); 
