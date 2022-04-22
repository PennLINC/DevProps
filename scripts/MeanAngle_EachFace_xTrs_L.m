% add paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

%%% Load in surface data
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
surfL = [SubjectsFolder '/surf/lh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
%%% Convert PGG to tangent-plane circular coordinates
% First get sphere coordinates in az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);

%%%% Calculate all PGG angles
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
% extract face-wise vector cartesian vector components
gPGx_L=gPGg_L(:,1);
gPGy_L=gPGg_L(:,2);
gPGz_L=gPGg_L(:,3);

% use native freesurfer command for mw mask indices
surfML = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
% face mask indices
fmwIndVec_l=find(F_MW_L);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);

% load in subjs list
Subjs=readtable('~/PWs/rs_subs.csv');
% initialize mean Dirs by face array for populating
FacemeanDirs=zeros(height(Subjs),5120);
% for each subj, read left opFlow
for s=1:height(Subjs)
	tic
	subj=table2array(Subjs{s,2});
	OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'];
	data=load(OpFlFp);
	% vector fields
	vfl=data.us.vf_left;
	% get length of OpFl pairs
	lenOpFl=length(data.us.vf_left);
	% array for Angles over Faces by TRs
	AngArr=zeros(5120,lenOpFl);
	% INITIALIZE SEP ARRAY FOR PGG DISTANCE SPECIFICALLY
	AngDistArr=zeros(5120,1);
	% convert to x y (az/el)
	for i=1:length(azd_L)
	    % get PGG angle for this face
	    gvs_L=cart2sphvec(double([gPGx_L(i);gPGy_L(i);gPGz_L(i)]),azd_L(i),eld_L(i));
	    gazes_L(i)=gvs_L(1);
	    gels_L(i)=gvs_L(2);
	    PGGAng=angle(gvs_L(1)+j*gvs_L(2));
	    % get angles from optical flow vector fields for this face
	    for fr=1:lenOpFl
                % current vector field
                relVf_L=vfl{fr};
                % xyz components
                xComp_L=relVf_L(i,1);
                yComp_L=relVf_L(i,2);
                zComp_L=relVf_L(i,3);
                % convert to spherical coord system
                vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(i),eld_L(i));
	        % get angle relative to pos. x-axis
	        AngArr(i,fr)=angle(vs_L(1)+j*vs_L(2));
            end
	    % get mean over TRs for each face
	    FacemeanDirs(s,i)=circ_mean(AngArr(i,:));
	    % GET IT RELATIVE TO PGG
	    AngDistArr(i)=circ_mean(AngArr(i,:))- PGGAng;
	end
	toc	
	% CONVERT SEP. ARRAY TO 0-(1/2 pi) RANGE: see Convert_MeanMeans_to90.m for more detail
        % REFLECTION 1
	AngDistArr=abs(AngDistArr);
	% 2
	AngDistArr=AngDistArr-pi;
	AngDistArr=abs(AngDistArr);
	AngDistArr=AngDistArr*-1;
	AngDistArr=AngDistArr+pi;
	% save out angle array
	AngArrFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj 'AngDist90_L.csv'];
	writetable(table(AngDistArr(g_noMW_combined_L)),AngArrFp);
	s
end
