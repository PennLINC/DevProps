% add paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

%%% Load in surface data
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage4';
surfR = [SubjectsFolder '/surf/rh.sphere'];
% surface topography
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_r = faces_r + 1;
% faces_L
F_R=faces_r;
% vertices V
V_R=vx_r;
% get incenters of triangles
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;
%%% Convert PGG to tangent-plane circular coordinates
% First get sphere coordinates in az/el/r
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

%%%% Calculate all PGG angles
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% calculate group PG gradient on sphere
gPGg_R = grad(F_R, V_R, gPG_RH);
% extract face-wise vector cartesian vector components
gPGx_R=gPGg_R(:,1);
gPGy_R=gPGg_R(:,2);
gPGz_R=gPGg_R(:,3);

% use native freesurfer command for mw mask indices
surfMR = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/label/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);
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
	vfr=data.us.vf_right;
	% get length of OpFl pairs
	lenOpFl=length(data.us.vf_right);
	% array for Angles over Faces by TRs
	AngArr=zeros(5120,lenOpFl);
	% INITIALIZE SEP ARRAY FOR PGG DISTANCE SPECIFICALLY
	AngDistArr=zeros(5120,1);
	% convert to x y (az/el)
	for i=1:length(azd_R)
	    % get PGG angle for this face
	    gvs_R=cart2sphvec(double([gPGx_R(i);gPGy_R(i);gPGz_R(i)]),azd_R(i),eld_R(i));
	    gazes_R(i)=gvs_R(1);
	    gels_R(i)=gvs_R(2);
	    PGGAng=angle(gvs_R(1)+j*gvs_R(2));
	    % get angles from optical flow vector fields for this face
	    for fr=1:lenOpFl
                % current vector field
                relVf_R=vfr{fr};
                % xyz components
                xComp_R=relVf_R(i,1);
                yComp_R=relVf_R(i,2);
                zComp_R=relVf_R(i,3);
                % convert to spherical coord system
                vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(i),eld_R(i));
	        % get angle relative to pos. x-axis
	        AngArr(i,fr)=angle(vs_R(1)+j*vs_R(2));
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
	AngArrFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj 'AngDist90_R.csv'];
	writetable(table(AngDistArr(g_noMW_combined_R)),AngArrFp);
	s
end
