function PGG_AngDistCalc(subj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load in optical flow data and calculate gradient gradient data, compare angular distance between them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

%%%% Load in fsav5 opflow calc
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs5.mat'];
data=load(OpFlFp)

%%% Load in surface data
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/7.2.0/subjects/fsaverage5';
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
% vector fields
vfl=data.us.vf_left;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
% vector fields
vfr=data.us.vf_right;
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;

% load in subject's PG
%%%%% DOOO IT

% calculate PG gradient on sphere
PGg_L = grad(F_L, V_L, PG_LH);
PGg_R = grad(F_R, V_R, PG_RH);

% extract face-wise vector cartesian vector components
PGx_L=PGg_L(:,1);
PGy_L=PGg_L(:,2);
PGz_L=PGg_L(:,3);
PGx_R=PGg_R(:,1);
PGy_R=PGg_R(:,2);
PGz_R=PGg_R(:,3);

% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));

% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

% translate xyz vector components at coordinates to az/el/r
azes_L=zeros(1,length(azd_L));
els_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    vs_L=cart2sphvec(double([PGx_L(i);PGy_L(i);PGz_L(i)]),azd_L(i),eld_L(i));
    azes_L(i)=vs_L(1);
    els_L(i)=vs_L(2);
    rvec=vs_L(3);
end
% right hemi
azes_R=zeros(1,length(azd_R));
els_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    vs_R=cart2sphvec(double([PGx_R(i);PGy_R(i);PGz_R(i)]),azd_R(i),eld_R(i));
    azes_R(i)=vs_R(1);
    els_R(i)=vs_R(2);
    rvec=vs_R(3);
end

% get length of OpFl pairs
lenOpFl=

% translate xyz vector fields from opfl to az/el/r
azesOpf_L=zeros(length(azd_L),lenOpFl);
elsOpf_L=zeros(length(azd_L),lenOpFl);
thetasOpf_L=zeros(length(azd_L),lenOpFl);
for i=1:length(azd_L)
    for fr=1:lenOpFl
	% current vector field
        relVf_L=vfl{fr};
	% xyz components
        xComp_L=relVf_L(i,1);
        yComp_L=relVf_L(i,2);
        zComp_L=relVf_L(i,3);
	% convert to spherical coord system
        vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(i),eld_L(i));
        % store in output vector (r is redundant across all vecs)
	azesOpf_L(i,fr)=vs_L(1);
        elsOpf_L(i,fr)=vs_L(2);
        rvec=vs_L(3);
    end
end
% right hemisphre
azesOpf_R=zeros(length(azd_R),lenOpFl);
elsOpf_R=zeros(length(azd_R),lenOpFl);
thetasOpf_R=zeros(length(azd_R),lenOpFl);
for i=1:length(azd_R)
    for fr=1:lenOpFl
	% current vector field
        relVf_R=vfr{fr};
	% xyz components
        xComp_R=relVf_R(i,1);
        yComp_R=relVf_R(i,2);
        zComp_R=relVf_R(i,3);
	% convert to spherical coord system
        vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(i),eld_R(i));
        % store in output vector (r is redundant across all vecs)
	azesOpf_R(i,fr)=vs_R(1);
        elsOpf_R(i,fr)=vs_R(2);
        rvec=vs_R(3);
    end
end

% convert pseudo y x to single angle (theta), az as x, el as y
PGg_thetas_L=atan2(els_L,azes_L);
% initialize
Opfl_thetas_L=zeros(length(azd_L),lenOpFl);
% loop over each vertex
for i=1:length(azd_L)
    for fr=1:lenOpFl
        Opfl_thetas_L(i,fr)=atan2(elsOpf_L(i,fr),azesOpf_L(i,fr));
    end
end
% right hemi
PGg_thetas_R=atan2(els_R,azes_R);
% initialize
Opfl_thetas_R=zeros(length(azd_R),lenOpFl);
% loop over each vertex
for i=1:length(azd_R)
    for fr=1:lenOpFl
        Opfl_thetas_R(i,fr)=atan2(elsOpf_R(i,fr),azesOpf_R(i,fr));
    end
end

% actually calculate Opflow vector vs PGG vector angular distance
angDist_L=zeros(lenOpFl,length(azd_L));
% for each vertex
for Vert=1:length(azd_L)
    % note azimuth elevation ordering for atan2d
    PGvec_L=[azes_L(Vert) els_L(Vert)];
    for fr=1:lenOpFl
        OpFlVec_L=[azesOpf_L(Vert,fr) elsOpf_L(Vert,fr)];
	%%%%%% TRIP CHECK THIS CALC, FIND LINK
        a = atan2d(PGvec_L(1)*OpFlVec_L(2)-PGvec_L(2)*OpFlVec_L(1),PGvec_L(1)*PGvec_L(2)+OpFlVec_L(1)*OpFlVec_L(2));
        angDist_L(fr,Vert) = a;
    end
end
% right hemi
angDist_R=zeros(lenOpFl,length(azd_R));
% for each vertex
for Vert=1:length(azd_R)
    % note azimuth elevation ordering for atan2d
    PGvec_L=[azes_R(Vert) els_R(Vert)];
    for fr=1:lenOpFl
        OpFlVec_R=[azesOpf_R(Vert,fr) elsOpf_R(Vert,fr)];
        a = atan2d(PGvec_R(1)*OpFlVec_R(2)-PGvec_R(2)*OpFlVec_R(1),PGvec_R(1)*PGvec_R(2)+OpFlVec_R(1)*OpFlVec_R(2));
        angDist_R(fr,Vert) = a;
    end
end

% save it out
AngDist=struct;
AngDist.Left=angDist_L;
AngDist.Rifht=angDist_R;
AngDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistMat.mat'];
save(AngDistFP,'AngDist')
