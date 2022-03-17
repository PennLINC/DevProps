function PGG_AngDistCalc(subj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load in optical flow data and calculate gradient gradient data, compare angular distance between them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add OFD toolbox to path
%addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

%%%% Load in fsav5 opflow calc
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs5_c.mat'];
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

% load in GROUP curvature map: please excuse maintained "PG" variable names
gpg=load('~/data/fs5curv.mat');
gPG_LH=gpg.gpg.gPG_LH;
gPG_RH=gpg.gpg.gPG_RH;

% AND FOR GROUP
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);

% extract face-wise vector cartesian vector components
gPGx_L=gPGg_L(:,1);
gPGy_L=gPGg_L(:,2);
gPGz_L=gPGg_L(:,3);
gPGx_R=gPGg_R(:,1);
gPGy_R=gPGg_R(:,2);
gPGz_R=gPGg_R(:,3);

% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));

% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

% translate xyz vector components at coordinates to az/el/r
gazes_L=zeros(1,length(azd_L));
gels_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    gvs_L=cart2sphvec(double([gPGx_L(i);gPGy_L(i);gPGz_L(i)]),azd_L(i),eld_L(i));
    gazes_L(i)=gvs_L(1);
    gels_L(i)=gvs_L(2);
	% drop the third vector, as each point is equidistant from the center of the sphere
end
% right hemi
gazes_R=zeros(1,length(azd_R));
gels_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    gvs_R=cart2sphvec(double([gPGx_R(i);gPGy_R(i);gPGz_R(i)]),azd_R(i),eld_R(i));
    gazes_R(i)=gvs_R(1);
    gels_R(i)=gvs_R(2);
end

% get length of OpFl pairs
lenOpFl=length(data.us.vf_left);

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

% actually calculate Opflow vector vs PGG vector angular distance
gangDist_L=zeros(lenOpFl,length(azd_L));

% for each vertex
for Vert=1:length(azd_L)
    % note azimuth elevation ordering for atan2d
    gPGvec_L=[gazes_L(Vert) gels_L(Vert)]; 
    % PG GROUP LOAD IN
    for fr=1:lenOpFl
        OpFlVec_L=[azesOpf_L(Vert,fr) elsOpf_L(Vert,fr)];
	% go with this top rec as primary https://stackoverflow.com/questions/40461268/calculate-angle-between-two-vectors-matlab
	% note it returns same value as two alternative methods below
	%dotUV = dot(PGvec_L,OpFlVec_L);
	%normU = norm(PGvec_L);
	%normV = norm(OpFlVec_L);
	%a_alt = acosd(dotUV/(normU * normV));
	%a_alt2 = atan2d(norm(cross([PGvec_L 0],[OpFlVec_L 0])),dot([PGvec_L 0],[OpFlVec_L 0]))
    	ga = acosd(min(1,max(-1, gPGvec_L(:).' *OpFlVec_L(:) / norm(gPGvec_L) / norm(OpFlVec_L) )));
	gangDist_L(fr,Vert) = ga;	
    end
end
% and for group
gangDist_R=zeros(lenOpFl,length(azd_R));

% for each vertex
for Vert=1:length(azd_R)
    % note azimuth elevation ordering for atan2d
    gPGvec_R=[gazes_R(Vert) gels_R(Vert)];
    for fr=1:lenOpFl
        OpFlVec_R=[azesOpf_R(Vert,fr) elsOpf_R(Vert,fr)];
        %dotUV = dot(PGvec_R,OpFlVec_R);
        %normU = norm(PGvec_R);
        %normV = norm(OpFlVec_R);
        %a_alt = acosd(dotUV/(normU * normV));
        %a_alt2 = atan2d(norm(cross([PGvec_R 0],[OpFlVec_R 0])),dot([PGvec_R 0],[OpFlVec_R 0]))
        ga = acosd(min(1,max(-1, gPGvec_R(:).' *OpFlVec_R(:) / norm(gPGvec_R) / norm(OpFlVec_R) )));
        gangDist_R(fr,Vert) = ga;
    end
end

% save it out
AngDist=struct;
AngDist.gLeft=gangDist_L;
AngDist.gRight=gangDist_R;
AngDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_curvAngDistMat_c.mat'];
save(AngDistFP,'AngDist')
