%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Save out the distribution of vector angles relative to arb. x=1 y=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% add subject name in manually
%subj='';
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'];
data=load(OpFlFp)
% vector fields
vfl=data.us.vf_left;
vfr=data.us.vf_right;

%%% Load in surface data
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
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert Op Flow angles to tangent-plane circular coordinates
% First get sphere coordinates in az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);
% get length of OpFl pairs
lenOpFl=length(data.us.vf_left);
% translate xyz vector fields from opfl to az/el/r
rel2x_L=zeros(lenOpFl,length(eld_L));
rel2x_R=zeros(lenOpFl,length(eld_R));
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
        % get angle relative to pos. x-axis
        rel2x_L(fr,i)=angle(vs_L(1)+j*vs_L(2));
    end
end
% right hemisphre
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
        rel2x_R(fr,i)=angle(vs_R(1)+j*vs_R(2));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% apply mw mask
rel2x_L=rel2x_L(:,g_noMW_combined_L);
rel2x_R=rel2x_R(:,g_noMW_combined_R);

% circular mean over TRs
rel2x_L_me=zeros(1,length(g_noMW_combined_L));
rel2x_R_me=zeros(1,length(g_noMW_combined_R));
for i=1:length(g_noMW_combined_L)
    rel2x_L_me(i)=circ_mean(rel2x_L(:,i));
end
% right hemisphre
for i=1:length(g_noMW_combined_R)
    rel2x_R_me(i)=circ_mean(rel2x_R(:,i));
end

% save for python
writetable(table(rel2x_L_me),'~/data/rel2xL_exSubj2.csv')
writetable(table(rel2x_R_me),'~/data/rel2xR_exSubj2.csv')
