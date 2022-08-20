function CircFC_low(subj)
%%%% Derive Circular Correlation Matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in subj data and freesurfer info
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
% Load in fsav4 opflow calc
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4_VeryLow.mat'];
data=load(OpFlFp)
% Load in surface data
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
% get temporal length
lenOpFl=length(vfl);
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output arrays
AnglesL=zeros(length(g_noMW_combined_L),lenOpFl);
AnglesR=zeros(length(g_noMW_combined_R),lenOpFl);
CFC_L=zeros(5120);
CFC_R=zeros(5120);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FOR ALL LEFT FACES %%%%%%
%Get az/el, derive circ corrs%
% I. Az/el
% note we are loopign over the non-mw-face indices, skipping medial wall faces but leaving 0's in their stead
for F=g_noMW_combined_L
        % translate xyz vector fields from opfl to az/el/r to polar angles
        % note, by not saving rho (just theta), we are discarding magnitude information at this point
        for fr=1:lenOpFl
                % current vector field
                relVf_L=vfl{fr};
                % xyz components
                xComp_L=relVf_L(F,1);
                yComp_L=relVf_L(F,2);
                zComp_L=relVf_L(F,3);
                % convert to spherical coord system
                vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
                % store in output vector (r is redundant across all vecs, only using az and el)
                AnglesL(F,fr)=cart2pol(vs_L(1),vs_L(2));
        end
end	
% II. derive circ corrs
% get triu, loop over triu coordinates
[i,j]=meshgrid(g_noMW_combined_L,g_noMW_combined_L);
UpperCoords_x=i(i<j);
UpperCoords_y=j(i<j);
% loop over upper triangle of adjacency matrix
for F=1:length(UpperCoords_x);
	F;
	% get x and y coord
	Xcoord=UpperCoords_x(F);
	Ycoord=UpperCoords_y(F);
	[CFC_L(Xcoord,Ycoord) pval]=circ_corrcc(AnglesL(Xcoord,:),AnglesL(Ycoord,:));
end	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FOR ALL RIGHT FACES %%%%%
%Get az/el, derive circ corrs%
% I. Az/el
% note we are loopign over the non-mw-face indices, skipping medial wall faces but leaving 0's in their stead
for F=g_noMW_combined_R
        % translate xyz vector fields from opfl to az/el/r to polar angles
        % note, by not saving rho (just theta), we are discarding magnitude information at this point
        for fr=1:lenOpFl
                % current vector field
                relVf_R=vfr{fr};
                % xyz components
                xComp_R=relVf_R(F,1);
                yComp_R=relVf_R(F,2);
                zComp_R=relVf_R(F,3);
                % convert to spherical coord system
                vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(F),eld_R(F));
                % store in output vector (r is redundant across all vecs, only using az and el)
                AnglesR(F,fr)=cart2pol(vs_R(1),vs_R(2));
        end
end
% II. derive circ corrs
% get triu, loop over triu coordinates
[i,j]=meshgrid(g_noMW_combined_R,g_noMW_combined_R);
UpperCoords_x=i(i<j);
UpperCoords_y=j(i<j);
% loop over upper triangle of adjacency matrix
for F=1:length(UpperCoords_x);
        % get x and y coord
        Xcoord=UpperCoords_x(F);
        Ycoord=UpperCoords_y(F);
        [CFC_R(Xcoord,Ycoord) pval]=circ_corrcc(AnglesR(Xcoord,:),AnglesR(Ycoord,:));
end

% mirror
CFC_L=triu(CFC_L)+triu(CFC_L,1)';
CFC_R=triu(CFC_R)+triu(CFC_R,1)';
% mask mw 0s out
CFC_L=CFC_L(g_noMW_combined_L,g_noMW_combined_L);
CFC_R=CFC_R(g_noMW_combined_R,g_noMW_combined_R);
% saveout
adjMats=struct('L',{CFC_L},'R',{CFC_R});
fn=['~/results/PWs/Proced/' subj '/' subj '_CircFC_VeryLow.mat'];
save(fn,'adjMats')
