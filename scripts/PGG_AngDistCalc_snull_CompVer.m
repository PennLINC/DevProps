function PGG_AngDistCalc_snull(subj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load in optical flow data and calculate gradient gradient data, compare angular distance between them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add OFD toolbox to path
%addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

%%%% Load in fsav5 opflow calc
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'];
data=load(OpFlFp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load in surface data %%%%%
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4';
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
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the same with the smoothed PG: but dont use the mask from it
% load in GROUP PG
gpg=load('~/data/gpg_fs4.mat');
gPG_LH=gpg.gpg.gPG_LH;
gPG_RH=gpg.gpg.gPG_RH;
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);
% cartesian components
gPGx_L=gPGg_L(:,1);
gPGy_L=gPGg_L(:,2);
gPGz_L=gPGg_L(:,3);
gPGx_R=gPGg_R(:,1);
gPGy_R=gPGg_R(:,2);
gPGz_R=gPGg_R(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% true gazes and gels for comparison dip
tgazes_L=zeros(1,length(azd_L));
tgels_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    gvs_L=cart2sphvec(double([gPGx_L(i);gPGy_L(i);gPGz_L(i)]),azd_L(i),eld_L(i));
    tgazes_L(i)=gvs_L(1);
    tgels_L(i)=gvs_L(2);
	% drop the third vector, as each point is equidistant from the center of the sphere
end
% right hemi
tgazes_R=zeros(1,length(azd_R));
tgels_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    gvs_R=cart2sphvec(double([gPGx_R(i);gPGy_R(i);gPGz_R(i)]),azd_R(i),eld_R(i));
    tgazes_R(i)=gvs_R(1);
    tgels_R(i)=gvs_R(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load in spun PGG vectors
sp_PGGs=load('/cbica/projects/pinesParcels/results/aggregated_data/PGGPermuts_fs4.mat');
gazes_L_all=squeeze(sp_PGGs.spGPGGs(1,:,1:5120));
gels_L_all=squeeze(sp_PGGs.spGPGGs(2,:,1:5120));
gazes_R_all=squeeze(sp_PGGs.spGPGGs(1,:,5121:10240));
gels_R_all=squeeze(sp_PGGs.spGPGGs(2,:,5121:10240));
%load in spun PGG scalar data, for masking (note difference between PG permutations and PGG permutations)
sp_fp='/cbica/projects/pinesParcels/results/aggregated_data/PGPermuts_fs4.mat';
sp=load(sp_fp);
% initialize facewise output structs to bear mask indices
sp_mw_L=struct;
sp_mw_R=struct;
% unfortunately we need to convert these vertices to faces for equiv masking. Easiest/most equiv way to do this below
for i=1:1000
	sp_gPG_LH=sp.bigrotl(i,:);
	sp_gPG_RH=sp.bigrotr(i,:);
	%%% find where verts are zero or 100: indicates medial wall location and displaced MW, respectively
	% use found verts to denote invalid faces: those that draw from one of the mw or displaced mw verts
	mwInds_L=find(sp_gPG_LH==100);
	mwInds_R=find(sp_gPG_RH==100);
	displ_mwInds_L=find(sp_gPG_LH==0);
	displ_mwInds_R=find(sp_gPG_RH==0);
	% combine them
	mw_L=union(mwInds_L,displ_mwInds_L);
	mw_R=union(mwInds_R,displ_mwInds_R);	
	% all faces that touch a medial wall vertex to be masked
	MW_f1_L=find(ismember(F_L(:,1),mw_L));
	MW_f2_L=find(ismember(F_L(:,2),mw_L));
	MW_f3_L=find(ismember(F_L(:,3),mw_L));
	% rh
	MW_f1_R=find(ismember(F_R(:,1),mw_R));
	MW_f2_R=find(ismember(F_R(:,2),mw_R));
	MW_f3_R=find(ismember(F_R(:,3),mw_R));
	% combine... again...
	MW_combined_L=union(MW_f1_L,MW_f2_L);
	MW_combined_L=union(MW_combined_L,MW_f3_L);
	% now for right hemisphere
	MW_combined_R=union(MW_f1_R,MW_f2_R);
	MW_combined_R=union(MW_combined_R,MW_f3_R);	
	sp_mw_L.inds{i}=MW_combined_L;
	sp_mw_R.inds{i}=MW_combined_R;	
end

% get length of OpFl pairs
lenOpFl=length(data.us.vf_left);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% translate xyz vector fields from opfl to az/el/r %%%%%%%%%%
azesOpf_L=zeros(length(azd_L),lenOpFl);
elsOpf_L=zeros(length(azd_L),lenOpFl);
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
disp('done converting opfl vectors from cartesian')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate Opflow vector vs PGG vector angular distance, diptest %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output vector for dip test results
DTres=zeros(1,1000);
% for each spin
for S=1:1000
	tic
	S
	% initialize ang dist over all vectors vector
	gangDistL=zeros(lenOpFl,length(azd_L));
	gangDistR=zeros(lenOpFl,length(azd_R));
	% this iterations gazels and gel
	gazes_L=gazes_L_all(S,:);
	gels_L=gels_L_all(S,:);
	gazes_R=gazes_R_all(S,:);
	gels_R=gels_R_all(S,:);
	% for each vertex
	for Vert=1:length(azd_L)
	    % note azimuth elevation ordering for atan2d
	    gPGvec_L=[gazes_L(Vert) gels_L(Vert)]; 
	    % PG GROUP LOAD IN
	    for fr=1:lenOpFl
	        OpFlVec_L=[azesOpf_L(Vert,fr) elsOpf_L(Vert,fr)];
    		ga = acosd(min(1,max(-1, gPGvec_L(:).' *OpFlVec_L(:) / norm(gPGvec_L) / norm(OpFlVec_L) )));
		gangDist_L(fr,Vert) = ga;	
    	     end
	end
	% right hemi
	% for each vertex
	for Vert=1:length(azd_R)
	    % note azimuth elevation ordering for atan2d
	    gPGvec_R=[gazes_R(Vert) gels_R(Vert)];
	    for fr=1:lenOpFl
	        OpFlVec_R=[azesOpf_R(Vert,fr) elsOpf_R(Vert,fr)];
        	ga = acosd(min(1,max(-1, gPGvec_R(:).' *OpFlVec_R(:) / norm(gPGvec_R) / norm(OpFlVec_R) )));
        	gangDist_R(fr,Vert) = ga;
    	    end
	end
	toc
	% get index of where they are 0 in all directions
	sp_gPGg_L0=find(all(sp_gPG_LH')==0);
	sp_gPGg_R0=find(all(sp_gPG_RH')==0);
	% get union of two medial wall masks
	MasterMaskL=union(sp_gPGg_L0,sp_mw_L.inds{S});
	MasterMaskR=union(sp_gPGg_R0,sp_mw_R.inds{S});
	% include OG MW mask to remove OpFl vectors there
	MasterMaskL=union(fmwIndVec_l,MasterMaskL);
        MasterMaskR=union(fmwIndVec_r,MasterMaskR);
	% extract data outside of these masks	
	OutOfMaskL=setdiff([1:5120],MasterMaskL);
	OutOfMaskR=setdiff([1:5120],MasterMaskR);
	% merge ang distances
	sp_AngDists=horzcat(gangDist_L(:,OutOfMaskL),gangDist_R(:,OutOfMaskR));
	% dip test
	[dip,xl,xu,ifault,gcm,lcm,mn,mj]=HartigansDipTest(sp_AngDists);
	DTres(S)=dip;
	% save out some example distributions
	%if (S<10)
	%	SpunAngDist_exFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_SpunDistr' num2str(S) '.csv'];
	%	writetable(table(sp_AngDists),SpunAngDist_exFP)	
	%end
end

% run once on actual PGG for comparison
True_gangDist_L=zeros(lenOpFl,length(azd_L));
True_gangDist_R=zeros(lenOpFl,length(azd_R));
% for each vertex
for Vert=1:length(azd_L)
    % note azimuth elevation ordering for atan2d
    gPGvec_L=[tgazes_L(Vert) tgels_L(Vert)]; 
    % PG GROUP LOAD IN
    for fr=1:lenOpFl
        OpFlVec_L=[azesOpf_L(Vert,fr) elsOpf_L(Vert,fr)];
    	ga = acosd(min(1,max(-1, gPGvec_L(:).' *OpFlVec_L(:) / norm(gPGvec_L) / norm(OpFlVec_L) )));
	True_gangDist_L(fr,Vert) = ga;	
    end
end
for Vert=1:length(azd_R)
    % note azimuth elevation ordering for atan2d
    gPGvec_R=[tgazes_R(Vert) tgels_R(Vert)];
    % PG GROUP LOAD IN
    for fr=1:lenOpFl
        OpFlVec_R=[azesOpf_R(Vert,fr) elsOpf_R(Vert,fr)];
        ga = acosd(min(1,max(-1, gPGvec_R(:).' *OpFlVec_R(:) / norm(gPGvec_R) / norm(OpFlVec_R) )));
        True_gangDist_R(fr,Vert) = ga;
    end
end

% make the last one "true" observed from masked data
Tr_AngDists=horzcat(True_gangDist_L(g_noMW_combined_L),True_gangDist_R(g_noMW_combined_R));
[dip,xl,xu,ifault,gcm,lcm,mn,mj]=HartigansDipTest(Tr_AngDists);
DTres(1001)=dip;

% save it out
SpunAngDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_SpunDips4.csv'];
writetable(table(DTres),SpunAngDistFP)
