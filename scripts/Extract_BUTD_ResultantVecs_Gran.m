function Extract_BUTD_ResultantVecs(OpFlFp,outFP_L,outFP_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each face. Measured relative to gPGG.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temp to evaluate compute time
tic

% load in gPPG ang dist: for thresholdingdd OFD toolbox to path
% addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% Load in fsav4 opflow calc
%OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4_run_' runNum '.mat'];
data=load(OpFlFp)
% Load in surface data
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
surfR = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use native freesurfer command for mw mask indices
surfML = '/oak/stanford/groups/leanew1/users/apines/surf/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/oak/stanford/groups/leanew1/users/apines/surf/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;

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
Lvertex_nomw=setdiff([1:2562],mwIndVec_l);
Rvertex_nomw=setdiff([1:2562],mwIndVec_l);

% load in GPG
gpg=load('/oak/stanford/groups/leanew1/users/apines/maps/gpg_fs4.mat');
gPG_LH=gpg.gpg.gPG_LH;
gPG_RH=gpg.gpg.gPG_RH;
% pg2
gpg2=load('/oak/stanford/groups/leanew1/users/apines/maps/gpg2_fs4.mat');
gPG2_LH=gpg2.gpg.gPG_LH;
gPG2_RH=gpg2.gpg.gPG_RH;
% pg3
gpg3=load('/oak/stanford/groups/leanew1/users/apines/maps/gpg3_fs4.mat');
gPG3_LH=gpg3.gpg.gPG_LH;
gPG3_RH=gpg3.gpg.gPG_RH;
% beta
gb=load('/oak/stanford/groups/leanew1/users/apines/maps/gB_fs4.mat');
B_LH=gb.gpg.gPG_LH;
B_RH=gb.gpg.gPG_RH;
% gamma 1
gg1=load('/oak/stanford/groups/leanew1/users/apines/maps/g1_fs4.mat');
g1_LH=gg1.gpg.gPG_LH;
g1_RH=gg1.gpg.gPG_RH;
% gamma 2
gg2=load('/oak/stanford/groups/leanew1/users/apines/maps/g2_fs4.mat');
g2_LH=gg2.gpg.gPG_LH;
g2_RH=gg2.gpg.gPG_RH;

% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);
gPG2g_L = grad(F_L, V_L, gPG2_LH);
gPG2g_R = grad(F_R, V_R, gPG2_RH);
gPG3g_L = grad(F_L, V_L, gPG3_LH);
gPG3g_R = grad(F_R, V_R, gPG3_RH);
BgL = grad(F_L, V_L, B_LH);
BgR = grad(F_R, V_R, B_RH);
G1gL = grad(F_L, V_L, g1_LH);
G1gR = grad(F_R, V_R, g1_RH);
G2gL = grad(F_L, V_L, g2_LH);
G2gR = grad(F_R, V_R, g2_RH);

% extract face-wise vector cartesian vector components
gPGx_L=gPGg_L(:,1);
gPGy_L=gPGg_L(:,2);
gPGz_L=gPGg_L(:,3);
gPGx_R=gPGg_R(:,1);
gPGy_R=gPGg_R(:,2);
gPGz_R=gPGg_R(:,3);
% pg2
gPG2x_L=gPG2g_L(:,1);
gPG2y_L=gPG2g_L(:,2);
gPG2z_L=gPG2g_L(:,3);
gPG2x_R=gPG2g_R(:,1);
gPG2y_R=gPG2g_R(:,2);
gPG2z_R=gPG2g_R(:,3);
% pg3
gPG3x_L=gPG3g_L(:,1);
gPG3y_L=gPG3g_L(:,2);
gPG3z_L=gPG3g_L(:,3);
gPG3x_R=gPG3g_R(:,1);
gPG3y_R=gPG3g_R(:,2);
gPG3z_R=gPG3g_R(:,3);
%%%% grad of freqs
% B
gBx_L=BgL(:,1);
gBy_L=BgL(:,2);
gBz_L=BgL(:,3);
gBx_R=BgR(:,1);
gBy_R=BgR(:,2);
gBz_R=BgR(:,3);
% G1
gG1x_L=G1gL(:,1);
gG1y_L=G1gL(:,2);
gG1z_L=G1gL(:,3);
gG1x_R=G1gR(:,1);
gG1y_R=G1gR(:,2);
gG1z_R=G1gR(:,3);
% G2
gG2x_L=G2gL(:,1);
gG2y_L=G2gL(:,2);
gG2z_L=G2gL(:,3);
gG2x_R=G2gR(:,1);
gG2y_R=G2gR(:,2);
gG2z_R=G2gR(:,3);

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
end
% right hemi
gazes_R=zeros(1,length(azd_R));
gels_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    gvs_R=cart2sphvec(double([gPGx_R(i);gPGy_R(i);gPGz_R(i)]),azd_R(i),eld_R(i));
    gazes_R(i)=gvs_R(1);
    gels_R(i)=gvs_R(2);
end
% PG2 translate xyz vector components at coordinates to az/el/r
pg2azes_L=zeros(1,length(azd_L));
pg2els_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    pg2vs_L=cart2sphvec(double([gPG2x_L(i);gPG2y_L(i);gPG2z_L(i)]),azd_L(i),eld_L(i));
    pg2azes_L(i)=pg2vs_L(1);
    pg2els_L(i)=pg2vs_L(2);
end
% right hemi
pg2azes_R=zeros(1,length(azd_R));
pg2els_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    pg2vs_R=cart2sphvec(double([gPG2x_R(i);gPG2y_R(i);gPG2z_R(i)]),azd_R(i),eld_R(i));
    pg2azes_R(i)=pg2vs_R(1);
    pg2els_R(i)=pg2vs_R(2);
end
% translate xyz vector components at coordinates to az/el/r
pg3azes_L=zeros(1,length(azd_L));
pg3els_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    pg3vs_L=cart2sphvec(double([gPG3x_L(i);gPG3y_L(i);gPG3z_L(i)]),azd_L(i),eld_L(i));
    pg3azes_L(i)=pg3vs_L(1);
    pg3els_L(i)=pg3vs_L(2);
end
% right hemi
pg3azes_R=zeros(1,length(azd_R));
pg3els_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    pg3vs_R=cart2sphvec(double([gPG3x_R(i);gPG3y_R(i);gPG3z_R(i)]),azd_R(i),eld_R(i));
    pg3azes_R(i)=pg3vs_R(1);
    pg3els_R(i)=pg3vs_R(2);
end


%%%% populate azes and els of freq maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BETA
% translate xyz vector components at coordinates to az/el/r
bazes_L=zeros(1,length(azd_L));
bels_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    bvs_L=cart2sphvec(double([gBx_L(i);gBy_L(i);gBz_L(i)]),azd_L(i),eld_L(i));
    bazes_L(i)=bvs_L(1);
    bels_L(i)=bvs_L(2);
end
% right hemi
bazes_R=zeros(1,length(azd_R));
bels_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    bvs_R=cart2sphvec(double([gBx_R(i);gBy_R(i);gBz_R(i)]),azd_R(i),eld_R(i));
    bazes_R(i)=bvs_R(1);
    bels_R(i)=bvs_R(2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%5
% GAMMA 1
% translate xyz vector components at coordinates to az/el/r
g1azes_L=zeros(1,length(azd_L));
g1els_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    g1vs_L=cart2sphvec(double([gG1x_L(i);gG1y_L(i);gG1z_L(i)]),azd_L(i),eld_L(i));
    g1azes_L(i)=g1vs_L(1);
    g1els_L(i)=g1vs_L(2);
end
% right hemi
g1azes_R=zeros(1,length(azd_R));
g1els_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    g1vs_R=cart2sphvec(double([gG1x_R(i);gG1y_R(i);gG1z_R(i)]),azd_R(i),eld_R(i));
    g1azes_R(i)=g1vs_R(1);
    g1els_R(i)=g1vs_R(2);
end


%%%%%%%%%%%%%%%%%%%%%5
% GAMMA 2
% translate xyz vector components at coordinates to az/el/r
g2azes_L=zeros(1,length(azd_L));
g2els_L=zeros(1,length(eld_L));
for i=1:length(azd_L)
    g2vs_L=cart2sphvec(double([gG2x_L(i);gG2y_L(i);gG2z_L(i)]),azd_L(i),eld_L(i));
    g2azes_L(i)=g2vs_L(1);
    g2els_L(i)=g2vs_L(2);
end
% right hemi
g2azes_R=zeros(1,length(azd_R));
g2els_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    g2vs_R=cart2sphvec(double([gG2x_R(i);gG2y_R(i);gG2z_R(i)]),azd_R(i),eld_R(i));
    g2azes_R(i)=g2vs_R(1);
    g2els_R(i)=g2vs_R(2);
end

%
%%%
%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output (bu_az,bu_el,td_az,td_el,propBu for each face)
% and last 2 columns are BU vec length, TD vec length
%%% TS COLUMNS: MAGNITUDE, THETA, AngD from PG, B, G1, G2, PG2, PG3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%%%
%

% initialize output scalars
OutDf_L=cell(length(F_L),7);
OutDf_R=cell(length(F_R),7);

% count Number of TRs once rather than iteratively
NumTRs=size(vfl);
NumTRs=NumTRs(2);
lenOpFl=NumTRs;

% time series
OutTs_L=cell(length(F_L),8,lenOpFl);
OutTs_R=cell(length(F_R),8,lenOpFl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Starting for the left hemisphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note we are loopign over the non-mw-face indices, skipping medial wall faces but leaving 0's in their stead
for F=g_noMW_combined_L
	% get angles into needed format
	% note azimuth elevation ordering for atan2d
    	gPGvec_L=[gazes_L(F) gels_L(F)];
	gPG2vec_L=[pg2azes_L(F) pg2els_L(F)];
	gPG3vec_L=[pg3azes_L(F) pg3els_L(F)];
    	Bvec_L=[bazes_L(F) bels_L(F)];
    	G1vec_L=[g1azes_L(F) g1els_L(F)];
    	G2vec_L=[g2azes_L(F) g2els_L(F)];
	% translate xyz vector fields from opfl to az/el/r to polar angles
	% note, by not saving rho (just theta), we are discarding magnitude information at this point
	Thetas_L=zeros(1,lenOpFl);
	Mags_L=zeros(1,lenOpFl);
	gangDist_L=zeros(1,lenOpFl);
	BangDist_L=zeros(1,lenOpFl);
	G1angDist_L=zeros(1,lenOpFl);
	G2angDist_L=zeros(1,lenOpFl);
	PG1angDist_L=zeros(1,lenOpFl);
        PG2angDist_L=zeros(1,lenOpFl);
	for fr=1:lenOpFl
		% current vector field
	        relVf_L=vfl{fr};
		% xyz components
	        xComp_L=relVf_L(F,1);
	        yComp_L=relVf_L(F,2);
	        zComp_L=relVf_L(F,3);
		% convert to spherical coord system
	        vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
		% convert to spherical coordinates
        	OpFlVec_L= [vs_L(1) vs_L(2)];
	        % store in output vector (r is redundant across all vecs, only using az and el)
		[Thetas_L(fr),Mags_L(fr)]=cart2pol(vs_L(1),vs_L(2));
		% store rho in magnitude series
		% store az and el
		%OutTs_L(F,3,fr) = {vs_L(1)};
		%OutTs_L(F,4,fr) = {vs_L(2)};	
		%%%% angular distance calculations
		% PG ang d
		OutTs_L(F,3,fr)={acosd(min(1,max(-1, gPGvec_L(:).' *OpFlVec_L(:) / norm(gPGvec_L) / norm(OpFlVec_L) )))};
        	% beta g ang d
		OutTs_L(F,4,fr) ={acosd(min(1,max(-1, Bvec_L(:).' *OpFlVec_L(:) / norm(Bvec_L) / norm(OpFlVec_L) )))};
        	% g1 ang d
		OutTs_L(F,5,fr) ={acosd(min(1,max(-1, G1vec_L(:).' *OpFlVec_L(:) / norm(G1vec_L) / norm(OpFlVec_L) )))};
        	% g2 ang d
		OutTs_L(F,6,fr) ={acosd(min(1,max(-1, G2vec_L(:).' *OpFlVec_L(:) / norm(G2vec_L) / norm(OpFlVec_L) )))};
		% PG2 ang d
                OutTs_L(F,7,fr) ={acosd(min(1,max(-1, gPG2vec_L(:).' *OpFlVec_L(:) / norm(gPG2vec_L) / norm(OpFlVec_L) )))};
                % PG3 ang d
                OutTs_L(F,8,fr) ={acosd(min(1,max(-1, gPG3vec_L(:).' *OpFlVec_L(:) / norm(gPG3vec_L) / norm(OpFlVec_L) )))};
	end
	% delineate this face's resultant vector angle
	L_CM=circ_mean(Thetas_L);
	% circular SD
	L_CSD=circ_std(Thetas_L);
	OutDf_L(F,9)=num2cell(L_CSD);
	% arbitrary but consistent circ_mean
	OutDf_L(F,10)=num2cell(L_CM);
	% get angular distance from gPGG
	PGGang_L=cart2pol(gazes_L(F),gels_L(F));
	OutDf_L(F,8)=num2cell(PGGang_L-L_CM);
	% add magnitude time series
	OutTs_L(F,2,:)=num2cell(Mags_L);
	% add theta time series
	OutTs_L(F,1,:)=num2cell(Thetas_L);
	% extract indices of bottom up (broadly, dist < 90)
        FaceAngDistPGG_L=OutTs_L{F,3,:};
        BU_Trs_L=find(FaceAngDistPGG_L<90);
        % get inverse indices for top down (broadly, dist > 90)
        TD_Trs_L=find(FaceAngDistPGG_L>90);
        % get proportion of TRs that are BU for this face
        propBU_L=(length(BU_Trs_L))/NumTRs;
        % plop prop in output df, 1st column
        OutDf_L(F,1)=num2cell(propBU_L);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Now for the right hemisphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note we are loopign over the non-mw-face indices, skipping medial wall faces but leaving 0's in their stead
for F=g_noMW_combined_R
        % get angles into needed format
        % note azimuth elevation ordering for atan2d
        gPGvec_R=[gazes_R(F) gels_R(F)];
        gPG2vec_R=[pg2azes_R(F) pg2els_R(F)];
        gPG3vec_R=[pg3azes_R(F) pg3els_R(F)];
	Bvec_R=[bazes_L(F) bels_R(F)];
        G1vec_R=[g1azes_R(F) g1els_R(F)];
        G2vec_R=[g2azes_R(F) g2els_R(F)];
	% translate xyz vector fields from opfl to az/el/r to polar angles
        % note, by not saving rho (just theta), we are discarding magnitude information at this point
        Thetas_R=zeros(1,lenOpFl);
        Mags_R=zeros(1,lenOpFl);
        gangDist_R=zeros(1,lenOpFl);
        BangDist_R=zeros(1,lenOpFl);
        G1angDist_R=zeros(1,lenOpFl);
        G2angDist_R=zeros(1,lenOpFl);
	PG1angDist_R=zeros(1,lenOpFl);
        PG2angDist_R=zeros(1,lenOpFl);
        for fr=1:lenOpFl
                % current vector field
                relVf_R=vfr{fr};
                % xyz components
                xComp_R=relVf_R(F,1);
                yComp_R=relVf_R(F,2);
                zComp_R=relVf_R(F,3);
                % convert to spherical coord system
                vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(F),eld_R(F));
                % convert to spherical coordinates
                OpFlVec_R= [vs_R(1) vs_R(2)];
                % store in output vector (r is redundant across all vecs, only using az and el)
                [Thetas_R(fr),Mags_R(fr)]=cart2pol(vs_R(1),vs_R(2));
                % store rho in magnitude series
                % store az and el
                %OutTs_R(F,3,fr) = {vs_R(1)};
                %OutTs_R(F,4,fr) = {vs_R(2)};
                %%%% angular distance calculations
                % PG ang d
                OutTs_R(F,3,fr)={acosd(min(1,max(-1, gPGvec_R(:).' *OpFlVec_R(:) / norm(gPGvec_R) / norm(OpFlVec_R) )))};
                % beta g ang d
                OutTs_R(F,4,fr) ={acosd(min(1,max(-1, Bvec_R(:).' *OpFlVec_R(:) / norm(Bvec_R) / norm(OpFlVec_R) )))};
                % g1 ang d
                OutTs_R(F,5,fr) ={acosd(min(1,max(-1, G1vec_R(:).' *OpFlVec_R(:) / norm(G1vec_R) / norm(OpFlVec_R) )))};
                % g2 ang d
                OutTs_R(F,6,fr) ={acosd(min(1,max(-1, G2vec_R(:).' *OpFlVec_R(:) / norm(G2vec_R) / norm(OpFlVec_R) )))};
        	% PG2 ang d
		OutTs_R(F,7,fr) ={acosd(min(1,max(-1, gPG2vec_R(:).' *OpFlVec_R(:) / norm(gPG2vec_R) / norm(OpFlVec_R) )))};
                % PG3 ang d
                OutTs_R(F,8,fr) ={acosd(min(1,max(-1, gPG3vec_R(:).' *OpFlVec_R(:) / norm(gPG3vec_R) / norm(OpFlVec_R) )))};
	end
        % delineate this face's resultant vector angle
        R_CM=circ_mean(Thetas_R);
        % circular SD
        R_CSD=circ_std(Thetas_R);
        OutDf_R(F,9)=num2cell(R_CSD);
        % arbitrary but consistent circ_mean
        OutDf_R(F,10)=num2cell(R_CM);
        % get angular distance from gPGG
        PGGang_R=cart2pol(gazes_R(F),gels_R(F));
        OutDf_R(F,8)=num2cell(PGGang_R-R_CM);
        % add magnitude time series
        OutTs_R(F,2,:)=num2cell(Mags_R);
        % add theta time series
        OutTs_R(F,1,:)=num2cell(Thetas_R);
	 % extract indices of bottom up (broadly, dist < 90)
        FaceAngDistPGG_R=OutTs_R{F,3,:};
        BU_Trs_R=find(FaceAngDistPGG_R<90);
        % get inverse indices for top down (broadly, dist > 90)
        TD_Trs_R=find(FaceAngDistPGG_R>90);
        % get proportion of TRs that are BU for this face
        propBU_R=(length(BU_Trs_R))/NumTRs;
        % plop prop in output df, 1st column
        OutDf_R(F,1)=num2cell(propBU_R);
end

%%% save output df
%outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_run_' runNum '_BUTD_L_resultantVecs.mat'];
%save(outFP_L,'OutDf_L')
%outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_run_' runNum '_BUTD_R_resultantVecs.mat'];
%save(outFP_R,'OutDf_R')
%toc
%%% MW MASK
%OutDf_L_masked=OutDf_L(g_noMW_combined_L,:);
%OutDf_R_masked=OutDf_R(g_noMW_combined_R,:);
% saveout
%outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_run_' runNum '_BUTD_L_resultantVecs_masked.mat'];
%save(outFP_L,'OutDf_L_masked')
%outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_run_' runNum '_BUTD_R_resultantVecs_masked.mat'];
%save(outFP_R,'OutDf_R_masked')

%%%% faces-to-vertices converter

% init vertices cell array with []s for each vertex as a cell
LvertTS=cell(length(vx_l),8,lenOpFl);
RvertTS=cell(length(vx_r),8,lenOpFl);

% for each face
for F=g_noMW_combined_L;
	% get F values
	Fvalues=OutTs_L(F,:,:);
	% extract vertices comrpising F
	curVert=faces_l(F,:);
	% for each timepoint
	for tp=1:lenOpFl
		% append vertex1 cell
		for m=1:8
			LvertTS{curVert(1),m,tp}=[LvertTS{curVert(1),m,tp} OutTs_L(F,m,tp)];
			% append vertex2 cell
			LvertTS{curVert(2),m,tp}=[LvertTS{curVert(2),m,tp} OutTs_L(F,m,tp)];
			% append vertex3 cell
			LvertTS{curVert(3),m,tp}=[LvertTS{curVert(3),m,tp} OutTs_L(F,m,tp)];
		end
	end
end
% now average over each vertex value per timepoints
for v=1:length(vx_l)
	for tp=1:lenOpFl
		% sep mean calc needed for theta
		for m=1
			if ~isempty(LvertTS{v,m,tp});
                                InterfacingFaces=LvertTS{v,m,tp};
                                LvertTS{v,m,tp}=circ_mean([InterfacingFaces{:}]);
                        end
		end
		for m=2:8
			if ~isempty(LvertTS{v,m,tp});
				InterfacingFaces=LvertTS{v,m,tp};
				LvertTS{v,m,tp}=mean([InterfacingFaces{:}]);
			end
		end
	end
end
% right
for F=g_noMW_combined_R;
        % get F values
        Fvalues=OutTs_R(F,:,:);
        % extract vertices comrpising F
        curVert=faces_r(F,:);
        % for each timepoint
        for tp=1:lenOpFl
                % append vertex1 cell
                for m=1:8
                        RvertTS{curVert(1),m,tp}=[RvertTS{curVert(1),m,tp} OutTs_R(F,m,tp)];
                        % append vertex2 cell
                        RvertTS{curVert(2),m,tp}=[RvertTS{curVert(2),m,tp} OutTs_R(F,m,tp)];
                        % append vertex3 cell
                        RvertTS{curVert(3),m,tp}=[RvertTS{curVert(3),m,tp} OutTs_R(F,m,tp)];
                end
        end
end

% now average over each vertex value per timepoints
for v=1:length(vx_r)
        for tp=1:lenOpFl
                % sep mean calc needed for theta
		for m=1
			if ~isempty(RvertTS{v,m,tp});
				InterfacingFaces=RvertTS{v,m,tp};
				RvertTS{v,m,tp}=circ_mean([InterfacingFaces{:}]);
			end
		end
		for m=2:8
                        if ~isempty(RvertTS{v,m,tp});
				InterfacingFaces=RvertTS{v,m,tp};
                        	RvertTS{v,m,tp}=mean([InterfacingFaces{:}]);
                	end
		end
        end
end

% saveout
%outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_run_' runNum '_BUTD_L_resultantVecs_masked_ts.mat'];
save(outFP_L,'LvertTS','-v7.3')
%outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_run_' runNum '_BUTD_R_resultantVecs_masked_ts.mat'];
save(outFP_R,'RvertTS','-v7.3')
