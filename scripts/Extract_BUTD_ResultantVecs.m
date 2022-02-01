function Extract_BUTD_ResultantVecs(subj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each face. Measured relative to gPGG.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temp to evaluate compute time
tic

% load in gPPG ang dist: for thresholdingdd OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% Load in fsav4 opflow calc
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'];
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

% load in GROUP PG
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_L_3k.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/hcp.gradients_R_3k.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);

%%% For masking out the medial wall
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);
% get index of where they are 0 in all directions
gPGg_L0=find(all(gPGg_L')==0);
gPGg_R0=find(all(gPGg_R')==0);
% get inverse for indexing : faces that ARE NOT touching mW verts
g_noMW_combined_L=setdiff([1:5120],gPGg_L0);
g_noMW_combined_R=setdiff([1:5120],gPGg_R0);

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
end
% right hemi
gazes_R=zeros(1,length(azd_R));
gels_R=zeros(1,length(eld_R));
for i=1:length(azd_R)
    gvs_R=cart2sphvec(double([gPGx_R(i);gPGy_R(i);gPGz_R(i)]),azd_R(i),eld_R(i));
    gazes_R(i)=gvs_R(1);
    gels_R(i)=gvs_R(2);
end

% load in gPGG angular distances for parsing into top-down and bottom-up in the loops
AngDistFP=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistMat4.mat'];
AngDist=load(AngDistFP);
% just pretend this line doesn't exist
AngDist=AngDist.AngDist;

% initialize output (bu_az,bu_el,td_az,td_el,propBu for each face)
% and last 2 columns are BU vec length, TD vec length
OutDf_L=cell(length(F_L),7);
OutDf_R=cell(length(F_R),7);

% count Number of TRs once rather than iteratively
NumTRs=size(AngDist.gLeft);
NumTRs=NumTRs(1);
lenOpFl=NumTRs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Starting for the left hemisphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note we are loopign over the non-mw-face indices, skipping medial wall faces but leaving 0's in their stead
for F=g_noMW_combined_L
	% extract indices of bottom up (broadly, dist < 90)
	FaceAngDistPGG_L=AngDist.gLeft(:,F);
	BU_Trs_L=find(FaceAngDistPGG_L<90);
	% get inverse indices for top down (broadly, dist > 90)
	TD_Trs_L=find(FaceAngDistPGG_L>90);

	% get proportion of TRs that are BU for this face
	propBU_L=(length(BU_Trs_L))/NumTRs;

	% plop prop in output df, 1st column
	OutDf_L(F,1)=num2cell(propBU_L);	
	
	% get angles into needed format

	% translate xyz vector fields from opfl to az/el/r to polar angles
	% note, by not saving rho (just theta), we are discarding magnitude information at this point
	Thetas_L=zeros(1,lenOpFl);
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
		[Thetas_L(fr),rho]=cart2pol(vs_L(1),vs_L(2));
	end

%%%%%%%%%%%% Bottom-up

% contingent on BU TRs existing. Note this is within the F loop despite indentation.
if length(BU_Trs_L)>0
	% unfortunately, will require splitting up left and right and making EACH a contingency (4 contg, BULBUR,TDLTDR)
	% get bottom up TRs (broadly)	
	BUAngs_L=Thetas_L(BU_Trs_L);

	% get resVec thetas
	BU_L_CM=circ_mean(BUAngs_L);

	% center on angular distance from gPGG
	PGGang_L=cart2pol(gazes_L(F),gels_L(F));
	BU_L_CM_rel=PGGang_L-BU_L_CM;

	% get BU resultant vector angle, convert back to cartesian in the proccess (1 as fill-in for rho, discards OpFl magn.)
	[BUHorzC_L,BUVertC_L]=pol2cart(BU_L_CM_rel,1);

        % get BU resultant vector length
        VL_L=circ_r(BUAngs_L);
        OutDf_L(F,6)=num2cell(VL_L);

	% if in valid face, scale x y to vector length and plop in output df
	if (std(FaceAngDistPGG_L)~=0)
		BUHorzC_L=BUHorzC_L*VL_L;
		BUVertC_L=BUVertC_L*VL_L;
		OutDf_L(F,2)=num2cell(BUHorzC_L);
		OutDf_L(F,3)=num2cell(BUVertC_L);
	end
end
% end "if bottom down TRs instances at this face exist" contingency

%%%%%%%%%% Top-down

% same contingency as above
if length(TD_Trs_L)>0
        % get topdown TRs (broadly)
        TDAngs_L=Thetas_L(TD_Trs_L);

         % get resVec thetas
        TD_L_CM=circ_mean(TDAngs_L);

        % center on angular distance from gPGG
        TD_L_CM_rel=PGGang_L-TD_L_CM;

        % get TD resultant vector angle
        [TDHorzC_L,TDVertC_L]=pol2cart(TD_L_CM_rel,1);

        % get TD resultant vector length
        VL_L=circ_r(TDAngs_L);
        OutDf_L(F,7)=num2cell(VL_L);

        % if in valid face, scale x y to vector length and plop in output df
        if (std(FaceAngDistPGG_L)~=0)
                TDHorzC_L=TDHorzC_L*VL_L;
                TDVertC_L=TDVertC_L*VL_L;
                OutDf_L(F,4)=num2cell(TDHorzC_L);
                OutDf_L(F,5)=num2cell(TDVertC_L);
        end
end
% end looping over every face
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% the right hemisphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for F=g_noMW_combined_R
	FaceAngDistPGG_R=AngDist.gRight(:,F);
	BU_Trs_R=find(FaceAngDistPGG_R<90);
	TD_Trs_R=find(FaceAngDistPGG_R>90);
	propBU_R=(length(BU_Trs_R))/NumTRs;
	OutDf_R(F,1)=num2cell(propBU_R);
	Thetas_R=zeros(1,lenOpFl);
        for fr=1:lenOpFl
                % current vector field
                relVf_R=vfr{fr};
                % xyz components
                xComp_R=relVf_R(F,1);
                yComp_R=relVf_R(F,2);
                zComp_R=relVf_R(F,3);
                % convert to spherical coord system
                vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(F),eld_R(F));
                % store in output vector (r is redundant across all vecs)
                [Thetas_R(fr),rho]=cart2pol(vs_R(1),vs_R(2));
        end
%%%%%%% Bottom up
if length(BU_Trs_R)>0
	BUAngs_R=Thetas_R(BU_Trs_R);
	BU_R_CM=circ_mean(BUAngs_R);
	PGGang_R=cart2pol(gazes_R(F),gels_R(F));
	BU_R_CM_rel=PGGang_R-BU_R_CM;	
	[BUHorzC_R,BUVertC_R]=pol2cart(BU_R_CM_rel,1);
	VL_R=circ_r(BUAngs_R);
	OutDf_R(F,6)=num2cell(VL_R);
	if (std(FaceAngDistPGG_R)~=0)
                BUHorzC_R=BUHorzC_R*VL_R;
                BUVertC_R=BUVertC_R*VL_R;
                OutDf_R(F,2)=num2cell(BUHorzC_R);
                OutDf_R(F,3)=num2cell(BUVertC_R);
        end
end

%%%%%% Top-down
if length(TD_Trs_R)>0
        TDAngs_R=Thetas_R(TD_Trs_R);
	TD_R_CM=circ_mean(TDAngs_R);
	TD_R_CM_rel=PGGang_R-TD_R_CM;
        [TDHorzC_R,TDVertC_R]=pol2cart(TD_R_CM_rel,1);
        VL_R=circ_r(TDAngs_R);
        OutDf_R(F,7)=num2cell(VL_R);	
	if (std(FaceAngDistPGG_R)~=0)
                TDHorzC_R=TDHorzC_R*VL_R;
                TDVertC_R=TDVertC_R*VL_R;
                OutDf_R(F,4)=num2cell(TDHorzC_R);
                OutDf_R(F,5)=num2cell(TDVertC_R);
        end	
end

% end looping over every face
end

%%% save output df
outFP_L=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_L_resultantVecs.mat'];
save(outFP_L,'OutDf_L')
outFP_R=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_BUTD_R_resultantVecs.mat'];
save(outFP_R,'OutDf_R')


toc
