function Bin_And_Aggregate_PGGDistributions(subj)
% read in angular distance from PGG: we want circular and linear angular distances aggregated across subjs

% to make this computationally feasible, read each subj and collapse across all TR across all faces into BINS

% load in gPPG ang dist: for thresholdingdd OFD toolbox to path
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% Load in fsav4 opflow calc
OpFlFp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_OpFl_fs4.mat'];
data=load(OpFlFp);
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

% initialize output (counts in angular bins: 37 for 0 - pi to-be circular hist, but to have distinct bars at 0 and pi)
OutDf=zeros(1,36);
Bins36=-180:10:180;
% initialize facewise vectors for 1-36 mode for 
faceModesL=zeros(1,length(g_noMW_combined_L));
faceModesR=zeros(1,length(g_noMW_combined_R));


%%% WILL NEED SEP SCRIPT TO AGGREGATE ACROSS SUBJS, CALC DIF HISTS IN ACCORDANCE WITH RS VS TASK, OLD VS. YOUNG, ETC



% count Number of TRs once rather than iteratively
NumTRs=size(AngDist.gLeft);
NumTRs=NumTRs(1);
lenOpFl=NumTRs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Starting for the left hemisphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note we are looping over the non-mw-face indices, skipping medial wall faces but leaving 0's in their stead
for F=g_noMW_combined_L
	% get angular distance from gPGG
	PGGang_L=[gazes_L(F) gels_L(F)];
	% initialize vector for angle at every TR pair
	FAngles=zeros(1,lenOpFl);
	for fr=1:lenOpFl
		% current vector field
	        relVf_L=vfl{fr};
		% xyz components
	        xComp_L=relVf_L(F,1);
	        yComp_L=relVf_L(F,2);
	        zComp_L=relVf_L(F,3);
		% convert to spherical coord system
	        vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
		% extract OpFl angle
		Ang_L=[vs_L(1) vs_L(2)];
		%https://de.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
		x1=PGGang_L(1);
		y1=PGGang_L(2);
		x2=Ang_L(1);
		y2=Ang_L(2);
		AngDist= atan2d(x1*y2-y1*x2,x1*x2+y1*y2);	
		% flag distance stepping out of -pi to pi range
		if AngDist > 180
		disp ('dist > 180')
		%	AngDist = AngDist - pi;
		end
		if AngDist < -180
		%	AngDist = AngDist + pi;
		disp('dist <-180')
		end
		FAngles(fr)=AngDist;
	end
	% discretize
	Disc_FAngles=histcounts(FAngles,Bins36);
	OutDf=OutDf+Disc_FAngles;
	% get modes for this face
	[M,I]=max(Disc_FAngles);
	faceModesL(F)=I;
end

% mask (face 5120 saved in as entry 5120 of facemodes, etc.)
faceModesL=faceModesL(g_noMW_combined_L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% now the right hemisphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note we are looping over the non-mw-face indices, skipping medial wall faces but leaving 0's in their stead
for F=g_noMW_combined_R
        % get angular distance from gPGG
        PGGang_R=[gazes_R(F) gels_R(F)];
        % initialize vector for angle at every TR pair
        FAngles=zeros(1,lenOpFl);
        for fr=1:lenOpFl
                % current vector field
                relVf_R=vfr{fr};
                % xyz components
                xComp_R=relVf_R(F,1);
                yComp_R=relVf_R(F,2);
                zComp_R=relVf_R(F,3);
                % convert to spherical coord system
                vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(F),eld_R(F));
                % extract OpFl angle
                Ang_R=[vs_R(1) vs_R(2)];
                %https://de.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
                x1=PGGang_R(1);
                y1=PGGang_R(2);
                x2=Ang_R(1);
                y2=Ang_R(2);
                AngDist= atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
                % flag distance stepping out of -pi to pi range
                if AngDist > 180
                disp ('dist > 180')
                %       AngDist = AngDist - pi;
                end
                if AngDist < -180
                %       AngDist = AngDist + pi;
                disp('dist <-180')
		end
		FAngles(fr)=AngDist;
        end
        % discretize
        Disc_FAngles=histcounts(FAngles,Bins36);
        OutDf=OutDf+Disc_FAngles;
        % get modes for this face
        [M,I]=max(Disc_FAngles);
        faceModesR(F)=I;
end

% mask (face 5120 saved in as entry 5120 of facemodes, etc.)
faceModesR=faceModesR(g_noMW_combined_R);

% save outdf
fn=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistHist.mat'];
save(fn,'OutDf')

fn=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_DistHist.csv'];
writetable(table(OutDf),fn)


% save facial modes
fn=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistFacial_360ModesL.csv'];
writetable(table(faceModesL),fn)
fn=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_AngDistFacial_360ModesR.csv'];
writetable(table(faceModesR),fn)
