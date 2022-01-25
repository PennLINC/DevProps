function ClassifyAngDistrs(subj)

% designed to operate at the pixel-level within subjects
% four outcomes as is: Random (fails FDR'ed HR test), Unimodal (fails hartigan's dip test), or bimodal (passes HR, passes hartigan's)
% unimodal further distinguished into BU or TD with circular mean: save into different struct

% add matlab path for used functions
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% convert to character vector
sname=char(subj);
childfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/'];
% in this case output dir is same as input dir
parentfp=childfp;
OpFlfn=[childfp sname '_OpFl_fs5.mat'];
% load in this subjs opFlow
OpFlmat=load(OpFlfn);

% vector fields
vfl=OpFlmat.us.vf_left;
vfr=OpFlmat.us.vf_right;

%%%%%%% pull in PGG masks
% Load in surface data
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
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
% now load in medial wall VERTICES
mw_v_l=read_medial_wall_label([SubjectsFolder '/label/lh.Medial_wall.label']);
mw_v_r=read_medial_wall_label([SubjectsFolder '/label/rh.Medial_wall.label']);
% all faces that touch a medial wall vertex to be masked
MW_f1_L=find(ismember(F_L(:,1),mw_v_l));
MW_f2_L=find(ismember(F_L(:,2),mw_v_l));
MW_f3_L=find(ismember(F_L(:,3),mw_v_l));
% rh
MW_f1_R=find(ismember(F_R(:,1),mw_v_r));
MW_f2_R=find(ismember(F_R(:,2),mw_v_r));
MW_f3_R=find(ismember(F_R(:,3),mw_v_r));
% inclusive - mask if face involves ANY mw vertices
MW_combined_L=union(MW_f1_L,MW_f2_L);
MW_combined_L=union(MW_combined_L,MW_f3_L);
% now for right hemisphere
MW_combined_R=union(MW_f1_R,MW_f2_R);
MW_combined_R=union(MW_combined_R,MW_f3_R);
% get inverse for indexing : faces that ARE NOT touching mW verts
noMW_combined_L=setdiff([1:20480],MW_combined_L);
noMW_combined_R=setdiff([1:20480],MW_combined_R);
% further mask derivation
% MASK WHERE PGG = 0: individ AND group
% load in GROUP PG
gLPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/Gradients.lh.fsaverage5.func.gii'];
gLPGf=gifti(gLPGfp);
gPG_LH=gLPGf.cdata(:,1);
% right hemi
gRPGfp=['/cbica/projects/pinesParcels/data/princ_gradients/Gradients.rh.fsaverage5.func.gii'];
gRPGf=gifti(gRPGfp);
gPG_RH=gRPGf.cdata(:,1);
% load in subject's PG
LPGfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/' sname '_PG_L_10k_rest.func.gii'];
LPGf=gifti(LPGfp);
PG_LH=LPGf.cdata(:,1);
% right hemi
RPGfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' subj '/' subj '_PG_R_10k_rest.func.gii'];
RPGf=gifti(RPGfp);
PG_RH=RPGf.cdata(:,1);
% calculate PG gradient on sphere
PGg_L = grad(F_L, V_L, PG_LH);
PGg_R = grad(F_R, V_R, PG_RH);
% calculate group PG gradient on sphere
gPGg_L = grad(F_L, V_L, gPG_LH);
gPGg_R = grad(F_R, V_R, gPG_RH);
% get index of where they are 0 in all directions
PGg_L0=find(all(PGg_L')==0);
gPGg_L0=find(all(gPGg_L')==0);
PGg_R0=find(all(PGg_R')==0);
gPGg_R0=find(all(gPGg_R')==0);
% continue to get unions
ind_MW_combined_L=union(MW_combined_L,PGg_L0);
gro_MW_combined_L=union(MW_combined_L,gPGg_L0);
% and right hemi
ind_MW_combined_R=union(MW_combined_R,PGg_R0);
gro_MW_combined_R=union(MW_combined_R,gPGg_R0);
% get inverse for indexing : faces that ARE NOT touching mW verts
i_noMW_combined_L=setdiff([1:20480],ind_MW_combined_L);
i_noMW_combined_R=setdiff([1:20480],ind_MW_combined_R);
g_noMW_combined_L=setdiff([1:20480],gro_MW_combined_L);
g_noMW_combined_R=setdiff([1:20480],gro_MW_combined_R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;
%%%% get spherical coordinates
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get length of angular time series
lenOpFl=length(OpFlmat.us.vf_left);
% get size of utilized sphere in faces
singleVecFieldsizeL=size(vfl);
singleVecFieldsizeR=size(vfr);
sphSizeL=singleVecFieldsizeL(1);
sphSizeR=singleVecFieldsizeR(1);

% translate xyz vector fields from opfl to az/el/r
azesOpf_L=zeros(sphSizeL,lenOpFl);
elsOpf_L=zeros(sphSizeL,lenOpFl);
for i=1:sphSizeL
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
azesOpf_R=zeros(sphSizeR,lenOpFl);
elsOpf_R=zeros(sphSizeR,lenOpFl);
for i=1:sphSizeR
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

% mask the angles
azs_L_masked=azesOpf_L(g_noMW_combined_L,:);
els_L_masked=elsOpf_L(g_noMW_combined_L,:);
azs_R_masked=azesOpf_R(g_noMW_combined_R,:);
els_R_masked=elsOpf_R(g_noMW_combined_R,:);

% initialize output map to store class membership of each pixel
FaceClassL=zeros(1,length(azs_L_masked));
FaceClassR=zeros(1,length(azs_R_masked));
% initialize output map for circular mean of unimodal/non-randomly distributed circular means
FaceMeanL=zeros(1,length(azs_L_masked));
FaceMeanR=zeros(1,length(azs_R_masked));

%%%%%%%%%% HR TEST FOR NON-NORMALITY OF CIRCULAR DISTRIBUTIONS
% initialize output arrays
% HR test - Left
pVecL=zeros(1,length(azs_L_masked));
%tVecL=zeros(1,length(Lrow));
% t-stat map = left
HR_t_L=zeros(1,length(azs_L_masked));
% HR test each face
for P=1:length(azs_L_masked);
	P
	% extract the angular distribution of this pixel over all frames and segments
	FaceDistrX=azs_L_masked(P,:);
	FaceDistrY=els_L_masked(P,:);
	% convert to radians
	FaceRads=cart2pol(FaceDistrX,FaceDistrY);
	% HR test
	[p,T]=hrtest(FaceRads(:),1000);
	pVecL(P)=p;
	% save out test stat purely for viz
	HR_t_L(P)=T;
end

% HR test - Right
pVecR=zeros(1,length(azs_R_masked));
%tVecL=zeros(1,length(Lrow));
% t-stat map = left
HR_t_R=zeros(1,length(azs_R_masked));
% HR test each face
for P=1:length(azs_R_masked);
        P
        % extract the angular distribution of this pixel over all frames and segments
        FaceDistrX=azs_R_masked(P,:);
        FaceDistrY=els_R_masked(P,:);
        % convert to radians
        FaceRads=cart2pol(FaceDistrX,FaceDistrY);
        % HR test
        [p,T]=hrtest(FaceRads(:),1000);
        pVecR(P)=p;
        % save out test stat purely for viz
        HR_t_R(P)=T;
end

% saveout
HRT_outFN_L=[parentfp sname '_HRT_L.mat'];
save(HRT_outFN_L,'HR_t_L')
HRT_outFN_R=[parentfp sname '_HRT_R.mat'];
save(HRT_outFN_R,'HR_t_R')
HRT_outFN_L=[parentfp sname '_pvecL.mat'];
save(P_outFN_L,'pVecL')
HRT_outFN_R=[parentfp sname '_pvecR.mat'];
save(P_outFN_R,'pVecR')


%%%%%%%%%%%%%%%%% not adapted past this point


%%%% FDR
% merge together pvecs for full-brain MC correction
%HR_mergedPvec=cat(2,pVecL,pVecR);
%FDRed=mafdr(HR_mergedPvec);

% split back into 2
%FDRed_L=FDRed(1:length(Lrow));
%FDRed_R=FDRed((length(Lrow)+1):(length(Lrow)+length(Rrow)));

% classify passed HR-FDR as 1 - Left
%NonRand_L=find(FDRed_L<0.05);
%PixClassL(NonRand_L)=1;
% update Lrow and Lcol with only non-randomly distributed pixels
%Lrow_NR=Lrow(NonRand_L);
%Lcol_NR=Lcol(NonRand_L);

% classify passed HR-FDR as 1 - Right
%NonRand_R=find(FDRed_R<0.05);
%PixClassR(NonRand_R)=1;
% update Lrow and Lcol with only non-randomly distributed pixels
%Rrow_NR=Rrow(NonRand_R);
%Rcol_NR=Rcol(NonRand_R);

%%%%%%%%%% END OF HR TESTS

% comment out later, but nice update in the meantime
%length(Lrow_NR)
%length(Rrow_NR)

%%%%%%%%% HARTIGAN'S DIP TEST FOR NON-UNIMODALITY

% subject remaining to hartigan's dip test - run on angular distance from gradient gradients
% load in in angular distance files
%AngDistFn=[childfp sname '_BU_angDist.csv'];
%AngDist=readtable(AngDistFn);

% initialize hartigan vectors
%hVecL_p=zeros(1,length(Lrow));
%hVecL_d=zeros(1,length(Lrow));
% set invalids to NaN to avoid confusion with low p-values
%hVecL_p(hVecL_p==0)=NaN;
%hVecL_d(hVecL_d==0)=NaN;

%for P=1:length(Lrow);
%	P
%	% if this distribution is non-random
%	if FDRed_L(P) < 0.05
%		FDRed_L(P)
%		% extract angular distance distribution from PG in this pixel
%		AngDistaDistr=table2array(AngDist(:,P));
%		% Hartigan's dip test - https://snl.salk.edu/~jude/waveform_public/HartigansDipSignifTest.m
%		[dip,p_value,xlow,xup]=HartigansDipSigniftest(AngDistaDistr,1000);
%		% pack the dip into the vector
%		hVecL_p(P)=p_value;
%		hVecL_d(P)=dip;
%	end
%end

% initialize hartigan vectors
%hVecR_p=zeros(1,length(Rrow_NR));
%hVecR_d=zeros(1,length(Rrow_NR));
% set invalids to NaN to avoid confusion with low p-values
%hVecR_p(hVecR_p==0)=NaN;
%hVecR_d(hVecR_d==0)=NaN;

%for P=1:length(Rrow);
%        P
%        % if this distribution is non-random
%	if FDRed_R(P) < 0.05
%		FDRed_R(P)
%		% extract angular distance distribution from PG in this pixel
 %       	AngDistaDistr=table2array(AngDist(:,P+length(Lrow)));
        	% Hartigan's dip test - https://snl.salk.edu/~jude/waveform_public/HartigansDipSignifTest.m
  %      	[dip,p_value,xlow,xup]=HartigansDipSigniftest(AngDistaDistr,1000);
        	% pack the dip into the vector
  %      	hVecR_p(P)=p_value;
   %     	hVecR_d(P)=dip;
%	end
%end


% merge dip test pvalus to one map for FDR
%HAR_mergedPvec=cat(2,hVecL_p,hVecR_p);
%HAR_FDRed=mafdr(HAR_mergedPvec);

% split back into 2
%HAR_FDRed_L=HAR_FDRed(1:length(Lrow));
%HAR_FDRed_R=HAR_FDRed((length(Lrow)+1):(length(Lrow)+length(Rrow)));

% change pixClass from 1 to 2 for non-unimodal 
%NonUni_L=find(HAR_FDRed_L<0.05);
%PixClassL(NonUni_L)=2;
% and for right hemi
%NonUni_R=find(HAR_FDRed_R<0.05);
%PixClassL(NonUni_R)=2;
%%%%%%%%%%%%%%%% END HARTIGAN'S

% finally, for unimodal distributions, calculate circular mean for viz
% use downsampled flat PG to reconstruct original ordering of this table

%%%%%%%%%% CIRCULAR MEANS
% intialize circMean vectors
%CM_VecL=zeros(1,length(Lrow));
% set invalids to NaN to avoid confusion with low p-values
%CM_VecL(CM_VecL==0)=NaN;

% loop over every pixel
%for P=1:length(Lrow);
%        P
        % if this distribution is non-random
%        if FDRed_L(P) < 0.05
		% AND if the distribution is not multi-modal
%                if HAR_FDRed_L(P) > 0.05
	               % get row and column of this pixel
%		        Row=Lrow(P);
%		        Col=Lcol(P);
		        % extract the angular distribution of this pixel over all frames and segments
%		        PixelDistrX=real(LeftVFs(Row,Col,:));
%		        PixelDistrY=imag(LeftVFs(Row,Col,:));
		        % convert to radians
%		        PixelRads=cart2pol(PixelDistrX,PixelDistrY);
			% calculate circular mean
%			[mu,ul,ll]=circ_mean(PixelRads(:));
	                % inset circular mean into CM vec
%			CM_VecL(P)=mu;
%		end
%        end
%end

% intialize circMean vectors
%CM_VecR=zeros(1,length(Rrow));
% set invalids to NaN to avoid confusion with low p-values
%CM_VecR(CM_VecR==0)=NaN;

% loop over every right-hemi pixel
%for P=1:length(Rrow);
%        P
        % if this distribution is non-random
%        if FDRed_R(P) < 0.05
                % AND if the distribution is not multi-modal
%                if HAR_FDRed_R(P) > 0.05
                       % get row and column of this pixel
%                        Row=Rrow(P);
%                        Col=Rcol(P);
                        % extract the angular distribution of this pixel over all frames and segments
%                        PixelDistrX=real(RightVFs(Row,Col,:));
%                        PixelDistrY=imag(RightVFs(Row,Col,:));
                        % convert to radians
%                        PixelRads=cart2pol(PixelDistrX,PixelDistrY);
                        % calculate circular mean
%                        [mu,ul,ll]=circ_mean(PixelRads(:));
                        % inset circular mean into CM vec
%                        CM_VecR(P)=mu;
%                end
%        end
%end


% further classify those failing the dip test on basis of circular means
% change pixClass from 1 to 2 for non-unimodal 
%UniCM_L=find(~isnan(CM_VecL));
%PixMeanL(UniCM_L)=CM_VecL(UniCM_L);
% and for right hemi
%UniCM_L=find(~isnan(CM_VecL));
%PixMeanL(UniCM_L)=CM_VecL(UniCM_L);

% ... maybe map back to real space before leaving this maze?

% use gradient gradient map for template
%circMeanMapL=PGBU_L_x;
%classMapL=PGBU_L_x;
% for left
%for i=1:length(Lrow);
%        Row=Lrow(i);
%        Col=Lcol(i);
%	circMeanMapL(Row,Col)=CM_VecL(i);
%	classMapL(Row,Col)=PixClassL(i);
%end

% use gradient gradient map for template
%circMeanMapR=PGBU_R_x;
%classMapR=PGBU_R_x;
% for Right
%for i=1:length(Rrow);
%        Row=Rrow(i);
%        Col=Rcol(i);
%	circMeanMapR(Row,Col)=CM_VecR(i);
%        classMapR(Row,Col)=PixClassR(i);
%end

% saveout left hemi
%circM_outFN_L=[parentfp sname '_circM_L.mat'];
%classM_outFN_L=[parentfp sname '_classM_L.mat'];
%save(circM_outFN_L,'circMeanMapL')
%save(classM_outFN_L,'classMapL')
%HRT_outFN_L=[parentfp sname '_HRT_L.mat'];
%save(HRT_outFN_L,'HR_t_L')

% and saveout right hemi
%circM_outFN_R=[parentfp sname '_circM_R.mat'];
%classM_outFN_R=[parentfp sname '_classM_R.mat'];
%save(circM_outFN_R,'circMeanMapR')
%save(classM_outFN_R,'classMapR')
%HRT_outFN_R=[parentfp sname '_HRT_R.mat'];
%save(HRT_outFN_R,'HR_t_R')



