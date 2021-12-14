function ClassifyAngDistrs(subj)

% designed to operate at the pixel-level within subjects
% four outcomes as is: Random (fails FDR'ed HR test), Unimodal (fails hartigan's dip test), or bimodal (passes HR, passes hartigan's)
% unimodal further distinguished into BU or TD with circular mean: save into different struct

% add matlab path for used functions
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/hermans-rasson'));

% convert to character vector
sname=char(subj);
childfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/'];

% load in this subjs opFlow
OpFlfn=['/cbica/projects/pinesParcels/results/OpFl_output/' sname '/OpFlowResults.mat'];
OpFlmat=load(OpFlfn);

% merge vector fields across runs (concatenate across time dimension)
LeftVFs=cat(3,OpFlmat.MegaStruct.Vf_Left{:});
RightVFs=cat(3,OpFlmat.MegaStruct.Vf_Right{:});

% use downsampled flat PG for reference
PGBU_L_x=load([childfp sname '_PG_GxL_BU.csv']);
PGBU_R_x=load([childfp sname '_PG_GxR_BU.csv']);
[Lrow,Lcol]=find(~isnan(PGBU_L_x));
[Rrow,Rcol]=find(~isnan(PGBU_R_x));

% initialize output map to store class membership of each pixel
PixClassL=zeros(1,length(Lrow));
PixClassR=zeros(1,length(Rrow));
% initialize output map for circular mean of unimodal/non-randomly distributed circular means
PixMeanL=zeros(1,length(Lrow));
PixMeanR=zeros(1,length(Rrow));

%%%%%%%%%% HR TEST FOR NON-NORMALITY OF CIRCULAR DISTRIBUTIONS
% initialize output arrays
% HR test - Left
pVecL=zeros(1,length(Lrow));
tVecL=zeros(1,length(Lrow));
% t-stat map = left
HR_t_L=PGBU_L_x;

% HR test each pixel
for P=1:length(Lrow);
	P
	% get row and column of this pixel
	Row=Lrow(P);
	Col=Lcol(P);
	% extract the angular distribution of this pixel over all frames and segments
	PixelDistrX=real(LeftVFs(Row,Col,:));
	PixelDistrY=imag(LeftVFs(Row,Col,:));
	% convert to radians
	PixelRads=cart2pol(PixelDistrX,PixelDistrY);
	% HR test
	[p,T]=hrtest(PixelRads(:),1000);
	pVecL(P)=p;
	tVecL(P)=T;
	% save out test stat purely for viz
	HR_t_L(Row,Col)=T;
end

% right
pVecR=zeros(1,length(Rrow));
tVecR=zeros(1,length(Rrow));
% t-stat map = left
HR_t_R=PGBU_R_x;
for P=1:length(Rrow);
        % get row and column of this pixel
        Row=Rrow(P);
        Col=Rcol(P);
        % extract the angular distribution of this pixel over all frames and segments
        PixelDistrX=real(RightVFs(Row,Col,:));
        PixelDistrY=imag(RightVFs(Row,Col,:));
        % convert to radians
        PixelRads=cart2pol(PixelDistrX,PixelDistrY);
        % HR test
        [p,T]=hrtest(PixelRads(:),1000);
        pVecR(P)=p;
        tVecR(P)=T;
	% save out test stat purely for viz
	HR_t_R(Row,Col)=T;
end

%%%% FDR
% merge together pvecs for full-brain MC correction
HR_mergedPvec=cat(2,pVecL,pVecR);
FDRed=mafdr(HR_mergedPvec);

% split back into 2
FDRed_L=FDRed(1:length(Lrow));
FDRed_R=FDRed((length(Lrow)+1):(length(Lrow)+length(Rrow)));

% classify passed HR-FDR as 1 - Left
NonRand_L=find(FDRed_L<0.05);
PixClassL(NonRand_L)=1;
% update Lrow and Lcol with only non-randomly distributed pixels
Lrow_NR=Lrow(NonRand_L);
Lcol_NR=Lcol(NonRand_L);

% classify passed HR-FDR as 1 - Right
NonRand_R=find(FDRed_R<0.05);
PixClassR(NonRand_R)=1;
% update Lrow and Lcol with only non-randomly distributed pixels
Rrow_NR=Rrow(NonRand_R);
Rcol_NR=Rcol(NonRand_R);

%%%%%%%%%% END OF HR TESTS

% comment out later, but nice update in the meantime
length(Lrow_NR)
length(Rrow_NR)

%%%%%%%%% HARTIGAN'S DIP TEST FOR NON-UNIMODALITY

% subject remaining to hartigan's dip test - run on angular distance from gradient gradients
% load in in angular distance files
AngDistFn=[childfp sname '_BU_angDist.csv'];
AngDist=readtable(AngDistFn);

% initialize hartigan vectors
hVecL_p=zeros(1,length(Lrow));
hVecL_d=zeros(1,length(Lrow));
% set invalids to NaN to avoid confusion with low p-values
hVecL_p(hVecL_p==0)=NaN;
hVecL_d(hVecL_d==0)=NaN;

for P=1:length(Lrow);
	P
	% if this distribution is non-random
	if FDRed_L(P) < 0.05
		FDRed_L(P)
		% extract angular distance distribution from PG in this pixel
		AngDistaDistr=table2array(AngDist(:,P));
		% Hartigan's dip test - https://snl.salk.edu/~jude/waveform_public/HartigansDipSignifTest.m
		[dip,p_value,xlow,xup]=HartigansDipSigniftest(AngDistaDistr,1000);
		% pack the dip into the vector
		hVecL_p(P)=p_value;
		hVecL_d(P)=dip;
	end
end

% initialize hartigan vectors
hVecR_p=zeros(1,length(Rrow_NR));
hVecR_d=zeros(1,length(Rrow_NR));
% set invalids to NaN to avoid confusion with low p-values
hVecR_p(hVecR_p==0)=NaN;
hVecR_d(hVecR_d==0)=NaN;

for P=1:length(Rrow);
        P
        % if this distribution is non-random
	if FDRed_R(P) < 0.05
		FDRed_R(P)
		% extract angular distance distribution from PG in this pixel
        	AngDistaDistr=table2array(AngDist(:,P+length(Lrow)));
        	% Hartigan's dip test - https://snl.salk.edu/~jude/waveform_public/HartigansDipSignifTest.m
        	[dip,p_value,xlow,xup]=HartigansDipSigniftest(AngDistaDistr,1000);
        	% pack the dip into the vector
        	hVecR_p(P)=p_value;
        	hVecR_d(P)=dip;
	end
end


% merge dip test pvalus to one map for FDR
HAR_mergedPvec=cat(2,hVecL_p,hVecR_p);
HAR_FDRed=mafdr(HAR_mergedPvec);

% split back into 2
HAR_FDRed_L=HAR_FDRed(1:length(Lrow));
HAR_FDRed_R=HAR_FDRed((length(Lrow)+1):(length(Lrow)+length(Rrow)));

% change pixClass from 1 to 2 for non-unimodal 
NonUni_L=find(HAR_FDRed_L<0.05);
PixClassL(NonUni_L)=2;
% and for right hemi
NonUni_R=find(HAR_FDRed_R<0.05);
PixClassL(NonUni_R)=2;
%%%%%%%%%%%%%%%% END HARTIGAN'S

% finally, for unimodal distributions, calculate circular mean for viz
% use downsampled flat PG to reconstruct original ordering of this table

%%%%%%%%%% CIRCULAR MEANS
% intialize circMean vectors
CM_VecL=zeros(1,length(Lrow));
% set invalids to NaN to avoid confusion with low p-values
CM_VecL(CM_VecL==0)=NaN;

% loop over every pixel
for P=1:length(Lrow);
        P
        % if this distribution is non-random
        if FDRed_L(P) < 0.05
		% AND if the distribution is not multi-modal
                if HAR_FDRed_L(P) > 0.05
	               % get row and column of this pixel
		        Row=Lrow(P);
		        Col=Lcol(P);
		        % extract the angular distribution of this pixel over all frames and segments
		        PixelDistrX=real(LeftVFs(Row,Col,:));
		        PixelDistrY=imag(LeftVFs(Row,Col,:));
		        % convert to radians
		        PixelRads=cart2pol(PixelDistrX,PixelDistrY);
			% calculate circular mean
			[mu,ul,ll]=circ_mean(PixelRads(:));
	                % inset circular mean into CM vec
			CM_VecL(P)=mu;
		end
        end
end

% intialize circMean vectors
CM_VecR=zeros(1,length(Rrow));
% set invalids to NaN to avoid confusion with low p-values
CM_VecR(CM_VecR==0)=NaN;

% loop over every right-hemi pixel
for P=1:length(Rrow);
        P
        % if this distribution is non-random
        if FDRed_R(P) < 0.05
                % AND if the distribution is not multi-modal
                if HAR_FDRed_R(P) > 0.05
                       % get row and column of this pixel
                        Row=Rrow(P);
                        Col=Rcol(P);
                        % extract the angular distribution of this pixel over all frames and segments
                        PixelDistrX=real(RightVFs(Row,Col,:));
                        PixelDistrY=imag(RightVFs(Row,Col,:));
                        % convert to radians
                        PixelRads=cart2pol(PixelDistrX,PixelDistrY);
                        % calculate circular mean
                        [mu,ul,ll]=circ_mean(PixelRads(:));
                        % inset circular mean into CM vec
                        CM_VecR(P)=mu;
                end
        end
end


% further classify those failing the dip test on basis of circular means
% change pixClass from 1 to 2 for non-unimodal 
UniCM_L=find(~isnan(CM_VecL));
PixMeanL(UniCM_L)=CM_VecL(UniCM_L);
% and for right hemi
UniCM_L=find(~isnan(CM_VecL));
PixMeanL(UniCM_L)=CM_VecL(UniCM_L);

% ... maybe map back to real space before leaving this maze?

% use gradient gradient map for template
circMeanMapL=PGBU_L_x;
classMapL=PGBU_L_x;
% for left
for i=1:length(Lrow);
	circMeanMapL(Lrow,Lcol)=CM_VecL(i);
	classMapL(Lrow,Lcol)=PixClassL(i);
end

% use gradient gradient map for template
circMeanMapR=PGBU_R_x;
classMapR=PGBU_R_x;
% for left
for i=1:length(Lrow);
        circMeanMapL(Lrow,Lcol)=CM_VecL(i);
        classMapL(Lrow,Lcol)=PixClassL(i);
end

% saveout left hemi
circM_outFN_L=[parentfp sname '_circM_L.mat'];
classM_outFN_L=[parentfp sname '_classM_L.mat'];
save(circM_outFN_L,'circMeanMapL')
save(classM_outFN_L,'classMapL')
HRT_outFN_L=[parentfp sname '_HRT_L.mat'];
save(HRT_outFN_L,'HR_t_L')

% and saveout right hemi
circM_outFN_R=[parentfp sname '_circM_R.mat'];
classM_outFN_R=[parentfp sname '_classM_R.mat'];
save(circM_outFN_R,'circMeanMapR')
save(classM_outFN_R,'classMapR')
HRT_outFN_R=[parentfp sname '_HRT_R.mat'];
save(HRT_outFN_R,'HR_t_R')



