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

% initialize output arrays
% HR test - Left
pVecL=zeros(1,length(Lrow));
tVecL=zeros(1,length(Lrow));
% t-stat map = left
HR_t_L=PGBU_L_x;

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
FDRed_R=FDRde((length(Lrow)+1):(length(Lrow)+length(Rrow)));

% classify passed HR-FDR as 1 - Left
NonRand_L=find(FDRed_L<0.05);
PixClassL(NonRand_L)=1;

% classify passed HR-FDR as 1 - Right
NonRand_R=find(FDRed_R<0.05);
PixClassR(NonRand_R)=1;

% subject remaining to hartigan's dip test - run on angular distance from gradient gradients
% load in in angular distance files
AngDistFn=[childfp sname '_BU_angDist.csv'];
AngDist=load(AngDistFn);

% use downsampled flat PG to reconstruct original ordering of this table

% hart dip each

% merge dip test pvalus to one map for FDR

% FDR
mafdr

% split back

% further classify those failing the dip test on basis of circular means

% circ mean

% save out circular mean of unimodal pixels : value of Bottomup is high, TD is low

% save out HR t stat maps

% save out classification maps

