function AngMods(subj)

% designed to operate at the pixel-level within subjects
% 24'ile angles into bins

% add matlab path for used functions
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/hermans-rasson'));

% convert to character vector
sname=char(subj);
childfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/'];
% in this case output dir is same as input dir
parentfp=childfp;
% load in this subjs opFlow
OpFlfn=['/cbica/projects/pinesParcels/results/OpFl_output/' sname '/OpFlowResults.mat'];
OpFlmat=load(OpFlfn);

% merge vector fields across runs (concatenate across time dimension)
LeftVFs=cat(3,OpFlmat.MegaStruct.Vf_Left{:});
RightVFs=cat(3,OpFlmat.MegaStruct.Vf_Right{:});

% use downsampled flat PG for reference
LHlrPG=load([childfp sname '_PG_lowResFlat_L.csv']);
RHlrPG=load([childfp sname '_PG_lowResFlat_R.csv']);
[Lrow,Lcol]=find(~isnan(LHlrPG));
[Rrow,Rcol]=find(~isnan(RHlrPG));

% initialize output map to store mode of each pixel
PixModeL=zeros(1,length(Lrow));
PixModeR=zeros(1,length(Rrow));
% and to save relative prominence of the mode
PixMode_promL=zeros(1,length(Lrow));
PixMode_promR=zeros(1,length(Rrow));

% Bin and count each pixel (left)
for P=1:length(Lrow);
	% get row and column of this pixel
	Row=Lrow(P);
	Col=Lcol(P);
	% extract the angular distribution of this pixel over all frames and segments
	PixelDistrX=real(LeftVFs(Row,Col,:));
	PixelDistrY=imag(LeftVFs(Row,Col,:));
	% convert to radians
	PixelRads=cart2pol(PixelDistrX,PixelDistrY);
	% bin with histc, discritize is giving unequal bin ranges at the edges
	binranges = -pi:.5236:3.14161;
	[bincounts] = histc(PixelRads(:),binranges);
	[maximum,index]=max(bincounts);
	PixModeL(P)=index;
	PixMode_promL(P)=maximum;
end

% Bin and count each pixel (right)
for P=1:length(Rrow);
        % get row and column of this pixel
        Row=Rrow(P);
        Col=Rcol(P);
        % extract the angular distribution of this pixel over all frames and segments
        PixelDistrX=real(RightVFs(Row,Col,:));
        PixelDistrY=imag(RightVFs(Row,Col,:));
        % convert to radians
        PixelRads=cart2pol(PixelDistrX,PixelDistrY);
        % bin with histc, discritize is giving unequal bin ranges at the edges
        binranges = -pi:.5236:3.14161;
        [bincounts] = histc(PixelRads(:),binranges);
        [maximum,index]=max(bincounts);
        PixModeR(P)=index;
	PixMode_promR(P)=maximum;
end

% insert back into map

% same loop to vectoriz big map
modeLmap=LHlrPG;
modePromLmap=LHlrPG;
for P=1:length(Lrow)
	Row=Lrow(P);
	Col=Lcol(P);
	modeLmap(Row,Col)=PixModeL(P);
	modePromLmap(Row,Col)=PixMode_promL(P);
end

% and for right
modeRmap=RHlrPG;
modePromRmap=RHlrPG;
for P=1:length(Rrow)
        Row=Rrow(P);
        Col=Rcol(P);
        modeRmap(Row,Col)=PixModeR(P);
        modePromRmap(Row,Col)=PixMode_promR(P);
end

% saveout left hemi
mode_outFN_L=[parentfp sname '_mode_L.mat'];
modeProm_outFN_L=[parentfp sname '_modeProm_L.mat'];
save(mode_outFN_L,'modeLmap')
save(modeProm_outFN_L,'modePromLmap')

% and saveout right hemi
mode_outFN_R=[parentfp sname '_mode_R.mat'];
modeProm_outFN_R=[parentfp sname '_modeProm_R.mat'];
save(mode_outFN_R,'modeLmap')
save(modeProm_outFN_R,'modePromRmap')


