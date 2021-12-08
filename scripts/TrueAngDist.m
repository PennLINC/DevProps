function TrueAngDist(subj)

%%% pull out PG grad, OpFl distribution, calculate angular distance per pixel

% char conversion and input file locations
sname=char(subj);
childfp=['/cbica/projects/pinesParcels/results/PWs/Proced/' sname '/'];

% extract PG gradients 
PGBU_L_x=load([childfp sname '_PG_GxL_BU.csv']);
PGBU_L_y=load([childfp sname '_PG_GyL_BU.csv']);
PGBU_R_x=load([childfp sname '_PG_GxR_BU.csv']);
PGBU_R_y=load([childfp sname '_PG_GyR_BU.csv']);
PGTD_L_x=load([childfp sname '_PG_GxL_TD.csv']);
PGTD_L_y=load([childfp sname '_PG_GyL_TD.csv']);
PGTD_R_x=load([childfp sname '_PG_GxR_TD.csv']);
PGTD_R_y=load([childfp sname '_PG_GyR_TD.csv']);

% extract OpFl results
OpFlResfn=['/cbica/projects/pinesParcels/results/OpFl_output/' sname '/OpFlowResults.mat'];
OpFlRes=load(OpFlResfn);

% pull in mask
%FlatMask_L=load([childfp sname '_FlatMask_L.csv']);
%FlatMask_R=load([childfp sname '_FlatMask_R.csv']);

% merge vector fields across runs (concatenate across time dimension)
LeftVFs=cat(3,OpFlRes.MegaStruct.Vf_Left{:});
RightVFs=cat(3,OpFlRes.MegaStruct.Vf_Right{:});

% extract available number of frames
sizeOfVfs=size(LeftVFs);
NumFrames=sizeOfVfs(3);

% get coordinates of each viable pixel
[Lrow,Lcol]=find(~isnan(PGBU_L_x));
[Rrow,Rcol]=find(~isnan(PGBU_R_x));

% initialize TD,BU,and angle-doubled outArray
BU_angDist=zeros(NumFrames,(length(Lrow)+length(Rrow)));
TD_angDist=zeros(NumFrames,(length(Lrow)+length(Rrow)));
angDoubDist=zeros(NumFrames,(length(Lrow)+length(Rrow)));

% for each Left pixel
for P = 1:length(Lrow)
	% get coordinates
	Row=Lrow(P);
	Col=Lcol(P);
	% PG Vectors
	PGVecBU=[PGBU_L_x(Row,Col) PGBU_L_y(Row,Col)];
	PGVecTD=[PGTD_L_x(Row,Col) PGTD_L_y(Row,Col)];
	% for each vector field (each TR)
	for V = 1:NumFrames
		% Vf Vectors
		XVec=real(LeftVFs(Row,Col,V));
		YVec=imag(LeftVFs(Row,Col,V));
		VFVec=[XVec YVec];
		% calculate distance between angles: 
		% mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
		BUCosTheta = max(min(dot(PGVecBU,VFVec)/(norm(PGVecBU)*norm(VFVec)),1),-1);
		BUThetaInDegrees = real(acosd(BUCosTheta));
		TDCosTheta = max(min(dot(PGVecTD,VFVec)/(norm(PGVecTD)*norm(VFVec)),1),-1);
                TDThetaInDegrees = real(acosd(TDCosTheta));
		%% convert vectors to angles for doubling, measure relative to positive y axis
	%	origin=[0 1];
		%mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
	%	Org_BUDegrees = atan2d(PGVecBU(1)*origin(2)-PGVecBU(2)*origin(1),PGVecBU(1)*origin(1)+PGVecBU(2)*origin(2));
	%	Org_VFDegrees = atan2d(VFVec(1)*origin(2)-VFVec(2)*origin(1),VFVec(1)*origin(1)+VFVec(2)*origin(2));	
	%	Org_BUDegrees=Org_BUDegrees+180;
	%	Org_VFDegrees=Org_VFDegrees+180;
		% ang-doub: webspace.ship.edu/pgmarr/Geo441/Lectures/Lec%2016%20-%20Directional%20Statistics.pdf
	%	DoubPGAng=2*Org_BUDegrees;
	%	if DoubPGAng > 360;
	%		DoubPGAng = DoubPGAng - 360;
	%	end
	%	DoubVFAng=2*Org_VFDegrees;
	%	if DoubVFAng > 360;
	%		DoubVFAng =  DoubVFAng - 360;
	%	end
	%	% Convert to radians for conversion to vectors
	%	DoubPGAngRad=deg2rad(DoubPGAng);
	%	DoubVFAngRad=deg2rad(DoubVFAng);
	%	% convert to vectors
	%	[xDoubPGVec,yDoubPGVec]=pol2cart(DoubPGAngRad,1);
	%	[xDoubVFVec,yDoubVFVec]=pol2cart(DoubVFAngRad,1);
	%	DoubPGVec=[xDoubPGVec yDoubPGVec];
	%	DoubVFVec=[xDoubVFVec yDoubVFVec];
	%	% once again find difference
	%	DoubCosTheta = max(min(dot(DoubPGVec,DoubVFVec)/(norm(DoubPGVec)*norm(DoubVFVec)),1),-1);
        %        DoubThetaInDegrees = real(acosd(DoubCosTheta));	
		% throw em in the ang Dist vectors. To be averaged over TRs
		BU_angDist(V,P)=BUThetaInDegrees;
		TD_angDist(V,P)=TDThetaInDegrees;
		%angDoubDist(V)=DoubThetaInDegrees;	
	end
% end for each pixel
end

% for each Left pixel
for P = 1:length(Rrow)
        % get coordinates
        Row=Rrow(P);
        Col=Rcol(P);
        % PG Vectors
        PGVecBU=[PGBU_R_x(Row,Col) PGBU_R_y(Row,Col)];
        PGVecTD=[PGTD_R_x(Row,Col) PGTD_R_y(Row,Col)];
        % for each vector field (each TR)
        for V = 1:NumFrames
                % Vf Vectors
                XVec=real(RightVFs(Row,Col,V));
                YVec=imag(RightVFs(Row,Col,V));
                VFVec=[XVec YVec];
                % calculate distance between angles: 
                % mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
                BUCosTheta = max(min(dot(PGVecBU,VFVec)/(norm(PGVecBU)*norm(VFVec)),1),-1);
                BUThetaInDegrees = real(acosd(BUCosTheta));
                TDCosTheta = max(min(dot(PGVecTD,VFVec)/(norm(PGVecTD)*norm(VFVec)),1),-1);
                TDThetaInDegrees = real(acosd(TDCosTheta));
                % throw em in the ang Dist vectors. To be averaged over TRs
                BU_angDist(V,P+length(Lrow))=BUThetaInDegrees;
                TD_angDist(V,P+length(Lrow))=TDThetaInDegrees;
        end
end

% save out files
writetable(table(BU_angDist),[childfp sname '_BU_angDist.csv']);
writetable(table(TD_angDist),[childfp sname '_TD_angDist.csv']);
