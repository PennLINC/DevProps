% add needed path
addpath(genpath('~/dropbox/cifti-matlab/'));
% set input fp for parcel
ProjectFolder = '/cbica/projects/hcpd/data/SSPs';

% load in y7 labels
y7=cifti_read('~/data/Yeo7_Yeo-Buckner-Choi-Raut.dlabel.nii');

% cortex only
y7lab=y7.cdata(1:59412);

% load in hard parcel
initName = [ProjectFolder '/SingleParcellation/fslrSurf_Initialization/init.mat'];
SP=load(initName);
% convert to hard parcel
[ ~ , hardParcel]=max(SP.initV,[],2);

% network name vector
networkName = {'visual', 'motor', 'dan', 'van', 'limbic', 'fp', 'dmn'};
% get correspondence: save to df if wanted later

%Correspondence_Yeo7Systems = zeros(1, 7);

% transmodality data
pg=cifti_read('~/data/hcp.gradients.dscalar.nii');
pgcort=pg.cdata(1:59412,1);

% network gradient loadings initialization
netGradLoads=zeros(1,17);
% vertices plotted as respective network's average hierarchy loading bector
netHLoads=zeros(59412,1);
for i = 1:17
    Index = find(hardParcel == i);
    % Yeo 7 systems
    Label_Yeo_System7 = y7lab(Index);
    % get breakdown of y7 labels in this network
    bd=tabulate(Label_Yeo_System7);
    % get maximum constituent network
    [maximum,indmax]=max(bd(:,3));
    % print out overlap network
    i
    maximum
    networkName(indmax)
    % print out transmodality
    mean(pgcort(Index))
    % save to vector
    netGradLoads(i)=mean(pgcort(Index));
    % save vertex-wise values out, average over this network
    netHLoads(Index)=mean(pgcort(Index));
end

% write over original PG data for output file (to use cifti template)
pg.cdata(1:59412,1)=netHLoads;
outputfile=['~/results/PWs_GParcel_HLoadings.dscalar.nii'];
cifti_write(pg,outputfile)

