% add needed path
addpath(genpath('/cbica/projects/hcpd/scripts/tools'));
% set input fp for parcel
ProjectFolder = '/cbica/projects/hcpd/data/SSPs';

% load in y7 labels
y7=read_cifti('~/data/Yeo7_Yeo-Buckner-Choi-Raut.dlabel.nii');

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
pg=read_cifti('~/data/hcp.gradients.dscalar.nii');
pgcort=pg.cdata(1:59412);

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
end
