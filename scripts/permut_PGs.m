function permut_PGs(subj)
% create the permutation (null) maps of hierarchy directionality

% spin the gradient 1000 times
spin_pg(subj)

% load in spun values

% initialize 4d array
SpunPG_Gs=zeros(1000,64,89,2);

% For each spin, project to high-res, get gradient gradient, insert into 4d array, delete high-res map. 
for s=1:1000
% convert data vector to brainmap (fsaverage5 space)
% project to high-res
% derive high rest gradient gradient
% resample to flatmap
% insert resample into 4d array
SpunPG_Gs(s,:,:,1)=x grad;
SpunPG_Gs(s,:,:,2)=y grad;
end

% save 4d array
