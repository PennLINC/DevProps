% load OFD for mask
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))

% load subjs list
%%%%%%% MODULAR: WILL NEED RS, RS MATCHED FROM _C, and _C SUBJ LIST
%%%% Age Tertiles as well
Subjs=readtable('~/PWs/rs_cmatch_subs.csv')
%initialize output array
% 36 bins for angular distances
OutDf=zeros(1,36);
% 36 bins for modes
OutDf_modes=zeros(1,36);

% loop over each subj
for s = 1:height(Subjs)
        Subj=table2array(Subjs(s,2));
        fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
        % if file exists, then
        if isfile([fp '/' Subj{:} '_AngDistHist.mat']);
                s
                % load in subj's distr
                Angs=load([fp '/' Subj{:} '_AngDistHist.mat']);
                distHist=Angs.OutDf;
		% add to OutDf
		OutDf=OutDf+distHist;
		% load in modal distributions
		ModD_L=readtable([fp '/' Subj{:} '_AngDistFacial_360ModesL.csv']);
		ModD_R=readtable([fp '/' Subj{:} '_AngDistFacial_360ModesR.csv']);
		% combined
		subjModes=tabulate(table2array(horzcat(ModD_L,ModD_R)));
		% discretize : get all TR pairs from non mw-indices
                OutDf_modes=OutDf_modes+subjModes(:,2);
	end
end

%saveout 
csvwrite('~/PWs/rs_cmatch_subs_AngDistHist360.csv',OutDf)
csvwrite('~/PWs/rs_cmatch_subs_AngDistHist360_modedOverTRs.csv',OutDf_modes)
