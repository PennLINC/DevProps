% aim - to simulate synthetic data with wave propagation
% author - saggar@stanford.edu
% date - Aug 23, 2022

addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'))
addpath(genpath('~/dropbox/cifti-matlab/'));

% vertices for the brain surface
g_L = gifti('~/data/Surfs/S900.L.inflated_MSMAll.32k_fs_LR.surf.gii');
g_R = gifti('~/data/Surfs/S900.R.inflated_MSMAll.32k_fs_LR.surf.gii');

% dummy data to begin with
cifti = cifti_read('~/data/exemplar_dtseries.dtseries.nii'); % just need this file for structure, take any HCP file 32k cifti

%% create wavepropagation only
wave_dir = 2; %x=1, y=2, z=3;
wave_width = 4000
wave_bin = normpdf(linspace(-wave_width,wave_width,2*wave_width+1),0,wave_width/3);
wave_amp = 10000;
wave_step = round(length(wave_bin)/10);% reduce wave_step for longer duration files

cifti2 = cifti;
cifti2.cdata =zeros(91282,1);

% generating wave data, ignoring non-cortical indices
g_L.vertices(setdiff([1:size(g_L.vertices,1)], cifti2.diminfo{1}.models{1}.vertlist+1),:) = [];
g_R.vertices(setdiff([1:size(g_R.vertices,1)], cifti2.diminfo{1}.models{2}.vertlist+1),:) = [];

[~, indx_L] = sort(g_L.vertices(:,wave_dir));
[~, indx_R] = sort(g_R.vertices(:,wave_dir));


k = 1;
lh_vertcount = cifti2.diminfo{1}.models{1}.count;
lh_vertlist = [1:lh_vertcount]';

rh_vertcount = cifti2.diminfo{1}.models{2}.count;
rh_vertlist = [1+lh_vertcount:lh_vertcount+rh_vertcount]';

for i = 1:wave_step:length(indx_L)-length(wave_bin)
 
    wave_vert_L = indx_L(i:i+length(wave_bin)-1);
  
    cifti2.cdata(1:lh_vertcount, k) = randn([1, lh_vertcount])';
    cifti2.cdata(lh_vertlist(wave_vert_L), k) = wave_bin*wave_amp;
    

    k = k+1
    
end
k = 1;
for i = 1:wave_step:length(indx_R)-length(wave_bin)

  wave_vert_R = indx_R(i:i+length(wave_bin)-1);
  cifti2.cdata(1+lh_vertcount:rh_vertcount+lh_vertcount, k) = randn([1, rh_vertcount])';
  cifti2.cdata(rh_vertlist(wave_vert_R), k) = wave_bin*wave_amp;
  k = k + 1;
end
    
cifti2.diminfo{2}.length = k-1;
cifti_write(cifti2, '~/data/ciftiout2.dtseries.nii');

%% adding realistic power spectra, adapted code from Tim Laumann's method to simulated power spectra from real data
% see Laumann et al. 2017 cer cortex paper
TRin = 1;  % randomly picked
TRout = 1;  % randomly picked

% run pwelch loop over regions
cifti_real = cifti_read('~/data/exemplar_dtseries.dtseries.nii');

parcelPS = [];
for reg = 1:1:91282
    reg
    [parcelPS(reg,:),freq] = pwelch(cifti_real.cdata(reg,:),[],[],[],1/TRin);
end

ps_mean = mean(parcelPS);
P_target = ps_mean;
timelength = cifti2.diminfo{2}.length;

wn = cifti2.cdata';

%%Normalize to unit variance
wn_std = std(wn);
wn_mean = mean(wn);
wn_unit = (wn-repmat(wn_mean,timelength,1))./repmat(wn_std,timelength,1);


%%Power spectral density resampling
if mod(timelength,2)
    midpoint = ceil(timelength./2);
else
    midpoint = timelength./2+1;
end

% Find frequency indices to interpolate
if mod(length(P_target),2)
    P_target_length = ceil(length(P_target)/2)+1;
else
    P_target_length = length(P_target)/2+1;
end

count = 1;
for k = 0:timelength/2;
    Vq(count) = (k*TRin*length(P_target))/(TRout*(timelength))+1;
    count = count + 1;
end

% Interpolate 
resamp_P_target = interp1(P_target(1:P_target_length),Vq);

% Flip over nyquist
resamp_P_target = [resamp_P_target fliplr(resamp_P_target(2:(midpoint-1)))]';

if mod(timelength,2) == 1
    resamp_P_target = [0; resamp_P_target];
end

%%Enforce power spectra
fft_timecourse_sim = fft(wn_unit);
resamp_P_target_rep = repmat(resamp_P_target,[1 size(fft_timecourse_sim,2)]);
F_sim = fft_timecourse_sim.*sqrt(resamp_P_target_rep);

timecourse_sim = ifft(F_sim,'symmetric');
timecourse_sim = (timecourse_sim-repmat(mean(timecourse_sim),size(timecourse_sim,1),1))./repmat(std(timecourse_sim),size(timecourse_sim,1),1);
cifti2.cdata = timecourse_sim';
cifti_write(cifti2,'~/data/ciftiout_AntPost_Sym.dtseries.nii');
