%add needed paths
addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));

% add surfaces (pial for equiv distances)
surfL=read_surf('/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/surf/lh.pial');
surfR=read_surf('/cbica/software/external/freesurfer/centos7/6.0.0/subjects/fsaverage4/surf/rh.pial');

% big diStance matrix
bdsml=zeros(length(surfL),length(surfL));
bdsmr=zeros(length(surfR),length(surfR));

% for all vertices
for V=1:length(bdsml);
		% vertex props (x,y,z coords)
		initVertL=surfL(V,:);
		xVL=initVertL(1);
		yVL=initVertL(2);
		zVL=initVertL(3);
                initVertR=surfR(V,:);
                xVR=initVertR(1);
                yVR=initVertR(2);
                zVR=initVertR(3);
        
		% search through all of them for eucl. dist. calc.
		for i=1:length(bdsml);
			xL=surfL(i,1);
			yL=surfL(i,2);
			zL=surfL(i,3);
			eucld_L=sqrt((xL-xVL)^2+(yL-yVL)^2+(zL-zVL)^2);
			xR=surfR(i,1);
                        yR=surfR(i,2);
                        zR=surfR(i,3);
                        eucld_R=sqrt((xR-xVR)^2+(yR-yVR)^2+(zR-zVR)^2);
			bdsml(V,i)=eucld_L;
			bdsmr(V,i)=eucld_R;
		end	
end
save('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_left_fsaverage4.mat','bdsml');
save('/cbica/projects/pinesParcels/data/aggregated_data/euclidean_distance_right_fsaverage4.mat','bdsmr');
