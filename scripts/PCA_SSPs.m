%%%%%%%%%%%%%%%%%%%%%%% Aggregate Subj Loadings
% get subj list
% initialize matrix
% for each subj
% load in mat
% end
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% convert to faces
% convert vertex ts to face ts
sizeTsl=size(TRs_l);
lengthTSl=sizeTsl(2);
sizeTsr=size(TRs_r);
lengthTSr=sizeTsr(2);
% initialize face TS : # of faces and # of TRs
faceTS_l=zeros(length(faces_l),sizeTsl(2));
faceTS_r=zeros(length(faces_r),sizeTsr(2));
% there is probably a more efficient way to do this than loop over every TR
for t=1:lengthTSl
        ts_this_tr_l=TRs_l(:,t);
        faceTS_l(:,t)=sum(ts_this_tr_l(faces_l),2)./3;
        ts_this_tr_r=TRs_r(:,t);
        faceTS_r(:,t)=sum(ts_this_tr_r(faces_r),2)./3;
end
%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Run PCA across subjs
%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%% print out some PCs
% subj PC values
% visualization of PCA loadings
%%%%%%%%%%%%%%%%%%%%%%%%
