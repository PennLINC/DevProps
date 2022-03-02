# import necces. libs
import os
import scipy.io as sio
import numpy as np
import time
import random
from sklearn import model_selection
from sklearn import linear_model
from sklearn import preprocessing
from joblib import Parallel, delayed
from sklearn.cross_decomposition import PLSCanonical
from sklearn.experimental import enable_iterative_imputer  
from sklearn.impute import IterativeImputer
import statsmodels.formula.api as sm
from scipy.io import loadmat

# thank you ZC for the template!

# load in PCA values
Vec_Mat_Path='/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjVecsPCA.mat'
SSP_Mat_Path='/cbica/projects/pinesParcels/results/PWs/FaceSpace_SubjTopogPCA.mat'
data = sio.loadmat(SSP_Mat_Path)
SSP_Matrix_str = data['SSP_PC_struct']
data = sio.loadmat(Vec_Mat_Path)
Vec_Matrix_str = data['Vecs_PC_struct']

# extract subject-wise PC scores for features
SSP_Matrix=SSP_Matrix_str[0,0]['score']
Vec_Matrix=Vec_Matrix_str[0,0]['score']

# Optimize number of PCs to use 
numPCs=np.arange(30,100)

# start with conditions in ZC paper
Fold_Quantity=2
Components_Number=30
RepeatTimes=101

# Get random splits of subjects for folds, save
Subjects_Quantity = np.shape(SSP_Matrix)[0]
EachFold_Size = np.int(np.fix(np.divide(Subjects_Quantity, Fold_Quantity)))
Remain = np.mod(Subjects_Quantity, Fold_Quantity)
RandIndex = np.arange(Subjects_Quantity)
RandIndex_Mat = {'RandIndex': RandIndex}
sio.savemat('/cbica/projects/pinesParcels/results/PWs/RandIndex.mat', RandIndex_Mat)

# save number of features for future reference
SFeatures_Quantity = np.shape(SSP_Matrix)[1]
VFeatures_Quantity = np.shape(Vec_Matrix)[1]
# initialize mean out-of-sample correlations
MeanCorrs = np.zeros((2, RepeatTimes))
# initialize a "selected PCs" vector
SelecPCsNum= np.zeros((RepeatTimes))

# initialize feature loading vectors
# actually do the thing you said you would do in the comment, adam

# for each repeat
for i in np.arange(RepeatTimes):
	np.random.shuffle(RandIndex)
	# for each repeat: split into 3rds
	first2_3rds=RandIndex[0:258]
	last3rd=RandIndex[259:388]
	# initialize to select optimal PC number in inner loop 
	numPCsCor=np.zeros((10,100,Components_Number))
	# for each possible number of PCs
	for p in numPCs:
		SSP_Matrix_smol=SSP_Matrix[:,0:p]
		Vec_Matrix_smol=Vec_Matrix[:,0:p]
		# get your 10 folds, right here
		Kf=model_selection.KFold(n_splits=10)
		iterator=0
		# 10 times: loop over possible numPCs in 1st 2 3rds: 1st is to get coefs for p PCs, 2nd to test p PCs
		for train_index, test_index in Kf.split(first2_3rds):
			innertrain_index=first2_3rds[train_index]
			innertest_index=first2_3rds[test_index]
			SSP_test = SSP_Matrix_smol[innertest_index, :]
			SSP_train = SSP_Matrix_smol[innertrain_index, :]
			Vec_test = Vec_Matrix_smol[innertest_index, :]
			Vec_train = Vec_Matrix_smol[innertrain_index, :]
			# normalize all features
			normalize = preprocessing.StandardScaler()
			SSP_train = normalize.fit_transform(SSP_train)
			SSP_test = normalize.fit_transform(SSP_test)
			Vec_train = normalize.fit_transform(Vec_train)
			Vec_test = normalize.fit_transform(Vec_test)
			# create plsca object
			plsca = PLSCanonical(n_components=Components_Number, algorithm='svd')
			# fit it (train it)
			plsca.fit(SSP_train, Vec_train)
			SSP_test_ca, Vec_test_ca = plsca.transform(SSP_test, Vec_test)
			# Correlation on testing data
			for k in np.arange(Components_Number):
				Fold_J_Corr_tmp = np.corrcoef(SSP_test_ca[:,k], Vec_test_ca[:,k])
				Fold_J_Corr_tmp = Fold_J_Corr_tmp[0,1]
				numPCsCor[iterator,p,k]=Fold_J_Corr_tmp	
			iterator=iterator+1
	# select number of PCs corresponding to max in-fold pred.
	comp1Pred=numPCsCor[:,:,0]
	OptimalPCnum=np.where(np.mean(comp1Pred,axis=0)==np.amax(np.mean(comp1Pred,axis=0)))
	# python feeling like matlab rn
	OptimalPCnum=np.int(np.round(np.mean(OptimalPCnum[0])))
	print(OptimalPCnum)
	SelecPCsNum[i]=OptimalPCnum
	# now train model coeffs on this numba
	SSP_test = SSP_Matrix[last3rd,0:OptimalPCnum]
	SSP_train = SSP_Matrix[first2_3rds,0:OptimalPCnum]
	Vec_test = Vec_Matrix[last3rd,0:OptimalPCnum]
	Vec_train = Vec_Matrix[first2_3rds,0:OptimalPCnum]
	plsca = PLSCanonical(n_components=Components_Number, algorithm='svd')
	plsca.fit(SSP_train, Vec_train)
	SSP_test_ca, Vec_test_ca = plsca.transform(SSP_test, Vec_test)
	c1cor=np.corrcoef(SSP_test_ca[:,0], Vec_test_ca[:,0])
	MeanCorrs[0,i]=c1cor[0,1]
	c2cor=np.corrcoef(SSP_test_ca[:,1], Vec_test_ca[:,1])
	MeanCorrs[1,i]=c2cor[0,1]
print('Mean out-of-samp cor: plsc1:',np.mean(MeanCorrs[0,:]))
print('Mean out-of-samp cor: plsc2:',np.mean(MeanCorrs[1,:]))
# save out performance
sio.savemat('/cbica/projects/pinesParcels/results/PWs/PLSCorrelations.mat',{'MeanCorrs':MeanCorrs})
# extract median PCs retained
medPCs=np.int(np.round(np.median(SelecPCsNum)))
print('Median PCs retained:',medPCs)
# fit on whole sample for feat weights
normalize = preprocessing.StandardScaler()
SSP_norm = normalize.fit_transform(SSP_Matrix[:,0:medPCs])
Vec_norm = normalize.fit_transform(Vec_Matrix[:,0:medPCs])
plsca = PLSCanonical(n_components=Components_Number, algorithm='svd')
plsca.fit(SSP_norm, Vec_norm)
# save out weights
SSP_Weights = plsca.x_weights_
Vec_Weights = plsca.y_weights_
sio.savemat('/cbica/projects/pinesParcels/results/PWs/SSP_Weights.mat',{'SSP_Weights':SSP_Weights})
sio.savemat('/cbica/projects/pinesParcels/results/PWs/Vec_Weights.mat',{'Vec_Weights':Vec_Weights})
