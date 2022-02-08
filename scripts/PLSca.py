# import necces. libs
import os
import scipy.io as sio
import numpy as np
import time
import random
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

# 388x150?
data = sio.loadmat(SSP_Mat_Path)
SSP_Matrix_str = data['SSP_PC_struct']
data = sio.loadmat(Vec_Mat_Path)
Vec_Matrix_str = data['Vecs_PC_struct']

# extract subject-wise PC scores for features
SSP_Matrix=SSP_Matrix_str[0,0]['score']
Vec_Matrix=Vec_Matrix_str[0,0]['score']

# truncate: 150 so training subj # is half the number of features on either side
SSP_Matrix=SSP_Matrix[:,0:100]
Vec_Matrix=Vec_Matrix[:,0:100]

# start with conditions in ZC paper
Fold_Quantity=2
Components_Number=30
CVRepeatTimes=101

# Get random splits of subjects for folds, save
Subjects_Quantity = np.shape(SSP_Matrix)[0]
EachFold_Size = np.int(np.fix(np.divide(Subjects_Quantity, Fold_Quantity)))
Remain = np.mod(Subjects_Quantity, Fold_Quantity)
RandIndex = np.arange(Subjects_Quantity)
np.random.shuffle(RandIndex)
RandIndex_Mat = {'RandIndex': RandIndex}
sio.savemat('/cbica/projects/pinesParcels/results/PWs/RandIndex.mat', RandIndex_Mat)

# save number of features for future reference
SFeatures_Quantity = np.shape(SSP_Matrix)[1]
VFeatures_Quantity = np.shape(Vec_Matrix)[1]
# initialize fold_corr
Fold_Corr = np.zeros((Fold_Quantity, Components_Number))

# for each fold
for j in np.arange(Fold_Quantity):
	# get indices corresponding to this fold: note python iterator starts w/ 0
	Fold_J_Index = RandIndex[EachFold_Size * j + np.arange(EachFold_Size)]
	# extract ssp and vec train and test data
	SSP_test = SSP_Matrix[Fold_J_Index, :]
	SSP_train = np.delete(SSP_Matrix, Fold_J_Index, axis=0)	
	Vec_test = Vec_Matrix[Fold_J_Index, :]
	Vec_train = np.delete(Vec_Matrix, Fold_J_Index, axis=0)
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
	# Correlation on training data
	Fold_J_Corr_Training = []
	for k in np.arange(Components_Number):
		Fold_J_Corr_tmp = np.corrcoef(plsca.x_scores_[:,k], plsca.y_scores_[:,k])
		Fold_J_Corr_tmp = Fold_J_Corr_tmp[0,1]
		Fold_J_Corr_Training.append(Fold_J_Corr_tmp)
	# Correlation on testing data
	Fold_J_Corr = []
	for k in np.arange(Components_Number):
		Fold_J_Corr_tmp = np.corrcoef(SSP_test_ca[:,k], Vec_test_ca[:,k])
		Fold_J_Corr_tmp = Fold_J_Corr_tmp[0,1]
		Fold_J_Corr.append(Fold_J_Corr_tmp)	

	# save out weights
	SSP_Weights = plsca.x_weights_
	Vec_Weights = plsca.y_weights_
	Fold_J_result = {'Index': Fold_J_Index, 'Fold_J_Corr': Fold_J_Corr, \
		'SSP_test_ca': SSP_test_ca, 'Vec_test_ca': Vec_test_ca, \
		'SSP_Weights': SSP_Weights, 'Vec_Weights': Vec_Weights, \
		'Fold_J_Corr_Training': Fold_J_Corr_Training }
	Fold_J_FileName = 'Fold_' + str(j) + '_Score.mat'
	ResultantFile = os.path.join('/cbica/projects/pinesParcels/results/PWs/', Fold_J_FileName)
	sio.savemat(ResultantFile, Fold_J_result)
	Fold_Corr[j,:] = Fold_J_Corr
	# record mean correlation
	Mean_Corr = np.mean(Fold_Corr, axis = 0)
	Res_NFold = {'Mean_Corr':Mean_Corr}
	ResultantFile = os.path.join('/cbica/projects/pinesParcels/results/PWs/', 'Res_NFold.mat')
	sio.savemat(ResultantFile, Res_NFold)
	print(Mean_Corr)

