# from local
import numpy as np
from scipy import special as special
from scipy.optimize import minimize as sp_min
import matplotlib.pyplot as pl
from scipy.io import loadmat

######

# Subj Data
testSubj_L=loadmat('C:/Users/adam/testSubjMat_L.mat', squeeze_me=True)
testSubj_R=loadmat('C:/Users/adam/testSubjMat_R.mat', squeeze_me=True)

# COL 0 = PROP. OF TRs THAT ARE BU
# COL 1 = BU VECTOR HORIZONTAL COMPONENT
# COL 2 = BU VECTOR VERTICAL COMPONENT
# COL 3 = TD VECTOR HORIZONTAL COMPONENT
# COL 4 = TD VECTOR VERTICAL COMPONENT
# COL 5 = BU VECTOR LENGTH
# COL 6 = TD VECTOR LENGTH

LeftHemi=testSubj_L['OutDf_L']
RightHemi=testSubj_R['OutDf_R']

# plot out normative BU and TD vector
#pl.scatter([LeftHemi[400,1],LeftHemi[400,3]],[LeftHemi[400,2],LeftHemi[400,4]])
#pl.xlim([-1,1])
#pl.ylim([-1,1])
#pl.show()

# scatter plot of resultant length of vectors: seems to indicate strong TD exitsts where strong TD exists
pl.scatter(LeftHemi[1:1500,3],LeftHemi[1:1500,4])
pl.show()

np.corrcoef(LeftHemi[1:1500,3].astype(float),LeftHemi[1:1500,4].astype(float))
# Normative PGG directions

# calculate distance in radians so as to produce 2pi of variation

# convert input data to radians

# will need mask indiceshttps://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.vonmises.html as well

#####

