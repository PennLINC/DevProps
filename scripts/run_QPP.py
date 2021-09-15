import time
# record start time
start_time = time.time()
import nibabel as nb
import numpy as np
from QPP_cpac import detect_qpp
data=nb.load('/cbica/projects/abcdfnets/dropbox/xxx.dscalar.nii').dataobj
for I in range(100):
	QPPout=detect_qpp(np.transpose(data),1,25,10,.1,I,convergence_iterations=1)
	print("Iteration" + str(I) + " took %s seconds ---" % (time.time() - start_time))
