import hcp_utils as hcp
import numpy as np
import nibabel as nb
import sys

# Subject is set to the passed argument
subj = sys.argv[1]

# load in subject parcel as template cifti
parcelLoc='/cbica/projects/abcdfnets/results/SingleParcel_1by1/' + str(subj) + '/' + str(subj) + '_Parcel.dscalar.nii'
parcel=nb.load(parcelLoc)
parcelCort=parcel.dataobj[:,hcp.struct.cortex]

# child (output) filepath
childfp='/cbica/projects/abcdfnets/results/wave_output/' + str(subj) + '/'

# load in gradient
subjPGfn=childfp+str(subj)+'_PG1.npy'
grads=np.load(subjPGfn)

# jump through cifti hoops to make a new dscalar to viz
newDat=parcel.get_fdata()
newDat[:,hcp.struct.cortex]=grads
new_img=nb.Cifti2Image(newDat,header=parcel.header,nifti_header=parcel.nifti_header)

# output filename
subjPG_viz_fn=childfp+str(subj)+'_PG1.dscalar.nii'
new_img.to_filename(subjPG_viz_fn)
