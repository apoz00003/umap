#!/usr/bin/env python
# This script is used to convert a CT image into a virtual U-map
import nibabel as nib
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ct', help="CT image to register")

parser.add_argument('u_bone', help="The output umap")
args = parser.parse_args()

ctt = nib.load(args.ct)

ct_data = ctt.get_data()

kt=np.array(ct_data)
nans= np.isnan(kt)
#kt [nans] = 0

# CONVERSION from HU to PET u values
# Conversion based on:
# Carney J P et al. (2006) Med. Phys. 33 976-83

# This is a test change for git tutorial

a = 0.0000564
b = 0.0408
BP = 1040.

low = kt < BP
kt[low] = 0.000096 * (1000. + kt[low])
high = kt>= BP
kt[high] = a * (1000. + kt[high]) + b

low_u =kt<0.0134
kt[low_u]=0.0134
# End of conversion

u_soft_fixed = 0.1
u_air = 0.
kt=10000. * kt


#umap = 1000. * (u_air + u_bone+(1-u_bone)* u_soft_fixed)
#umap = 10000. * (u_air *np.array(air.get_data()) + u_bone + tissu* u_soft_fixed)

#umap = (u_air * np.array(air.get_data()) + u_bone * np.array(bones.get_data()) + (1 - np.array(bones.get_data())) * (1 - np.array(air.get_data())) * u_soft_fixed)
nans=np.isnan(kt)
kt[nans] = 0
save_im = nib.Nifti1Image(kt, affine=ctt.affine,header=ctt.header)
nib.save(save_im, args.u_bone)

