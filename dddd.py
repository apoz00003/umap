#!/usr/bin/env python
# This script is used to convert a CT image into a virtual U-map
import nibabel as nib
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ct', help="CT image to register")
parser.add_argument('bone', help="bone mask image to register")

parser.add_argument('umap', help="The output umap")
args = parser.parse_args()

ctt = nib.load(args.ct)
bone = nib.load(args.bone)

bone_data = bone.get_data()
bone_mask = np.array(bone_data)


ct_data = ctt.get_data()


ct= 1000. * np.log(np.array(ct_data))/0.36
nans= np.isnan(ct)
ct[nans] = 0

# CONVERSION from HU to PET u values
# Conversion based on:
# Carney J P et al. (2006) Med. Phys. 33 976-83

a = 0.0000564
b = 0.0408
BP = 1030.

low = ct < BP
ct[low] = 0.000096 * (1000. + ct[low])
high = ct >= BP
ct[high] = a * (1000. + ct[high]) + b

low_u = ct < 0.1134
ct[low_u] =0.1134
# End of conversion

u_soft_fixed = 0.1
u_air = 0.
tissue_mask= 1-(bone_mask)
#umap = 10000. * (u_air *0.+ ct * bone_mask  + tissue_mask * u_soft_fixed)
umap = 10000. * (u_air + ct + tissue_mask )* u_soft_fixed


nans=np.isnan(umap)
umap[nans] = 0
save_im = nib.Nifti1Image(umap, affine=ctt.affine,header=ctt.header)
nib.save(save_im, args.umap)




#!/usr/bin/env python
# This script is used to convert a CT image into a virtual U-map
import nibabel as nib
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ct', help="CT image to register")
parser.add_argument('umap', help="The output umap")
args = parser.parse_args()

ctt = nib.load(args.ct)
ct_data = ctt.get_data()



ct=1000. * np.log(np.array(ct_data))/ 0.36
#ct= 1000. * np.log(np.array(ct_data)) / 2.39
nans = np.isnan(ct)
ct[nans] = 0

#ct= 0.000001351 * (ct_map**3) - 0.003617 * (ct_map**2) + 3.841 * ct_map - 19.46
# CONVERSION from HU to PET u values
# Conversion based on:
# Carney J P et al. (2006) Med. Phys. 33 976-83

a = 0.0000564
b = 0.0408
BP = 1030.

low = ct < BP
ct[low] = 0.000096 * (1000. + ct[low])
high = ct >= BP
ct[high] = a * (1000. + ct[high]) + b

low_u = ct < 0.
ct[low_u] =0.
# End of conversion

u_soft_fixed = 0.1
u_air = 0.

umap = 10000. * (u_air + ct )* u_soft_fixed

nans = np.isnan(umap)
umap[nans] = 0
save_im = nib.Nifti1Image(umap, affine=ctt.affine,header=ctt.header)
nib.save(save_im, args.umap)


