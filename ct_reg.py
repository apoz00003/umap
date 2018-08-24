#!/usr/bin/env python
# This script is used to convert a CT image into a virtual U-map
import nibabel as nib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ct', help="CT image to register")
parser.add_argument('umap', help="The output umap")
args = parser.parse_args()

ct = nib.load(args.ct)
# CONVERSION from HU to PET u values
# Conversion based on:
# Carney J P et al. (2006) Med. Phys. 33 976-83

a = 0.000051
b = 0.0471
BP = 1047.

low = ct < BP
ct[low] = 0.000096 * (1000. + ct[low])
high = ct >= BP
ct[high] = a * (1000. + ct[high]) + b

low_u = ct < 0.1134
ct[low_u] = 0.1134
# End of conversion

u_soft_fixed = 0.1
u_air = 0.

umap = 10000. * (u_air * np.array(air.get_data()) +
                 ct * np.array(bones.get_data()) +
                 (1 - np.array(bones.get_data())) *
                 (1 - np.array(air.get_data())) * u_soft_fixed)

nans = np.isnan(umap)
umap[nans] = 0
save_im = nib.Nifti1Image(umap, affine=ct.affine)
nib.save(save_im, args.umap)
