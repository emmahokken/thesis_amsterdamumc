import matplotlib.pyplot as plt 
import nibabel as nib 
import numpy as np 
import os
from scipy import ndimage 


# load in r2* image 
subj = 8
r2_star_path = f'../sharpness/data/{subj:03}_r2star_gt.nii'
r2_star_full_path = f'../../../../../data/projects/ahead/raw_gdata/Subcortex_{subj:04}_006_R02/nii/r2scorr.nii'
r2_star = nib.load(r2_star_path).get_fdata()
r2_star_full = nib.load(r2_star_full_path).get_fdata()
print(r2_star.shape)

# plt.imshow(ndimage.rotate(r2_star_full[:,:,123],90))
# plt.show()

plt.imshow(ndimage.rotate(r2_star_full[80:170,100:170,121],90), cmap='gray')
plt.axis('off')
plt.savefig('ventricle.pdf')