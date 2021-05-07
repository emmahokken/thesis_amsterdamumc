import matplotlib.pyplot as plt 
import nibabel as nib 
import numpy as np 
import os
from scipy import ndimage 


# load in r2* image 
subj = 8
r2_star_path = f'../sharpness/data/{subj:03}_r2star_gt.nii'
r2_star_full_path = f'../../../../../data/projects/ahead/raw_gdata/Subcortex_{subj:04}_006_R02/nii/r1corr.nii'
r2_star = nib.load(r2_star_path).get_fdata()
r2_star_full = nib.load(r2_star_full_path).get_fdata()

# plt.imshow(ndimage.rotate(r2_star_full[:,:,123],90))
# plt.show()
s = 124
fig, ax = plt.subplots(figsize=(15,15))
ax.imshow(ndimage.rotate(r2_star_full[80:170,100:170,s],90), cmap='gray')
ax.axis('off')
plt.savefig('ventricle.pdf')
plt.show()

plt.imshow(ndimage.rotate(r2_star_full[:,:,s],90), cmap='gray')
plt.axis('off')
plt.show()
plt.savefig('whole_brain.pdf')