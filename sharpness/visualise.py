import nibabel as nib 
import matplotlib.pyplot as plt 
import numpy as np

striatum = 'data/segm/sub-mask-str_hem-l_lvlreg-uncorr_def-img.nii'
vent4 = 'data/segm/sub-mask-vent_hem-4_lvlreg-uncorr_def-img.nii'
tha = 'data/segm/sub-mask-tha_hem-r_lvlreg-corr_def-img.nii'
ants = 'data/segm/sub-t1corr_cruise-gwb_ants-def0.nii' 
ants2 = 'data/segm/sub-t1uncorr_cruise-gwb_ants-def0.nii'
r1 = 'data/r1corr.nii'
test_nii = '../../data/recon/test_all_3_ssim/Subcortex_0005_axial_121.nii'
f = nib.load(test_nii).get_fdata()

print(f.shape)

plt.imshow(f[:,:,0], cmap='gray')
plt.show()