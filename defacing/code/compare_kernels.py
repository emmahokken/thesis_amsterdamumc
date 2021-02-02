import cv2 
import matplotlib.patches as patches 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres
import numpy as np 
import os 
from scipy import ndimage

from defacing import deface

# decalre root path and scan type
root_path = '../../../../../data/projects/ahead/raw_gdata'
scan_type = 'inv2_te1_m_corr'
brain_mask_type = 'mask_inv2_te2_m_corr'
save_path = '../defacing/dilation'
subject = 'Subcortex_0004_002_R02'

kernel_five = np.ones((5,5), dtype=np.uint8)
kernel_twenty = np.ones((61,61), dtype=np.uint8)

brain_mask_path = f'{root_path}/{subject}/nii/{brain_mask_type}.nii'
scan_path = f'{root_path}/{subject}/nii/{scan_type}.nii'

# load in brain mask (in two steps for saving later)
brain_mask_nii = nib.load(brain_mask_path)
brain_mask = brain_mask_nii.get_fdata()
inv2_nii = nib.load(scan_path)
inv2 = inv2_nii.get_fdata()

# dilate brain mask and fill in with (second inversion) scan 
dilation_five = cv2.dilate(brain_mask, kernel_five, iterations=15)
filled_in_five = dilation_five*inv2

# dilate with 20x20 kernel
dilation_twenty = cv2.dilate(brain_mask, kernel_twenty, iterations=1)
filled_in_twenty = dilation_twenty*inv2
difference = dilation_twenty - dilation_five

plt.imshow(ndimage.rotate(difference[:,:,63], 90), cmap='gray')
plt.show()

fig, ax = plt.subplots(1,2)
ax[0].imshow(ndimage.rotate(dilation_five[:,:,63], 90), cmap='gray')
ax[0].axis('off')
ax[0].set_title('5x5 kernel, 15 iterations')
ax[1].imshow(ndimage.rotate(dilation_twenty[:,:,63], 90), cmap='gray')
ax[1].axis('off')
ax[1].set_title('50x50 kernel, 1 iterations')
plt.show()

exit()

# save new scan 
if not os.path.exists(f'{save_path}/{d}'):
    os.makedirs(f'{save_path}/{d}')
save_dir = f'{save_path}/{d}/dilated_{scan_type}.nii'

nifit_image = nib.Nifti1Image(dataobj=filled_in, header=inv2_nii.header, affine=inv2_nii.affine)
nib.save(img=nifit_image, filename=save_dir)
