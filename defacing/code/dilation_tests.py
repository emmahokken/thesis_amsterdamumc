import cv2 
import matplotlib.patches as patches 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres
import numpy as np 
import os 
from scipy import ndimage

from defacing import deface

# load face files 
inv1_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/inv1_te1_m_corr.nii'
inv2_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/inv2_te1_m_corr.nii'
t1corr_path = '../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/t1corr.nii'
face_allm_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/Subcortex_0004_002_R02_allm.nii'
inv2_n4itk_path = '../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/n4itk_inv2_te1_m_corr.nii'

t1corr_removed_background_path = '../defacing/t1corr_removed_background.nii'
t1corr_scaled_path = '../defacing/t1corr_scaled.nii'
t1corr_scaled_removed_background_path = '../defacing/t1corr_scaled_removed_background.nii'
t1corr_n4itk_removed_background_path = '../defacing/t1corr_n4itk_removed_background.nii'
t1corr_n4itk_removed_background_dilated_path = '../defacing/t1corr_n4itk_removed_background_dilated.nii'

no_face_inv1 = f"../defacing/{inv1_path.split('/')[-1]}"
no_face_inv2 = f"../defacing/{inv2_path.split('/')[-1]}"

brain_mask_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/mask_inv2_te2_m_corr.nii'
test_path = '../defacing/001.nii'
t1_mask_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/minv2_te2_m_corr.nii'
c1_te2_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/inv1_te1_ph_corr.nii'
t1_mask = nib.load(t1_mask_path).get_fdata()

#  load in relevant files
brain_mask_nii = nib.load(brain_mask_path)
brain_mask = brain_mask_nii.get_fdata()
inv1 = nib.load(inv1_path).get_fdata()
inv2_d = nib.load(inv2_path)
inv2 = inv2_d.get_fdata()
t1corr = nib.load(t1corr_path)
t1corr_data = t1corr.get_fdata()
c1_te2 = nib.load(c1_te2_path).get_fdata()
t1corr_n4itk = nib.load(inv2_n4itk_path).get_fdata()

# declare kernel as a 5x5 matrix 
kernel = np.ones((5,5), dtype=np.uint8)


# ''' FILLED IN BRAIN MASK '''


# fig, ax = plt.subplots(1,5)
# ax[0].imshow(ndimage.rotate(brain_mask[:,:,64], 90), cmap='gray')
# ax[0].axis('off')
# ax[0].set_title('Original brain mask')
# for i, iteration in enumerate([5, 10, 15, 20]):
#     # dilate brain mask
#     dilation = cv2.dilate(brain_mask, kernel, iterations=iteration)

#     # plot
#     ax[i+1].imshow(ndimage.rotate(dilation[:,:,64], 90), cmap='gray')
#     ax[i+1].axis('off')
#     ax[i+1].set_title(f'Dilation at {iteration} iterations')

# plt.savefig('../images/dilation/dilated_brain_mask.pdf')
# # plt.show()

# ''' FILLED IN BRAIN MASK '''

# filled_in = brain_mask*inv2

# fig, ax = plt.subplots(1,5)
# ax[0].imshow(ndimage.rotate(filled_in[:,:,64], 90), cmap='gray')
# ax[0].axis('off')
# ax[0].set_title('Original brain mask')
# for i, iteration in enumerate([5, 10, 15, 20]):
#     # dilate brain mask
#     dilation = cv2.dilate(brain_mask, kernel, iterations=iteration)
    
#     # fill in dilated brain mask with second inversion
#     filled_in = dilation*inv2

#     # plot
#     ax[i+1].imshow(ndimage.rotate(filled_in[:,:,64], 90), cmap='gray')
#     ax[i+1].axis('off')
#     ax[i+1].set_title(f'Dilation at {iteration} iterations')

# plt.savefig('../images/dilation/dilated_inv2.pdf')
# plt.show()

horizontal = np.array(
            [[0,0,0,0,0], 
            [0,0,0,0,0],
            [1,1,1,1,1],
            [0,0,0,0,0],
            [0,0,0,0,0]], dtype=np.uint8)
cross = np.array(
            [[0,0,1,0,0], 
            [0,0,1,0,0],
            [1,1,1,1,1],
            [0,0,1,0,0],
            [0,0,1,0,0]], dtype=np.uint8)
vertical = np.array(
            [[0,0,1,0,0], 
            [0,0,1,0,0],
            [0,0,1,0,0],
            [0,0,1,0,0],
            [0,0,1,0,0]], dtype=np.uint8)
dilation_k = cv2.dilate(brain_mask, kernel, iterations=10)
dilation_c = cv2.dilate(brain_mask, cross, iterations=15)
dilation_h = cv2.dilate(brain_mask, horizontal, iterations=15)
dilation_v = cv2.dilate(brain_mask, vertical, iterations=15)
# filled_in = dilation*inv2

t1 = dilation_k * t1corr_data
inv2 = dilation_k * inv2

fig, ax = plt.subplots(1,3)
ax[0].imshow(ndimage.rotate(dilation_k[:,:,64], 90), cmap='gray')
ax[0].axis('off')
ax[0].set_title(f'Brain mask')
ax[1].imshow(ndimage.rotate(inv2[:,:,64], 90), cmap='gray')
ax[1].axis('off')
ax[1].set_title(f'inv2_te1')
ax[2].imshow(ndimage.rotate(t1[:,:,64], 90), cmap='gray')
ax[2].axis('off')
ax[2].set_title(f't1corr')

plt.show()