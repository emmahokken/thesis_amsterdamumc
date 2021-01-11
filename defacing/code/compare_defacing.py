import cv2 
import matplotlib.patches as patches 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres
import numpy as np 
import os 
from scipy import ndimage
from skimage import filters 
 
from defacing import deface

# load defaced images (for subject 0004)
fs_inv2 = nib.load('../defacing/freesurfer_defaced_inv2.nii').get_fdata()
pd_inv2 = nib.load('../defacing/pydeface_inv2.nii').get_fdata()

# declare image paths 
t1corr_path = '../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/t1corr.nii'
inv2_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/inv2_te1_m_corr.nii'
brain_mask_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/mask_inv2_te2_m_corr.nii'

# load images 
t1corr = nib.load(t1corr_path).get_fdata()
inv2 = nib.load(inv2_path).get_fdata()

# deface both t1corr and inv2
# t1corr_defaced_freesurfer, fs_t1_mask = deface(t1corr_path, brain_mask_path, '../defacing/freesurfer_t1corr.nii', 'freesurfer')
# inv2_defaced_freesurfer, fs_inv2_mask = deface(inv2_path, brain_mask_path, '../defacing/freesurfer_inv2.nii', 'freesurfer')

# t1corr_defaced_pydeface, pd_t1_mask = deface(t1corr_path, brain_mask_path, '../defacing/pydeface_t1corr.nii', 'pydeface')
# inv2_defaced_pydeface, pd_inv2_mask = deface(inv2_path, brain_mask_path, '../defacing/pydeface_inv2_new.nii', 'pydeface')


fs_t1corr = nib.load('../defacing/freesurfer_t1corr.nii').get_fdata()
fs_inv2 = nib.load('../defacing/freesurfer_inv2.nii').get_fdata()

pd_t1corr = nib.load('../defacing/pydeface_t1corr.nii').get_fdata()
pd_inv2 = nib.load('../defacing/pydeface_inv2_new.nii').get_fdata()

# fs_inv2[fs_inv2 > 0] = 1 
# pd_inv2[pd_inv2 > 0] = 1 

# fs_inv2 = t1corr * fs_inv2
# pd_inv2 = t1corr * pd_inv2

fig, ax= plt.subplots(2,2)
ax[0,0].imshow(ndimage.rotate(fs_t1corr[:,:,63], 90), cmap='gray')
ax[0,0].axis('off')
ax[0,0].set_ylabel('t1corr')
ax[0,0].set_title('Freesurfer defacing')
ax[1,0].imshow(ndimage.rotate(fs_inv2[:,:,63], 90), cmap='gray')
ax[1,0].axis('off')
ax[0,1].set_title('Pydeface defacing')
ax[0,1].imshow(ndimage.rotate(pd_t1corr[:,:,63], 90), cmap='gray')
ax[0,1].axis('off')
ax[0,1].set_ylabel('inv2')
ax[1,1].imshow(ndimage.rotate(pd_inv2[:,:,63], 90), cmap='gray')
ax[1,1].axis('off')
plt.savefig('../defacing/freesurfer_vs_pydeface_start.pdf')
plt.show()