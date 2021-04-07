import cv2 
import matplotlib.patches as patches 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres
import numpy as np 
import os 
from scipy import ndimage
# from skimage import filters 
 
from defacing import deface

# declare image paths 
prefix = '../../../../../..'
t1corr_path = prefix + '/data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/t1corr.nii'
inv2_path = prefix + '/data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/inv2_te1_m_corr.nii'
brain_mask_path = prefix + '/data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/mask_inv2_te2_m_corr.nii'

# load images 
# t1corr = nib.load(t1corr_path).get_fdata()
# inv2 = nib.load(inv2_path).get_fdata()

t1corr_defaced_freesurfer_path_04 = '../results/freesurfer/freesurfer_t1corr_0004.nii'
inv2_defaced_freesurfer_path_04 = '../results/freesurfer/freesurfer_inv2_0004.nii' 
t1corr_defaced_freesurfer_path_64 = '../results/freesurfer/freesurfer_t1corr_0064.nii'
inv2_defaced_freesurfer_path_64 = '../results/freesurfer/freesurfer_inv2_0064.nii' 

t1corr_defaced_pydeface_path = '../results/pydeface/pydeface_t1corr.nii'
inv2_defaced_pydeface_path = '../results/pydeface/pydeface_inv2_new.nii'

# deface both t1corr and inv2
if not os.path.exists(t1corr_defaced_freesurfer_path_04):
    t1corr_defaced_freesurfer, fs_t1_mask = deface(t1corr_path, brain_mask_path, t1corr_defaced_freesurfer_path_04, 'freesurfer')
    inv2_defaced_freesurfer, fs_inv2_mask = deface(inv2_path, brain_mask_path, inv2_defaced_freesurfer_path_04, 'freesurfer')

    # t1corr_defaced_pydeface, pd_t1_mask = deface(t1corr_path, brain_mask_path, t1corr_defaced_pydeface_path, 'pydeface')
    # inv2_defaced_pydeface, pd_inv2_mask = deface(inv2_path, brain_mask_path, inv2_defaced_pydeface_path, 'pydeface')

# load original image 
t1corr = nib.load(t1corr_path).get_fdata()
inv2 = nib.load(inv2_path).get_fdata()

# load defaced images (for subject 0004)
fs_t1corr_04 = nib.load(t1corr_defaced_freesurfer_path_04).get_fdata()
fs_inv2_04 = nib.load(inv2_defaced_freesurfer_path_04).get_fdata()
fs_t1corr_64 = nib.load(t1corr_defaced_freesurfer_path_64).get_fdata()
fs_inv2_64 = nib.load(inv2_defaced_freesurfer_path_64).get_fdata()

pd_t1corr = nib.load(t1corr_defaced_pydeface_path).get_fdata()
pd_inv2 = nib.load(inv2_defaced_pydeface_path).get_fdata()

# fs_inv2[fs_inv2 > 0] = 1 
# pd_inv2[pd_inv2 > 0] = 1 

# fs_inv2 = t1corr * fs_inv2
# pd_inv2 = t1corr * pd_inv2

fig, ax= plt.subplots(2,2)
ax[0,0].imshow(ndimage.rotate(fs_t1corr_04[:,:,63], 90), cmap='gray')
ax[0,0].axis('off')
ax[0,0].set_ylabel('t1corr')
ax[0,0].set_title('Freesurfer defacing')
ax[1,0].imshow(ndimage.rotate(fs_inv2_04[:,:,63], 90), cmap='gray')
ax[1,0].axis('off')
ax[0,1].set_title('Pydeface defacing')
ax[0,1].imshow(ndimage.rotate(pd_t1corr[:,:,63], 90), cmap='gray')
ax[0,1].axis('off')
ax[0,1].set_ylabel('inv2')
ax[1,1].imshow(ndimage.rotate(pd_inv2[:,:,63], 90), cmap='gray')
ax[1,1].axis('off')
plt.savefig('../results/figures/freesurfer_vs_pydeface_start.pdf')
# plt.show()


fig, ax= plt.subplots(2,2)
ax[0,0].imshow(ndimage.rotate(t1corr[:,:,63], 90), cmap='gray')
ax[0,0].set_ylabel('Original', rotation=0, labelpad=30)
ax[0,0].set_title('T1 corrected')
ax[1,0].imshow(ndimage.rotate(fs_t1corr_04[:,:,63], 90), cmap='gray')
ax[1,0].set_ylabel('Defaced', rotation=0, labelpad=30)
ax[0,1].set_title('Second inversion')
ax[0,1].imshow(ndimage.rotate(inv2[:,:,63], 90), cmap='gray')
ax[1,1].imshow(ndimage.rotate(fs_inv2_04[:,:,63], 90), cmap='gray')

ax[0,0].set_xticks([])
ax[1,0].set_xticks([])
ax[0,1].set_xticks([])
ax[1,1].set_xticks([])

ax[0,0].set_yticks([])
ax[1,0].set_yticks([])
ax[0,1].set_yticks([])
ax[1,1].set_yticks([])
plt.tight_layout()
# fig.subplots_adjust(wspace=0,hspace=0)
plt.savefig('../results/figures/freesurfer_comparison.pdf')
plt.show()
