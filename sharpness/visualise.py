import nibabel as nib 
import matplotlib.pyplot as plt 
import numpy as np
import os
from tqdm import tqdm

subj = 64
base = '../../../../..'
segm_path = f'/data/projects/ahead/segmentations/mapped/'
inv2_path = f'/data/projects/ahead/raw_gdata/Subcortex_00{subj}_0{subj}_R02/nii/inv2_te1_m_corr.nii'
r1corr_path = f'/data/projects/ahead/raw_gdata/Subcortex_00{subj}_0{subj}_R02/nii/r1corr.nii'
recon_nii_path = f'../../data/recon/test_all_3_ssim/R2star_map_rim/results_Subcortex_00{subj}_axial_loss_mse.nii'
brain_mask_path = f'/data/projects/ahead/raw_gdata/Subcortex_00{subj}_0{subj}_R02/nii/mask_inv2_te2_m_corr.nii'

# load in main images 
inv2 = nib.load(inv2_path).get_fdata()
r1corr = nib.load(r1corr_path).get_fdata()

# iterate over all segmentations
for subdir, dirs, files in os.walk(base+segm_path):
    for file in tqdm(files):
        if f'sub-0{subj}_mask' in file:
            s_nii = nib.load(base+segm_path+file)
            s = s_nii.get_fdata()

            # create mask of segmentation
            s[s > 0] = 0
            s[s < 0] = 1

            plt.imshow(s[:,:,140])
            plt.plot(np.linspace(0,200))
            plt.show()

            exit()

plt.imshow(f[:,:,0], cmap='gray')
plt.show()