import cv2 
import nibabel as nib 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import numpy as np
import os
from scipy import ndimage
from tqdm import tqdm

subj = 5
subj2 = 3
base = '../../../../..'
segm_path = f'/data/projects/ahead/segmentations/mapped/'
inv2_path = f'/data/projects/ahead/raw_gdata/Subcortex_{subj:04}_0{subj2:03}_R02/nii/inv2_te1_m_corr.nii'
r1corr_path = f'/data/projects/ahead/raw_gdata/Subcortex_{subj:04}_0{subj2:03}_R02/nii/r1corr.nii'
recon_nii_path = f'../../data/recon/test_all_3_ssim/R2star_map_rim/results_Subcortex_{subj:04}_axial_loss_mse.nii'
brain_mask_path = f'/data/projects/ahead/raw_gdata/Subcortex_{subj:04}_{subj2:03}_R02/nii/mask_inv2_te2_m_corr.nii'
coil_path = '/data/projects/ahead/raw_gdata/Subcortex_0005_002_R02/Subcortex_0005_002_R02_inv2_2_gdataCorrected.nii.gz'

# load in main images 
inv2 = nib.load(inv2_path).get_fdata()
r1corr = nib.load(r1corr_path).get_fdata()
coil = nib.load(coil_path).get_fdata()
print(coil.shape)
print(inv2.shape)

exit()


kernel = np.array([[1,1,1,1,1], [1,1,1,1,1], [1,1,1,1,1], [1,1,1,1,1], [1,1,1,1,1]])

fig, ax = plt.subplots(1,2)
z = 149


# iterate over all ROIs
for subdir, dirs, files in os.walk(base+segm_path):

    # create color range for scatter 
    files = [f for f in files if f'sub-0{subj}_mask' in f]
    colours = cm.rainbow(np.linspace(0, 1, len(files)))

    for i, file in tqdm(enumerate(files)):
        s_nii = nib.load(base+segm_path+file)
        s = s_nii.get_fdata()

        # create mask of segmentation
        s[s > 0] = 0
        s[s < 0] = 1
        
        # calculate boundary locations for mask
        distance = ndimage.distance_transform_cdt(s[:,:,z], 'taxicab') == 1
        boundary_idx = np.vstack(np.where(distance == 1))

        # zip x and y axes 
        # boundary = list(zip(boundary_idx[0], boundary_idx[1]))

        # create label for ROI
        names = file.split('-')
        roi, _ = names[2].split('_')
        hem = names[3][0]
        side = 'left' if hem == 'l' else 'right' if hem == 'r' else 'fourth'
        
        ax[0].scatter(boundary_idx[0], boundary_idx[1], color=colours[i], s=0.3, marker='o', label=f'{roi}, {side}')
        ax[1].scatter(boundary_idx[1], boundary_idx[0], color=colours[i], s=0.3, marker='o', label=f'{roi}, {side}')

    break 
       
ax[0].imshow(inv2[:,:,z], cmap='gray')
ax[1].imshow(ndimage.rotate(inv2[:,:,z], 90), cmap='gray')
# plt.axis('off')
plt.legend(loc='upper left', bbox_to_anchor=(1,1))
plt.show()
