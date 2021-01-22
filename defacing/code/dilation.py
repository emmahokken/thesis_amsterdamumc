import cv2 
import matplotlib.patches as patches 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres
import numpy as np 
import os 
from scipy import ndimage
from tqdm import tqdm

from defacing import deface


# decalre root path and scan type
root_path = '../../../../../data/projects/ahead/raw_gdata'
scan_type = 'inv2_te1_m_corr'
brain_mask_type = 'mask_inv2_te2_m_corr'
save_path = '../results/dilation'

# dilation parameters 
kernel = np.ones((5,5), dtype=np.uint8)
iters = 13

# iterate over all subjects
for subdir, dirs, files in os.walk(root_path):
    for d in tqdm(dirs):
        # get brain mask and scan from within directory
        brain_mask_path = f'{root_path}/{d}/nii/{brain_mask_type}.nii'
        path = f'{root_path}/{d}/nii/{scan_type}.nii'

        # load in brain mask (in two steps for saving later)
        # try:
        brain_mask_nii = nib.load(brain_mask_path)
        brain_mask = brain_mask_nii.get_fdata()
        inv2_nii = nib.load(path)
        inv2 = inv2_nii.get_fdata()
    
        # small = inv2[:40,:40,64]
        # print('mean', np.mean(small), 'std:', np.std(small))
        # # plt.imshow(small)
        # # plt.show()
        # exit()

        # dilate brain mask
        dilation = cv2.dilate(brain_mask, kernel, iterations=iters)

        # smooth edges
        gaus_dilation = cv2.GaussianBlur(dilation, (5,5), 0)

        # create inverse mask of dilation
        inv_dilation = abs(dilation - 1)

        # add random noise to background 
        small = inv2[:40,:40,64]
        noise = np.random.normal(np.mean(small), np.std(small),inv_dilation.shape)
        noise_mask = noise * inv_dilation

        # fill in with (second inversion) scan and add noise to background 
        filled_in = (inv2*gaus_dilation) + noise_mask
       
        # save new scan 
        if not os.path.exists(f'{save_path}/{d}'):
            os.makedirs(f'{save_path}/{d}')
        save_dir = f'{save_path}/{d}/dilated_{scan_type}.nii'

        nifit_image = nib.Nifti1Image(dataobj=filled_in, header=inv2_nii.header, affine=inv2_nii.affine)
        nib.save(img=nifit_image, filename=save_dir)

    break 
