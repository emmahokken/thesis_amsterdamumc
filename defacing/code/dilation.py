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
wrong_files = []

test_scan = 'Subcortex_0054_054_R02'
# dilation parameters 
kernel = np.ones((5,5), dtype=np.uint8)
iters = 13
wrong = 0
f = nib.load(f'{root_path}/{test_scan}/nii/{brain_mask_type}.nii')

# ah = nib.load("../../../../../data/projects/ahead/raw_gdata/Subcortex_0054_054_R02/nii/mask_inv2_te2_m_corr.nii")
# exit()
# iterate over all subjects
for subdir, dirs, files in os.walk(root_path):
    for d in tqdm(dirs):
        # get brain mask and scan from within directory
        brain_mask_path = f'{root_path}/{d}/nii/{brain_mask_type}.nii'
        path = f'{root_path}/{d}/nii/{scan_type}.nii'

        # load in brain mask (in two steps for saving later)
        try:
            brain_mask_nii = nib.load(brain_mask_path)
            brain_mask = brain_mask_nii.get_fdata()
            inv2_nii = nib.load(path)
            inv2 = inv2_nii.get_fdata()
        except:
        # except nib.filebasedimages.ImageFileError:
            # print('things went wrong')
            wrong += 1
            wrong_files.append(d)

        # # dilate brain mask and fill in with (second inversion) scan 
        # dilation = cv2.dilate(brain_mask, kernel, iterations=iters)
        # filled_in = dilation*inv2
        
        # # save new scan 
        # if not os.path.exists(f'{save_path}/{d}'):
        #     os.makedirs(f'{save_path}/{d}')
        # save_dir = f'{save_path}/{d}/dilated_{scan_type}.nii'

        # nifit_image = nib.Nifti1Image(dataobj=filled_in, header=inv2_nii.header, affine=inv2_nii.affine)
        # nib.save(img=nifit_image, filename=save_dir)

    break 


print(f'things went wrong {wrong} times')