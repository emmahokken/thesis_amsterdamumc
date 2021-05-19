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


def dilate(file_path, save_path='', save_mask=True, kw=5, kh=5, kd=5, iters=10, brain_mask_type='mask_inv2_te2_m_corr'):
    ''' 
    Dilate a given brain mask. 

    Args: 
        file_path: file path to the brain mask 
        save_path: file path to save the dilated brain mask to. default is current directory.
        save_mask: boolean value whether to save the mask
        kw, kh, kd: dimensions (width x height x depth) of dilation kernel
        iters: the number of iterations dilation should be applied for. default is 15.
        brain_mask_type: filename of brain mask. default is 'mask_inv2_te2_m_corr'.

    Returns: 
        gaus_dilation: smoothed and dilated brain mask
        inv_dilation: inverse of the smoothed and dilated brain mask 
    '''
    
    # dilation parameters 
    kernel = np.ones((kw,kh,kd), dtype=np.uint8)

    # get brain mask and scan from within directory
    brain_mask_path = f'{file_path}/nii/{brain_mask_type}.nii'

    # load in brain mask (in two steps for saving later)
    brain_mask_nii = nib.load(brain_mask_path)
    brain_mask = brain_mask_nii.get_fdata()

    # dilate brain mask and correct to 1 (was 0.999...)
    dilation = ndimage.binary_dilation(brain_mask, kernel, iterations=iters)
    dilation[dilation > 0] = 1

    # plt.imshow(dilation[:,130,:])
    # plt.show()
    
    # smooth edges
    gaus_dilation = cv2.GaussianBlur(np.float32(dilation), (5, 5), 0)

    # create inverse mask of dilatio
    inv_dilation = abs(gaus_dilation - 1)
    inv_dilation[inv_dilation < 1] = 0
    
    gaus_dilation = np.moveaxis(gaus_dilation, 2, 0)
    inv_dilation = np.moveaxis(inv_dilation, 2, 0)

    gaus_dilation = np.fliplr(gaus_dilation)
    inv_dilation = np.fliplr(inv_dilation)

    if save_mask:
        nifti_image = nib.Nifti1Image(gaus_dilation, header=brain_mask_nii.header, affine=brain_mask_nii.affine)
        nib.save(nifti_image, f'{save_path}/{brain_mask_type}_dilated.nii')
    
    return gaus_dilation, inv_dilation

