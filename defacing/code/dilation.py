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


def dilate(file_path, kw=5, kh=5, iters=13, brain_mask_type='mask_inv2_te2_m_corr'):

    # dilation parameters 
    kernel = np.ones((kw,kh), dtype=np.uint8)
  
    # get brain mask and scan from within directory
    brain_mask_path = f'{file_path}/nii/{brain_mask_type}.nii'

    # load in brain mask (in two steps for saving later)
    brain_mask_nii = nib.load(brain_mask_path)
    brain_mask = brain_mask_nii.get_fdata()

    # dilate brain mask
    dilation = cv2.dilate(brain_mask, kernel, iterations=iters)

    # smooth edges
    gaus_dilation = cv2.GaussianBlur(dilation, (kw,kh), 0)

    # create inverse mask of dilatio
    inv_dilation = abs(dilation - 1)
    inv_dilation[inv_dilation < 1] = 0
    
    return gaus_dilation, inv_dilation

   
    # fill in with (second inversion) scan and add noise to background 
    filled_in = (inv2*gaus_dilation)
    filled_in_noise = (inv2*gaus_dilation) + noise_mask

