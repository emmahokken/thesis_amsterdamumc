import cv2
import matplotlib.pyplot as plt 
import matplotlib.patches as patches 
import nibabel as nib 
import nighres 
import numpy as np 
import os 
from scipy import ndimage
from scipy.ndimage import gaussian_filter
import sys 

def deface(image_path, brain_mask_path, defaced_path, method='freesurfer'):
    ''' 
    Function that performs the defacing and smoothing procedure for a given image and brain mask path. 

    Args:
        image_path: path to the input (NIFTI) image. 
        brain_mask_path: path to the corresponding (NIFTI) brain mask to the image. Must be same dimensionality as input image. 
    
    Returns:
        gaus_img: the Gaussian smoothed defaced image. 
     '''
   
    ''' DEFACING '''
    if method == 'freesurfer':
        # perform defacing (this should be done in the command line, hence the os.system call)
        os.system(f'mri_deface {image_path} ../defacing/talairach_mixed_with_skull.gca ../defacing/face.gca {defaced_path}')
    elif method == 'pydeface':
        os.system(f'pydeface {image_path} --outfile {defaced_path}')
    else:
        sys.exit("Please provide a defacing method. Must be either \'freesurfer\' or \'pydeface\'.")

    # load in needed images 
    img = nib.load(image_path).get_fdata()
    nf_img = nib.load(defaced_path).get_fdata()
    brain_mask = nib.load(brain_mask_path).get_fdata()
      
    # correct for removal of brain tissue
    nf_img[brain_mask > 0] = img[brain_mask > 0]
    
    ''' SMOOTHING '''
    # create (inverse) mask of defaced region
    face_mask = np.array((img - nf_img) == 0, dtype=float)

    # perform a guassian filter 
    gaus_mask = cv2.GaussianBlur(face_mask, (5,5), 0)
    gaus_img = img*(gaus_mask)

    return gaus_img, gaus_mask