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

# ''' LOADING IMAGES '''
# # get filenames 
# face = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/inv2_te1_m_corr.nii'

# face_allm = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/Subcortex_0004_002_R02_allm.nii'
# filenm = face.split('/')
# no_face = f'../defacing/fs1{filenm[-1]}'
# brain_mask = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/mask_inv2_te2_m_corr.nii'
# test = '../defacing/001.nii'

# # perform defacing (this should be done in the command line, hence the os.system call)
# # os.system(f'mri_deface {face} ../defacing/talairach_mixed_with_skull.gca ../defacing/face.gca {no_face}')
# # os.system(f'pydeface {face} --outfile {no_face} --force')

# # create image of original scan 
# img = nib.load(face).get_fdata()
# conv = nib.load(test).get_fdata()

# # load and fiddle with alternate images 
# alt_img = nib.load(face_allm).get_fdata()
# alt_img = np.transpose(alt_img,(1,2,0))
# alt_img = alt_img[::-1,:,::-1]

# # create image of brain mask
# brain_mask = nib.load(brain_mask).get_fdata()

# # create image of defaced scan
# nf_img = nib.load(no_face).get_fdata()

# # correct for removal of brain tissue
# nf_img[brain_mask > 0] = img[brain_mask > 0]

# ''' SMOOTHING '''
# # create (inverse) mask of defaced region
# face_mask = np.array((img - nf_img) == 0, dtype=float)

# # perform a guassian filter 
# gaus_mask = gaussian_filter(face_mask, sigma=5)
# gaus_img = img*(gaus_mask)

# alt = alt_img*(gaus_mask)


# # fig, ax = plt.subplots(1,2)
# # ax[0].imshow(ndimage.rotate(gaus_img[60],90))
# # ax[1].imshow(ndimage.rotate(alt[60],90))
# # plt.show()

# ''' PLOTTING '''
# # declare x,y,z slices
# x, y, z = 58, 68, 62

# fig, ax = plt.subplots(3,3)
# # plt.style.use('dark_background')

# # with face
# ax[0,0].imshow(ndimage.rotate(img[x],180), cmap='gray')
# ax[0,1].imshow(ndimage.rotate(img[:,y],180), cmap='gray')
# ax[0,2].imshow(ndimage.rotate(img[:,:,z],90), cmap='gray')
# ax[0,1].set_title('With face',color='w')

# # wihout face
# ax[1,0].imshow(ndimage.rotate(nf_img[x],180), cmap='gray')
# ax[1,1].imshow(ndimage.rotate(nf_img[:,y],180), cmap='gray')
# ax[1,2].imshow(ndimage.rotate(nf_img[:,:,z],90), cmap='gray')
# ax[1,1].set_title('Without face',color='w')

# # difference
# ax[2,0].imshow(ndimage.rotate(gaus_img[x],180), cmap='gray')
# ax[2,1].imshow(ndimage.rotate(gaus_img[:,y],180), cmap='gray')
# ax[2,2].imshow(ndimage.rotate(gaus_img[:,:,z],90), cmap='gray')
# ax[2,1].set_title('Difference', color='w')

# # mask
# # ax[2,0].imshow(brain_mask[x], cmap='gray')
# # ax[2,1].imshow(brain_mask[:,y], cmap='gray')
# # ax[2,2].imshow(brain_mask[:,:,brain_mask.shape[0]//2], cmap='gray')
# # ax[2,1].set_title('Mask')


# # adjust space between plots
# plt.subplots_adjust(hspace=0.3,wspace=0.1)
# # remove tickmarks on axes
# for a in ax:
#     for b in a:
#         b.set_xticks([])
#         b.set_yticks([])

# plt.savefig(f'../defacing_{filenm}.png')
# plt.close()

# # # plot smoothing difference (if any)
# # fig, ax = plt.subplots(1,2)
# # old = nib.load(no_face).get_fdata()
# # old = ndimage.rotate(old, 90)
# # ax[0].imshow(old[:,:,old.shape[0]//2], cmap='gray')
# # ax[0].set_title('without smoothing')
# # ax[1].imshow(nf_img[:,:,nf_img.shape[0]//2], cmap='gray')
# # ax[1].set_title('with smoothing')
# # plt.savefig('../smoothing.png')



# ''' CORRECT DEFACING ''' 
# # check whether brain volume is removed by comparing defaced scan to mask

# # if mask has white value but defaced has black, put value back. 


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