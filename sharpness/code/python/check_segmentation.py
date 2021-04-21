import cv2 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres 
import numpy as np 
import os 
import sys 
from tqdm import tqdm 
from scipy import ndimage

prefix = '../../../../../../..'
infile_path = prefix + '/data/projects/ahead/segmentations2/qMRI-Recomputed/'

infile = f'{infile_path}sub-005/ses-1/anat/wb/qmri/sub-005_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
coreg_mat_file = '../../data/coregistration/sub-005_outfile_mat.nii.gz'
coreg_r2star_file = '../../data/coregistration/sub-005_outfile_r2star.nii.gz'
mat_path = '../../../../data/converted_mat_echo1/Subcortex_0005_002_R02/Subcortex_0005_002_R02_allm.nii'
coreg_r1corr_file = f'../../data/coregistration/sub-005_outfile_r1corr.nii.gz'
r1corr_file = f'../../data/coregistration/sub-005_r1corr_std.nii.gz'
bin_map_file = '../../data/coregistration/sub-005bin_map.nii'

coreg_mat = nib.load(coreg_mat_file).get_fdata()
coreg_r2star = nib.load(coreg_r2star_file).get_fdata()
mat = nib.load(mat_path).get_fdata()
inf = nib.load(infile).get_fdata()
r1corr = nib.load(r1corr_file).get_fdata()
coreg_r1corr = nib.load(coreg_r1corr_file).get_fdata()
# bin_map = nib.load(bin_map_file).get_fdata()

print(coreg_r1corr.shape)
print(r1corr.shape)
print(inf.shape)
print(bin_map.shape)

plt.subplot(221)
plt.imshow(coreg_r1corr[:,:,146], cmap='gray')
plt.title('Coregisterd allm file')
plt.axis('off')
# plt.subplot(222)
# plt.imshow(bin_map[:,:,146], cmap='gray')
# plt.title('Binary segmentation map')
# plt.axis('off')
plt.subplot(223)
plt.imshow(r1corr[:,:,146], cmap='gray')
plt.axis('off')
plt.title('Original allm file')
plt.subplot(224)
plt.imshow(inf[:,:,146], cmap='gray')
plt.axis('off')
plt.title('Segmentation aligned file')
plt.show()