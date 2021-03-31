'''Coregister MRI images, done in python file for convenience'''

import matplotlib.pyplot as plt
import nibabel as nib 
import numpy as np
import os 
import sys 

subj = 5
prefix = '../../../../..'
infile = f'/data/projects/ahead/segmentations2/qMRI-Recomputed/sub-{subj:03}/ses-1/anat/wb/qmri/sub-{subj:03}_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
reffile = f'/data/projects/ahead/raw_gdata/Subcortex_{subj:04}_002_R02/nii/r1corr.nii'
outfile = f'data/coregistration/sub-{subj:03}_outfile.nii'
trans_mat = f'data/coregistration/sub-{subj:03}_transform_mat.txt'

segm2_file = f'/data/projects/ahead/segmentations2/Automated-Parcellation/qmri2/sub-{subj:03}_atlas10_31struct_qmri2_massp-label.nii.gz'
os.system('echo Hello')

# os.system(f'flirt -in {prefix + infile} -ref {prefix + reffile} -out {outfile} -omat {trans_mat}')


transm = np.loadtxt(trans_mat)

segm2 = nib.load(segm2_file).get_fdata()
print(segm2.shape)
plt.imshow(segm2[:,100,:])
plt.show()

