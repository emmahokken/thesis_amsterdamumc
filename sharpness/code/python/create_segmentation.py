import cv2 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres 
import numpy as np 
import os 
import sys 
from tqdm import tqdm 


# s = '../../data/segm/sub-000_mask-vent_hem-l_lvlreg-corr_def-img.nii'
# s = nib.load(s).get_fdata()
# print(s.shape)
# plt.imshow(s[:,:,190], cmap='gray')
# plt.show()
# exit()

# declare relevant prefixes and directories 
prefix = '../../../../../../..'
ref_path = prefix + '/data/projects/ahead/raw_gdata/'
infile_path = prefix + '/data/projects/ahead/segmentations2/qMRI-Recomputed/'
segm_path = prefix + '/data/projects/ahead/segmentations2/Automated-Parcellation/qmri2/'

# gather all files 
infile_directories = [d[1] for d in os.walk(infile_path)][0]
subjects = ['0'+d[4:] for d in infile_directories]
ref_directories_all = [d[1] for d in os.walk(ref_path)][0]

ref_directories = []
for r in ref_directories_all:
    if  any(x in r for x in subjects):
        ref_directories.append(r)

segm_files = [f[2] for f in os.walk(segm_path)][0]
segm_files.sort()
infile_directories.sort()
ref_directories.sort()

# iterate over all 
for s, i, r in zip(segm_files, infile_directories, ref_directories):

    infile = f'{infile_path}{i}/ses-1/anat/wb/qmri/{i}_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
    reffile = f'{ref_path}{r}/nii/r1corr.nii'
    outfile = f'../../data/coregistration/{i}_outfile.nii'
    transfile = f'../../data/coregistration/{i}_transform_mat.txt'

    # perform coregistration if this has not been done yet
    if not os.path.exists(transfile):
        os.system(f'flirt -in {infile} -ref {reffile} -out {outfile} -omat {transfile}')

    trans_mat = np.loadtxt(transfile)

    # load in binary map
    bin_map_file = nib.load(segm_path + s)
    bin_map = bin_map_file.get_fdata()
    regions = np.unique(bin_map)


    # transform segmented file if this has not been done yet
    segm_outfile = f'../../data/coregistration/{i}_transformed_{s}'
    if not os.path.exists(segm_outfile):
        os.system(f'flirt -in {segm_path + s} -ref {reffile} -out {segm_outfile} -init {transfile} -applyxfm')

    trans_bin_map_file = nib.load(segm_outfile)
    trans_bin_map = trans_bin_map_file.get_fdata()
    trans_bin_map = np.moveaxis(trans_bin_map, 2,0)
    trans_bin_map = trans_bin_map.astype(int)
    trans_regions = np.unique(trans_bin_map)

    # plt.plot(trans_regions)
    # plt.show()
    # print(trans_regions)
    # print(regions)
    # print(np.min(trans_bin_map), np.max(trans_bin_map))
    # print(np.min(bin_map), np.max(bin_map))
    # plt.subplot(121)
    # plt.imshow(trans_bin_map[130,144:155,144:155], vmin=0, vmax=32)

    # plt.subplot(122)
    # plt.imshow(bin_map[130,144:155,144:155], vmin=0,vmax=32)
    # plt.show()


    # exit()
    # iterate over numbers (aka regions) 
    for r in trans_regions[1:]:

        # because the mask is altered, the orignal mask needs to be copied
        working_map = np.copy(trans_bin_map)

        # find indices that do not beling to region and set those to 0
        region = np.where(trans_bin_map != r)
        working_map[region] = 0

        nifit_image = nib.Nifti1Image(dataobj=working_map, header=trans_bin_map_file.header, affine=trans_bin_map_file.affine)
        # nib.save(img=nifit_image, filename=save_dir)

        # compute distance map from that region 
        levelset = nighres.surface.probability_to_levelset(nifit_image, save_data=True)

        levelset_file = levelset['result']
        levelset = levelset_file.get_fdata()


        
    exit()
