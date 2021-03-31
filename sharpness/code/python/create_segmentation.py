import cv2 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres 
import numpy as np 
import os 
import sys 
from tqdm import tqdm 
from scipy import ndimage


# s = '../../data/segm/sub-000_mask-vent_hem-l_lvlreg-corr_def-img.nii'
# s = nib.load(s).get_fdata()
# print(s.shape)
# plt.imshow(s[:,:,190], cmap='gray')
# plt.show()
# exit()

# declare relevant prefixes and directories 
prefix = '../../../../../../..'
r1corr_path = prefix + '/data/projects/ahead/raw_gdata/'
infile_path = prefix + '/data/projects/ahead/segmentations2/qMRI-Recomputed/'
segm_path = prefix + '/data/projects/ahead/segmentations2/Automated-Parcellation/qmri2/'
r2star_path = '../../../data/recon/test_all_3_ssim/R2star_map_gt/'
recon_path = '../../../data/recon/test_all_3_ssim/R2star_map_rim/'

# gather all files 
segm_files = [f[2] for f in os.walk(segm_path)][0]
rim_files = [f[2] for f in os.walk(recon_path)][0]
r2star_files = [f[2] for f in os.walk(r2star_path)][0]
infile_directories = [d[1] for d in os.walk(infile_path)][0]
subjects = ['0'+d[4:] for d in infile_directories]
r1corr_directories_all = [d[1] for d in os.walk(r1corr_path)][0]

# filter out any directories that do not belong to the subjects we will investigiate 
r1corr_directories = []
for r in r1corr_directories_all:
    if  any(x in r for x in subjects):
        r1corr_directories.append(r)

# sort files 
segm_files.sort()
infile_directories.sort()
r1corr_directories.sort()
rim_files.sort()
r2star_files.sort()


# (adjusted) labels (in order) for regions 
labels_structures = ['','str_hem-l','str_hem-r', 'stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r', 'rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
                    'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r', '3v','vent_hem-4','amg_hem-l','amg_hem-r','ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem','pag-l','pag_hem-r',
                    'ppn_hem-l','ppn_hem-r','cl_hem-l','cl_hem-r']

# iterate over all 
for seg, i, r1, rim, r2 in zip(segm_files, infile_directories, r1corr_directories, rim_files, r2star_files):

    infile = f'{infile_path}{i}/ses-1/anat/wb/qmri/{i}_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
    r1corr = f'{r1corr_path}{r1}/nii/r1corr.nii'
    rim_file = f'{recon_path}{rim}'
    r2star_file = f'{r2star_path}{r2}'
    
    outfile_r1corr = f'../../data/coregistration/{i}_outfile_r1cor.nii'
    outfile_rim = f'../../data/coregistration/{i}_outfile_rim.nii'
    outfile_r2star = f'../../data/coregistration/{i}_outfile_r2star.nii'
    transfile_r1corr = f'../../data/coregistration/{i}_transform_mat_r1corr.txt'
    transfile_rim = f'../../data/coregistration/{i}_transform_mat_rim.txt'
    transfile_r2star = f'../../data/coregistration/{i}_transform_mat_r2star.txt'

    # perform coregistration if this has not been done yet
    if not os.path.exists(transfile_r1corr):
        os.system(f'flirt -in {infile} -ref {r1corr} -out {outfile_r1corr} -omat {transfile_r1corr}')

        os.system(f'flirt -in {infile} -ref {rim_file} -out {outfile_rim} -omat {transfile_rim}')
        os.system(f'flirt -in {infile} -ref {r2star_file} -out {outfile_r2star} -omat {transfile_r2star}')

    # transform segmented file if this has not been done yet
    segm_outfile = f'../../data/coregistration/{i}_transformed_{seg}'
    if not os.path.exists(segm_outfile):
        os.system(f'flirt -in {segm_path + seg} -ref {r1corr} -out {segm_outfile} -init {transfile_r1corr} -applyxfm')

    segm_outfile_r2star= f'../../data/coregistration/{i}_transformed_r2star_{seg}'
    segm_outfile_rim = f'../../data/coregistration/{i}_transformed_rim_{seg}'


    # os.system(f'flirt -in {segm_path + seg} -ref {r2star_file} -out {segm_outfile_r2star} -init {transfile_r2star} -applyxfm')
    # os.system(f'flirt -in {segm_path + seg} -ref {rim_file} -out {segm_outfile_rim} -init {transfile_rim} -applyxfm')
    trans_mat = np.loadtxt(transfile_rim)

    # load in binary map
    bin_map_file = nib.load(segm_path + seg)
    bin_map = bin_map_file.get_fdata()
    regions = np.unique(bin_map)

    trans_rim_map_file = nib.load(segm_outfile_rim)
    trans_rim_map = trans_rim_map_file.get_fdata()
    trans_r2star_map_file = nib.load(segm_outfile_r2star)
    trans_r2star_map = trans_r2star_map_file.get_fdata()

    trans_rim_map = np.moveaxis(trans_rim_map, 0,2)
    trans_rim_map = trans_rim_map.astype(int)
    trans_rim_regions = np.unique(trans_rim_map)

    trans_r2star_map = np.moveaxis(trans_r2star_map, 0,2)
    trans_r2star_map = trans_r2star_map.astype(int)
    trans_r2star_regions = np.unique(trans_r2star_map)

    plt.subplot(131)
    plt.imshow(ndimage.rotate(bin_map[:,bin_map.shape[1]//2,:],90))
    plt.subplot(132)
    plt.imshow(trans_rim_map[:,trans_rim_map.shape[1]//2,:])
    plt.subplot(133)
    plt.imshow(trans_r2star_map[:,trans_r2star_map.shape[1]//2,:])
    plt.show()

    print(bin_map.shape)

    print(trans_rim_map.shape)
    print(trans_r2star_map.shape)
    exit()
    # iterate over numbers (aka regions) 
    for r in trans_rim_regions[1:]:

        # because the mask is altered, the orignal mask needs to be copied
        working_map_rim = np.copy(trans_rim_map)
        working_map_r2star = np.copy(trans_r2star_map)
        # find indices that do not beling to region and set those to 0
        region_rim = np.where(trans_rim_map != r)
        working_map_rim[region_rim] = 0
    
        region_r2star = np.where(trans_r2star_map != r)
        working_map_r2star[region_r2star] = 0

        nifit_image_rim = nib.Nifti1Image(dataobj=working_map_rim, header=trans_rim_map_file.header, affine=trans_rim_map_file.affine)
        nifit_image_r2star = nib.Nifti1Image(dataobj=working_map_r2star, header=trans_rim_map_file.header, affine=trans_rim_map_file.affine)
        # nib.save(img=nifit_image, filename=save_dir)

        region_name = labels_structures[r]
        
        r2star_levelset_outfile = f'{i}_mask-{region_name}_lvlreg-gt_def-img.nii'
        rim_levelset_outfile = f'{i}_mask-{region_name}_lvlreg-rim_def-img.nii'
        output_dir = '../../data/coregistration/'
        # compute distance map from that region 
        levelset = nighres.surface.probability_to_levelset(nifit_image_rim, save_data=True, output_dir=output_dir, file_name=rim_levelset_outfile)
        levelset = nighres.surface.probability_to_levelset(nifit_image_r2star, save_data=True, output_dir=output_dir, file_name=r2star_levelset_outfile)

        # levelset_file = levelset['result']
        # levelset = levelset_file.get_fdata()


        
    exit()
