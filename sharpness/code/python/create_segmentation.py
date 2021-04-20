import cv2 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres 
import numpy as np 
import os 
import sys 
from tqdm import tqdm 
from scipy import ndimage


# declare relevant prefixes and directories 
prefix = '../../../../../../..'
infile_path = prefix + '/data/projects/ahead/segmentations2/qMRI-Recomputed/'
segm_path = prefix + '/data/projects/ahead/segmentations2/Automated-Parcellation/qmri2/'
r2star_path = '../../../../data/recon/test_all_3_ssim/R2star_map_gt/'
mat_path = '../../../../data/converted_mat_echo1/'

# gather all files 
segm_files = [f[2] for f in os.walk(segm_path)][0]
r2star_files = [f[2] for f in os.walk(r2star_path)][0]
infile_directories = [d[1] for d in os.walk(infile_path)][0]
subjects = ['0'+d[4:] for d in infile_directories]
mat_directoriess = [f[1] for f in os.walk(mat_path)][0]

# sort files 
segm_files.sort()
infile_directories.sort()
r2star_files.sort()

# (adjusted) labels (in order) for regions 
labels_structures = ['','str_hem-l','str_hem-r', 'stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r', 'rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
                    'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r', '3v','vent_hem-4','amg_hem-l','amg_hem-r','ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem','pag-l','pag_hem-r',
                    'ppn_hem-l','ppn_hem-r','cl_hem-l','cl_hem-r']

# iterate over new segmentations, infiles, r2star maps
for seg, i, r2, mat in zip(segm_files, infile_directories, r2star_files, mat_directoriess):

    print('Working on', i)
    
    infile = f'{infile_path}{i}/ses-1/anat/wb/qmri/{i}_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
    r2star_file = f'{r2star_path}{r2}'
    mat_file = f'{mat_path}{mat}/{mat}_allm.nii'
    outfile_r2star = f'../../data/coregistration/{i}_outfile_mat.nii'
    transfile = f'../../data/coregistration/{i}_transform_mat_mat.txt'
    

    # perform coregistration if this has not been done yet
    if not os.path.exists(outfile_r2star):
        print('coregistering with mat file')
        os.system(f'flirt -in {infile} -ref {mat_file} -out {outfile_r2star} -omat {transfile}')

    # load in binary map
    bin_map_file = nib.load(segm_path + seg)
    bin_map = bin_map_file.get_fdata()
    regions = np.unique(bin_map)

    # iterate over numbers (aka regions) 
    for r in regions[1:]:

        region_name = labels_structures[int(r)]

        # because the mask is altered, the orignal mask needs to be copied
        working_map = np.copy(bin_map)
       
        # find indices that do not beling to region and set those to 0
        inv_region = np.where(bin_map != r)
        working_map[inv_region] = 0
        region = np.where(bin_map == r)
        working_map[region] = 1

        # change back into nifti image (needed for nighres)
        nifti_image = nib.Nifti1Image(dataobj=working_map, header=bin_map_file.header, affine=bin_map_file.affine)
        
        levelset_outfile = f'{i}_mask-{region_name}_lvlreg-gt_def-img_2.nii'
        output_dir = '../../data/coregistration/'

        # compute distance map from that region 
        levelset = nighres.surface.probability_to_levelset(nifti_image, save_data=True, output_dir=output_dir, file_name=levelset_outfile)

        # transform segmented file
        os.system(f'flirt -in {levelset["result"]} -ref {r2star_file} -out {output_dir+levelset_outfile} -init {transfile} -applyxfm')
        
        # remove file created with nighres
        os.remove(levelset['result'])

        # load created levelset
        levelset_file = nib.load(f'{output_dir}{levelset_outfile}.gz')
        levelset = levelset_file.get_fdata()

        # sometimes gzip doesn't work, so remove gzp alltogether 
        os.remove(f'{output_dir}{levelset_outfile}.gz')

        # flip any axes?
        # flipped = np.flip(levelset, axis=(2))
        flipped = levelset

        # save flipped image
        nifti_image = nib.Nifti1Image(dataobj=flipped, header=levelset_file.header, affine=levelset_file.affine)
        nib.save(nifti_image, '../../data/segm/'+levelset_outfile)

