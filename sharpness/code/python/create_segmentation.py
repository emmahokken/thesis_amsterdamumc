import cv2 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres 
import numpy as np 
import os 
import sys 
from tqdm import tqdm 
from scipy import ndimage



def create_segmentation(create_gp=False):
    '''
    Creates levelmap segmentations for 32 new regions. 
    '''

    # declare relevant prefixes and directories 
    prefix = '../../../../../../..'
    infile_path = prefix + '/data/projects/ahead/segmentations2/qMRI-Recomputed/'
    segm_path = prefix + '/data/projects/ahead/segmentations2/Automated-Parcellation/qmri2/'
    r1corr_path = prefix + '/data/projects/ahead/raw_gdata/'
    output_dir = '../../data/coregistration/'

    # load files
    segm_files, infile_directories, r1corr_directories = load_files(infile_path, segm_path, r1corr_path)

    # (adjusted) labels (in order) for regions 
    labels_structures = ['','str_hem-l','str_hem-r', 'stn_hem-l','stn_hem-r','sn_hem-l','sn_hem-r', 'rn_hem-l','rn_hem-r','gpi_hem-l','gpi_hem-r','gpe_hem-l','gpe_hem-r',
                        'tha_hem-l','tha_hem-r','vent_hem-l','vent_hem-r', '3v','vent_hem-4','amg_hem-l','amg_hem-r','ic_hem-l','ic_hem-r','vta_hem-l','vta_hem-r','fx_hem','pag-l','pag_hem-r',
                        'ppn_hem-l','ppn_hem-r','cl_hem-l','cl_hem-r']

    gp_labels = ['gp_hem-l','gp_hem-r']

    # iterate over new segmentations, infiles, r2star maps
    for seg, i, r1 in zip(segm_files, infile_directories, r1corr_directories):

        print('Working on', i)

        r1corr_std = f'{output_dir}{i}/{i}_r1corr_std.nii'       
        transfile = f'{output_dir}{i}/{i}_transform_mat.txt'

        # perform coregistration on files 
        coregister(segm_path, infile_path, r1corr_path, seg, i, r1, output_dir, r1corr_std, transfile)

        # load in binary map
        bin_map_file = nib.load(segm_path + seg)
        bin_map = bin_map_file.get_fdata()
        regions = np.unique(bin_map)
        
        # z = bin_map.shape[2]
        count = 0 
        
        # iterate over numbers (aka regions) 
        for r in regions[1:]:

            region_name = labels_structures[int(r)]
            # if 'str' not in region_name:
            #     continue 
            
            if create_gp:
                if 'gp' not in region_name:
                    continue 

                count += 1
                if count > 2:
                    break

            # because the mask is altered, the orignal mask needs to be copied
            working_map = np.copy(bin_map)
        
            # find indices that do not beling to region and set those to 0
            inv_region = np.where(bin_map != r)
            working_map[inv_region] = 0
            region = np.where(bin_map == r)
            working_map[region] = 1
                
            # add globus pallidus externa to region as well
            if create_gp:
                region = np.where(bin_map == r+2)
                working_map[region] = 1
            
            # change back into nifti image (needed for nighres)
            nifti_image = nib.Nifti1Image(dataobj=working_map, header=bin_map_file.header, affine=bin_map_file.affine)

            # choose name for correct levelset outfile
            if create_gp:
                levelset_outfile = f'{i}_mask-{gp_labels[count-1]}_lvlreg-gt_def-img.nii'
            else:            
                levelset_outfile = f'{i}_mask-{region_name}_lvlreg-gt_def-img.nii'

            # compute distance map from that region 
            levelset = nighres.surface.probability_to_levelset(nifti_image, save_data=True, output_dir=output_dir, file_name=levelset_outfile)

            # transform segmented file
            os.system(f'flirt -in {levelset["result"]} -ref {r1corr_std} -out {output_dir+levelset_outfile} -init {transfile} -applyxfm')
            
            # remove file created with nighres
            os.remove(levelset['result'])

            # load created levelset
            levelset_file = nib.load(f'{output_dir}{levelset_outfile}.gz')
            levelset = levelset_file.get_fdata()

            # sometimes gzip doesn't work, so remove gzp alltogether 
            os.remove(f'{output_dir}{levelset_outfile}.gz')

            # flip first and last axes and crop image
            # flipped = np.flip(levelset, axis=(2))
            flipped = levelset 
            z = flipped.shape[2]
            flipped = flipped[:,:,z//2-25:z//2+25]

            # save cropped and flipped image
            nifti_image = nib.Nifti1Image(dataobj=flipped, header=levelset_file.header, affine=levelset_file.affine)
            nib.save(nifti_image, '../../data/segm/'+levelset_outfile)
            nib.save(nifti_image, f'{output_dir}{i}/segm/{levelset_outfile}')


def load_files(infile_path, segm_path, r1corr_path):
    ''' 
    Load files needed for creation of levelmaps.

    Args:
        infile_path: 
        segm_path:
        r1corr_path:
        r2star_path:
    
    Returns:
        infile_directories: loaded scans needed for coregistration
        segm_files: loaded files with binary segmentations
        r1corr_directories: 

    '''
   
    # gather all files 
    segm_files = [f[2] for f in os.walk(segm_path)][0]
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

    return segm_files, infile_directories, r1corr_directories


def coregister(segm_path, infile_path, r1corr_path, seg, i, r1, output_dir, r1corr_std, transfile):
    '''
    Perform coregistration of two scans. Needed to transform levelset later on.
    
    Args:
        segm_path: path to binary segmentation files
        infile_path: path to scans used for coregistration
        r1corr_path: path to scans to coregister to
        seg: current segmentation file
        i: current input file for coregistration
        r1: current file to coregister to
        output_dir: directory where to save coregistered files 
        r1corr_std: filename to save R1 corrected reorented file     
    '''

    infile = f'{infile_path}{i}/ses-1/anat/wb/qmri/{i}_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
    r1corr = f'{r1corr_path}{r1}/nii/r1corr.nii'
    outfile_r1corr = f'{output_dir}{i}/{i}_r1corr_coregistered.nii'

    if not os.path.exists(f'{output_dir}{i}/segm/'):
        os.makedirs(f'{output_dir}{i}/segm/')

    # orientate infile to be in standard orientation
    os.system(f'fslreorient2std {r1corr} {r1corr_std}')

    # perform coregistration if this has not been done yet
    if not os.path.exists(outfile_r1corr+'.gz'):
        os.system(f'flirt -in {infile} -ref {r1corr_std} -out {outfile_r1corr} -omat {transfile}')

    os.system(f'flirt -in {segm_path + seg} -ref {r1corr_std} -out {output_dir}{i}/{i}_binary_map.nii -init {transfile} -applyxfm')

    coreg_bin_map_file = nib.load(f'{output_dir}{i}/{i}_binary_map.nii.gz')
    coreg_bin_map = coreg_bin_map_file.get_fdata()
    coreg_bin_map = coreg_bin_map.astype(np.int64)
    
    z = coreg_bin_map.shape[2]

    # flip first and last axes and crop image
    flipped = np.flip(coreg_bin_map, axis=(2))
    flipped = flipped[:,:,z//2-25:z//2+25]

    # save cropped and flipped image
    nifti_image = nib.Nifti1Image(dataobj=flipped, header=coreg_bin_map_file.header, affine=coreg_bin_map_file.affine)
    nib.save(nifti_image, f'{output_dir}{i}/{i}_binary_map.nii.gz')


if __name__ == '__main__':
    create_segmentation(create_gp=True)