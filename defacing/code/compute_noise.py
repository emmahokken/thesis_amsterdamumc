import matplotlib.pyplot as plt
import nibabel as nib 
import numpy as np 
import os
import sys 

# load relevant files
subj = 5
prefix = '../../../../../..'
root_path = prefix + '/data/projects/ahead/raw_gdata/'
basefile = f'/data/projects/ahead/raw_gdata/Subcortex_{subj:04}_002_R02/'
# infile = f'/data/projects/ahead/segmentations2/qMRI-Recomputed/sub-{subj:03}/ses-1/anat/wb/qmri/sub-{subj:03}_ses-1_acq-wb2_mod-r1hz_orient-std_brain.nii.gz'
# reffile = f'/data/projects/ahead/raw_gdata/Subcortex_{subj:04}_002_R02/nii/r1corr.nii'
# outfile = f'data/coregistration/sub-{subj:03}_outfile.nii'
# trans_mat = f'data/coregistration/sub-{subj:03}_transform_mat.txt'

# gather all directories 
directories = [d[1] for d in  os.walk(root_path)][0]

for d in directories: 
    # iterate over the four echo times
    for echo in range(1,5):
        file = f'{root_path}{d}/{d}_inv2_{echo}_gdataCorrected.nii.gz'

        # load specific echo file 
        data_load = nib.load(file)
        data = data_load.get_fdata(dtype=np.complex64)

        # transform to kspace, normalize by scale of 1/sqrt(n)
        kspace = np.fft.fft(data, norm='ortho')
    
        noise_per_coil = []

        # iterate over coils 
        for coil in range(kspace.shape[3]):

            coil_data = kspace[:,:,:,coil]

            cube_front = 15
            cube_back = 16

            # select 15*15*15 cube from the 8 corners of k-space
            cube1 = coil_data[0:cube_front,0:cube_front,0:cube_front]       # 000
            cube2 = coil_data[0:cube_front,0:cube_front,-cube_back:-1]      # 001
            cube3 = coil_data[0:cube_front,-cube_back:-1:,-cube_back:-1]    # 011
            cube4 = coil_data[-cube_back:-1:,-cube_back:-1,-cube_back:-1]   # 111
            cube6 = coil_data[-cube_back:-1,-cube_back:-1,0:cube_front]     # 110
            cube5 = coil_data[-cube_back:-1,0:cube_front,0:cube_front]      # 100
            cube7 = coil_data[0:cube_front,-cube_back:-1,0:cube_front]      # 010
            cube8 = coil_data[-cube_back:-1,0:cube_front,-cube_back:-1]     # 101

            # concatinate cubes 
            cubes = np.concatenate((cube1, cube2, cube3, cube4, cube5, cube6, cube7, cube8))
            
            # calculate mean and std 
            mean_real = np.mean(cubes.real)
            std_real = np.std(cubes.real)

            mean_imag = np.mean(cubes.imag)
            std_imag = np.std(cubes.imag)

            # generate new Gaussian noise from that distribution 
            noise_real = np.random.normal(mean_real, std_real, coil_data.shape)
            noise_imag = np.random.normal(mean_imag, std_imag, coil_data.shape)

            noise_per_coil.append([noise_real, noise_imag])
        # somehow combine 32 coils 
        # add noise to defaced data 



        exit()