import matplotlib.pyplot as plt
import nibabel as nib 
import numpy as np 
import os
import sys 
from tqdm import tqdm 

from dilation import dilate

# declare relevant prefixes and directories 
prefix = '../../../../../..'
root_path = prefix + '/data/projects/ahead/raw_gdata/'
save_path = '../results/dilation/'

# gather all directories 
directories = [d[1] for d in  os.walk(root_path)][0]

for d in tqdm(directories[20:]): 

    with open(f'../data/{d}/noise_stats.txt', 'r') as f:
        # skip over first line 
        f.readline()

        # read file, column wise
        rows = [[x for x in line.split(',')] for line in f]
        cols = np.array([np.array(col, dtype=float) for col in zip(*rows)])

        mean_real, mean_imag, std_real, std_imag = [c for c in cols]

    # get dilated brain mask for defacing and reorient axes
    dilation, inv_dilation = dilate(f'{root_path}{d}')
    dilation = np.moveaxis(dilation, 2, 0)
    inv_dilation = np.moveaxis(inv_dilation, 2, 0)

    # repeat array n of coils times so we don't have to iterate over each coil seperately
    dilation = np.repeat(dilation[:,:,:,np.newaxis], 32, axis=3)
    inv_dilation = np.repeat(inv_dilation[:,:,:,np.newaxis], 32, axis=3)
    # print(dilation.shape, inv_dilation.shape)

    # iterate over the four echo times
    for echo in range(1,5):
        file = f'{root_path}{d}/{d}_inv2_{echo}_gdataCorrected.nii.gz'
       
        # load specific echo file 
        data_load = nib.load(file)
        data = data_load.get_fdata(dtype=np.complex64)

        # defaced_image = []

        # # iterate over coils 
        # for coil in range(data.shape[3]):

        #     coil_data = data[:,:,:,coil]

        # deface coil image 
        dilated_data = data*dilation

        # generate new Gaussian noise from that distribution and make space for brain in noise
        noise_real =  np.moveaxis([np.random.normal(mean_real[coil], std_real[coil], data.shape[:-1]) for coil in range(data.shape[3])], 0, 3) * inv_dilation
        # print(noise_real.shape)
        noise_imag = np.moveaxis([np.random.normal(mean_imag[coil], std_imag[coil], data.shape[:-1]) for coil in range(data.shape[3])], 0, 3) * inv_dilation

        # combine noise and coil image
        final_image_real = dilated_data.real + noise_real
        final_image_imag = dilated_data.imag + noise_imag
        
        final_image = final_image_real + final_image_imag
        # print(final_image.shape)

        # plt.imshow(final_image[:,:,140,3].real, cmap='gray')
        # plt.show()

        # defaced_image = np.moveaxis(np.array(defaced_image), 0, 3)
       
        # save new scan 
        if not os.path.exists(f'{save_path}{d}'):
            os.makedirs(f'{save_path}/{d}')
        save_dir = f'{save_path}/{d}/dilated_{d}_inv2_{echo}_gdataCorrected.nii.gz'

        nifit_image = nib.Nifti1Image(dataobj=final_image, header=data_load.header, affine=data_load.affine)
        nib.save(img=nifit_image, filename=save_dir)


