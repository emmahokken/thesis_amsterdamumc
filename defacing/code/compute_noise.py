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

for d in tqdm(directories): 

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

    # iterate over the four echo times
    for echo in range(1,5):
        file = f'{root_path}{d}/{d}_inv2_{echo}_gdataCorrected.nii.gz'
       
        # load specific echo file 
        data_load = nib.load(file)
        data = data_load.get_fdata(dtype=np.complex64)

        defaced_image = []

        # iterate over coils 
        for coil in range(data.shape[3]):

            coil_data = data[:,:,:,coil]

            print(np.mean(coil_data[:15,:15,:15].real))
            print(mean_real[coil])

            # deface coil image 
            dilated_data = coil_data*dilation

            # generate new Gaussian noise from that distribution 
            noise_real = np.random.normal(mean_real[coil], std_real[coil], coil_data.shape)
            noise_imag = np.random.normal(mean_imag[coil], std_imag[coil], coil_data.shape)

            # convert back to image space 
            # imspace_noise_real = np.fft.ifft(noise_real, norm='ortho')
            # imspace_noise_imag = np.fft.ifft(noise_imag, norm='ortho')

            # make space for brain in noise
            noise_real = noise_real * inv_dilation
            noise_imag = noise_imag * inv_dilation

            # combine noise and coil image
            final_coil_image_real = dilated_data.real + noise_real
            final_coil_image_imag = dilated_data.imag + noise_imag
            
            defaced_image.append(final_coil_image_real + final_coil_image_imag)

            fig, ax = plt.subplots(1,2)

            ax[0].imshow((final_coil_image_real + final_coil_image_imag)[:15,:15,140].real, cmap='gray')
            ax[0].set_title('added noise')
            ax[0].axis('off')
            ax[1].imshow(coil_data[:15,:15,140].real, cmap='gray')
            ax[1].set_title('original noise')
            ax[1].axis('off')
            plt.show()
            plt.show()


        defaced_image = np.array(defaced_image)
        print('done', defaced_image.shape)
        fig, ax = plt.subplots(1,2)

        ax[0].imshow(data[:,:,140,2].real, cmap='gray')
        ax[0].set_title('defaced image')
        ax[0].axis('off')
        ax[1].imshow(defaced_image[:,:,140,2].real, cmap='gray')
        ax[1].set_title('defaced image')
        ax[1].axis('off')
        plt.show()

        exit()

        # save new scan 
        if not os.path.exists(f'{save_path}{d}'):
            os.makedirs(f'{save_path}/{d}')
        save_dir = f'{save_path}/{d}/dilated_{d}_inv2_{echo}_gdataCorrected.nii.gz'

        nifit_image = nib.Nifti1Image(dataobj=defaced_image, header=data_load.header, affine=data_load.affine)
        nib.save(img=nifit_image, filename=save_dir)


