import matplotlib.pyplot as plt
import nibabel as nib 
import numpy as np 
import os
import sys 

from dilation import dilate

# declare relevant prefixes and directories 
prefix = '../../../../../..'
root_path = prefix + '/data/projects/ahead/raw_gdata/'

# gather all directories 
directories = [d[1] for d in  os.walk(root_path)][0]

for d in directories: 

    with open(f'../data/{d}/noise_stats.txt', 'r') as f:
        # skip over first line 
        f.readline()

        # read file, column wise
        rows = [[x for x in line.split(',')] for line in f]
        cols = np.array([np.array(col, dtype=float) for col in zip(*rows)])

        mean_real, mean_imag, std_real, std_imag = [c for c in cols]

    dilation = dilate(f'{root_path}{d}')

    print(dilation.shape)

    # iterate over the four echo times
    for echo in range(1,5):
        file = f'{root_path}{d}/{d}_inv2_{echo}_gdataCorrected.nii.gz'
       
        # load specific echo file 
        data_load = nib.load(file)
        data = data_load.get_fdata(dtype=np.complex64)

        # transform to kspace, normalize by scale of 1/sqrt(n)
        # data_ifft = np.fft.ifftshift(data)
        # kspace = np.fft.fft2(data, norm='ortho')
        # kspace_ifft = np.fft.fftshift(np.fft.fft2(data_ifft, norm='ortho'))
        # kspace_n = np.fft.fftshift(np.fft.ifftn(data_ifft, norm='ortho'))

        # kn = np.fft.fftshift(np.fft.fftn(data, axes=(0,1,2), norm='ortho'))

        # plt.imshow(kn[:,:,100,0].real, cmap='gray')
        # plt.title('fftn')

        # plt.show()

        # exit()

        # iterate over coils 
        for coil in range(data.shape[3]):

            coil_data = data[:,:,:,coil]
            print(coil_data.shape)

            # deface coil image 
            dilated_data = coil_data*dilation

            print(dilation.shape)
            # generate new Gaussian noise from that distribution 
            noise_real = np.random.normal(mean_real[coil], std_real[coil], coil_data.shape)
            noise_imag = np.random.normal(mean_imag[coil], std_imag[coil], coil_data.shape)

            # convert back to image space 
            imspace_noise_real = np.fft.ifft(noise_real, norm='ortho')
            imspace_noise_imag = np.fft.ifft(noise_imag, norm='ortho')

            noise_per_coil.append([imspace_noise_real, imspace_noise_imag])

        # somehow combine 32 coils 
        # add noise to defaced data 



            exit()
    # save new scan 
    if not os.path.exists(f'{save_path}/{d}'):
        os.makedirs(f'{save_path}/{d}')
    save_dir = f'{save_path}/{d}/dilated_{scan_type}.nii'

    nifit_image = nib.Nifti1Image(dataobj=filled_in, header=inv2_nii.header, affine=inv2_nii.affine)
    nib.save(img=nifit_image, filename=save_dir)


