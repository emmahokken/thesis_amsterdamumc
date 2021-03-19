import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import nibabel as nib 
import numpy as np 
import os
from scipy import ndimage
import sys 
from tqdm import tqdm 

def check_noise(plot=False):
    print(plot)

    prefix = '../../../../../..'
    root_path = prefix + '/data/projects/ahead/raw_gdata/'
    save_path = '../results/dilation/'

    subjects = ['0004', '0016', '0083', '0067']

    subj = 'Subcortex_0083_085_R02'
    directories = [d[1] for d in  os.walk(save_path)][0]

    plot = False 

    # caclulate for different coils
    coils = [7,15,23,31]
    x = 50
    y = 50
    z = 50
    s = 140

    for d in directories:
        print(d)
        original = nib.load(f'{root_path}{d}/{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
        defaced = nib.load(save_path + d + f'/dilated_{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
        no_noise = nib.load(save_path + d + f'/dilated_no_noise_{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)

        for coil in coils:
            print('Coil:', coil + 1)
            # compute std for both imageinary and real channels from original and defaced 
            std_real_orig = np.std(original[:x,:y,:z, coil].real)
            std_imag_orig = np.std(original[:x,:y,:z, coil].imag)

            std_real_def = np.std(defaced[:x,:y,:z, coil].real)
            std_imag_def = np.std(defaced[:x,:y,:z, coil].imag)

            print('Original:')
            print('Real:', std_real_orig, 'Imaginary:', std_imag_orig)
            print('Defaced:')
            print('Real:', std_real_def, 'Imaginary:', std_imag_def)


            if plot:
                plt.subplot(131)
                plt.imshow(ndimage.rotate(original[149,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
                plt.title('Original image')
                plt.subplot(132)
                plt.imshow(ndimage.rotate(no_noise[140,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
                plt.title('Defaced, without noise')
                plt.subplot(133)
                plt.imshow(ndimage.rotate(defaced[140,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
                plt.title('Defaced, with noise')

                plt.colorbar(shrink=0.5)
                plt.show()

if __name__ == '__main__':
    check_noise(plot=True)