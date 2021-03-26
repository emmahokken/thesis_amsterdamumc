import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import nibabel as nib 
import numpy as np 
import os
from scipy import ndimage
import sys 
from tqdm import tqdm 

def check_noise(plot=False):
    prefix = '../../../../../..'
    root_path = prefix + '/data/projects/ahead/raw_gdata/'
    save_path = '../results/dilation/'
    # save_path = '../../../../Documents/defacing/results/'

    subjects = ['0004', '0016', '0083', '0067']

    subj = 'Subcortex_0083_085_R02'
    directories = [d[1] for d in  os.walk(save_path)][0]
    # print(directories)

    # caclulate for different coils
    coils = [7,15,23,31]
    n = 50
    s = 140

    for d in directories:
        print(d)
        original = nib.load(f'{root_path}{d}/{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
        defaced = nib.load(save_path + d + f'/dilated_{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
        no_noise = nib.load(save_path + d + f'/dilated_no_noise_{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)

        x, y, z = original.shape[:3]

        # fig, ax = plt.subplots()
        # ax.imshow(original[s,:,:,7].real, cmap='gray')
        # rect = patches.Rectangle((y-n,z-n), n, n, linewidth=1, edgecolor='lime', facecolor='none')
        # ax.add_patch(rect)
        # ax.axis('off')
        # plt.show()

        difference_real = []
        difference_imag = []
        percent_real = []
        percent_imag = []

        for coil in range(original.shape[3]):
            print('Coil:', coil + 1)
            # compute std for both imageinary and real channels from original and defaced 
            std_real_orig = np.std(original[x-n:,y-n:,z-n:, coil].real)
            std_imag_orig = np.std(original[x-n:,y-n:,z-n:, coil].imag)

            std_real_def = np.std(defaced[x-n:,y-n:,z-n:, coil].real)
            std_imag_def = np.std(defaced[x-n:,y-n:,z-n:, coil].imag)

            print('Original:')
            print('Real:', std_real_orig, 'Imaginary:', std_imag_orig)
            print('Defaced:')
            print('Real:', std_real_def, 'Imaginary:', std_imag_def)

            difference_real.append(std_real_def - std_real_orig)
            difference_imag.append(std_imag_def - std_imag_orig)

            p_real = ((std_real_def - std_real_orig) / ((std_real_def + std_real_orig) / 2)) * 100
            p_imag = ((std_imag_def - std_imag_orig) / ((std_imag_def + std_imag_orig) / 2)) * 100
            percent_real.append(p_real)
            percent_imag.append(p_imag)

            if plot:
                plt.subplot(131)
                plt.imshow(ndimage.rotate(original[s,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
                plt.title('Original image')
                plt.subplot(132)
                plt.imshow(ndimage.rotate(no_noise[s,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
                plt.title('Defaced, without noise')
                plt.subplot(133)
                plt.imshow(ndimage.rotate(defaced[s,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
                plt.title('Defaced, with noise')

                plt.colorbar(shrink=0.5)
                plt.show()
            
        # plt.hist(difference_real, bins=len(difference_real))
        # plt.show()

        x = np.arange(len(difference_real))
        w = 0.35

        fig, ax = plt.subplots()
        rec1 = ax.bar(x - w/2, difference_real, w, label='Real STD')
        rec2 = ax.bar(x + w/2, difference_imag, w, label='Imaginary STD')

        ax.set_xticks(x)
        ax.set_ylabel('STD Difference')
        ax.set_xlabel('Coils')

        ax.legend()
        plt.show()

        x = np.arange(len(percent_real))
        w = 0.35

        fig, ax = plt.subplots()
        rec1 = ax.bar(x - w/2, percent_real, w, label='Real')
        rec2 = ax.bar(x + w/2, percent_imag, w, label='Imaginary')

        ax.set_xticks(x)
        ax.set_ylabel('% difference')
        ax.set_xlabel('Coils')

        ax.legend()
        plt.show()

        print(np.mean(percent_real), np.mean(percent_imag))

if __name__ == '__main__':
    check_noise(plot=False)