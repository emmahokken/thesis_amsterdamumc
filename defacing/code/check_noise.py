import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import nibabel as nib 
import numpy as np 
import os
from scipy import ndimage
import sys 
from tqdm import tqdm 

import pandas as pd 

def check_noise(plot=False):
    ''' Iterates over all defaced files and compute the background noise of a specific region. 

    Args: 
        plot: boolean determining whether to plot the noise box. 
    '''
    prefix = '../../../../../..'
    root_path = prefix + '/data/projects/ahead/raw_gdata/'
    save_path = prefix + '/data/projects/ahead/defaced/'
   
    # caclulate for different coils
    n = 50
    s = 140
    df = pd.DataFrame(columns=['subj', 'echo', 'coil', 'diff_real', 'diff_imag','perc_real','perc_imag'])

    for d in tqdm(directories):
     
        for echo in range(1,5):
            data = {'subj': np.full(32,int(d.split('_')[1])), 'echo':np.full(32,echo), 'coil': range(1,33)}
            
            try: 
                original = nib.load(f'{root_path}{d}/{d}_inv2_{echo}_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
                defaced = nib.load(f'{save_path}{d}/{d}_inv2_{echo}_gdataCorrected_defaced.nii.gz').get_fdata(dtype=np.complex64)
                no_noise = nib.load(f'{save_path}{d}/{d}_inv2_{echo}_gdataCorrected_no_noise.nii.gz').get_fdata(dtype=np.complex64)
            except:
                print(d)
                continue

            x, y, z = original.shape[:3]       

            difference_real = []
            difference_imag = []
            percent_real = []
            percent_imag = []
            all_std_real_orig = []
            all_std_imag_orig = []
            all_std_real_def = []
            all_std_imag_def = []
            

            for coil in range(original.shape[3]):
                # print('Coil:', coil + 1)

                # compute std for both imageinary and real channels from original and defaced 
                std_real_orig = np.std(original[0:n,0:n,z-n:, coil].real)
                std_imag_orig = np.std(original[0:n,0:n,z-n:, coil].imag)

                std_real_def = np.std(defaced[0:n,0:n,z-n:, coil].real)
                std_imag_def = np.std(defaced[0:n,0:n,z-n:, coil].imag)

                all_std_real_orig.append(std_real_orig)
                all_std_imag_orig.append(std_imag_orig)
                all_std_real_def.append(std_real_def)
                all_std_imag_def.append(std_imag_def)

                difference_real.append(std_real_def - std_real_orig)
                difference_imag.append(std_imag_def - std_imag_orig)

                p_real = ((std_real_def - std_real_orig) / ((std_real_def + std_real_orig) / 2)) * 100
                p_imag = ((std_imag_def - std_imag_orig) / ((std_imag_def + std_imag_orig) / 2)) * 100
                percent_real.append(p_real)
                percent_imag.append(p_imag)


            data['diff_real'] = difference_real
            data['diff_imag'] = difference_imag
            data['perc_real'] = percent_real
            data['perc_imag'] = percent_imag

            data['std_real_orig'] = all_std_real_orig
            data['std_imag_orig'] = all_std_imag_orig
            data['std_real_def'] = all_std_real_def
            data['std_imag_def'] = all_std_imag_def


            df = df.append(pd.DataFrame(data),ignore_index=True)

    df.to_csv('defacing_background_noise_adjusted_box.csv')


def plot_defacing_comparison(original, no_noise, defaced, s, coil=7):
    ''' Plots the original, defaced, and defaced with no noise images. '''

    plt.subplot(131)
    plt.imshow(ndimage.rotate(original[s,:,:,coil].real, 90), vmin=-100, vmax=100, cmap='gray')
    plt.title('Original image')
    plt.subplot(132)
    plt.imshow(ndimage.rotate(no_noise[s,:,:,coil].real, 90), vmin=-100, vmax=100, cmap='gray')
    plt.title('Defaced, without noise')
    plt.subplot(133)
    plt.imshow(ndimage.rotate(defaced[s,:,:,coil].real, 90), vmin=-100, vmax=100, cmap='gray')
    plt.title('Defaced, with noise')

    plt.colorbar(shrink=0.5)
    plt.show()


def plot_noise_box(original, n, coil=7):
    ''' Plots a box from which the noise in compared. '''

    x, y, z = original.shape[:3]

    fig, ax = plt.subplots()
    ax.imshow(original[s,:,:,coil].real, cmap='gray')
    rect = patches.Rectangle((y-n,z-n), n, n, linewidth=1, edgecolor='lime', facecolor='none')
    ax.add_patch(rect)
    ax.axis('off')
    plt.show()


if __name__ == '__main__':
    check_noise(plot=False)