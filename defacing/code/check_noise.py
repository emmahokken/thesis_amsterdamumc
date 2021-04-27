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
    save_path = prefix + '/data/projects/ahead/defaced/'
    
    subjects = ['0004', '0016', '0083', '0067']

    subj = 'Subcortex_0083_085_R02'
    directories = [d[1] for d in  os.walk(save_path)][0]

    # caclulate for different coils
    coils = [7,15,23,31]
    n = 50
    s = 140

    overall_difference_real = {1:[], 2:[], 3:[], 4:[]}
    overall_difference_imag = {1:[], 2:[], 3:[], 4:[]}
    overall_percent_real = {1:[], 2:[], 3:[], 4:[]}
    overall_percent_imag = {1:[], 2:[], 3:[], 4:[]}

    for d in tqdm(directories):
        # print(d)
        for echo in range(1,5):
            original = nib.load(f'{root_path}{d}/{d}_inv2_{echo}_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
            defaced = nib.load(f'{save_path}{d}/{d}_inv2_{echo}_gdataCorrected_defaced.nii.gz').get_fdata(dtype=np.complex64)
            # no_noise = nib.load(f'{save_path}{d}/{d}_inv2_{echo}_gdataCorrected_no_noise.nii.gz').get_fdata(dtype=np.complex64)

            x, y, z = original.shape[:3]       

            difference_real = []
            difference_imag = []
            percent_real = []
            percent_imag = []

            for coil in range(original.shape[3]):
                # print('Coil:', coil + 1)

                # compute std for both imageinary and real channels from original and defaced 
                std_real_orig = np.std(original[x-n:,y-n:,z-n:, coil].real)
                std_imag_orig = np.std(original[x-n:,y-n:,z-n:, coil].imag)

                std_real_def = np.std(defaced[x-n:,y-n:,z-n:, coil].real)
                std_imag_def = np.std(defaced[x-n:,y-n:,z-n:, coil].imag)

                # print('Original:')
                # print('Real:', std_real_orig, 'Imaginary:', std_imag_orig)
                # print('Defaced:')
                # print('Real:', std_real_def, 'Imaginary:', std_imag_def)

                difference_real.append(std_real_def - std_real_orig)
                difference_imag.append(std_imag_def - std_imag_orig)

                overall_difference_real[echo].append(std_real_def - std_real_orig)
                overall_difference_imag[echo].append(std_imag_def - std_imag_orig)

                p_real = ((std_real_def - std_real_orig) / ((std_real_def + std_real_orig) / 2)) * 100
                p_imag = ((std_imag_def - std_imag_orig) / ((std_imag_def + std_imag_orig) / 2)) * 100
                percent_real.append(p_real)
                percent_imag.append(p_imag)

                overall_percent_real[echo].append(p_real)
                overall_percent_imag[echo].append(p_imag)

            # plot_difference_bar(overall_difference_real[1], overall_difference_imag[1], echo)

       
        # print(np.mean(percent_real), np.mean(percent_imag))

    for echo in range(1,5):
        plot_difference_bar(overall_difference_real[echo], overall_difference_imag[echo], echo)
        plot_difference_bar(overall_percent_real[echo], overall_percent_imag[echo], echo, percent=True)

        print('Overall difference mean')
        print('Real:', np.mean(overall_percent_real[echo]), 'Imaginary:', np.mean(overall_percent_imag[echo]))


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

def plot_difference_bar(difference_real, difference_imag, echo, percent=False):
    ''' Plots a bar graph showing the difference in background noise for the orignal and the defaced image. '''

    x = np.arange(len(difference_real))
    w = 0.5

    fig, ax = plt.subplots()
    fig.set_figwidth(10)
    rec1 = ax.bar(x - w/2, difference_real, w, label='Real', color='rebeccapurple')
    rec2 = ax.bar(x + w/2, difference_imag, w, label='Imaginary', color='orange')

    ax.set_xticks(x)
    ax.set_xlabel('Coils')

    ax.legend()
    if percent:
        ax.set_ylabel('% difference')
        plt.savefig(f'../results/figures/background_noise_percent_difference_{echo}.pdf')
    else:
        ax.set_ylabel('STD Difference')
        plt.savefig(f'../results/figures/background_noise_std_difference_{echo}.pdf')

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