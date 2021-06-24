import nibabel as nib 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
from scipy import ndimage, stats

def main():

    # plot_noise_boundary()
    
    # load in data
    df = pd.read_csv('../results/files/defacing_background_noise.csv')
    print(df.columns)
    print('Mean difference', df.perc_real.mean(),df.perc_imag.mean())
    print('Significant difference between real and imaginary?', stats.ttest_ind(df.perc_real, df.perc_imag))
    print('Significant difference between original and defaced', stats.ttest_ind(df.std_real_orig, df.std_real_def))
    print('Significant difference between original and defaced', stats.ttest_ind(df.std_imag_orig, df.std_imag_def))
    print('Levene', stats.levene(df.std_real_orig, df.std_real_def))
    print('Levene part 2', stats.levene(df.std_imag_orig, df.std_imag_def))
   
    # iterate over coils
    for e in df.echo.unique():
        differences = pd.DataFrame()
        for c in df.coil.unique():
            data = df[(df.coil == c) & (df.echo == e)]
            num_subj = data.subj.unique().shape[0]
            diff_real = data.diff_real.sum() / num_subj
            diff_imag = data.diff_imag.sum() / num_subj
            perc_real = data.perc_real.sum() / num_subj
            perc_imag = data.perc_imag.sum() / num_subj

            differences = differences.append({'diff_real':diff_real, 'diff_imag':diff_imag, 'perc_real':perc_real, 'perc_imag': perc_imag}, ignore_index=True)

        # plot the two versions of the differences, percentual and real valued 
        plot_difference_bar(differences.diff_real, differences.diff_imag, e, percent=False)
        plot_difference_bar(differences.perc_real, differences.perc_imag, e, percent=True)


def plot_difference_bar(difference_real, difference_imag, echo, percent=False):
    ''' Plots a bar graph showing the difference in background noise for the orignal and the defaced image. '''

    x = np.arange(1, difference_real.shape[0]+1)
    w = 0.35

    fig, ax = plt.subplots()
    fig.set_figwidth(10)
    rec1 = ax.bar(x - w/2, difference_real, w, label='Real', color='rebeccapurple', edgecolor='k', linewidth=0.4)
    rec2 = ax.bar(x + w/2, difference_imag, w, label='Imaginary', color='orange', edgecolor='k', linewidth=0.4)

    me = pd.concat([difference_real, difference_imag])
    ax.axhline(me.mean(), label='Combined mean', color='cornflowerblue')

    ax.set_xlabel('Coils')
    ax.legend()
    if percent:
        # create line and ytick for mean 
        ax2 = ax.twinx()
        ax2.set_yticks([me.mean()])
        ax2.set_ylim(0,100)

        ax.set_xticks(x)
        ax.set_title(f'Difference in noise (in percent) for echo time {echo}')
        ax.set_yticks(range(0,101,10))
        ax.set_ylabel('STD difference in percent')
        plt.savefig(f'../results/figures/background_noise_percent_difference_{echo}.pdf')
    else:
        ax.set_title(f'Difference in noise for echo time {echo}')
        ax.set_ylabel('STD Difference')
        plt.savefig(f'../results/figures/background_noise_std_difference_{echo}.pdf')

    plt.show()

def plot_noise_boundary():
    prefix = '../../../../../..'
    root_path = prefix + '/data/projects/ahead/raw_gdata/'
    save_path = prefix + '/data/projects/ahead/defaced/'
    
    subj = 'Subcortex_0083_085_R02'
    for echo in range(1,5):
        original = nib.load(f'{root_path}{subj}/{subj}_inv2_{echo}_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
        defaced = nib.load(f'{save_path}{subj}/{subj}_inv2_{echo}_gdataCorrected_defaced.nii.gz').get_fdata(dtype=np.complex64)
        no_noise = nib.load(f'{save_path}{subj}/{subj}_inv2_{echo}_gdataCorrected_no_noise.nii.gz').get_fdata(dtype=np.complex64)

        fig, ax = plt.subplots(1,3)
        ax[0].imshow(ndimage.rotate(original[140,:,:,2].real,90), cmap='gray', vmin=-10000, vmax=0)
        ax[0].set_title('Origiinal image')
        ax[0].axis('off')
        ax[1].imshow(ndimage.rotate(defaced[140,:,:,2].real,90), cmap='gray', vmin=-10000, vmax=0)
        ax[1].set_title('Defaced image \nwith added noise')
        ax[1].axis('off')
        ax[2].imshow(ndimage.rotate(no_noise[140,:,:,2].real,90), cmap='gray', vmin=-10000, vmax=0)
        ax[2].set_title('Defaced image \nwithout added noise')
        ax[2].axis('off')
        plt.tight_layout()
        plt.savefig(f'../results/figures/noise_boundary_comparison_{echo}.pdf')
        # plt.show()

        # breakpoint()

if __name__ == '__main__':
    main()

