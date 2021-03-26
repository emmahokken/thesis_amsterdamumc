import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import nibabel as nib 
import numpy as np 
import os
from scipy import ndimage
import sys 
from tqdm import tqdm 

from dilation import dilate


def add_noise(r_outer, r_inner, plot=False):
        
    # declare relevant prefixes and directories 
    prefix = '../../../../../..'
    root_path = prefix + '/data/projects/ahead/raw_gdata/'
    save_path = '../results/dilation/'

    # gather all directories 
    directories = [d[1] for d in  os.walk(root_path)][0]
    subjects = ['0004', '0016', '0067', '0083']
    for d in tqdm(directories): 
        print(d)
        # if not any(x in d for x in subjects):
        #     continue

        # get dilated brain mask for defacing and reorient axes
        dilation, inv_dilation = dilate(f'{root_path}{d}')
        dilation = np.moveaxis(dilation, 2, 0)
        inv_dilation = np.moveaxis(inv_dilation, 2, 0)

        dilation = np.fliplr(dilation)
        inv_dilation = np.fliplr(inv_dilation)

        # iterate over the four echo times
        for echo in range(1,5):
            print('Echo:', echo)
            file = f'{root_path}{d}/{d}_inv2_{echo}_gdataCorrected.nii.gz'
        
            # load specific echo file 
            data_load = nib.load(file)
            data = data_load.get_fdata(dtype=np.complex64)

            kspace = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(data), axes=(0,1,2)))

            defaced_image = []
            no_noise = []

            s = 291
            indices = create_kspace_mask(data, kspace, r_outer, r_inner, s, plot=plot)
    
            # iterate over coils 
            for coil in range(data.shape[3]):
                # print('Coil:', coil)

                coil_data = data[:,:,:,coil]

                # deface coil image 
                dilated_data = coil_data*dilation
                without_noise = dilated_data.copy()
                no_noise.append(without_noise)

                # imspace = np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(kspace[x_mid-25:x_mid+25,y_mid-25:y_mid+25,290,coil]), norm='ortho'))
                # calculate std for noise generation 
                std_real = np.std(kspace[indices[0], indices[1], s, coil].real)
                std_imag = np.std(kspace[indices[0], indices[1], s, coil].imag)
                
                # print('STD real:', std_real)
                # print('STD imag:', std_imag)
                # generate new Gaussian noise from that distribution and make space for brain in noise
                noise_real = np.random.normal(0, std_real, coil_data.shape) * inv_dilation
                noise_imag = np.random.normal(0, std_imag, coil_data.shape) * inv_dilation

                # combine noise and coil image
                dilated_data.real += noise_real
                dilated_data.imag += noise_imag
                
                defaced_image.append(dilated_data)
            
            # print('Done with coils.')

            defaced_image = np.moveaxis(np.array(defaced_image), 0, 3)
            no_noise = np.moveaxis(np.array(no_noise), 0, 3)

            # save new scan 
            if not os.path.exists(f'{save_path}{d}'):
                os.makedirs(f'{save_path}/{d}')
            save_dir = f'{save_path}/{d}/dilated_{d}_inv2_{echo}_gdataCorrected.nii.gz'

            nifit_image = nib.Nifti1Image(dataobj=defaced_image, header=data_load.header, affine=data_load.affine)
            nib.save(img=nifit_image, filename=save_dir)

            save_dir = f'{save_path}/{d}/dilated_no_noise_{d}_inv2_{echo}_gdataCorrected.nii.gz'

            nifit_image = nib.Nifti1Image(dataobj=no_noise, header=data_load.header, affine=data_load.affine)
            nib.save(img=nifit_image, filename=save_dir)



def create_kspace_mask(data, kspace, r_outer, r_inner, s, plot=False):
    ''' Create the mask used to gather noise statistics from kspace. 
    
    Args: 
        data: MRI data
        outer_h, outer_w, inner_h, inner_w: measurements of the box that is used to calculate STD
        s: slice to plot
        plot: whether to plot the box, default set to False
    
    Returns: 
        indices: indices of the kspace mask  ''' 

    # determine middle for noise generation 
    x_mid = data.shape[0] // 2
    y_mid = data.shape[1] // 2

    # create mask for ring of kspace for noise 
    h, w = data.shape[:2]
    y, x = np.ogrid[:h, :w]

    circ_mask = np.sqrt((x-y_mid)**2 + (y-x_mid)**2) < r_outer
    circ_mask_inner = np.sqrt((x- y_mid)**2 + (y- x_mid)**2) < r_inner
    circ_mask[circ_mask_inner] = 0 

    indices = np.where(circ_mask > 0)

    if plot:
        # plot box around region
        fig, ax = plt.subplots()
        ax.imshow(ndimage.rotate(kspace[:,:,s,8].real, 90), cmap='gray')
        circ_outer = patches.Circle((x_mid,y_mid), radius=r_outer, linewidth=1, edgecolor='lime', facecolor='none')
        ax.add_patch(circ_outer)
        circ_inner = patches.Circle((x_mid,y_mid), radius=r_inner, linewidth=1, edgecolor='c', facecolor='none')
        ax.add_patch(circ_inner)
        ax.axis('off')
        plt.show()

    return indices


if __name__ == '__main__':

 
    r_outer = 117
    r_inner = 90

    add_noise(r_outer, r_inner, plot=False)