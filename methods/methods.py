import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt 
from scipy import ndimage

def main():
    ''' Create overview of available AHEAD data. '''


    inv1_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/inv1_te1_m_corr.nii'
    inv2_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/inv2_te1_m_corr.nii'
    t1corr_path = '../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/t1corr.nii'
    face_allm_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/Subcortex_0004_002_R02_allm.nii'

    no_face_inv1 = f"../defacing/{inv1_path.split('/')[-1]}"
    no_face_inv2 = f"../defacing/{inv2_path.split('/')[-1]}"

    brain_mask_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/mask_inv2_te2_m_corr.nii'
    test_path = '../defacing/001.nii'
    t1_mask_path = '../../../../data/projects/ahead/raw_gdata/Subcortex_0004_002_R02/nii/minv2_te2_m_corr.nii'

    brain_mask = nib.load(brain_mask_path).get_fdata()
    inv1 = nib.load(inv1_path).get_fdata()
    inv2 = nib.load(inv2_path).get_fdata()
    t1corr = nib.load(t1corr_path)
    t1corr_data = t1corr.get_fdata()

    img = brain_mask

    fig, ax = plt.subplots(1,4)
    ax[0].imshow(ndimage.rotate(inv1[:,:,64],90), cmap='gray')
    ax[0].axis('off')
    ax[0].set_title('First inversion')
    ax[1].imshow(ndimage.rotate(inv2[:,:,64],90), cmap='gray')
    ax[1].axis('off')
    ax[1].set_title('Second inversion')
    ax[2].imshow(ndimage.rotate(t1corr_data[:,:,64],90), cmap='gray')
    ax[2].axis('off')
    ax[2].set_title('T1-weigthed image')
    ax[3].imshow(ndimage.rotate(brain_mask[:,:,64],90), cmap='gray')
    ax[3].axis('off')
    ax[3].set_title('Brain mask')

    plt.show()

if __name__ == '__main__':
    main()