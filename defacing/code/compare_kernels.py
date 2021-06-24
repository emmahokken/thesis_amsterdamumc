import cv2 
import matplotlib.patches as patches 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres
import numpy as np 
import os 
from scipy import ndimage

from defacing import deface
from dilation import dilate 

def main():
    '''
    Plot various figures to investigate dilation results.
    '''

    # decalre root path and scan type
    root_path = '../../../../../data/projects/ahead/raw_gdata'
    # scan_type = 'inv1_te1_m_corr'
    scan_type = 't1corr'
    brain_mask_type = 'mask_inv2_te2_m_corr'
    save_path = '../defacing/dilation'
    subject_5 = 'Subcortex_0005_002_R02'
    subject_8 = 'Subcortex_0008_006_R02'
    subject_64 = 'Subcortex_0064_064_R02'

    kernel_five = np.ones((5,5), dtype=np.uint8)
    kernel_fifty = np.ones((50,50), dtype=np.uint8)

    brain_mask_path_5 = f'{root_path}/{subject_5}/nii/{brain_mask_type}.nii'
    scan_path_5 = f'{root_path}/{subject_5}/nii/{scan_type}.nii'
    brain_mask_path_8 = f'{root_path}/{subject_8}/nii/{brain_mask_type}.nii'
    scan_path_8 = f'{root_path}/{subject_8}/nii/{scan_type}.nii'
    brain_mask_path_64 = f'{root_path}/{subject_64}/nii/{brain_mask_type}.nii'
    scan_path_64 = f'{root_path}/{subject_64}/nii/{scan_type}.nii'

    # load in brain mask (in two steps for saving later)
    brain_mask_nii = nib.load(brain_mask_path_5)
    brain_mask = brain_mask_nii.get_fdata()
    inv2_nii = nib.load(scan_path_5)
    inv2_5 = inv2_nii.get_fdata()
    inv2_nii = nib.load(scan_path_8)
    inv2_8 = inv2_nii.get_fdata()
    inv2_nii = nib.load(scan_path_64)
    inv2_64 = inv2_nii.get_fdata()
    # inv2 = np.moveaxis(inv2, 2, 0)
    # brain_mask = np.moveaxis(brain_mask, 2, 0)

    # dilate brain mask and fill in with (second inversion) scan 
    dilation_five, inv_five = dilate(f'{root_path}/{subject_5}', save_mask=False, iters=5)
    dilation_ten_5, inv_ten = dilate(f'{root_path}/{subject_5}', save_mask=False,iters=10)
    dilation_fifteen, inv_fifteen = dilate(f'{root_path}/{subject_5}', save_mask=False, iters=15)
    dilation_twenty, inv_twenty = dilate(f'{root_path}/{subject_5}', save_mask=False, iters=20)

    dilation_ten_8, inv_ten = dilate(f'{root_path}/{subject_8}', save_mask=False,iters=10)
    dilation_ten_64, inv_ten = dilate(f'{root_path}/{subject_64}', save_mask=False,iters=10)

    # plot_different_planes(dilation_ten_5, brain_mask, '../results/figures/dilation_different_planes_scipy_mask.pdf')
    # plot_different_planes(dilation_ten_5, inv2_5, '../results/figures/dilation_different_planes_scipy_scan.pdf')
    # plot_different_iterations_sagital(brain_mask, dilation_five, dilation_ten_5, dilation_fifteen, dilation_twenty, inv2_5, '../results/figures/kernel_iters=5-20_sagital.pdf')
    # plot_different_iterations_coronal(brain_mask, dilation_five, dilation_ten_5, dilation_fifteen, dilation_twenty, inv2_5, '../results/figures/kernel_iters=5-20_coronal.pdf')
    plot_different_subjects(dilation_ten_5, dilation_ten_8, dilation_ten_64, inv2_5, inv2_8, inv2_64, '../results/figures/dilation_different_subjects.pdf')

    # compare_kernels(dilation_five, dilation_twenty, '../results/figures/dilation_kernel_comparison.pdf')

def compare_kernels(kernel1, kernel2, save_path):
    '''
    Plots and saves a comparison of a bigger kernel and a smaller kernel with more iterations.

    Args:
        kernel1: smaller kernel with more iterations
        kernel2: bigger kernel with less iterations 
        save_path: file path to save the plot to 
    '''
    s = 63

    fig, ax = plt.subplots(1,2)
    ax[0].imshow(ndimage.rotate(kernel1[:,:,s], 90), cmap='gray')
    ax[0].axis('off')
    ax[0].set_title('5x5 kernel, 15 iterations')
    ax[1].imshow(ndimage.rotate(kernel2[:,:,s], 90), cmap='gray')
    ax[1].axis('off')
    ax[1].set_title('50x50 kernel, 1 iteration')
    plt.savefig(save_path)
    plt.show()


def plot_different_planes(dilation, brain, save_path):
    ''' 
    Plots and saves scan or mask in different planal views. 
    
    Args:
        dilation: dilated brain mask
        brain: original MRI scan or brain mask   
        save_path: file path to save the plot to 
    '''
    x = 49
    y = 61
    z = 74

    plt.subplot(231)
    plt.title('Sagital view')
    plt.imshow(ndimage.rotate((brain)[:,:,z],90), cmap='gray', aspect='auto')
    plt.xticks([])
    plt.yticks([])
    plt.ylabel('Original', rotation=0, labelpad=30)
    plt.subplot(234)
    plt.imshow(ndimage.rotate((dilation*brain)[:,:,z],90), cmap='gray', aspect='auto')
    plt.xticks([])
    plt.yticks([])
    plt.ylabel('Defaced', rotation=0, labelpad=30)
    plt.subplot(232)
    plt.title('Coronal view')
    plt.imshow(ndimage.rotate((brain)[x,:,:],180), cmap='gray', aspect='auto')
    plt.axis('off')
    plt.subplot(235)
    plt.imshow(ndimage.rotate((dilation*brain)[x,:,:],180), cmap='gray', aspect='auto')
    plt.axis('off')
    plt.subplot(233)
    plt.title('Axial view')
    plt.imshow(ndimage.rotate((brain)[:,y,:],180), cmap='gray', aspect='auto')
    plt.axis('off')
    plt.subplot(236)
    plt.imshow(ndimage.rotate((dilation*brain)[:,y,:],180), cmap='gray', aspect='auto')
    plt.axis('off')
    plt.savefig(save_path)
    plt.show()


def plot_different_iterations_coronal(brain_mask, dilation1, dilation2, dilation3, dilation4, brain, save_path):
    ''' 
    Plots and saves a comparison of different amount of iterations for the dilation algorithm.

    Args:
        brain_mask: original brain mask 
        dilation1, dilation2, dilation3, dilation4: four dilations, each with a different amount of iterations performed
        brain: scan to apploy the dilated mask to
        save_path: file path to save the plot to 
    '''
    s = 49

    plt.subplot(231)
    plt.imshow(ndimage.rotate((brain_mask*brain)[s,:,:], 180), cmap='gray')
    plt.axis('off')
    plt.title('Brain mask')
    plt.subplot(232)
    plt.imshow(ndimage.rotate((dilation1*brain)[s,:,:], 180), cmap='gray')
    plt.axis('off')
    plt.title('5 iterations')
    plt.subplot(233)
    plt.imshow(ndimage.rotate((dilation2*brain)[s,:,:], 180), cmap='gray')
    plt.axis('off')
    plt.title('10 iterations')
    plt.subplot(234)
    plt.imshow(ndimage.rotate((dilation3*brain)[s,:,:], 180), cmap='gray')
    plt.axis('off')
    plt.title('15 iterations')
    plt.subplot(235)
    plt.imshow(ndimage.rotate((dilation4*brain)[s,:,:], 180), cmap='gray')
    plt.axis('off')
    plt.title('20 iterations')
    plt.subplot(236)
    plt.imshow(ndimage.rotate((brain)[s,:,:], 180), cmap='gray')
    plt.axis('off')
    plt.title('Full image')
    plt.savefig(save_path)
    plt.show()

def plot_different_iterations_sagital(brain_mask, dilation1, dilation2, dilation3, dilation4, brain, save_path):
    ''' 
    Plots and saves a comparison of different amount of iterations for the dilation algorithm.

    Args:
        brain_mask: original brain mask 
        dilation1, dilation2, dilation3, dilation4: four dilations, each with a different amount of iterations performed
        brain: scan to apploy the dilated mask to
        save_path: file path to save the plot to 
    '''
    s = 74

    plt.subplot(231)
    plt.imshow(ndimage.rotate((brain_mask*brain)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    plt.title('Brain mask')
    plt.subplot(232)
    plt.imshow(ndimage.rotate((dilation1*brain)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    plt.title('5 iterations')
    plt.subplot(233)
    plt.imshow(ndimage.rotate((dilation2*brain)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    plt.title('10 iterations')
    plt.subplot(234)
    plt.imshow(ndimage.rotate((dilation3*brain)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    plt.title('15 iterations')
    plt.subplot(235)
    plt.imshow(ndimage.rotate((dilation4*brain)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    plt.title('20 iterations')
    plt.subplot(236)
    plt.imshow(ndimage.rotate((brain)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    plt.title('Full image')
    plt.savefig(save_path)
    plt.show()

def plot_different_subjects(dilation_ten_5, dilation_ten_8, dilation_ten_64, inv2_5, inv2_8, inv2_64, save_path):
    s = 140

    plt.subplot(131)
    plt.imshow(ndimage.rotate((inv2_5*dilation_ten_5)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    # plt.title('Brain mask')
    plt.subplot(132)
    plt.imshow(ndimage.rotate((inv2_8*dilation_ten_8)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    # plt.title('5 iterations')
    plt.subplot(133)
    plt.imshow(ndimage.rotate((inv2_64*dilation_ten_64)[:,:,s], 90), cmap='gray')
    plt.axis('off')
    # plt.title('10 iterations')
    plt.savefig(save_path)
    plt.show()

if __name__ == '__main__':
    main()