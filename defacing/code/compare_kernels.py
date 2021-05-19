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
    scan_type = 'inv1_te1_m_corr'
    # scan_type = 't1corr'
    brain_mask_type = 'mask_inv2_te2_m_corr'
    save_path = '../defacing/dilation'
    subject = 'Subcortex_0005_002_R02'

    kernel_five = np.ones((5,5), dtype=np.uint8)
    kernel_fifty = np.ones((50,50), dtype=np.uint8)

    brain_mask_path = f'{root_path}/{subject}/nii/{brain_mask_type}.nii'
    scan_path = f'{root_path}/{subject}/nii/{scan_type}.nii'

    # load in brain mask (in two steps for saving later)
    brain_mask_nii = nib.load(brain_mask_path)
    brain_mask = brain_mask_nii.get_fdata()
    inv2_nii = nib.load(scan_path)
    inv2 = inv2_nii.get_fdata()
    # inv2 = np.moveaxis(inv2, 2, 0)
    # brain_mask = np.moveaxis(brain_mask, 2, 0)

    # dilate brain mask and fill in with (second inversion) scan 
    dilation_five, inv_five = dilate(f'{root_path}/{subject}', save_mask=False, iters=5)
    dilation_ten, inv_ten = dilate(f'{root_path}/{subject}', save_mask=False,iters=10)
    dilation_fifteen, inv_fifteen = dilate(f'{root_path}/{subject}', save_mask=False, iters=15)
    dilation_twenty, inv_twenty = dilate(f'{root_path}/{subject}', save_mask=False, iters=20)

    plot_different_planes(dilation_ten, brain_mask, '../results/figures/dilation_different_planes_scipy_mask.pdf')
    plot_different_planes(dilation_ten, inv2, '../results/figures/dilation_different_planes_scipy_scan.pdf')
    plot_different_iterations(brain_mask, dilation_five, dilation_ten, dilation_fifteen, dilation_twenty, inv2, '../results/figures/kernel_size_kernels=15-60.pdf')

    compare_kernels(dilation_five, dilation_twenty, '../results/figures/dilation_kernel_comparison.pdf')

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
    plt.ylabel('Dilated', rotation=0, labelpad=30)
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

def plot_different_iterations(brain_mask, dilation1, dilation2, dilation3, dilation4, brain, save_path):
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

if __name__ == '__main__':
    main()