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

# dil = inv2*dilation_fifteen
# plt.imshow(dil[:,:,120])
# plt.show()
# nif = nib.Nifti1Image(dil, affine=inv2_nii.affine, header=inv2_nii.header)
# nib.save(nif,'triedilated_in12.nii')

# # dilate with 20x20 kernel
# dilation_mask_five = cv2.dilate(brain_mask, kernel_five, iterations=13)
# filled_in_five = dilation_twenty*inv2
# dilation_mask_fifty = cv2.dilate(brain_mask, kernel_fifty, iterations=1)
# filled_in_fifty = dilation_twenty*inv2
# difference = dilation_mask_fifty - dilation_mask_five

# plt.imshow(ndimage.rotate(difference[:,:,63], 90), cmap='gray')
# plt.show()

# fig, ax = plt.subplots(1,2)
# ax[0].imshow(ndimage.rotate(dilation_five[:,:,63], 90), cmap='gray')
# ax[0].axis('off')
# ax[0].set_title('5x5 kernel, 15 iterations')
# ax[1].imshow(ndimage.rotate(dilation_twenty[:,:,63], 90), cmap='gray')
# ax[1].axis('off')
# ax[1].set_title('50x50 kernel, 1 iteration')
# plt.savefig('../results/figures/dilation_kernel_comparison.pdf')
# plt.show()
s = 74


''' PLOT DIFFERENT ITERATIONS '''
plt.subplot(231)
plt.imshow(ndimage.rotate((brain_mask*inv2)[:,:,s], 90), cmap='gray')
plt.axis('off')
plt.title('Brain mask')
plt.subplot(232)
plt.imshow(ndimage.rotate((dilation_five*inv2)[:,:,s], 90), cmap='gray')
plt.axis('off')
plt.title('5 iterations')
plt.subplot(233)
plt.imshow(ndimage.rotate((dilation_ten*inv2)[:,:,s], 90), cmap='gray')
plt.axis('off')
plt.title('10 iterations')
plt.subplot(234)
plt.imshow(ndimage.rotate((dilation_fifteen*inv2)[:,:,s], 90), cmap='gray')
plt.axis('off')
plt.title('15 iterations')
plt.subplot(235)
plt.imshow(ndimage.rotate((dilation_twenty*inv2)[:,:,s], 90), cmap='gray')
plt.axis('off')
plt.title('20 iterations')
plt.subplot(236)
plt.imshow(ndimage.rotate((inv2)[:,:,s], 90), cmap='gray')
plt.axis('off')
plt.title('Full image')
plt.savefig('../results/figures/kernel_size_kernels=15-60.pdf')
plt.show()

''' PLOT DIFFERENT PLANES MASK '''
plt.subplot(231)
plt.title('Sagital view')
plt.imshow(ndimage.rotate((brain_mask)[:,:,74],90), cmap='gray', aspect='auto')
plt.xticks([])
plt.yticks([])
plt.ylabel('Original', rotation=0, labelpad=30)
plt.subplot(234)
plt.imshow(ndimage.rotate((dilation_ten)[:,:,74],90), cmap='gray', aspect='auto')
plt.xticks([])
plt.yticks([])
plt.ylabel('Dilated', rotation=0, labelpad=30)
plt.subplot(232)
plt.title('Coronal view')
plt.imshow(ndimage.rotate((brain_mask)[49,:,:],180), cmap='gray', aspect='auto')
plt.axis('off')
plt.subplot(235)
plt.imshow(ndimage.rotate((dilation_ten)[49,:,:],180), cmap='gray', aspect='auto')
plt.axis('off')
plt.subplot(233)
plt.title('Axial view')
plt.imshow(ndimage.rotate((brain_mask)[:,61,:],180), cmap='gray', aspect='auto')
plt.axis('off')
plt.subplot(236)
plt.imshow(ndimage.rotate((dilation_ten)[:,61,:],180), cmap='gray', aspect='auto')
plt.axis('off')
plt.savefig('../results/figures/dilation_different_planes_scipy_mask.pdf')
plt.show()

''' PLOT DIFFERENT PLANES SCAN '''
plt.subplot(231)
plt.title('Sagital view')
plt.imshow(ndimage.rotate((inv2)[:,:,74],90), cmap='gray', aspect='auto')
plt.xticks([])
plt.yticks([])
plt.ylabel('Original', rotation=0, labelpad=30)
plt.subplot(234)
plt.imshow(ndimage.rotate((dilation_ten*inv2)[:,:,74],90), cmap='gray', aspect='auto')
plt.xticks([])
plt.yticks([])
plt.ylabel('Dilated', rotation=0, labelpad=30)
plt.subplot(232)
plt.title('Coronal view')
plt.imshow(ndimage.rotate((inv2)[49,:,:],180), cmap='gray', aspect='auto')
plt.axis('off')
plt.subplot(235)
plt.imshow(ndimage.rotate((dilation_ten*inv2)[49,:,:],180), cmap='gray', aspect='auto')
plt.axis('off')
plt.subplot(233)
plt.title('Axial view')
plt.imshow(ndimage.rotate((inv2)[:,61,:],180), cmap='gray', aspect='auto')
plt.axis('off')
plt.subplot(236)
plt.imshow(ndimage.rotate((dilation_ten*inv2)[:,61,:],180), cmap='gray', aspect='auto')
plt.axis('off')
plt.savefig('../results/figures/dilation_different_planes_scipy_scan.pdf')
plt.show()