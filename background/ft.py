'''Create an image that contains both a scan in kspace and a scan in image space'''

import matplotlib.pyplot as plt 
import nibabel as nib 
import numpy as np
from scipy import ndimage


prefix = '../../../../..'
inv1_path = prefix + '/data/projects/ahead/raw_gdata/Subcortex_0005_002_R02/nii/inv2_te2_m_corr.nii'
inv1 = nib.load(inv1_path).get_fdata()
kspace = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(inv1), axes=(0,1,2), norm='ortho'))

kspace_real = kspace.real 

plt.subplot(121)
plt.imshow(ndimage.rotate(kspace_real[:,:,63],90), cmap='gray', vmin=0, vmax=150000)
plt.title('k-space')
plt.axis('off')
plt.subplot(122)
plt.imshow(ndimage.rotate(inv1[:,:,63],90), cmap='gray', vmin=0, vmax=2500000)
plt.title('image space')
plt.axis('off')
plt.savefig('kspace_imspace_comparison.pdf')
plt.show()