import h5py
import nibabel as nib 
import numpy as np
import matplotlib.pyplot as plt 

path = '../../../../../data/projects/recon/data/qMRI/Brain_MEGRE/test_all_3_ssim/'
file = 'Subcortex_0005_axial_121.h5'

f = h5py.File(path + file, 'r')
print(f.keys())
gt = f['qMRI_GT']
rim = f['qMRI_RIM']
print(rim)
gt = np.array(gt).squeeze()
rim = np.array(rim).squeeze()

# gt = gt.squeeze() 
fig, ax = plt.subplots(2,2)
ax[0,0].imshow(gt[0,:,:], cmap='gray')
ax[0,1].imshow(gt[1,:,:], cmap='gray')
ax[1,0].imshow(gt[2,:,:], cmap='gray')
ax[1,1].imshow(gt[3,:,:], cmap='gray')
plt.show()
print(gt.shape)