import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import nibabel as nib 
import numpy as np 
import os
from scipy import ndimage
import sys 
from tqdm import tqdm 

prefix = '../../../../../..'
root_path = prefix + '/data/projects/ahead/raw_gdata/'
save_path = '../results/dilation/'

subjects = ['0004', '0016', '0083', '0067']

subj = 'Subcortex_0083_085_R02'
directories = [d[1] for d in  os.walk(save_path)][0]


# original = nib.load(f'{root_path}{subj}/{subj}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
# defaced = nib.load(save_path + subj + f'/dilated_{subj}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
# no_noise = nib.load(save_path + subj + f'/dilated_no_noise_{subj}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)



# print(original.shape)
# print(defaced.shape)
# print(no_noise.shape)

# print(np.std(original[:,:,:,7]))
# print(np.std(original[:,:,:,15]))
# print(np.std(original[:50,:50,:50,7].real))
# print(np.std(original[:50,:50,:50,15].real))

# print(original[140:150,140:150,140,2].real)

# print(defaced[140:150,140:150,140,2].real)

# print(no_noise[140:150,140:150,140,2].real)

# fig, axes = plt.subplots(1,3)
# scans = [original, no_noise, defaced]
# files = [{'Original': 0}, {'Defaced, without noise': 1}, {'Defaced, with noise': 2}]
# for i, ax in enumerate(axes.flat):
#     scan = files[i].values()
#     title = files[i].keys()
#     print(scan)
#     im = ax.imshow(scan[:,:,140,2].real, cmap='gray')
#     ax.set_title(title)

# fig.subplots_adjust(right=0.8)
# fig.colorbar(im,ax=axes.ravel().tolist())
# print(np.min(original.real), np.max(original.real))

# print(np.min(defaced.real), np.max(defaced.real))

# plt.subplot(131)
# plt.imshow(ndimage.rotate(original[149,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
# plt.title('Original image')
# plt.subplot(132)
# plt.imshow(ndimage.rotate(no_noise[140,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
# plt.title('Defaced, without noise')
# plt.subplot(133)
# plt.imshow(ndimage.rotate(defaced[140,:,:,2].real, 90), vmin=-100, vmax=100, cmap='gray')
# plt.title('Defaced, with noise')

# plt.colorbar(shrink=0.5)
# plt.show()


# exit()

# find std of corner region 
x = 50
y = 50
z = 50 
s = 140

# plot box around region
# fig, ax = plt.subplots()
# ax.imshow(ndimage.rotate(original[s,:,:,8].real, 90), cmap='gray')
# rect = patches.Rectangle((original.shape[1]-x,0), x, y, linewidth=1, edgecolor='r', facecolor='none')
# ax.add_patch(rect)
# ax.axis('off')
# plt.show()

# caclulate for different coils
coils = [7,15,23,31]

for d in directories:
    print(d)
    original = nib.load(f'{root_path}{d}/{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
    defaced = nib.load(save_path + d + f'/dilated_{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)
    no_noise = nib.load(save_path + d + f'/dilated_no_noise_{d}_inv2_1_gdataCorrected.nii.gz').get_fdata(dtype=np.complex64)

    for coil in coils:
        print('Coil:', coil)
        # compute std for both imageinary and real channels from original and defaced 
        std_real_orig = np.std(original[:x,:y,:z, coil].real)
        std_imag_orig = np.std(original[:x,:y,:z, coil].imag)

        std_real_def = np.std(defaced[:x,:y,:z, coil].real)
        std_imag_def = np.std(defaced[:x,:y,:z, coil].imag)

        print('Original:', std_real_orig, std_imag_orig)
        print('Defaced:', std_real_def, std_imag_def)




    