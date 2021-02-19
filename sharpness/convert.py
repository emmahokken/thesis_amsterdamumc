import h5py
import nibabel as nib 
import numpy as np
import matplotlib.pyplot as plt 
import os 
from tqdm import tqdm

base = '../../../../..'
subj = 64
segm_path = f'/data/projects/ahead/segmentations/mapped/'
spec_segm = f'sub-0{subj}_mask-tha_hem-r_lvlreg-corr_def-img.nii.gz'
recon_h5_path = f'/data/projects/recon/data/qMRI/Brain_MEGRE/results/test_all_3_ssim/results_Subcortex_00{subj}_axial_loss_mse.h5'
recon_nii_path = f'../../data/recon/test_all_3_ssim/R2star_map_rim/results_Subcortex_00{subj}_axial_loss_mse.nii'
brain_mask_path = f'/data/projects/ahead/raw_gdata/Subcortex_00{subj}_0{subj}_R02/nii/mask_inv2_te2_m_corr.nii'
sharpness_path = f'data/segm/sub-mask-str_hem-l_lvlreg-uncorr_def-img.nii'
t1corr_cruise_path = f'/data/projects/ahead/segmentations/cortex/sub-0{subj}_t1corr_cruise-gwb_ants-def0.nii'
t1uncorr_cruise_path = f'/data/projects/ahead/segmentations/cortex/sub-0{subj}_t1uncorr_cruise-gwb_ants-def0.nii'
save_path = 'data/segm/cropped'

# segm = nib.load(base+segm_path+spec_segm).get_fdata()
# recon_nii = nib.load(recon_nii_path)
# r_nii = recon_nii.get_fdata()
# recon_h5 = h5py.File(base+recon_h5_path, 'r')
# recon_rim = recon_h5['R2star_map_gt']
# recon_mask = recon_h5['mask_brain']
# brain_mask = nib.load(base+brain_mask_path).get_fdata()
# sharp = nib.load(sharpness_path).get_fdata()
# sta = nib.load('data/segm/sub-t1corr_cruise-gwb_ants-def0.nii').get_fdata()

# # crop grey white barrier t1corr map 
# t1corr_cruise_nii = nib.load(base+t1corr_cruise_path)
# t1uncorr_cruise_nii = nib.load(base+t1uncorr_cruise_path)
# t1corr_cruise = t1corr_cruise_nii.get_fdata()
# t1uncorr_cruise = t1uncorr_cruise_nii.get_fdata()

reconscan = nib.load('data/064_r1corr.nii').get_fdata()
fig, ax = plt.subplots(3,4)
ax[0,0].imshow(reconscan[:,:,50], cmap='gray')
ax[0,0].axis('off')
ax[0,1].imshow(reconscan[:,:,100], cmap='gray')
ax[0,1].axis('off')
ax[0,2].imshow(reconscan[:,:,150], cmap='gray')
ax[0,2].axis('off')
ax[0,3].imshow(reconscan[:,:,200], cmap='gray')
ax[0,3].axis('off')
ax[1,0].imshow(reconscan[:,10,:], cmap='gray')
ax[1,0].axis('off')
ax[1,1].imshow(reconscan[:,20,:], cmap='gray')
ax[1,1].axis('off')
ax[1,2].imshow(reconscan[:,30,:], cmap='gray')
ax[1,2].axis('off')
ax[1,3].imshow(reconscan[:,40,:], cmap='gray')
ax[1,3].axis('off')
ax[2,0].imshow(reconscan[50,:,:], cmap='gray')
ax[2,0].axis('off')
ax[2,1].imshow(reconscan[100,:,:], cmap='gray')
ax[2,1].axis('off')
ax[2,2].imshow(reconscan[150,:,:], cmap='gray')
ax[2,2].axis('off')
ax[2,3].imshow(reconscan[200,:,:], cmap='gray')
ax[2,3].axis('off')
plt.show()


start = 121
stop = 171

t1corr_cruise = t1corr_cruise[:,start:stop,:]
t1uncorr_cruise = t1uncorr_cruise[:,start:stop,:]

# # save t1corr and t1uncorr maps 
# nifti_image = nib.Nifti1Image(dataobj=t1corr_cruise, header=t1corr_cruise_nii.header, affine=t1corr_cruise_nii.affine)
# nib.save(img=nifti_image, filename=save_path+'/'+t1corr_cruise_path.split('/')[-1])
# nifti_image = nib.Nifti1Image(dataobj=t1uncorr_cruise, header=t1uncorr_cruise_nii.header, affine=t1uncorr_cruise_nii.affine)
# nib.save(img=nifti_image, filename=save_path+'/'+t1uncorr_cruise_path.split('/')[-1])

# # rearange segmented image to be in same line as recon rim image.
# r_nii = np.moveaxis(r_nii, 0, -1)
# nifti_image = nib.Nifti1Image(dataobj=r_nii, header=recon_nii.header, affine=recon_nii.affine)
# nib.save(img=nifti_image, filename='data/rim_r1corr.nii')


for subdir, dirs, files in os.walk(segm_path):
    for file in tqdm(files):
        if f'sub-0{subj}_mask' in file:
            s_nii = nib.load(segm_path+file)
            s = s_nii.get_fdata()
            s = s[:,start:stop,:]
            
            # print(s[:20,:20,:20])
            s = np.flip(s,axis=0)
            
            # save new scan 
            if not os.path.exists(f'{save_path}'):
                os.makedirs(f'{save_path}')
            save_dir = f'{save_path}/{file[:-3]}'
           
            nifti_image = nib.Nifti1Image(dataobj=s, header=s_nii.header, affine=s_nii.affine)
            nib.save(img=nifti_image, filename=save_dir)

segm[segm > 0] = 0
segm[segm < 0] = 1
plt.imshow(segm[:,:,148])
plt.colorbar()
plt.show()


fig, ax = plt.subplots(2,2)
ax[0,0].imshow(segm[:,:,64], cmap='gray')
ax[0,0].set_title('Segmentation')
ax[0,0].axis('off')
ax[0,1].imshow(recon_nii[:,:,20], cmap='gray')
ax[0,1].set_title('Reconstructed scan')
ax[0,1].axis('off')
ax[1,0].imshow(test, cmap='gray')
ax[1,0].set_title('Segmented out')
ax[1,0].axis('off')
ax[1,1].imshow(recon_mask[20,:,:], cmap='gray')
ax[1,1].set_title('Recon mask')
ax[1,1].axis('off')
plt.show()

# try to see if brain mask is identical in some part to the mask from the original brain. 
# although we do have the ground truth
# no but for the segmented image you need the exact intersection of the recon image and the original because you need to crop the segmented image to that region
# i don't really understand these segmentations if i'm honest
