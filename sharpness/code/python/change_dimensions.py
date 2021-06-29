import cv2 
import matplotlib.pyplot as plt 
import nibabel as nib 
import nighres 
import numpy as np 
import os 
import sys 
from tqdm import tqdm 
from scipy import ndimage

def main():
    ''' Adjust dimensions for reconstructed R2* maps. '''

    # declare relevant prefixes and directories 
    gt_path = '../../../data/recon/test_all_3_ssim/R2star_map_gt/'
    rim_path = '../../../data/recon/test_all_3_ssim/R2star_map_rim/'

    # gather all files 
    gt_files = [f[2] for f in os.walk(gt_path)][0]
    rim_files = [f[2] for f in os.walk(rim_path)][0]

    gt_files.sort()
    rim_files.sort()

    for g, r in zip(gt_files, rim_files):

        gt_file = nib.load(gt_path + g)
        gt = gt_file.get_fdata()
        
        rim_file = nib.load(rim_path + r)
        rim = rim_file.get_fdata()


        print(gt.shape)
        print(rim.shape)

if __name__ == '__main__':
    main()