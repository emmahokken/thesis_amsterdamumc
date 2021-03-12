import matplotlib.pyplot as plt 
import numpy as np
import scipy.stats as sp

def add_noise(im):
    im_slice = im[:,:,64] # (particular slice through the brain)
    s = 5 # noise level (NB actual Rician stdev depends on signal, see ricestat)
    im_g = im_slice + 5 * np.random.normal(0, 1, im_slice.shape) # *Add* Gaussian noise
    r = sp.rice.rvs(s, size=im_slice.shape)
    im_r = ricernd(im, s) # "Add" Rician noise (make Rician distributed)
    # Compute ranges
    min_o = round(min(im_slice))   
    max_o = round(max(im))
    min_g = round(min(im_g))
    max_g = round(max(im_g))
    min_r = round(min(im_r))
    max_r = round(max(im_r))

    # Show each image with the same color scaling limits
    clim = [min_g max(max_g, max_r)]
    
    # colormap(data.map)
    fix, ax = plt.subplots(1,3)
    ax[0].imshow(im_slice, cmap='gray')

    ax[0].set_title('Original')
    ax[1].imshow(im_g, cmap='gray')
    ax[1].set_title('Gaussian noise')
    
    # imagesc(im_r, clim)
    ax[2].imshow(im_r, cmap='gray')
    ax[2].set_title('Rician noise')
