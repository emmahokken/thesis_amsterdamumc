import numpy as np
import matplotlib.pyplot as plt 
import cv2 as cv

kspace_mask = np.zeros((234, 290), dtype=int)
x_mid = kspace_mask.shape[0] // 2
y_mid = kspace_mask.shape[1] // 2

r_outer = 110
r_inner = 80
h, w = kspace_mask.shape[:2]
y, x = np.ogrid[:h, :w]
circ_mask = np.sqrt((x- y_mid)**2 + (y- x_mid)**2) < r_outer

circ_mask2 = np.sqrt((x- y_mid)**2 + (y- x_mid)**2) < r_inner

circ_mask[circ_mask2] = 0

# el = cv.ellipse(kspace_mask)

# circular_mask = (kspace_mask.shape[0]-x_mid)**2 + (kspace_mask.shape[1]-y_mid)**2 < r**2
# plt.imshow(el)
# plt.show()
x = 50
kspace_mask[:x,:x] = 1 
plt.imshow(kspace_mask)
plt.show()