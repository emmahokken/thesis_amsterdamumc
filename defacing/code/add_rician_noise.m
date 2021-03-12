path = '../../../../../data/projects/ahead/raw_gdata/Subcortex_0064_064_R02/nii/inv2_te1_m_corr.nii';
dilated_path = '../results/dilation/Subcortex_0064_064_R02/dilated_inv2_te1_m_corr.nii';
data = load_untouch_nii(dilated_path);
disp(size(data.img))
im = double(data.img(:,:,64)); % (particular slice through the brain)
s=5; % noise level (NB actual Rician stdev depends on signal, see ricestat)
im_g = im + 5 * randn(size(im)); % *Add* Gaussian noise
im_r = ricernd(im, s); % "Add" Rician noise (make Rician distributed)
% Compute ranges
min_o = round(min(im(:)));   max_o = round(max(im(:)));
min_g = round(min(im_g(:))); max_g = round(max(im_g(:)));
min_r = round(min(im_r(:))); max_r = round(max(im_r(:)));
% Show each image with the same color scaling limits
clim = [min_g max(max_g, max_r)];
figure('Position', [30 500 800 300]); 
colormap('gray');
subplot(1,3,1); imagesc(imrotate(im,90), clim); axis image;
title('Original');
xlabel(['Range: (' num2str(min_o) ', ' num2str(max_o) ')'])
subplot(1,3,2); imagesc(imrotate(im_g,90), clim); axis image;
title('Gaussian noise');
xlabel(['Range: (' num2str(min_g) ', ' num2str(max_g) ')'])
subplot(1,3,3); imagesc(imrotate(im_r,90), clim); axis image;
title('Rician noise')
xlabel(['Range: (' num2str(min_r) ', ' num2str(max_r) ')'])
saveas(gcf,'noise.png')

