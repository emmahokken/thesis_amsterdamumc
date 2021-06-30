# MSc AI thesis: Towards Open-Source MRI Data: Defacing and Sharpness Estimation  

Magnetic Resonance Imaging (MRI) is a widely used tool in both research and clinical settings, with image quality increasing every year. However, the acquisition time of a high-quality MRI can be very long. Accelerated MRI is a method in which this acquisition time is reduced by undersampling the available data points.
Undersampling causes aliasing, i.e. overlap of regions, and adds noise to the image. A variety of methods are in place to reconstruct the image, with recent years showing a spike in reconstruction methods that utilise Deep Learning. One such method is the Recurrent Inference Machine (RIM).

This work comes in two parts. The first part aims to stimulate the research into Deep Learning-based accelerated MRI reconstruction by releasing the Amsterdam Ultra-high field adult lifespan database (AHEAD), a 7T dataset. To comply with privacy regulations such as GDPR, these images are to be defaced by the dilation of the brain mask. To simulate a real scan as much as possible, noise was added to the background. While the defacing method was able to accurately and precisely deface scans, the additional noise proved to be overestimated when compared to an original scan.

The second part focuses on creating a content-wise evaluation method for reconstructed MRI images. Currently, Deep Learning reconstruction methods are evaluated using quantitative measures. None of these evaluations take the image content quality into account. This work proposes a method that estimates the edge sharpness across tissue boundaries within an MR image using a two-fold k-means clustering approach. This method showed a decrease in sharpness for an increase in acceleration factor, showing that it is a valid method to evaluate image content quality and can provide valuable insights in Deep Learning-based accelerated MRI reconstruction.

## Code base

The code is devided by topic. 

### Defacing
Defacing algorithms and results. 

`code/`: contains the code for the defacing algorithm and experiments. This folder contains:

* `resources/` contain necessary resources for Freesurfer defacing. 
* `check_noise.py` opens all (dilation) defaced images and calculates background noise from cube. has the option to plot the cube. 
* `compare_defacing.py` compares Freesurfer and Pydeface. 
* `compare_kernels.py` compares differently sized kernels and plots dilation from different planes. 
* `compare_noise.py` reads in the saved file with noise measurements (from `results/files/`), computes differences and creates figures. 
* `compute_noise.m` computes noise using MRecon. this method is not used in the final project. 
* `defacing.py` applies Freesurfer defacing.
* `defacing_pipeline.py` pipeline for defacing using dilation. applies dilation of brain mask and computes STD of noise, adds it to the background. 
* `dilation.py` performs dilation of the brain mask. 

`results/`: contains results. 

* `figures/` contains resulting figures
* `files/` contains the file with information on added Gaussian background noise measurements (STD)


### Sharpness
Sharpness Estimation code and results. 

`mainSharpnessFatNav.m` contains the main sharpness estimation pipeline. 

`code/`: contains the code for the sharpness estimation. This folder contains:

* `python/` contains Python code with the following files:
	* `change_dimensions.py` small script to change dimensions of R2* maps. 
	* `check_segmentation.py` checks whether the segmentation and coregistration were successful. 
	* `create_segmentation.py` coregister scans and create levelmaps from binary segmentation maps.
	* `move.py` copy files in one directory to another.
	* `plot_fwhm.py` plots various figures based on the FWHM results. 
	* `statistics_fwhm.py` perform various statistical analyses, including paired t-test. 
* `errorFunctionFATNAV.m` contains code that calculates the error function used in the estimation. 
* `getdataSharpnessFatNav.m` loads the levelmaps and MRI data. 
* `getmotionFatNat.m` calculates displacement between ground truth and accelerated image.  
* `getParsFatNav.m` saves the parameters used in the estimation. 
* `getsigmaclusterFatNav.m` performs clustering and main sigma calculation
* `plotROIsFatNav.m` plots the ROIs found in `plots_saved/ROIs_anat/`.
* `statsFatNav.m` calculates satistics (including FWHM) and saves CSV to `results/`. 

`data/`: main folder contains R2* maps usd in sharpness estmation. This folder also contains:

* `coregistration` temporary save location for coregistering pipeline. 
* `segm/` contains levelmaps of regions. 
* `variables/` saves clustering regions after initial runs to ensure stable ground truth FWHM. 

`plots_saved/`: contains plots saved 

* `ROIs_anat` contains images of the ROIs used in the sharpness estimation. 

`results/`: contains Matlab-generated CSV files used for statistical analysis in Python. 

