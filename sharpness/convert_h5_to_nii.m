files = {'test_all_3_ssim/','test_all_6_ssim/','test_all_9_ssim/','test_all_12_ssim/','test_all_15_ssim/'};
datadir = '../../../../data/projects/recon/data/qMRI/Brain_MEGRE/';
scanname = 'Subcortex_0005_axial_121.h5';
niidir = '../../data/recon/';
datasets = {'qMRI_GT', 'qMRI_RIM'};

% iterate over files 
for i=1:length(files)
    % select acceleration factor file 
    acc = files{i};

    % create new directory to save files
    mkdir(strcat(niidir,files{i}));
    accdir = strcat(datadir, acc);
    
    % get all files in acceleration directory 
    scans = dir(fullfile(accdir, '*.h5'));
    
    for j=1:length(scans)
        scan = scans(j).name;
        file = strcat(accdir, scan);
  
        % convert for both ground truth and reconstructed image 
        for k=1:length(datasets)
            r = h5read(file, strcat('/',datasets{k}));
            [p, name] = fileparts(file);

            final_niidir = strcat(niidir, acc);
            outfile = fullfile(final_niidir,[strcat(datasets{k},'_',name),'.nii']);
            spacing = [0.7,0.7,0.7];
            nii = make_nii(r, spacing);
            save_nii(nii, outfile);
        end 
    end 
end
    