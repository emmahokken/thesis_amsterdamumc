acceleration_factors = {'test_all_3_ssim/','test_all_6_ssim/','test_all_9_ssim/','test_all_12_ssim/','test_all_15_ssim/'};
datadir = '../../../../data/projects/recon/data/qMRI/Brain_MEGRE/results/';
datasets = {'/R2star_map_gt', '/R2star_map_recon', '/R2star_map_rim'};
niidir = '../../data/recon/';

% iterate over files 
for i=1:length(acceleration_factors)
    % select acceleration factor file 
    acc = acceleration_factors{i};

    % create new directory to save files
    mkdir(strcat(niidir,acceleration_factors{i}, '/',datasets{1}));
    mkdir(strcat(niidir,acceleration_factors{i}, '/',datasets{2}));
    mkdir(strcat(niidir,acceleration_factors{i}, '/',datasets{3}));

    accdir = strcat(datadir, acc);
    
    % get all files in acceleration directory 
    scans = dir(fullfile(accdir, '*.h5'));
    
    for j=1:length(scans)
        scan = scans(j).name;
        file = strcat(accdir, scan);

        % convert for both ground truth and reconstructed images 
        for k=1:length(datasets)
            r = h5read(file, datasets{k});
            [p, name] = fileparts(file);
            r = permute(r,[2 3 1]);
            disp(size(r))
            final_niidir = strcat(niidir, acc);
            outfile = fullfile(final_niidir,[strcat(datasets{k},'/',name),'.nii']);
            spacing = [0.7,0.7,0.7];
            nii = make_nii(r, spacing);
            save_nii(nii, outfile);
        end 
    end 
end
    