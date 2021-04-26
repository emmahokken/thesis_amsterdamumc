
datadir = '../../../../data/projects/ahead/raw_gdata/';

niidir = '../../data/converted_mat_echo1/';

% get all files in directory
S = dir(fullfile(datadir, '*'));
scans = setdiff({S([S.isdir]).name},{'.','..'});
scans = scans(56:end);

subj = {'98', '105'};
% '5', '8', '18', '25', '31', '64', '70', '80',
% iterate over files 
for i=1:length(scans)
    filepath = dir(fullfile(strcat(datadir,'/',scans{i}),'*01_allm.mat'));
    disp(filepath.folder)
    disp(filepath.name)
    subject = 0;
    for j=1:length(subj)
        filepath.name
% %         filepath
        subj{j}
%         contains(filepath.name(1:14), subj{j})
        if contains(filepath.name(1:14), subj{j})
            subject = 1
            break
        end
    end
%     
    if subject == 0
        continue
    end
        
    filepath.name

    file = fullfile(filepath.folder, filepath.name);
    image = load(file);
    image = abs(image.allm{1,1}.all_ims_corrected);
    scalingfactor = max(max(max(image, [], 3),[],2),[],1);
    image = int16(255*image/scalingfactor);

    mkdir(strcat(niidir,scans{i},'/'));
    
    nii = make_nii(image);
    save_nii(nii, strcat(niidir,scans{i},'/',filepath.name(1:end-4),'.nii'));
end 
    