addpath('../../../../../opt/amc/matlab/toolbox/MRecon-4.3.1');
% set bart env
setenv('TOOLBOX_PATH', '/opt/amc/bart-0.5.00/bin')
addpath(getenv('/opt/amc/bart-0.5.00/matlab'));

% function [mean_real, mean_imag, std_real, std_imag] = calc_noise_stats(file)

file1 = '../../../../../data/projects/ahead/raw_data/Subcortex_0005_002_R02/2018_06_21/Su_2919/su_21062018_1041475_14_2_wip_fatnavref220_senseV4.raw';
file2 = '../../../../../data/projects/ahead/raw_data/Subcortex_0064_064_R02/2018_07_16/Su_4040/su_16072018_1121033_6_2_wip_wb_inv2_senseV4.raw';
file3_inv2 = '../../../../../data/projects/ahead/raw_data/Subcortex_0005_002_R02/2018_06_21/Su_2919/su_21062018_0952593_6_2_wip_wb_inv2_senseV4.raw';
file3_inv1 = '../../../../../data/projects/ahead/raw_data/Subcortex_0005_002_R02/2018_06_21/Su_2919/su_21062018_0953202_7_2_wip_wb_inv1_senseV4.raw';

data_dir = '../../../../../data/projects/ahead/raw_data/';
files = dir(fullfile(data_dir, '**/**/*.raw'));

save_dir = '../data/';

for i=1:length(files)

    file = fullfile(files(i).folder, files(i).name);
    file_parts = split(file,'/');
    subj = file_parts(6);
    file_name = file_parts(9);
    save_file = fullfile(save_dir, subj{1});

    % continue to next iteration if filename is incorrect 
    if (contains(files(i).name,'wip_wb_inv2_senseV4') == 0) 
        continue 
    % also continue to next iteration if this subject has been handled 
    elseif (exist(save_file, 'dir') == 7)
        continue 
    end
    disp(file)
    
    % read in data 
    r = MRecon(file);
    
    [typ_size, second_dim] = size(r.Parameter.Parameter2Read.typ);
    
    if typ_size == 0
        continue 
    end 
    
    
    % select specific data type 
    % 1 is standard data, 3 is phase correction data, 5 is noise 
    r.Parameter.Parameter2Read.typ = 5;
    % r.Parameter.Parameter2Read.Update;
    r.ReadData;

    % perform corrections 
    r.DcOffsetCorrection;
    r.PDACorrection;
    r.RandomPhaseCorrection;
    r.MeasPhaseCorrection;
    r.SortData;
    r.GridData;

    d = r.Data;
     
    s = size(d);
        
    % translate to image space, but normalize
    % scale by 1/sqrt(N) for fft, by sqrt(N) for ifft (
    d = ifft(d) * sqrt(s(1));
    
    mean_real = cell(s(4),1);
    mean_imag = cell(s(4),1);
    std_real = cell(s(4),1);
    std_imag = cell(s(4),1);
    
    for coil=1:s(4)
        mean_real{coil} = mean(real(d(:,1,1,coil)));
        mean_imag{coil} = mean(imag(d(:,1,1,coil)));
        std_real{coil} = std(real(d(:,1,1,coil)));
        std_imag{coil} = std(imag(d(:,1,1,coil)));
    end
    tab = table(mean_real, mean_imag, std_real, std_imag);
    
    mkdir(save_file)
    save_file = fullfile(save_file, 'noise_stats.txt');
    
    % save table
    writetable(tab,save_file);
    disp('done saving')
end 
