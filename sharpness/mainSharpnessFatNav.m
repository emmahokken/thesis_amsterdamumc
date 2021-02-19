% Main prototype script for sharpness quantification in motion corrected 
% data using Fat Navigators.
%
% Demonstration for R1-maps in 1 subject. 
% By default, only ventricles are processed for short processing time, 
% change ROIs below if needed
%
% A number of plots will be generated in a *plots_saved* sub-directory.
%
% Requirements: 
%    Segmented ROIs using NighRes software with level-set distance maps
%
% Dependencies: 
%    MATLAB fitting and image processing toolboxes
%    load_untouch_nii.m for loading nifti files
%
% Matthan Caan, 2020
% Amsterdam UMC

%% setup
clear all
close all

if isempty(which('load_untouch_nii'))
  error('Install load_untouch_nii')
end

maindir=fileparts(which('mainSharpnessFatNav.m'));
subjdir=fullfile(maindir,'data');
codedir=fullfile(maindir,'code');
addpath(codedir)

plotsavedir=fullfile(maindir,'plots_saved'); % main dir for saving ploits
mkdir(plotsavedir)

%% Get data

% define some subject ID
FileID.uIDs={'064'};
FileID.type={'r2star'};

% define ROIs to process
% FileID.uROIs = {'vent', 'tha', 'str','gwb'};
FileID.uROIs = {'vent'};


FileID.uHEMs = {'l', 'r','4'};
[map_corrall, t1_corrall, map_uncorrall, t1_uncorrall, MPos, voxRes] = ...
  getdataSharpnessFatNav(FileID, subjdir, 1);
fields = fieldnames(map_corrall);
% parameters
run getParsFatNav

%% plot ROIs and save to pngs
run plotROIsFatNav

%% initialize
nsubj=1;
doSubj=1;

better_signed = cell(1,length(fields));
worse_signed = cell(1,length(fields));
AR = cell(1,length(fields));

%% run
subj_ii = find(~cellfun(@isempty,strfind(fields,FileID.uIDs{1})))';
nUsedROI = 0;
% loop over ROIs
for ii=subj_ii
  % get data
  field_name = fields{ii};
  map_corr = map_corrall.(fields{ii}).img;
  data_corr = t1_corrall.(fields{ii});

  map_uncorr = map_uncorrall.(fields{ii}).img;
  data_uncorr = t1_uncorrall.(fields{ii});
  data_uncorr=permute(data_uncorr,[2 3 1]);

  % do fitting
  [bs,ws,bc,bi,ar_list] = ...
    getsigmaclusterFatNav(map_corr,data_corr,data_uncorr,pars,field_name,plotsavedir);
    % contine to next iteration if ROI is outside of scan
    if strcmp(bs, 'outside')
        continue
    end
    disp('right')
    nUsedROI = nUsedROI + 1;
    [better_signed{ii},worse_signed{ii},brd_crds{ii},brd_ind{ii},AR{ii}] = ...
        deal(bs,ws,bc,bi,ar_list);

  % get motion parameters
  numcl(ii) = length(better_signed{ii});
  mot_mean{ii} = getmotionFatNav(map_corrall.(fields{ii}).img,brd_crds{ii},MPos.(fields{ii}),voxRes);
end
disp(length(AR))
nUsedROI
if length(AR) ~= nUsedROI
    disp(size(better_signed))
    AR(nUsedROI+1:end) = [];
end 
%% compute and print statistics
run statsFatNav

cd('../../')