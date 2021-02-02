function [ map_recon, subj_data_recon, map_gt, subj_data_gt, MPos,voxRes] = getdataSharpnessFatNav( FileID, subjdir, varargin)
% Get data for FileID. Output: label map (distance transform), t1_corr
% data, r1/t1_uncorr data and MPos motion parameters.
%
% Requires load_untouch_nii.m
%
% Matthan Caan, 2020
% Amsterdam UMC

doR1=0; % do R1 instead of T1
if nargin==3
  doR1=varargin{1};
end
if doR1
  filetag='rim_r1';    
else
  filetag='t1';
end
disp(filetag)
cortag='gwb';

lbldir = fullfile(subjdir,'segm/cropped');

datadir = subjdir;

% Load motion parameters - what do these do?
cd(subjdir)
mparsname = dir('mpars.mat');
mpars = load(char(mparsname.name));
ID = strcat(char('sub-'),FileID.uIDs{1});
ID = char('sub-064_');
id = strcat('s',FileID.uIDs{1});
disp(id)
disp('fieldid')
disp(FileID.uROIs)
% Iterate over the ROIs
for iROI = 1:length(FileID.uROIs)  
    if strcmp(FileID.uROIs{iROI},cortag)
      ROI=FileID.uROIs{iROI};
    else
      ROI = strcat(char('mask-'),FileID.uROIs{iROI});
    end
    disp(ROI)
    roi = FileID.uROIs{iROI};
    % Iterate over hemispheres 
    for iHEM = 1:length(FileID.uHEMs)
      if iHEM<3, doit=1; else doit=0; end
      if strcmp(FileID.uROIs{iROI},cortag) && iHEM==2
        doit=0;
      end
      % for ventricle, also load 4th ventricle data
      if strcmp(FileID.uHEMs{iHEM},'4') && strcmp(FileID.uROIs{iROI},'vent')
        doit=1;
      end
      % If currently looking at left or right hemisphere, or at the 4th
      % ventricle, do the following
      if doit
        HEM = strcat(char('hem-'),FileID.uHEMs{iHEM});
        hem = FileID.uHEMs{iHEM};
        nm = strcat(id,roi,hem);
        disp(nm)
        % Load uncorrected and corrected data
        cd(datadir) 
        % Change this to the reconstructed images instead of the motion
        % corrected ones 
        % Currently, filename is r1corr.nii
        filenamerecon = [filetag,'corr.nii'];
        tmp=load_untouch_nii(filenamerecon);
        voxRes=double(tmp.hdr.dime.dim(2:4));
        subj_data_recon.(nm)=tmp.img;
        filenamegt = [filetag,'uncorr.nii'];
        tmp=load_untouch_nii(filenamegt);
        subj_data_gt.(nm)=tmp.img;

        % Load label map (assuming corrected and uncorrected maps identical)
        if strcmp(FileID.uROIs{iROI},cortag)
          ROI=FileID.uROIs{iROI};
          cd(lbldir)
          % Change to correct filenems for reconstruction
          mapname=[ID 't1corr_cruise-' ROI '_ants-def0.nii'];
          map_recon.(nm) = load_untouch_nii(mapname);
          mapname=[ID 't1uncorr_cruise-' ROI '_ants-def0.nii'];
          map_gt.(nm) = load_untouch_nii(mapname);
        else
          % If no label map is available? 
          cd(lbldir)
          IDROI = strjoin({ID, ROI},'');
          mapname = strjoin({IDROI, HEM, char('lvlreg-corr_def-img.nii')},'_');
          map_recon.(nm) = load_untouch_nii(mapname);
          mapname = strjoin({IDROI, HEM, char('lvlreg-uncorr_def-img.nii')},'_');
          map_gt.(nm) = load_untouch_nii(mapname);
        end

        % Alternative way to store MPos data, now with same indexing as
        % other data
        MPos.(nm) = mpars.MPos;
      end
    end
end

end
