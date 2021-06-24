function [ map_rim, subj_data_rim, map_gt, subj_data_gt, MPos,voxRes] = getdata( FileID, subjdir, varargin)
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
  filetag='r1';    
else
  filetag='t1';
end
disp(filetag)
cortag='gwb';

lbldir = fullfile(subjdir,'segm');

datadir = subjdir;

% Load motion parameters
cd(subjdir)
mparsname = dir('mpars.mat');
mpars = load(char(mparsname.name));
ID = strcat(char('sub-'),FileID.uIDs{1});
id = strcat('s',FileID.uIDs{1});

for iROI = 1:length(FileID.uROIs)
    if strcmp(FileID.uROIs{iROI},cortag)
      ROI=FileID.uROIs{iROI};
    else
      ROI = strcat(char('mask-'),FileID.uROIs{iROI});
    end
    roi = FileID.uROIs{iROI};
    for iHEM = 1:length(FileID.uHEMs)
      if iHEM<3, doit=1; else doit=0; end
      if strcmp(FileID.uROIs{iROI},cortag) && iHEM==2
        doit=0;
      end
      % for ventricle, also load 4th ventricle data
      if strcmp(FileID.uHEMs{iHEM},'4') && strcmp(FileID.uROIs{iROI},'vent')
        doit=1;
      end
      if doit
        HEM = strcat(char('hem-'),FileID.uHEMs{iHEM});
        hem = FileID.uHEMs{iHEM};
        nm = strcat(id,roi,hem);
        
        % Load uncorrected and corrected data
        cd(datadir)  
        filenamerim = [FileID.uIDs{1}, '_', FileID.type{1},'_rim_', FileID.accFactor{1},'.nii'];
        tmp=load_untouch_nii(filenamerim);
        voxRes=double(tmp.hdr.dime.dim(2:4));
        subj_data_rim.(nm)=flip(tmp.img,2);
       
        filenamegt = [FileID.uIDs{1},'_',FileID.type{1},'_gt.nii'];
        tmp=load_untouch_nii(filenamegt);
        subj_data_gt.(nm)=flip(tmp.img,2);
        
        % Load label map (assuming corrected and uncorrected maps identical)
        if strcmp(FileID.uROIs{iROI},cortag)
          ROI=FileID.uROIs{iROI};
          cd(lbldir)
          mapname=[ID '_t1corr_cruise-' ROI '_ants-def0.nii'];
          map_rim.(nm) = load_untouch_nii(mapname);
          mapname=[ID '_t1uncorr_cruise-' ROI '_ants-def0.nii'];
          map_gt.(nm) = load_untouch_nii(mapname);
        else
          cd(lbldir)
          mapname = strjoin({ID, ROI, HEM, char('lvlreg-gt_def-img.nii')},'_');
%           mapname = strjoin({ID, ROI, HEM, char('lvlreg-corr_def-img.nii')},'_');
          map_rim.(nm) = load_untouch_nii(mapname);
          map_rim.(nm).img = permute(map_rim.(nm).img, [2 3 1]);
          map_rim.(nm).img = map_rim.(nm).img(:,end:-1:1,end:-1:1);
          mapname = strjoin({ID, ROI, HEM, char('lvlreg-gt_def-img.nii')},'_');
%           mapname = strjoin({ID, ROI, HEM, char('lvlreg-uncorr_def-img.nii')},'_');
          map_gt.(nm) = load_untouch_nii(mapname);
          map_gt.(nm).img = permute(map_gt.(nm).img, [2 3 1]);
          map_gt.(nm).img = map_gt.(nm).img(:,end:-1:1,end:-1:1);
        end
        % Alternative way to store MPos data, now with same indexing as
        % other data
        MPos.(nm) = mpars.MPos;
      end
    end
end

end
