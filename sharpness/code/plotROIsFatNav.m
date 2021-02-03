% plot ROIs over parameter maps and save to disk
% plotsavedir must be defined and data loaded

nsubj=length(FileID.uIDs);

p_thresh=1.1; % intensity threshold for images

savedir=fullfile(plotsavedir,'ROIs_anat');
mkdir(savedir)

for iSubj=1:nsubj
    disp(iSubj)
    disp(FileID.uIDs)
    
    disp(fields)    
   subj = find(~cellfun(@isempty,strfind(fields,FileID.uIDs{iSubj})))';
   cl=hsv(length(subj));
   idx=0;
   
   for iROI=subj
     idx=idx+1;
     map_recon = map_reconall.(fields{iROI}).img;
     p=(map_recon<0) .* (map_recon>-2);
     p=p>0;
     if idx==1
       ROImap=double(p);
     else
       ROImap(p)=idx;
     end
   end
   
   cd('../')
   
   sz=size(ROImap);
   ROImap=ROImap(:);
   d = t1_reconall.(fields{iROI});
   d(d>p_thresh)=p_thresh;d=d/p_thresh; % normalize and clip
   tmp=repmat(d(:),[1 3]);
   tmp(ROImap>0,:)=tmp(ROImap>0,:)+cl(ROImap(ROImap>0),:);
   tmp=reshape(tmp,[sz 3]);
   tmp=tmp/max(tmp(:));
%    make color image x,y,z,3
   im=im2mat(arr(rgb(tmp(:,:,50:10:200,:)))); 
%    im=im2mat(rgb(tmp(:,:,50:10:200,:))); 

   savename=fullfile(savedir,[FileID.uIDs{iSubj} '_ROIs.png']);
   imwrite(im,savename);
end
