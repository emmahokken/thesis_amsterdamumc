% plot ROIs over parameter maps and save to disk
% plotsavedir must be defined and data loaded

nsubj=length(FileID.uIDs);

p_thresh=100; % intensity threshold for images (1 for subj 000, 100 for others)

savedir=fullfile(plotsavedir,'ROIs_anat');
mkdir(savedir)

for iSubj=1:nsubj
   subj = find(~cellfun(@isempty,strfind(fields,FileID.uIDs{iSubj})))';
   cl=hsv(length(subj));
   idx=0;
   for iROI=subj
     idx=idx+1;
     map_corr = map_corrall.(fields{iROI}).img;
     p=(map_corr<0) .* (map_corr>-2);
     p=p>0;
     if idx==1
       ROImap=double(p);
     else
       ROImap(p)=idx;
     end
   end
   sz=size(ROImap);
   ROImap=ROImap(:);
   d = t1_corrall.(fields{iROI});
   d(d>p_thresh)=p_thresh;d=d/p_thresh; % normalize and clip
   tmp=repmat(d(:),[1 3]);
   tmp(ROImap>0,:)=tmp(ROImap>0,:)+cl(ROImap(ROImap>0),:);
   tmp=reshape(tmp,[sz 3]);
   tmp=tmp/max(tmp(:));
   % permute to show reconsttucted images properly 
   tmp=permute(tmp,[3 1 2 4]);

   %ims=im2mat(arr(rgb(tmp(50:10:200,:,:,:)))); 
   %ima=im2mat(arr(rgb(tmp(:,:,:,:)))); 
   imsmall=im2mat(arr(rgb(tmp(100:170,100:180,20,:)))); 
   imbig=im2mat(arr(rgb(tmp(:,:,1:4:40,:)))); 

   savename_small=fullfile(savedir,[FileID.uIDs{iSubj} '_ROIs_small.png']);
   savename_big=fullfile(savedir,[FileID.uIDs{iSubj} '_ROIs_big.png']);

   imwrite(imsmall,savename_small);
   imwrite(imbig,savename_big);
end
