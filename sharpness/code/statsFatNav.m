FWHM=0.7*2.355; % conversion sigma to FWHM, given imaging resolution of 0.7 mm


%% calculate some pre-stats
CRthresh=pars.mrelcb;


% first corr, then uncorr
for ii=1:length(AR)
  if ~isempty(AR{ii})
    % compute CRLB
    cbrim=2*AR{ii}.errrim./AR{ii}.sigmarim;
    cbgt=2*AR{ii}.errgt./AR{ii}.sigmagt;
    medcb(ii,:)=[nanmedian(cbgt) nanmedian(cbrim)];

    % update aggregate test score (binary)
    tmp=abs(AR{ii}.aggregates{1});
    tmp4=tmp>=4;
    tmp=mod(tmp,4);
    test=abs(cbrim)>CRthresh;
    agg1=mod(tmp,2)-2*test-4*tmp4;
    tmp=abs(AR{ii}.aggregates{2}); 
    tmp4=tmp>=4;
    tmp=mod(tmp,4);
    test=abs(cbgt)>CRthresh;
    agg2=mod(tmp,2)-2*test-4*tmp4;
    validboth{ii}=and(~agg1,~agg2);

    sr = AR{ii}.sigmarim;
    er = AR{ii}.errrim;
    sg = AR{ii}.sigmagt;
    eg = AR{ii}.errgt;
    unc = sqrt(abs(1./(sg.^2-sr.^2).*(sg.^2.*eg.^2+sr.^2.*er.^2)));

    sigma_diff_abs = validboth{ii}.*sqrt(abs(sg.^2-sr.^2));
    signed_indication = 2.*((sg.^2-sr.^2)>0)-1; % 1 if better, -1 if worse
    sigma_diff_signed{ii} = signed_indication.*sigma_diff_abs;

    better_signed{ii} = validboth{ii}.*sigma_diff_signed{ii}.*((sg-eg)>(sr+er));
    worse_signed{ii} = validboth{ii}.*sigma_diff_signed{ii}.*((sr-er)>(sg+eg));
    uncertainty{ii} = validboth{ii}.*unc;
  end
end


nROI=length(mot_mean)/nsubj;

% relative numberr of valid clusters
relvalid=cellfun(@sum,validboth)./cellfun(@length,validboth);
relvalid=reshape(relvalid,nROI,[]);
% absolute number of valid clusters 
validClusters=cellfun(@sum,validboth);
validClusters=reshape(validClusters,nROI,[]);

% number of clusters
ln=ones(length(AR),1); ln(:)=nan;

for ii=1:length(AR)
  try
  ln(ii)=length(AR{ii}.sigmarim);
  end
end

ln=reshape(ln,nROI,[]);


allsigma = []; allsigmarim=[]; allsigmagt=[];
allmotion = []; allu=[];
subj_name = cell(1,length(doSubj));
clear csigma cu
plotIDs=[];
i=iSubj;
  subj = find(~cellfun(@isempty,strfind(fields,FileID.uIDs{i})));
  plot_sigma = []; plot_sigma_rim=[]; plot_sigma_gt=[];
  plot_motion = []; plot_u=[];
  for j = 1:nROI
      inds = j;
      fieldinds = subj(inds);
      sigm_plot = []; sigm_mean_plot=[]; 
      sigm_rim=[]; sigm_gt=[]; 
        
      mot_plot = []; mot_mean_plot=[];
      k = 1;
      fieldind = fieldinds(k);
      sigm_data = better_signed{fieldind}+worse_signed{fieldind};
      signif = better_signed{fieldind}|worse_signed{fieldind};
      tmp=sigma_diff_signed{fieldind}(validboth{fieldind}>0);
      u=uncertainty{fieldind}(validboth{fieldind}>0);
      try
        csigma{i}=cat(2,csigma{i},FWHM*tmp);
        cu{i}=cat(2,cu{i},FWHM*u);
      catch
        csigma{i}=FWHM*tmp;
        cu{i}=FWHM*u;
      end
      sigm_mean_plot(k)=FWHM*median(tmp);
      uplot(k)=FWHM*nanmedian(u);
      tmp=AR{fieldind}.sigmarim(validboth{fieldind}>0);     % sigma after
      sigm_rim(k)=FWHM*median(tmp);
      tmp=AR{fieldind}.sigmagt(validboth{fieldind}>0);       % sigma after
      sigm_gt(k)=FWHM*median(tmp);
      mot_mean_plot(k) = (mot_mean{fieldind});
      mot_add = ones(sum(signif),1)*mot_mean{fieldind}; 
      mot_plot = cat(2,mot_plot,mot_add');
      plot_sigma = cat(2,plot_sigma,sigm_mean_plot);
    
      plot_u=cat(2,plot_u,uplot);
      plot_motion = cat(2,plot_motion,mot_mean_plot);
      plot_sigma_rim = cat(2,plot_sigma_rim,sigm_rim);
      plot_sigma_gt = cat(2,plot_sigma_gt,sigm_gt);
%       hold off
%       plot(plot_sigma_rim);
%       hold on
%       plot(plot_sigma_gt);
%       legend('rim', 'gt');
%       hold off
  end
  
  allsigma = cat(2,allsigma,plot_sigma);
  allu=cat(2,allu,plot_u);
  allsigmarim = cat(2,allsigmarim,plot_sigma_rim);
  allsigmagt = cat(2,allsigmagt,plot_sigma_gt);
  allmotion = cat(2,allmotion,plot_motion);
  
csigma=csigma(doSubj); cu=cu(doSubj);

%% Plot stats
plotfields = {};
for i =1:length(fields)
    plotfields{i} = fields{i}(5:end);
end

barColors = {[102 51 153] / 255, [255 165 0] / 255};
figure('visible', 'off');
b = bar(categorical(plotfields),[allsigmagt;allsigmarim].');
set(b,{'DisplayName'},{'Ground truth','RIM'}')
legend()
for i =1:2
    set(b(i),'CData',barColors{i});
    set(b(i),'FaceColor','flat');
end
title(strcat('FWHM of different structures for', {' '},FileID.accFactor{1}, 'x acceleration'));
ylabel('Sharpness in FWHM')
barFileName = strcat('../../plots_saved/', FileID.uIDs{1},'_',FileID.accFactor{1},'_FWHM_barchart.png');
saveas(gcf,barFileName);

%% Display stats

disp('ROIs:')
disp(fields')

disp('FWHM ground truth:')
disp(num2str(allsigmagt))
disp('FWHM reconstructed (rim):')
disp(num2str(allsigmarim))
disp('FWHM difference:')
disp(num2str(allsigma))

disp('Number of total clusters:')
disp(num2str(ln'))
disp('Number of valid clusters:')
disp(num2str(validClusters'))

disp('Relative nr of valid clusters:')
disp(num2str(relvalid'))

%% Create table for satistical analysis 

% Write ground truth data once to file 

sigma_rim = reshape(allsigmarim,[],1);
sigma_gt = reshape(allsigmagt,[],1);

acc_factor = zeros(size(sigma_rim));
acc_factor(acc_factor==0) = str2num(FileID.accFactor{1});
subj_id = zeros(size(sigma_rim));
subj_id(subj_id==0) = str2num(FileID.uIDs{1});

tabl = table(subj_id,fields, sigma_rim, sigma_gt, acc_factor)
writetable(tabl,strcat('../../results/',FileID.uIDs{1},'_',FileID.accFactor{1},'_FWHM.csv'),'WriteRowNames',true)
