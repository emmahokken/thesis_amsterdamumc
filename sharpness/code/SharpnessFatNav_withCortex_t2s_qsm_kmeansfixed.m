%clear all
%close all
function SharpnessFatNav_withCortex_t2s_qsm

addpath('/packages/matlab/toolbox/spm8/r6313')
addpath(genpath('/data1/projects/FatNav/Hannah/matlab/SharpnessFatNav'))


% set r2s (or qsm) here, also in plotsavedir
%r2s or qsm_inv2_te2_m_
%doR1=1;

doR1='r2s'
plotsavedir='/data1/projects/FatNav/Hannah/Plots_saved_t2s' % main dir for saving ploits

doR1='qsm'
plotsavedir='/data1/projects/FatNav/Hannah/Plots_saved_qsm' % main dir for saving ploits

doR1='mt1'
plotsavedir='/data1/projects/FatNav/Hannah/Plots_saved_mt1' % main dir for saving ploits

mkdir(plotsavedir)
mkdir(fullfile(plotsavedir,'Profiles'))

doDebug=0 % set to 1 to run only 2 subjects

skipSubjects=[9 15 20 26] % subjects to skip from analyses because of erroneous ROIs

%% Get all data

[FileID.uIDs, subjdir] = getIDs();
%
FileID.uROIs = {'vent', 'tha', 'str','gwb'};
%FileID.uROIs = {'gwb'};
%FileID.uHEMs = {'lr4', 'lr','lr',''};
FileID.uHEMs = {'l', 'r','4'};


% useSubj= 6
% FileID.uIDs = {FileID.uIDs{useSubj}};
% subjdir=subjdir(useSubj);

[map_corrall, t1_corrall, map_uncorrall, t1_uncorrall, MPos] = ...
  getdata_cortex(FileID, subjdir, doR1);
fields = fieldnames(map_corrall);

%% Get parameters

pars = getpars();

pars.mrelcb=.15 % FIXME: lower CRLB for cortex (maybe for all?)
pars.mrelcb=.5
%pars.mrelcb=1

%% plot ROIs and save to pngs
if 0
  run FatNav_plot_ROIs
end

%% initialize
if doDebug, 
  % 6 ROIs, two subjects
  %idx=1:12; 
  doSubj=1:2;
else
  %idx=1:length(fields); 
  %doSubj=1:length(FileID.uIDs)
  
  % define subjects to process
  nsubj=length(FileID.uIDs);
  doSubj=ones(1,nsubj);
  doSubj(skipSubjects)=0;
  doSubj=find(doSubj);
end 

better_signed = cell(1,length(fields));
worse_signed = cell(1,length(fields));
AR = cell(1,length(fields));

%% run
%for iSubj=10%doSubj
savedir='/data1/projects/FatNav/matthan';
saveID=['SharpessFatNav_' doR1 '_' randomtime];
matfile=fullfile(savedir,saveID)

for iSubj=doSubj%[21:25 27 28]
  subj_ii = find(~cellfun(@isempty,strfind(fields,FileID.uIDs{iSubj})))';
  iSubj
  for ii=subj_ii
    try
      ii
      
      uID=['s' FileID.uIDs{iSubj}];
      % get data
      field_name = fields{ii};
      map_corr = map_corrall.(fields{ii}).img;
      %data_corr = t1_corrall.(fields{ii});
      data_corr = t1_corrall.(uID);

      map_uncorr = map_uncorrall.(fields{ii}).img;
      %data_uncorr = t1_uncorrall.(fields{ii});
      data_uncorr = t1_uncorrall.(uID);

      % do fitting
      %if length(brd_crds)<=ii || isempty(brd_crds{ii})
      %MC: add subj ID for saving kmeans
      pars.uID=fields{ii}
      pars.savedir=savedir
      [better_signed{ii},worse_signed{ii},brd_crds{ii},brd_ind{ii},AR{ii}] = getsigmacluster(map_corr,data_corr,data_uncorr,pars,field_name,plotsavedir);
      %end
      % get motion parameters
      numcl(ii) = length(better_signed{ii});
      mot_mean_new{ii} = getmotionvector_editMatthan(map_corrall.(fields{ii}).img,brd_crds{ii},brd_ind{ii},MPos.(fields{ii}),numcl(ii));

      arr = AR{ii}.sigmadiff_signed.*AR{ii}.validboth;
      brd_sigma{ii} = arr(brd_ind{ii});
      arr2 = mot_mean_new{ii}.*AR{ii}.validboth;
      brd_motion{ii} = arr2(brd_ind{ii});

    catch
      disp(['Error for subject ' FileID.uIDs{iSubj} ' ii ' num2str(ii)])
      err=lasterror
      allerr{ii}=err;
    end
  end
  try
    disp('saving...')
    save(matfile,'better_signed','worse_signed','brd_crds','brd_ind','AR','mot_mean_new')
    disp('done.')
  catch
    disp('error saving')
  end
end



%%
% mot_mean = cell(1,length(idx));
% for ii = idx
%     mot_mean{ii} = getmotionvector(map_corrall.(fields{ii}).img,brd_crds{ii},brd_ind{ii},MPos.(fields{ii}),numcl(ii)); 
% end
% %%
% for ii = idx
%     mot_mean_new{ii} = getmotionvector_editMatthan(map_corrall.(fields{ii}).img,brd_crds{ii},brd_ind{ii},MPos.(fields{ii}),numcl(ii));
% end
% %%
% for ii = idx
%     arr = AR{ii}.sigmadiff_signed.*AR{ii}.validboth;
%     brd_sigma{ii} = arr(brd_ind{ii});
%     arr2 = mot_mean{ii}.*AR{ii}.validboth;
%     brd_motion{ii} = arr2(brd_ind{ii});
% end
%%
plotSubj=1; % define subject for plotting

plot_ii = find(~cellfun(@isempty,strfind(fields,FileID.uIDs{plotSubj})))'; % ROI indices

%run Clustersplot % Overview of clustering for all ROIs: intensity clusters, subclusters and validity indication
%run Sigmasplot % Detailed plot with errorbars of sigma for every cluster before and after correction. Only the left ventricle.
%run ROIsigmaplot % Boxplot of sigma before and after correction for all ROIs

%%
% run Wilcoxonplot % Boxplot of improvement/deterioration and wilcoxon test plot. 
% %run Motionvsimprovementplot % Plot motion vs improvement using old motion method
% run Motionvsimprovementplot_new % Plot motion vs improvement using new motion method
% run Motionseparate % Plots motion vs improvement for all subjects separately. Uses old motion method data


% TODO: create separate plots save dir, and doSave=1 in script below
Motionvsimprovementplot_cortex_t2s


error('stop here')

%% re-evaluate CRLB-validity
doSave=0
doSkipGMB=0 % skip gray/WM boundary
CRthresh=0.5;

% first corr, then uncorr
for ii=1:length(AR)
  if ~isempty(AR{ii})
    % compute CRLB
    cbcorr=2*AR{ii}.errcorr./AR{ii}.sigmacorr;
    cbuncorr=2*AR{ii}.erruncorr./AR{ii}.sigmauncorr;
    
    medcb(ii,:)=[nanmedian(cbuncorr) nanmedian(cbcorr)];

    % code from getsigmacluster.m
    % update aggregate test score (binary)
    tmp=abs(AR{ii}.aggregates{1}); 
    tmp4=tmp>=4;
    tmp=mod(tmp,4);
    test=abs(cbcorr)>CRthresh;
    agg1=mod(tmp,2)-2*test-4*tmp4;
    tmp=abs(AR{ii}.aggregates{2}); 
    tmp4=tmp>=4;
    tmp=mod(tmp,4);
    test=abs(cbuncorr)>CRthresh;
    agg2=mod(tmp,2)-2*test-4*tmp4;
    validboth{ii}=and(~agg1,~agg2);
    
    sc = AR{ii}.sigmacorr;%sigma_comp{1};
    ec = AR{ii}.errcorr;%err_comp{1};
    su = AR{ii}.sigmauncorr;%sigma_comp{2};
    eu = AR{ii}.erruncorr;%err_comp{2};
    unc = sqrt(abs(1./(su.^2-sc.^2).*(su.^2.*eu.^2+sc.^2.*ec.^2)));

    %sigma_diff_abs=AR{ii}.sigmadiff;

    sigma_diff_abs = validboth{ii}.*sqrt(abs(su.^2-sc.^2));
    signed_indication = 2.*((su.^2-sc.^2)>0)-1; % 1 if better, -1 if worse
    sigma_diff_signed{ii} = signed_indication.*sigma_diff_abs;

    better_signed{ii} = validboth{ii}.*sigma_diff_signed{ii}.*((su-eu)>(sc+ec));
    worse_signed{ii} = validboth{ii}.*sigma_diff_signed{ii}.*((sc-ec)>(su+eu));
    uncertainty{ii} = validboth{ii}.*unc;
  end
end
%cb_rel = 2*err(i)./sigma(i);

% relative nr of valid clusters
relvalid=cellfun(@sum,validboth)./cellfun(@length,validboth);
relvalid=reshape(relvalid,length(validboth),[])
nanmedian(relvalid,2)
nanstd(relvalid,[],2)

% number of clusters
ln=ones(length(AR),1); ln(:)=nan;
for ii=1:length(AR)
  try
  ln(ii)=length(AR{ii}.sigmacorr);
  end
end
ln=reshape(ln,length(validboth),[]);
nanmean(ln,2)
nanstd(ln,[],2)


FWHM=0.7*2.355; % conversion sigma to FWHM
% 
% if doSkipGMB
%   nROI=length(subj)-1
% else
%   nROI=length(subj)
% end
nROI=1

%figure
%colors = colormap(hsv(length(FileID.uIDs)));
colors = colormap(hsv(max(doSubj)));
allsigma = []; allsigmacorr=[]; allsigmauncorr=[];
allmotion = []; allu=[];
subj_name = cell(1,length(doSubj));
clear csigma cu
%ii=1;
%for i = 1:length(FileID.uIDs)
plotIDs=[];
for i = doSubj%1:length(FileID.uIDs)
  %subj = find(contains(fields,FileID.uIDs{i}));
  % MATLAB R2014a
  subj = find(~cellfun(@isempty,strfind(fields,FileID.uIDs{i})));
  color = colors(i,:);
  plot_sigma = []; plot_sigma_corr=[]; plot_sigma_uncorr=[];
  plot_motion = []; plot_u=[];
  %sigm_mean_plot = zeros(1,length(FileID.uROIs));
  %mot_mean_plot = zeros(1,length(FileID.uROIs));
 %for j = 1:length(FileID.uROIs)
  for j = 1:nROI
      inds = j;%[2*j-1,2*j];
      fieldinds = subj(inds);
      sigm_plot = []; sigm_mean_plot=[]; 
      sigm_corr=[]; sigm_uncorr=[];
      mot_plot = []; mot_mean_plot=[];
      k = 1;%:2 %left/right
              fieldind = fieldinds(k);
              sigm_data = better_signed{fieldind}+worse_signed{fieldind};
              signif = better_signed{fieldind}|worse_signed{fieldind};

              %sigm_mean_plot(ii) = mean(sigm_data(sigm_data>0));
              %tmp=AR{fieldind}.sigmadiff_signed;
              %tmp=AR{fieldind}.sigmadiff_signed(AR{fieldind}.validboth>0);% sigma difference
              %tmp=AR{fieldind}.sigmadiff_signed(validboth{fieldind}>0);% sigma difference
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
              %sigm_mean_plot(ii)=median(tmp);
              %tmp=AR{fieldind}.sigmacorr(AR{fieldind}.validboth>0);       % sigma after
              tmp=AR{fieldind}.sigmacorr(validboth{fieldind}>0);       % sigma after
              sigm_corr(k)=FWHM*median(tmp);
              %tmp=AR{fieldind}.sigmauncorr(AR{fieldind}.validboth>0);       % sigma after
              tmp=AR{fieldind}.sigmauncorr(validboth{fieldind}>0);       % sigma after
              sigm_uncorr(k)=FWHM*median(tmp);

              %sigm_plot = cat(2,sigm_plot,sigm_data(signif));
              mot_mean_plot(k) = (mot_mean_new{fieldind});
              mot_add = ones(sum(signif),1)*mot_mean_new{fieldind}; 
              mot_plot = cat(2,mot_plot,mot_add');
              %ii = ii+1;
       %plot_sigma = cat(2,plot_sigma,sigm_plot);
      %plot_motion = cat(2,plot_motion,mot_plot);
      plot_sigma = cat(2,plot_sigma,sigm_mean_plot);
      plot_u=cat(2,plot_u,uplot);
      plot_motion = cat(2,plot_motion,mot_mean_plot);
      plot_sigma_corr = cat(2,plot_sigma_corr,sigm_corr);
      plot_sigma_uncorr = cat(2,plot_sigma_uncorr,sigm_uncorr);
  end
  subj_name{i} = strcat('sub-',FileID.uIDs{i});
  allsigma = cat(2,allsigma,plot_sigma);
  allu=cat(2,allu,plot_u);
  allsigmacorr = cat(2,allsigmacorr,plot_sigma_corr);
  allsigmauncorr = cat(2,allsigmauncorr,plot_sigma_uncorr);
  allmotion = cat(2,allmotion,plot_motion);
  plotIDs=cat(2,plotIDs,ones(size(plot_sigma))*i);
end
csigma=csigma(doSubj); cu=cu(doSubj);

nSubj=length(doSubj);
idx=(1:nROI)'*ones(1,nSubj)


allmotion=reshape(allmotion,nROI,[]);
allsigmauncorr=reshape(allsigmauncorr,nROI,[]);
allsigmacorr=reshape(allsigmacorr,nROI,[]);
allsigma=reshape(allsigma,nROI,[]);

lgd=regexprep(fields(1:nROI),'s004','');

figure
subplot(1,3,1)
xx=median(allmotion,1);
yy=median(reshape(allsigmauncorr,nROI,[]),1);
%scatter(xx,yy,30,[0 0 0],'filled')
plot(xx,yy,'sk')
hold on
plot([0 2],[0.7 0.7],':','color',[.5 .5 .5])
%set(gca,'ylim',[0.5 1.8])
set(gca,'fontsize',14)
xlabel('motion (mm)')
ylabel('FWHM_{uncorr} (mm)')
title('uncorrected')
  [rho,pval] = corr(xx',yy','Type', 'Spearman')
  sign=''
  if pval<.05
    sign='*';
  end
  t=sprintf('%.2f',rho)
  text(0.5,1.5,['\rho = ' t sign],'interpreter','tex','fontsize',14)

subplot(1,3,2)
yy=median(reshape(allsigmacorr,nROI,[]),1);
%scatter(xx,yy,30,[0 0 0],'filled')
plot(xx,yy,'ok')
hold on
plot([0 2],[0.7 0.7],':','color',[.5 .5 .5])
%set(gca,'ylim',[0.5 1.8])
set(gca,'fontsize',14)
xlabel('motion (mm)')
ylabel('FWHM_{corr} (mm)')
title('corrected')
  [rho,pval] = corr(xx',yy','Type', 'Spearman')
  sign=''
  if pval<.05
    sign='*';
  end
  t=sprintf('%.2f',rho)
  text(0.5,1.5,['\rho = ' t sign],'interpreter','tex','fontsize',14)

subplot(1,3,3)
yy=median(reshape(allsigma,nROI,[]),1)
%hnd(2)=scatter(xx,yy,30,[0 0 0],'filled')
plot(xx,yy,'dk','markersize',7)
hold on
%plot([0 1.8],[0 0],':','color',[.5 .5 .5])
set(gca,'fontsize',14)
xlabel('motion (mm)')
ylabel('FWHM_{diff} (mm)')
%legend(hnd,{'all ROIs','subject median','negative ROIs'},'location','southeast')
title('improvement')
  [rho,pval] = corr(xx',yy','Type', 'Spearman')
  sign=''
  if pval<.05
    sign='*';
  end
  t=sprintf('%.2f',rho)
  text(0.5,1.15,['\rho = ' t sign],'interpreter','tex','fontsize',14)
% 
% if doSave
%   set(gcf,'papersize',[11 8.5])
%   p=get(gcf,'paperposition')
%   p(3:4)=[10 3]
%   set(gcf,'paperposition',p)
%   print(gcf,'fatnav_scatter.pdf','-dpdf')
% end


prctile(allsigmauncorr,50)
[prctile(allsigmauncorr,75)-prctile(allsigmauncorr,25)]


prctile(allsigmacorr,50)
[prctile(allsigmacorr,75)-prctile(allsigmacorr,25)]
