%% Open filters from all patients and assess PSF & CF characteristics
%==========================================================================

clc
clear all
cd('/home/aswinchari/Documents/NNAIS_0506/');
addpath('/home/aswinchari/Documents/NNAIS_0506/Code/FilterCode');
addpath('/home/aswinchari/Documents/GitHub/intersections');
addpath('/home/aswinchari/Documents/GitHub/edge-centric_demo-main/fcn/');

%% Create list of all neurons

filelist = dir('NNAIS*'); 

PSF = [];

for a = 1:length(filelist)
    
    filelist2 = dir(strcat(filelist(a).name,'/Spikes/*filters.mat'));
    
    for b = 1:length(filelist2)
        
        load(strcat(filelist2(b).folder,'/',filelist2(b).name));
    
        PSF = [PSF allclusters];
        
    end
end

clear a b allclusters

%% Filter out ones which are not single units & plot those properties

cols = cbrewer2('qual','Set2',3);

PSF = PSF([PSF.refrac] < 0.02');

for a = 1:length(PSF)
    PSF(a).cvISI = std(PSF(a).ISI)/mean(PSF(a).ISI);
end

clear a

allfreq = [PSF.freq];
allcvISI = [PSF.cvISI];
allspikeHW = [PSF.spikeHW]./30;

subplot(3,3,[1 2 3 4 5 6])
scatter3(allfreq(ismember([PSF.location],'cortex')),allcvISI(ismember([PSF.location],'cortex')),allspikeHW(ismember([PSF.location],'cortex')),20,cols(1,:),'filled');
hold on
scatter3(allfreq(ismember([PSF.location],'tuber')),allcvISI(ismember([PSF.location],'tuber')),allspikeHW(ismember([PSF.location],'tuber')),20,cols(2,:),'filled');
scatter3(allfreq(ismember([PSF.location],'tail')),allcvISI(ismember([PSF.location],'tail')),allspikeHW(ismember([PSF.location],'tail')),20,cols(3,:),'filled');
 xlabel('Frequency (Hz)');
 ylabel('Coeff Var ISI');
 zlabel('Spike Half Width (ms)');
 set(gca,'xscale','log')
 title('3D Plot of Spike Properties')
 legend({'Cortex', 'Tuber', 'Tail'})
 
subplot(3,3,7)
scatter(allfreq(ismember([PSF.location],'cortex')),allcvISI(ismember([PSF.location],'cortex')),20,cols(1,:),'filled');
hold on
scatter(allfreq(ismember([PSF.location],'tuber')),allcvISI(ismember([PSF.location],'tuber')),20,cols(2,:),'filled');
scatter(allfreq(ismember([PSF.location],'tail')),allcvISI(ismember([PSF.location],'tail')),20,cols(3,:),'filled');
 xlabel('Frequency (Hz)');
 ylabel('Coeff Var ISI');
 set(gca,'xscale','log')
 title('Frequency vs Coeff Var ISI')
 
subplot(3,3,8)
scatter(allfreq(ismember([PSF.location],'cortex')),allspikeHW(ismember([PSF.location],'cortex')),20,cols(1,:),'filled');
hold on
scatter(allfreq(ismember([PSF.location],'tuber')),allspikeHW(ismember([PSF.location],'tuber')),20,cols(2,:),'filled');
scatter(allfreq(ismember([PSF.location],'tail')),allspikeHW(ismember([PSF.location],'tail')),20,cols(3,:),'filled');
 xlabel('Frequency (Hz)');
 ylabel('Spike Half Width (ms)');
 set(gca,'xscale','log')
 title('Frequency vs Spike Half Width')
 
subplot(3,3,9)
scatter(allcvISI(ismember([PSF.location],'cortex')),allspikeHW(ismember([PSF.location],'cortex')),20,cols(1,:),'filled');
hold on
scatter(allcvISI(ismember([PSF.location],'tuber')),allspikeHW(ismember([PSF.location],'tuber')),20,cols(2,:),'filled');
scatter(allcvISI(ismember([PSF.location],'tail')),allspikeHW(ismember([PSF.location],'tail')),20,cols(3,:),'filled');
 xlabel('Coeff Var ISI');
 ylabel('Spike Half Width (ms)');
 title('Coeff Var ISI vs Spike Half Width')
 
set(gcf,'position',[20 20 1000 800]);
saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/SpikeProperties.png');

% T tests for this
% [p,tbl,stats] = anova1(allfreq,[PSF.location]);
% [p,tbl,stats] = anova1(allspikeHW,[PSF.location]);
% [p,tbl,stats] = anova1(allcvISI,[PSF.location]);

%% Plot PSF by location
  
% calculate AUC and correct the AUC to 1

for a = 1:length(PSF)
    
    PSF(a).area = sum(abs(cumtrapz(PSF(a).PSF)));
    tempPSF = PSF(a).PSF;
    PSF(a).normPSF = tempPSF./PSF(a).area;
    PSF(a).normarea = sum(abs(cumtrapz(PSF(a).normPSF)));
    clear tempPSF
    
end

% can show that AUC inversely proportional to freq
% scatter([PSF.freq],[PSF.area])
% xlabel('Freq');
% ylabel('PSF AUC');

% plot the PSFs

allPSF = [PSF.normPSF];

subplot(3,1,1)
plot(allPSF(:,ismember([PSF.location],'cortex')'));
title('Cortical PSFs')
ylim([-0.00002 0.00006])
xlabel('time (ms)')
ylabel('normalised gain')

subplot(3,1,2)
plot(allPSF(:,ismember([PSF.location],'tuber')'));
title('Tuber PSFs')
ylim([-0.00002 0.00006])
xlabel('time (ms)')
ylabel('normalised gain')

subplot(3,1,3)
plot(allPSF(:,ismember([PSF.location],'tail')'));
title('Tail PSFs')
ylim([-0.00002 0.00006])
xlabel('time (ms)')
ylabel('normalised gain')
  
set(gcf,'position',[20 20 1000 800]);
% saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/normalisedPSFs.png');

%% Are there differences in the PSFs? Need to do PCAs for this

[coeff, score, ~,~,explained,~] = pca(allPSF','Centered',false);

PC1{1} = score(ismember([PSF.location],'cortex'),1);
PC1{2} = score(ismember([PSF.location],'tuber'),1);
PC1{3} = score(ismember([PSF.location],'tail'),1);

PC2{1} = score(ismember([PSF.location],'cortex'),2);
PC2{2} = score(ismember([PSF.location],'tuber'),2);
PC2{3} = score(ismember([PSF.location],'tail'),2);

PC3{1} = score(ismember([PSF.location],'cortex'),3);
PC3{2} = score(ismember([PSF.location],'tuber'),3);
PC3{3} = score(ismember([PSF.location],'tail'),3);

PC4{1} = score(ismember([PSF.location],'cortex'),4);
PC4{2} = score(ismember([PSF.location],'tuber'),4);
PC4{3} = score(ismember([PSF.location],'tail'),4);

% set up colour vector
for a = 1:length(PSF)
    if PSF(a).location == string('cortex')      % green
        colvec(a,:) = cols(1,:);
        groups{a} = 'cortex';
    elseif PSF(a).location == string('tail')   % orange
        colvec(a,:) = cols(2,:);
        groups{a} = 'tail';
    else 
        colvec(a,:) = cols(3,:);                         % purple
        groups{a} = 'tuber';
    end
end

figure

subplot(3,6,[1 2 7 8 13 14])
scatter3(score(:,1),score(:,3),score(:,2),60,colvec,'filled')
xlabel('PC1')
ylabel('PC3')
zlabel('PC2')
title('Plot of first 3 PCs by group')

subplot(3,6,[3 9])
violin(PC1,'medc',[],'facecolor',cols,'facealpha',0.6);
scatter(ones(length(PC1{1}),1),PC1{1},[],cols(1,:),'filled')
scatter(ones(length(PC1{2}),1).*2,PC1{2},[],cols(2,:),'filled')
scatter(ones(length(PC1{3}),1).*3,PC1{3},[],cols(3,:),'filled')
legend('off')
title('Scores of PC1')
xticklabels({'cortex', 'tuber', 'tail'});
ylabel('PC1')

subplot(3,6,[4 10])
violin(PC2,'medc',[],'facecolor',cols,'facealpha',0.6);
scatter(ones(length(PC2{1}),1),PC2{1},[],cols(1,:),'filled')
scatter(ones(length(PC2{2}),1).*2,PC2{2},[],cols(2,:),'filled')
scatter(ones(length(PC2{3}),1).*3,PC2{3},[],cols(3,:),'filled')
legend('off');
title('Scores of PC2');
xticklabels({'cortex', 'tuber', 'tail'});
ylabel('PC2')

subplot(3,6,[5 11])
violin(PC3,'medc',[],'facecolor',cols,'facealpha',0.6);
scatter(ones(length(PC3{1}),1),PC3{1},[],cols(1,:),'filled')
scatter(ones(length(PC3{2}),1).*2,PC3{2},[],cols(2,:),'filled')
scatter(ones(length(PC3{3}),1).*3,PC3{3},[],cols(3,:),'filled')
legend('off');
title('Scores of PC3');
xticklabels({'cortex', 'tuber', 'tail'});
ylabel('PC3')

subplot(3,6,[6 12])
violin(PC4,'medc',[],'facecolor',cols,'facealpha',0.6);
scatter(ones(length(PC4{1}),1),PC4{1},[],cols(1,:),'filled')
scatter(ones(length(PC4{2}),1).*2,PC4{2},[],cols(2,:),'filled')
scatter(ones(length(PC4{3}),1).*3,PC4{3},[],cols(3,:),'filled')
legend('off');
title('Scores of PC4');
xticklabels({'cortex', 'tuber', 'tail'});
ylabel('PC4')

subplot(3,6,15)
plot(coeff(:,1),'k','LineWidth',3)
xlabel('time (ms)')
ylabel('gain')
yline(0)
title('Coefficients of PC1');

subplot(3,6,16)
plot(coeff(:,2),'k','LineWidth',3)
xlabel('time (ms)')
ylabel('gain')
yline(0)
title('Coefficients of PC2');

subplot(3,6,17)
plot(coeff(:,3),'k','LineWidth',3)
xlabel('time (ms)')
ylabel('gain')
yline(0)
title('Coefficients of PC3');

subplot(3,6,18)
plot(coeff(:,4),'k','LineWidth',3)
xlabel('time (ms)')
ylabel('gain')
yline(0)
title('Coefficients of PC4');

set(gcf,'position',[0,0,2000,1000])
saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/PCAuncentred.png');

% Statistical testing of this

% [p,tbl,stats] = anova1(score(:,1),[PSF.location]);
% [p,tbl,stats] = anova1(score(:,2),[PSF.location]);
% [p,tbl,stats] = anova1(score(:,3),[PSF.location]);
% [p,tbl,stats] = anova1(score(:,4),[PSF.location]);
%% Make node connectivity matrix & assess correlation of PSFs between locations

% correlation matrix from normalised dot product of PSFs
PSFFC = fcn_edgets2edgecorr(allPSF);

% for a = 1:length(PSF)
%     for b = 1:length(PSF)
%         PSFFC(a,b) = corr(PSF(a).normPSF,PSF(b).normPSF);
%     end
% end

% Do certain combos have different correlations from the rest?

cols = cbrewer2('qual','Dark2',6);

temp = PSFFC(ismember([PSF.location],'cortex'),ismember([PSF.location],'cortex'));
PSFcorr{1} = temp(triu(true(size(temp)),1));
clear temp

temp = PSFFC(ismember([PSF.location],'tuber'),ismember([PSF.location],'tuber'));
PSFcorr{2} = temp(triu(true(size(temp)),1));
clear temp

temp = PSFFC(ismember([PSF.location],'tail'),ismember([PSF.location],'tail'));
PSFcorr{3} = temp(triu(true(size(temp)),1));
clear temp

temp = PSFFC(ismember([PSF.location],'cortex'),ismember([PSF.location],'tuber'));
PSFcorr{4} = temp(:);
clear temp

temp = PSFFC(ismember([PSF.location],'cortex'),ismember([PSF.location],'tail'));
PSFcorr{5} = temp(:);
clear temp

temp = PSFFC(ismember([PSF.location],'tuber'),ismember([PSF.location],'tail'));
PSFcorr{6} = temp(:);
clear temp

% set up for ANOVA
KW = vertcat(PSFcorr{:});
groups = [repmat({'Cortex'},length(PSFcorr{1}),1); repmat({'Tuber'},length(PSFcorr{2}),1); repmat({'Tail'},length(PSFcorr{3}),1); repmat({'C-Tu'},length(PSFcorr{4}),1); repmat({'C-Ta'},length(PSFcorr{5}),1); repmat({'Tu-Ta'},length(PSFcorr{6}),1)];

[p, tbl, stats] = anova1(KW,groups);
[c,m] = multcompare(stats);

for a = 1:15
    psfsim(c(a,2),c(a,1)) = c(a,6);
end

psfsim = -log10(psfsim);
psfsim(isinf(psfsim))=0;
psfsim(:,6) = 0;

% plot this including the mean and SE

subplot(1,3,[1 2])
violin(PSFcorr,'medc',[],'facecolor',cols,'facealpha',0.6);
legend('off')
xticks([1:6])
xticklabels({''})
xticklabels({'Cortex-Cortex','Tuber-Tuber','Tail-Tail','Cortex-Tuber','Cortex-Tail','Tuber-Tail',});
title('PSF Similarity Within & Between Regions')
ylabel('Correlation')
ylim([-1.2 1.2])
xlim([0.5 6.5])
set(gca,'FontSize',10)

subplot(1,3,3)
imagesc(tril(psfsim))
colorbar
xlim([0.5 6.5])
ylim([0.5 6.5])
xticklabels({'Cortex-Cortex','Tuber-Tuber','Tail-Tail','Cortex-Tuber','Cortex-Tail','Tuber-Tail',});
xtickangle(90)
yticklabels({'Cortex-Cortex','Tuber-Tuber','Tail-Tail','Cortex-Tuber','Cortex-Tail','Tuber-Tail',});
xticks([1:6])
yticks([1:6])
title('Pairwise Comparisons (-log_1_0 of p-value)')
set(gca,'FontSize',10)

% subplot(3,1,3)
% errorbar([1:6],m(:,1),m(:,2).*1.96,'_');
% xlim([0.5 6.5])
% ylim([0.5 0.8])
% xticks([1:6])
% ylabel('Mean Correlation')
% xticklabels({'Cortex-Cortex','Tuber-Tuber','Tail-Tail','Cortex-Tuber','Cortex-Tail','Tuber-Tail',});
% xlabel('Regions')
% scatter(ones(length(PSFcorr{1}),1),PSFcorr{1},[],cols(1,:),'filled')
% scatter(ones(length(PSFcorr{2}),1).*2,PSFcorr{2},[],cols(2,:),'filled')
% scatter(ones(length(PSFcorr{3}),1).*3,PSFcorr{3},[],cols(3,:),'filled')
% scatter(ones(length(PSFcorr{4}),1).*4,PSFcorr{4},[],cols(4,:),'filled')
% scatter(ones(length(PSFcorr{5}),1).*5,PSFcorr{5},[],cols(5,:),'filled')
% scatter(ones(length(PSFcorr{6}),1).*6,PSFcorr{6},[],cols(6,:),'filled')

set(gcf,'position',[0,0,1600,400])
saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/FCCorr.png');

%% Ordering the PSFs to plot 

clearvars -except PSF allPSF

% order them by the area of the 100-130ms integral

for a = 1:length(PSF)
    
    PSF(a).thetaarea = sum(cumtrapz(PSF(a).normPSF(70:160)));
    
end

allPSF = [PSF.normPSF];
alllocation = [PSF.location];
allthetaarea = [PSF.thetaarea];

[~,order] = sort(allthetaarea,'descend');

orderedlocation = alllocation(order);
orderedPSF = allPSF(:,order);

% generate average PSFs by location

cortexPSF = orderedPSF(:,ismember(orderedlocation,'cortex'));
tuberPSF = orderedPSF(:,ismember(orderedlocation,'tuber'));
tailPSF = orderedPSF(:,ismember(orderedlocation,'tail'));

subplot(3,1,1)
imagesc(cortexPSF');
colorbar
caxis([-0.00005 0.00005])

subplot(3,1,2)
imagesc(tuberPSF');
colorbar
caxis([-0.00005 0.00005])

subplot(3,1,3)
imagesc(tailPSF');
colorbar
caxis([-0.00005 0.00005])

% Average histogram
% cent = linspace(-0.000005,0.000015,81);
% for a = 1:size(cortexPSF,1)
%     avcortexPSF(a,:) = hist(cortexPSF(a,:),cent);
% end
% imagesc(flipud(avcortexPSF'))
% caxis([0 25])

% Mean PSF +- CI
cols = cbrewer2('qual','Set2',3);

hold on
%cortex
patch([1:662 fliplr(1:662)], [(mean(cortexPSF,2)+std(cortexPSF,0,2))' (flipud(mean(cortexPSF,2)-std(cortexPSF,0,2)))'],cols(1,:),'FaceAlpha',0.2,'EdgeColor','none')
plot(mean(cortexPSF,2),'Color',cols(1,:),'LineWidth',2)
%tuber
patch([1:662 fliplr(1:662)], [(mean(tuberPSF,2)+std(tuberPSF,0,2))' (flipud(mean(tuberPSF,2)-std(tuberPSF,0,2)))'],cols(2,:),'FaceAlpha',0.2,'EdgeColor','none')
plot(mean(tuberPSF,2),'Color',cols(2,:),'LineWidth',2)
%tail
patch([1:662 fliplr(1:662)], [(mean(tailPSF,2)+std(tailPSF,0,2))' (flipud(mean(tailPSF,2)-std(tailPSF,0,2)))'],cols(3,:),'FaceAlpha',0.2,'EdgeColor','none')
plot(mean(tailPSF,2),'Color',cols(3,:),'LineWidth',2)