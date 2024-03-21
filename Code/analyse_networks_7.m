%% Analyse networks from CFs and PSFs
%==========================================================================

clc
clear all
cd('/home/aswinchari/Documents/NNAIS_0506/Networks');
load('networks.mat');

%% Calculate some nodal graph metrics for the networks

for a = 1:length(networks)
    
    disp(strcat('Network:',string(a),'of 7'));
    
    networks(a).PSFstrength         = strengths_und(networks(a).PSFnetwork);
    [networks(a).PSFcc,~,~,~]       = clustering_coef_wu_sign(networks(a).PSFnetwork,3);
    networks(a).PSFeigencent        = eigenvector_centrality_und(networks(a).PSFnetwork);
    
    networks(a).CFstrength         = strengths_und(networks(a).CFnetwork);
    [networks(a).CFcc,~,~,~]       = clustering_coef_wu_sign(networks(a).CFnetwork,3);
    networks(a).CFeigencent        = eigenvector_centrality_und(networks(a).CFnetwork);
    
    networks(a).CFstrength         = networks(a).CFstrength'; 
    networks(a).PSFstrength        = networks(a).PSFstrength';
end

save('networks.mat','networks');

%% Open networks and assess differences in graph metrics

clc
clear all
cd('/home/aswinchari/Documents/NNAIS_0506/Networks');
load('networks.mat');

%% Calculate percentile ranks of metrics

for a = 1:length(networks)
    
    % strength
    [~, temp] = sort(networks(a).PSFstrength);
    rankV(temp) = 1:numel(networks(a).PSFstrength);
    networks(a).PSFstrengthrank = (rankV./length(rankV))';
    clear temp rankV
    
    % cc
    [~, temp] = sort(networks(a).PSFcc);
    rankV(temp) = 1:numel(networks(a).PSFcc);
    networks(a).PSFccrank = (rankV./length(rankV))';
    clear temp rankV
    
    % eigencent
    [~, temp] = sort(networks(a).PSFeigencent);
    rankV(temp) = 1:numel(networks(a).PSFeigencent);
    networks(a).PSFeigencentrank = (rankV./length(rankV))';
    clear temp rankV
    
    
    % CFstrength
    [~, temp] = sort(networks(a).CFstrength);
    rankV(temp) = 1:numel(networks(a).CFstrength);
    networks(a).CFstrengthrank = (rankV./length(rankV))';
    clear temp rankV
    
    % CFcc
    [~, temp] = sort(networks(a).CFcc);
    rankV(temp) = 1:numel(networks(a).CFcc);
    networks(a).CFccrank = (rankV./length(rankV))';
    clear temp rankV
    
    % CFeigencent
    [~, temp] = sort(networks(a).CFeigencent);
    rankV(temp) = 1:numel(networks(a).CFeigencent);
    networks(a).CFeigencentrank = (rankV./length(rankV))';
    clear temp rankV
end
    
clear a

%% PSF graph metrics

allPSFstrength = vertcat(networks.PSFstrengthrank);
allPSFcc = vertcat(networks.PSFccrank);
allPSFeigencent = vertcat(networks.PSFeigencentrank);
allPSFlocation = vertcat(networks.PSFlocation);

[p,tbl,stats] = anova1(allPSFstrength,allPSFlocation);
[p,tbl,stats] = anova1(allPSFcc,allPSFlocation);
[p,tbl,stats] = anova1(allPSFeigencent,allPSFlocation);


%% CF graph metrics

allCFstrength = vertcat(networks.CFstrengthrank);
allCFcc = vertcat(networks.CFccrank);
allCFeigencent = vertcat(networks.CFeigencentrank);
allCFlocation = vertcat(networks.CFlocation);

[p,tbl,stats1] = anova1(allCFstrength,allCFlocation);
[p,tbl,stats2] = anova1(allCFcc,allCFlocation);
[p,tbl,stats3] = anova1(allCFeigencent,allCFlocation);
c1 = multcompare(stats1);
c2 = multcompare(stats2);
c3 = multcompare(stats3);

for a = 1:36
    strengthcompare(c1(a,2),c1(a,1)) = c1(a,6);
    cccompare(c2(a,2),c2(a,1)) = c2(a,6);
    eigencentcompare(c3(a,2),c3(a,1)) = c3(a,6);
end

% prepare p-value matrices for plots

strengthcompare = -log10(strengthcompare);
strengthcompare(isinf(strengthcompare))=0;
strengthcompare(:,9) = 0;

cccompare = -log10(cccompare);
cccompare(isinf(cccompare))=0;
cccompare(:,9) = 0;

eigencentcompare = -log10(eigencentcompare);
eigencentcompare(isinf(eigencentcompare))=0;
eigencentcompare(:,9) = 0;


%% Graph the CF metrics

STR{1} = allCFstrength(ismember([allCFlocation],'cortex-cortex'));
STR{2} = allCFstrength(ismember([allCFlocation],'cortex-tuber'));
STR{3} = allCFstrength(ismember([allCFlocation],'cortex-tail'));
STR{4} = allCFstrength(ismember([allCFlocation],'tuber-cortex'));
STR{5} = allCFstrength(ismember([allCFlocation],'tuber-tuber'));
STR{6} = allCFstrength(ismember([allCFlocation],'tuber-tail'));
STR{7} = allCFstrength(ismember([allCFlocation],'tail-cortex'));
STR{8} = allCFstrength(ismember([allCFlocation],'tail-tuber'));
STR{9} = allCFstrength(ismember([allCFlocation],'tail-tail'));

CC{1} = allCFcc(ismember([allCFlocation],'cortex-cortex'));
CC{2} = allCFcc(ismember([allCFlocation],'cortex-tuber'));
CC{3} = allCFcc(ismember([allCFlocation],'cortex-tail'));
CC{4} = allCFcc(ismember([allCFlocation],'tuber-cortex'));
CC{5} = allCFcc(ismember([allCFlocation],'tuber-tuber'));
CC{6} = allCFcc(ismember([allCFlocation],'tuber-tail'));
CC{7} = allCFcc(ismember([allCFlocation],'tail-cortex'));
CC{8} = allCFcc(ismember([allCFlocation],'tail-tuber'));
CC{9} = allCFcc(ismember([allCFlocation],'tail-tail'));

EIG{1} = allCFeigencent(ismember([allCFlocation],'cortex-cortex'));
EIG{2} = allCFeigencent(ismember([allCFlocation],'cortex-tuber'));
EIG{3} = allCFeigencent(ismember([allCFlocation],'cortex-tail'));
EIG{4} = allCFeigencent(ismember([allCFlocation],'tuber-cortex'));
EIG{5} = allCFeigencent(ismember([allCFlocation],'tuber-tuber'));
EIG{6} = allCFeigencent(ismember([allCFlocation],'tuber-tail'));
EIG{7} = allCFeigencent(ismember([allCFlocation],'tail-cortex'));
EIG{8} = allCFeigencent(ismember([allCFlocation],'tail-tuber'));
EIG{9} = allCFeigencent(ismember([allCFlocation],'tail-tail'));

cols = cbrewer2('seq','YlGnBu',9); 

subplot(3,3,[1 2])
violin(STR,'medc',[],'facecolor',cols,'facealpha',0.6);
legend('off')
title('Strength')
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xticks(1:9)
xtickangle(90)

subplot(3,3,[4 5])
violin(CC,'medc',[],'facecolor',cols,'facealpha',0.6);
legend('off')
title('Clustering Coefficient')
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xticks(1:9)
xtickangle(90)

subplot(3,3,[7 8])
violin(EIG,'medc',[],'facecolor',cols,'facealpha',0.6);
legend('off')
title('Eigenvector Centrality')
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xticks(1:9)
xtickangle(90)

subplot(3,3,3)
imagesc(tril(strengthcompare))
colorbar
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
yticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xtickangle(90)
xticks([1:9])
yticks([1:9])
title('Pairwise Comparisons (-log_1_0 of p-value)')

subplot(3,3,6)
imagesc(tril(cccompare))
colorbar
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
yticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xtickangle(90)
xticks([1:9])
yticks([1:9])
title('Pairwise Comparisons (-log_1_0 of p-value)')

subplot(3,3,9)
imagesc(tril(eigencentcompare))
colorbar
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
yticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xtickangle(90)
xticks([1:9])
yticks([1:9])
title('Pairwise Comparisons (-log_1_0 of p-value)')

set(gcf,'position',[0,0,1200,1000])
saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/CFgraphmetrics.png');