%% Open filters from all patients and assess CF characteristics
%==========================================================================

clc
clear all
cd('/home/aswinchari/Documents/NNAIS_0506/');
addpath('/home/aswinchari/Documents/NNAIS_0506/Code/FilterCode');
addpath('/home/aswinchari/Documents/GitHub/intersections');
addpath('/home/aswinchari/Documents/GitHub/edge-centric_demo-main/fcn/');

%% Create list of all eligible CFs

filelist = dir('NNAIS*'); 

CF = [];
CFlocation = [];

for a = 1:length(filelist)
    
    filelist2 = dir(strcat(filelist(a).name,'/Spikes/*filters.mat'));
    
    for b = 1:length(filelist2)
        
        load(strcat(filelist2(b).folder,'/',filelist2(b).name));
        
        include = [allclusters.refrac] < 0.02';
        allclusters = allclusters(include);
        
            for c = 1:length(allclusters)
        
                % normalise CFs

                allclusters(c).CF = allclusters(c).CF(:,include);

                for d = 1:size(allclusters(1).CF,2)     
                    temparea = sum(abs(cumtrapz(allclusters(c).CF(:,d))));
                    tempCF = allclusters(c).CF(:,d);
                    allclusters(c).normCF(:,d) = tempCF./temparea;
                    clear temparea tempCF
                end

            end   
    
        CF = [CF [allclusters.normCF]];
        
        from = repmat([allclusters.location]',length(allclusters),1);
        to = repelem([allclusters.location]',length(allclusters));
    
        for b = 1:length(from)
            loc(b,1) = strcat(from(b),'-',to(b));
        end
        
        CFlocation = [CFlocation; loc];
        
        clear c d include allclusters loc from to
        
    end
end

clear filelist filelist2 a b


%% Are there differences in the CF properties? Need to do PCAs for this

[coeff, score, ~,~,explained,~] = pca(CF','Centered',false);

PC1{1} = score(ismember([CFlocation],'cortex-cortex'),1);
PC1{2} = score(ismember([CFlocation],'cortex-tuber'),1);
PC1{3} = score(ismember([CFlocation],'cortex-tail'),1);
PC1{4} = score(ismember([CFlocation],'tuber-cortex'),1);
PC1{5} = score(ismember([CFlocation],'tuber-tuber'),1);
PC1{6} = score(ismember([CFlocation],'tuber-tail'),1);
PC1{7} = score(ismember([CFlocation],'tail-cortex'),1);
PC1{8} = score(ismember([CFlocation],'tail-tuber'),1);
PC1{9} = score(ismember([CFlocation],'tail-tail'),1);

PC2{1} = score(ismember([CFlocation],'cortex-cortex'),2);
PC2{2} = score(ismember([CFlocation],'cortex-tuber'),2);
PC2{3} = score(ismember([CFlocation],'cortex-tail'),2);
PC2{4} = score(ismember([CFlocation],'tuber-cortex'),2);
PC2{5} = score(ismember([CFlocation],'tuber-tuber'),2);
PC2{6} = score(ismember([CFlocation],'tuber-tail'),2);
PC2{7} = score(ismember([CFlocation],'tail-cortex'),2);
PC2{8} = score(ismember([CFlocation],'tail-tuber'),2);
PC2{9} = score(ismember([CFlocation],'tail-tail'),2);

PC3{1} = score(ismember([CFlocation],'cortex-cortex'),3);
PC3{2} = score(ismember([CFlocation],'cortex-tuber'),3);
PC3{3} = score(ismember([CFlocation],'cortex-tail'),3);
PC3{4} = score(ismember([CFlocation],'tuber-cortex'),3);
PC3{5} = score(ismember([CFlocation],'tuber-tuber'),3);
PC3{6} = score(ismember([CFlocation],'tuber-tail'),3);
PC3{7} = score(ismember([CFlocation],'tail-cortex'),3);
PC3{8} = score(ismember([CFlocation],'tail-tuber'),3);
PC3{9} = score(ismember([CFlocation],'tail-tail'),3);

PC4{1} = score(ismember([CFlocation],'cortex-cortex'),4);
PC4{2} = score(ismember([CFlocation],'cortex-tuber'),4);
PC4{3} = score(ismember([CFlocation],'cortex-tail'),4);
PC4{4} = score(ismember([CFlocation],'tuber-cortex'),4);
PC4{5} = score(ismember([CFlocation],'tuber-tuber'),4);
PC4{6} = score(ismember([CFlocation],'tuber-tail'),4);
PC4{7} = score(ismember([CFlocation],'tail-cortex'),4);
PC4{8} = score(ismember([CFlocation],'tail-tuber'),4);
PC4{9} = score(ismember([CFlocation],'tail-tail'),4);

% set up colour vector

cols = cbrewer2('seq','YlGnBu',9);

for a = 1:length(CF)
    if CFlocation(a) == string('cortex-cortex')      
        colvec(a,:) = cols(1,:);
    elseif CFlocation(a) == string('cortex-tuber')  
        colvec(a,:) = cols(2,:);
    elseif CFlocation(a) == string('cortex-tail')  
        colvec(a,:) = cols(3,:);
    elseif CFlocation(a) == string('tuber-cortex')  
        colvec(a,:) = cols(4,:);
    elseif CFlocation(a) == string('tuber-tuber')  
        colvec(a,:) = cols(5,:);
    elseif CFlocation(a) == string('tuber-tail')  
        colvec(a,:) = cols(6,:);
    elseif CFlocation(a) == string('tail-cortex')  
        colvec(a,:) = cols(7,:);
    elseif CFlocation(a) == string('tail-tuber')  
        colvec(a,:) = cols(8,:);
    elseif CFlocation(a) == string('tail-tail')  
        colvec(a,:) = cols(9,:);                      
    end
end

% Plot some stuff

figure

subplot(3,6,[1 2 7 8 13 14])
scatter3(score(:,1),score(:,3),score(:,2),60,colvec,'filled')
xlabel('PC1')
ylabel('PC3')
zlabel('PC2')
title('Plot of first 3 PCs by group')

subplot(3,6,[3 9])
violin(PC1,'medc',[],'facecolor',cols,'facealpha',0.6);
legend('off')
title('Scores of PC1')
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xticks(1:9)
xtickangle(90)
ylabel('PC1')

subplot(3,6,[4 10])
violin(PC2,'medc',[],'facecolor',cols,'facealpha',0.6);
legend('off')
title('Scores of PC2')
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xticks(1:9)
xtickangle(90)
ylabel('PC2')

subplot(3,6,[5 11])
violin(PC3,'medc',[],'facecolor',cols,'facealpha',0.6);
legend('off')
title('Scores of PC3')
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xticks(1:9)
xtickangle(90)
ylabel('PC3')

subplot(3,6,[6 12])
violin(PC4,'medc',[],'facecolor',cols,'facealpha',0.6);
legend('off')
title('Scores of PC4')
xticklabels({'co-co', 'co-tu', 'co-ta','tu-co', 'tu-tu', 'tu-ta','ta-co', 'ta-tu', 'ta-ta'});
xticks(1:9)
xtickangle(90)
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
saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/CFPCAuncentred.png');

%% Statistical testing of this

[p,tbl,stats] = anova1(score(:,1),CFlocation);
[p,tbl,stats] = anova1(score(:,2),CFlocation);
[p,tbl,stats] = anova1(score(:,3),CFlocation);
[p,tbl,stats] = anova1(score(:,4),CFlocation);

