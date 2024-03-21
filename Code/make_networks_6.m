%% Make networks from CFs and PSFs
%==========================================================================

clc
clear all
cd('/home/aswinchari/Documents/NNAIS_0506/');
addpath('/home/aswinchari/Documents/NNAIS_0506/Code/FilterCode');
addpath('/home/aswinchari/Documents/GitHub/edge-centric_demo-main/fcn/');

%% Copy all spike files to Networks folder and go into that

filelist = dir('NNAIS*'); 

mkdir('Networks');

for a = 1:length(filelist)
    
    filelist2 = dir(strcat(filelist(a).name,'/Spikes/*filters.mat'));
    
    for b = 1:length(filelist2)
        
       copyfile(strcat(filelist2(b).folder,'/',filelist2(b).name),strcat('Networks/',filelist(a).name,filelist2(b).name))
        
    end
end

%% Create PSF and CF networks

clc
clear all
cd('/home/aswinchari/Documents/NNAIS_0506/Networks');

filelist = dir('NNAIS*');

for a = 1:length(filelist)
    
    % load file
    
    load(filelist(a).name);
    
    include = [allclusters.refrac] < 0.02';
    allclusters = allclusters(include);
    
    for b = 1:length(allclusters)
        
        % normalise PSFs
        
        temparea = sum(abs(cumtrapz(allclusters(b).PSF)));
        tempPSF = allclusters(b).PSF;
        allclusters(b).normPSF = tempPSF./temparea;
        clear tempPSF temparea
    
        % normalise CFs
        
        allclusters(b).CF = allclusters(b).CF(:,include);
        
        for c = 1:size(allclusters(1).CF,2)     
            temparea = sum(abs(cumtrapz(allclusters(b).CF(:,c))));
            tempCF = allclusters(b).CF(:,c);
            allclusters(b).normCF(:,c) = tempCF./temparea;
            clear temparea tempCF
        end
        
    end
    
    clear b c include
    
    % make PSF and CF networks
    
    allPSF = [allclusters.normPSF];
    networks(a).PSFnetwork = fcn_edgets2edgecorr(allPSF);
    allCF = [allclusters.normCF];
    networks(a).CFnetwork = fcn_edgets2edgecorr(allCF);
    
    % locations for each network
    
    networks(a).PSFlocation = [allclusters.location]';
    
    from = repmat(networks(a).PSFlocation,length(allclusters),1);
    to = repelem(networks(a).PSFlocation,length(allclusters));
    
    for b = 1:length(from)
        networks(a).CFlocation(b,1) = strcat(from(b),'-',to(b));
    end

    clear allPSF allCF b from to allclusters
    
end

clear a filelist
    
save('networks.mat','networks');
    