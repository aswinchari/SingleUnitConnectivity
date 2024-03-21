%% Open filters files and run the spike properties for each one
% Need to run separately for each patient
%==========================================================================

clc
clear all
cd('/home/aswinchari/Documents/NNAIS_0506/NNAIS05/Spikes');
addpath('/home/aswinchari/Documents/NNAIS_0506/Code/FilterCode');
addpath('/home/aswinchari/Documents/GitHub/intersections');


%% Create file list and open them

filelist = dir('*filters.mat');
filelist2 = dir('*location.mat');

%% Open filters and add locations to them

for a = 1:length(filelist)
    
    % Load the file
    
    load(strcat(filelist(a).folder,'/',filelist(a).name));
    
    % add locations
    
    load(strcat(filelist2(a).folder,'/',filelist2(a).name));
    
    for b = 1:length(allclusters)
        allclusters(b).location = location(b);
    end
    
    % add some other data eg Hz, spike HW and peak ISI
    
    for b = 1:length(allclusters)
        allclusters(b).freq     = length(allclusters(b).spiketimes)/1800;        % can change depending on recording length
        allclusters(b).peakISI  = max(allclusters(b).ISI);
        if max(mean(allclusters(b).spikes)) > abs(min(mean(allclusters(b).spikes)));
        allclusters(b).spikedir = 'pos';
        else
        allclusters(b).spikedir = 'neg';  
        end
        
        if allclusters(b).spikedir == 'pos'
        [xint, yint] = intersections([-20:43],mean(allclusters(b).spikes),[-20:43],repelem(max(mean(allclusters(b).spikes)/2),64));
        allclusters(b).spikeHW  = xint(2) - xint(1);
        clear xint yint
        else
        [xint, yint] = intersections([-20:43],mean(allclusters(b).spikes),[-20:43],repelem(min(mean(allclusters(b).spikes)/2),64));
        allclusters(b).spikeHW  = xint(2) - xint(1);
        clear xint yint  
        end
    end
    
    save(filelist(a).name,'allclusters')
    
    clear allclusters location
    
end
    
