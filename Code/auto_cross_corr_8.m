%% Create some cross and auto correlation plots
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

PSF = PSF([PSF.refrac] < 0.02');

clear a b allclusters filelist filelist2

%% Plot some Autocorrs

for a = 1:10 %length(PSF)
    [counts, centres] = hist(PSF(a).spiketimes,[0.5:1:round(max(PSF(a).spiketimes))+1.5]);
    PSF(a).bincounts = counts
    [c, lags] = xcorr
    clear counts centres
end


