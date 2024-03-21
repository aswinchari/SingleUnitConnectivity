%% Open spikes files and run the post-spike and coupling filters for each one
%==========================================================================

clc
clear all
cd('/home/aswinchari/Documents/NNAIS_0506/NNAIS05/Spikes');
addpath('/home/aswinchari/Documents/NNAIS_0506/Code/FilterCode');

%% Create file list and open them

filelist = dir('*spikes.mat');

%% For each session of spikes, calculate PSF and CFs for each neuron

% fluff that's needed for other things
dt = 1e-3;
ihprs.ncols = 10;
ihprs.hpeaks = [.01, 0.5];  
ihprs.b = 0.5;
[iht, ihbasis] = makeBasis_PostSpike(ihprs, dt);
clear dt ihprs

for a = 1:length(filelist)
    
    disp(strcat('File',string(a)));
    
    % Load the file
    
    load(strcat(filelist(a).folder,'/',filelist(a).name));
    
    % Exclude some neurons because of high refractory period
    % Defined as <2ms interval, should be <2%
    
    for c = 1:length(allclusters)
        allclusters(c).ISI = diff(allclusters(c).spiketimes);
        allclusters(c).refrac = sum(allclusters(c).ISI<2)/length(allclusters(c).ISI);  
    end
    
    % allclusters = allclusters([allclusters.refrac]<0.02);
    
    % Create the cell array
    
    for c = 1:length(allclusters)
        neurons{c} = allclusters(c).spiketimes./1000;
    end
    
    for c = 1:length(allclusters)
        
        disp(strcat('Cluster',string(c)));
        
        % Create the PSF
        
        PSFoutput = glm_single_neuron_fit(neurons(c), [0 1800]);
        allclusters(c).PSF = ihbasis * PSFoutput.full.x(PSFoutput.basis_labels=="Postspike");
        clear PSFoutput

        % Create the CF 
        
        CFoutput = glm_ensemble_fit(neurons, [0 1800],0,c);
        allclusters(c).CF(:,c) = ihbasis * CFoutput.x(CFoutput.basis_labels=="PSF");
        for d = 1:length(allclusters)
            if d~=c
            allclusters(c).CF(:,d) = ihbasis * CFoutput.x(CFoutput.basis_labels==strcat("CPF",string(d)));   
            end
        end
        clear CFoutput
            
    end
    
    save(strcat(filelist(a).folder,'/',filelist(a).name,'_filters.mat'),'allclusters');
    clear allclusters
end
