%% Open EDF File and save spikes into separate folder using FT
%==========================================================================

cd('/home/aswinchari/Documents/NNAIS_0506/NNAIS07')
genpath('/home/aswinchari/Documents/GitHub/wave_clus-master');
addpath(ans);
clear ans;
savepath

%% Create file list and open them

filelist = dir('Data/*.edf');
load('channels.mat')

%% For each file, read EDF and save micro channels as mat file and run through spike detection

for a = 4:length(filelist)

    disp(strcat('Current File:',filelist(a).name));
    mkdir(strcat('Data/WaveClus_',filelist(a).name(1:end-4)));

    % Open one channel at a time

    for b = 1:length(channels{a})
    
    	disp(strcat('   Channel:',string(b)))
    
        cfg                         = [];
        cfg.dataset                 = strcat(filelist(a).folder,'/',filelist(a).name);
        cfg.channel                 = channels{a}(b)
        cfg.continuous              = 'yes';
        %cfg.dftfilter               = 'yes';
        cfg.trl                     = [1 54000000 0 ] %54000000 is 30 min

        eeg                        = ft_preprocessing(cfg);
        
        % isolate only that channel (needed for ns5 data only)
        % eeg.trial{1,1} = eeg.trial{1,1}(channels(b),:)
        
        % add in the dftfilter here
        
        cfg                         = [];
        cfg.dftfilter               = 'yes'; 
        
        eeg = ft_preprocessing(cfg,eeg);
        
        % prepare and save files for waveclus 
        
        data = eeg.trial{1,1}(1,:);
        sr = eeg.fsample;
        save(strcat('Data/WaveClus_',filelist(a).name(1:end-4),'/',filelist(a).name(1:end-4),'_channel_',string(channels{a}(b)),'.mat'),'data','sr');
        clear data eeg cfg;
    end 

% Prepare for WaveClus batch processing

cd(strcat('Data/WaveClus_',filelist(a).name(1:end-4)));

% Detect Spikes

disp('Starting Spike Detection')

filelist2 = dir('*.mat');

    for b = 1:length(filelist2)
        tosort{b,1} = strcat(filelist2(b).folder,'/',filelist2(b).name);
    end

writecell(tosort,'tosort.txt');
clear tosort

param.stdmin = 5;
param.sr = 30000;
param.detection = 'both';       % can be 'pos', 'neg' or 'both'
param.min_clus = 500;
param.segments_length = 30;

Get_spikes('tosort.txt','parallel',true,'par',param);


% Detect Clusters

disp('Starting Clustering')

filelist3 = dir('*spikes.mat');

    for c = 1:length(filelist3)
        tosort{c,1} = strcat(filelist3(c).folder,'/',filelist3(c).name);
    end

writecell(tosort,'tosort.txt');
clear tosort

Do_clustering('tosort.txt','parallel',true,'par',param,'make_plots',false);



cd ../..

end

clear all

disp('Over to you to refine the clusters!')

%% Manually load GUI and save clusters for each channel in each folder

% cd into the foldeer

wave_clus

% at this stage, can visualise unit using ViewSingleUnit.m or plot them all
% using PlotAllUnits.m

% cd back out of the folder

%% For each spike within each file, load it into a structure


filelist3 = dir('Data/WaveClus*')

for c = 1:length(filelist3)

    filelist4 = dir(strcat(filelist3(c).folder,'/',filelist3(c).name,'/times_*.mat'));
    
    allclusters = struct([]);
    
    for d = 1:length(filelist4)
        
        load(strcat(filelist4(d).folder,'/',filelist4(d).name));
        if max(cluster_class(:,1))~=0;
            for e = 1:max(cluster_class(:,1));
                cluster(e).file = filelist3(c).name;
                cluster(e).channel = filelist4(d).name(1:end-4);
                cluster(e).cluster = e;
                cluster(e).spikes = spikes(cluster_class(:,1)==e,:);
                cluster(e).spiketimes = cluster_class(cluster_class(:,1)==e,2); 
            end
            allclusters = [allclusters,cluster];
            clear cluster
        end
    end

mkdir('Spikes')
save(strcat('Spikes/',filelist3(c).name(10:end),'_spikes.mat'),'allclusters');

end