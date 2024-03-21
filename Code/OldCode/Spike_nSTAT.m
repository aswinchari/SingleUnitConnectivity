%% For each file, takes the spikes files
%==========================================================================

%% Load spikes

filelist = dir('SpikesLFP/*_spikes.mat');

fs = 30000;                 % sampling rate
tw = [-19:44] * (1/fs);     % time line for waveshapes (defined in waveclus parameters)

for a = 1:length(filelist)

    load(strcat(filelist(a).folder,'/',filelist(a).name));

    % define spike train properties

    for b=1:length(allclusters)

        allclusters(b).ISI = diff(allclusters(b).spiketimes);
        allclusters(b).ts = 0.5:max(allclusters(b).ISI);
        allclusters(b).scount =  hist(allclusters(b).ISI,  allclusters(b).ts);
        allclusters(b).aw = linspace(min(allclusters(b).spikes(:)),max(allclusters(b).spikes(:)), 100); 
        allclusters(b).Tw = repmat(tw,[size(allclusters(b).spikes,1),1]);
        allclusters(b).WSH = hist3([allclusters(b).Tw(:), allclusters(b).spikes(:)],'Ctrs',{tw allclusters(b).aw});
        allclusters(b).refrac = sum(allclusters(b).scount(1:3))/length(allclusters(b).spikes);
        
        subplot(length(allclusters),1,b)
        bar(allclusters(b).ts,allclusters(b).scount)
        title(strcat('Refractory =',num2str(round(allclusters(b).refrac*100,1)),'%'))
        
    end
   
    % create a raster
    
        
end