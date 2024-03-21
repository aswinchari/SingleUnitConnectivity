
%% View all units from single times_*.mat file

CurFile = uigetfile();
    
load(CurFile);

% find the cluster labels and remove 0 (noise cluster)
cluslabels = setdiff(unique(cluster_class(:,1)),0);

% make plot for every cluster

color = [0.6,0.6,0.6]; % gray

% prepare the waveshape plots:
fs = 30000; % sampling rate
tw = [-20:43] * (1/fs); % time line for waveshapes (defined in waveclus parameters)
   
    for c = 1:length(cluslabels)

        % collect spike times for current cluster
        ids = find(cluster_class(:,1) == c);

        % get spike times in ms
        Tspikes = cluster_class(ids,2);

        % compute ISIs
        % Tspikes contains the spikes ordered in time: 
        % ISI is the time difference between a spike and the next spike
        ISIs = diff(Tspikes);

        % prepare the ISI histogram
        ts = min(ISIs):max(ISIs); % prepare a time line with 1 ms timestep
        scount = hist(ISIs, ts);

        % compute waveshape density
        Waveshapes = spikes(ids,:);
        Tw = repmat(tw,[length(ids),1]);
        aw = linspace(min(Waveshapes(:)),max(Waveshapes(:)), 100); % break the waveshape amplitudes into bins
        wcount = hist3([Tw(:), Waveshapes(:)],'Ctrs',{tw aw});

        % plot
        figure; 
        subplot(2,3,1);
        imagesc(tw*1e3,aw, wcount');
        axis xy
        xlabel('Time (ms)')
        ylabel('Amplitude (\muV)')
        title('Waveshape histogram')

        subplot(2,3,2:3);
        bar(ts,scount, 'facecolor', color)
        xlabel('ISI (ms)')
        title('Inter-spike interval')
        text(max(ts)*0.6,max(scount)*0.8,['N_{spikes} = ',num2str(length(Tspikes))])

        subplot(2,3,4:6)
        scatter(Tspikes/1e3, rand(size(Tspikes)), 5,color,'filled')
        xlim([0,1800])
        xlabel('Time (s)');
        title('Rastergram')
        % seif ~isempty(cluslabels)t(gca,'ytick',[])

        % save plot here
        % saveas(gcf,join([(CurFile(1:end-9)),'unit',string(c),'.jpg']))

    end
