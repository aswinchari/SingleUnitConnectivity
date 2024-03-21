%% For each file, takes the spike and LFP files and displays phase preferences in a number of freq domains
%==========================================================================

%% Define freq bands

eegbands(1).name = 'delta';
eegbands(2).name = 'theta';
eegbands(3).name = 'alpha';
eegbands(4).name = 'beta';
eegbands(5).name = 'gamma';
eegbands(6).name = 'highgamma';
eegbands(7).name = 'ripple';
eegbands(8).name = 'fastripple';

eegbands(1).freq = [1 4];
eegbands(2).freq = [4 8];
eegbands(3).freq = [8 12];
eegbands(4).freq = [12 20];
eegbands(5).freq = [20 45];
eegbands(6).freq = [65 80];
eegbands(7).freq = [80 200];
eegbands(8).freq = [200 500];

for a = 1:length(eegbands)
    eegbands(a).midpoint = mean(eegbands(a).freq);
end 

%% Load EEG and Spike files

filelist = dir('SpikesLFP/*_LFP.mat');

for a = 1:length(filelist)
    
    load(strcat(filelist(a).folder,'/',filelist(a).name));
    load(strcat(filelist(a).folder,'/',filelist(a).name(1:end-7),'spikes.mat'));
    
    % fourier transform EEG files
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.output      = 'pow';
    cfg.taper       = 'hanning';
    cfg.foi         = [0.5:1:200];                      % 1-200Hz
    cfg.channel     = 'all';
    % cfg.toi         = 0:1:150                           % 0-150s of the clip
    % cfg.t_ftimwin   = ones(length(cfg.foi),1).*2;       % 2s time windows 
    
    [freq] = ft_freqanalysis(cfg,eeg);
    
    % Polt some PSDs 
    
    figure
    hold on
    for b = 1:length(freq.label)
        plot(freq.freq,freq.powspctrm(b,:))
    end
    legend(freq.label)
    title('Power in the Microelectrode Channels')
    xlabel('Frequency (Hz)')
    ylabel('Power (uV)')
    ylim([0.000001 1000])
    set(gca,'YScale','log')
    
    % saveas(gcf,'Figures/PSD.png');
    
    % Plot some time-frequency plots: time-freq calculation
    
    cfg = [];
    cfg.method      = 'mtmconvol';
    cfg.output      = 'pow';
    cfg.taper       = 'hanning';
    cfg.foi         = [0.5:1:30];                      % 1-200Hz
    cfg.channel     = 'all';
    cfg.toi         = 0:1:150                           % 0-150s of the clip
    cfg.t_ftimwin   = ones(length(cfg.foi),1).*2;       % 2s time windows 
    
    [freq] = ft_freqanalysis(cfg,eeg);
    
    % prepare a layout
    
    cfg             = [];
    cfg.columns     = 3;
    cfg.layout      = 'ordered';
    cfg.direction   = 'TBLR';
    
    layout = ft_prepare_layout(cfg,eeg);
    
    
    cfg = [];
    cfg.parameter   = 'powspctrm';
    cfg.maskstyle   = 'saturation';
    cfg.layout      = layout;
    cfg.channel     = 'all';
    cfg.showlabels  = 'yes';
    cfg.zlim        = [0 100];
    
    ft_multiplotTFR(cfg,freq)
    
    % saveas(gcf,'Figures/TFPlot.png');

    % Try to make ISI histograms & waveshape for each of the clusters
    
    fs = 30000;                 % sampling rate
    tw = [-19:44] * (1/fs);     % time line for waveshapes (defined in waveclus parameters)

    for b = 1:length(allclusters)
        allclusters(b).ISI = diff(allclusters(b).spiketimes);
        allclusters(b).ts = 0.5:max(allclusters(b).ISI);
        allclusters(b).scount =  hist(allclusters(b).ISI,  allclusters(b).ts);
        allclusters(b).aw = linspace(min(allclusters(b).spikes(:)),max(allclusters(b).spikes(:)), 100); 
        allclusters(b).Tw = repmat(tw,[size(allclusters(b).spikes,1),1]);
        allclusters(b).WSH = hist3([allclusters(b).Tw(:), allclusters(b).spikes(:)],'Ctrs',{tw allclusters(b).aw});
        allclusters(b).refrac = sum(allclusters(b).scount(1:3))/length(allclusters(b).spikes);
    
        % ft_spike_select them
        
        spike = ft_read_spike('WaveClus_NNAIS04_1/times_NNAIS04_1_channel_26.mat')
        
        cfg = [];
        cfg.spike
        
    end
    
    
    
    
end



