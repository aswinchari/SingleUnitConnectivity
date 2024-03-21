%% For each file, takes the LFP files and calculates power and entropy in different frequencies over time
%==========================================================================

%% Load EEG

filelist = dir('SpikesLFP/*_LFP.mat');

data = [];
freqs = [1:1:30];

for a = 1:length(filelist)
    
    data(a).filename = filelist(a).name(1:9);
    
    load(strcat(filelist(a).folder,'/',filelist(a).name));
     
    % Calculate power spectrum
            
        cfg = [];
        cfg.method      = 'mtmfft';
        cfg.output      = 'pow';
        cfg.taper       = 'hanning';
        cfg.foi         = freqs;
        cfg.channel     = 'all';

        [freq] = ft_freqanalysis(cfg,eeg);
       
        data(a).powspctrm = freq.powspctrm;
        
    end
    
    % FOOOF analysis
    
    settings = struct();
    f_range =  [2,30];
    
    fooof_results = fooof(freqs,data(a).powspctrm(1,:),f_range,settings);