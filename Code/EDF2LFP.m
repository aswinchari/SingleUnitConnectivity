%% Open EDF File and save LFP into separate folder using FT
%==========================================================================

cd('/home/aswinchari/Documents/NNAIS_0506/NNAIS06')
genpath('/home/aswinchari/Documents/GitHub/wave_clus-master');
addpath(ans);
clear ans;
savepath

%% Create file list and open them

filelist = dir('Data/');
load('channels.mat');
mkdir('LFP')
%%
for a = 3:length(filelist)

    disp(strcat('Current File:',filelist(a).name));

    % Open one channel at a time

    for b = 1:length(channels)
    
    	disp(strcat('   Channel:',string(b)))
        
        % First half
    
        cfg                         = [];
        cfg.dataset                 = strcat(filelist(a).folder,'/',filelist(a).name);
        cfg.channel                 = channels(b);
        cfg.continuous              = 'yes';
        %cfg.dftfilter               = 'yes';
        cfg.trl                     = [1 12000000 0 ]; %54000000 is 30 min, change based on file

        eeg                        = ft_preprocessing(cfg);
        
        % isolate only that channel (needed for ns5 data only)
        
        eeg.trial{1,1} = eeg.trial{1,1}(channels(b),:);
        
        % LP filter
       
        cfg                         = [];  
        cfg.lpfilter                = 'yes';
        cfg.dftfilter               = 'yes';
        cfg.lpfreq                  = 1000;
        
        eeg = ft_preprocessing(cfg,eeg);
        
        first = eeg.trial{1,1};
        
        clear eeg
        
        % Second half
    
        cfg                         = [];
        cfg.dataset                 = strcat(filelist(a).folder,'/',filelist(a).name);
        cfg.channel                 = channels(b);
        cfg.continuous              = 'yes';
        %cfg.dftfilter               = 'yes';
        cfg.trl                     = [12000001 24000000 0 ]; %54000000 is 30 min, change based on file

        eeg                        = ft_preprocessing(cfg);
        
        % isolate only that channel (needed for ns5 data only)
        
        eeg.trial{1,1} = eeg.trial{1,1}(channels(b),:);
        
        % LP filter
       
        cfg                         = [];  
        cfg.lpfilter                = 'yes';
        cfg.lpfreq                  = 1000;
        cfg.dftfilter               = 'yes';
        
        eeg = ft_preprocessing(cfg,eeg);
        
        second = eeg.trial{1,1};
        
        clear eeg
        
         % third half
    
        cfg                         = [];
        cfg.dataset                 = strcat(filelist(a).folder,'/',filelist(a).name);
        cfg.channel                 = channels(b);
        cfg.continuous              = 'yes';
        %cfg.dftfilter               = 'yes';
        cfg.trl                     = [24000001 36000000 0 ]; %54000000 is 30 min, change based on file

        eeg                        = ft_preprocessing(cfg);
        
        % isolate only that channel (needed for ns5 data only)
        
        eeg.trial{1,1} = eeg.trial{1,1}(channels(b),:);
        
        % LP filter
       
        cfg                         = [];  
        cfg.lpfilter                = 'yes';
        cfg.lpfreq                  = 1000;
        cfg.dftfilter               = 'yes';
        
        eeg = ft_preprocessing(cfg,eeg);
        
        third = eeg.trial{1,1};
        
        clear eeg
        
        % combine and save
        
        data = [first second third];
        
        save(strcat('LFP/',filelist(a).name(1:end-4),'_channel_',string(channels(b)),'_LFP','.mat'),'data');
        
        clear data first second
        
    end
end

