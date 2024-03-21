%% Open spikes files and run the post-spike and coupling filters for each one
%==========================================================================

cd('/home/aswinchari/Documents/NNAIS_0506/NNAIS07/Spikes');
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
    
%% Load all files
 
day1 = load(strcat(filelist(1).folder,'/',filelist(1).name));
day2 = load(strcat(filelist(2).folder,'/',filelist(2).name));
day3 = load(strcat(filelist(3).folder,'/',filelist(3).name)); 

%% Plot clusters by day

subplot(3,3,[1 2 3 4 5 6])
scatter3([day1.allclusters.freq],[day1.allclusters.peakISI],[day1.allclusters.spikeHW],20,'filled');
hold on
scatter3([day2.allclusters.freq],[day2.allclusters.peakISI],[day2.allclusters.spikeHW],20,'filled');
scatter3([day3.allclusters.freq],[day3.allclusters.peakISI],[day3.allclusters.spikeHW],20,'filled');
 xlabel('Frequency (Hz)');
 ylabel('Peak ISI (ms)');
 zlabel('Spike Half Width (ms)');
 set(gca,'xscale','log')
 set(gca,'yscale','log')
 title('3D plot')
 legend({'Day 1', 'Day 2', 'Day 3'})
 
subplot(3,3,7)
scatter([day1.allclusters.freq],[day1.allclusters.peakISI],20,'filled');
hold on
scatter([day2.allclusters.freq],[day2.allclusters.peakISI],20,'filled');
scatter([day3.allclusters.freq],[day3.allclusters.peakISI],20,'filled');
 xlabel('Frequency (Hz)');
 ylabel('Peak ISI (ms)');
 set(gca,'xscale','log')
 set(gca,'yscale','log')
 title('Frequency vs Peak ISI')
 
subplot(3,3,8)
scatter([day1.allclusters.freq],[day1.allclusters.spikeHW],20,'filled');
hold on
scatter([day2.allclusters.freq],[day2.allclusters.spikeHW],20,'filled');
scatter([day3.allclusters.freq],[day3.allclusters.spikeHW],20,'filled');
 xlabel('Frequency (Hz)');
 ylabel('Spike Half Width (ms)');
 set(gca,'xscale','log')
 title('Frequency vs Spike Half Width')
 
subplot(3,3,9)
scatter([day1.allclusters.peakISI],[day1.allclusters.spikeHW],20,'filled');
hold on
scatter([day2.allclusters.peakISI],[day2.allclusters.spikeHW],20,'filled');
scatter([day3.allclusters.peakISI],[day3.allclusters.spikeHW],20,'filled');
 xlabel('PeakISI (ms)');
 ylabel('Spike Half Width (ms)');
 set(gca,'xscale','log')
 title('Peak ISI vs Spike Half Width')
 
 set(gcf,'position',[20 20 1000 800]);
 saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/ClusterClusters05.png');
 
%% Plot clusters based on location

allclusters = [day1.allclusters]; %day2.allclusters day3.allclusters];
allfreq = [allclusters.freq];
allpeakISI = [allclusters.peakISI];
allspikeHW = [allclusters.spikeHW];

subplot(3,3,[1 2 3 4 5 6])
scatter3(allfreq(ismember([allclusters.location],'cortex')),allpeakISI(ismember([allclusters.location],'cortex')),allspikeHW(ismember([allclusters.location],'cortex')),20,'filled');
hold on
scatter3(allfreq(ismember([allclusters.location],'tuber')),allpeakISI(ismember([allclusters.location],'tuber')),allspikeHW(ismember([allclusters.location],'tuber')),20,'filled');
scatter3(allfreq(ismember([allclusters.location],'tail')),allpeakISI(ismember([allclusters.location],'tail')),allspikeHW(ismember([allclusters.location],'tail')),20,'filled');
 xlabel('Frequency (Hz)');
 ylabel('Peak ISI (ms)');
 zlabel('Spike Half Width (ms)');
 set(gca,'xscale','log')
 set(gca,'yscale','log')
 title('3D plot')
 legend({'Cortex', 'Tuber', 'Tail'})
 
subplot(3,3,7)
scatter(allfreq(ismember([allclusters.location],'cortex')),allpeakISI(ismember([allclusters.location],'cortex')),20,'filled');
hold on
scatter(allfreq(ismember([allclusters.location],'tuber')),allpeakISI(ismember([allclusters.location],'tuber')),20,'filled');
scatter(allfreq(ismember([allclusters.location],'tail')),allpeakISI(ismember([allclusters.location],'tail')),20,'filled');
 xlabel('Frequency (Hz)');
 ylabel('Peak ISI (ms)');
 set(gca,'xscale','log')
 set(gca,'yscale','log')
 title('Frequency vs Peak ISI')
 
subplot(3,3,8)
scatter(allfreq(ismember([allclusters.location],'cortex')),allspikeHW(ismember([allclusters.location],'cortex')),20,'filled');
hold on
scatter(allfreq(ismember([allclusters.location],'tuber')),allspikeHW(ismember([allclusters.location],'tuber')),20,'filled');
scatter(allfreq(ismember([allclusters.location],'tail')),allspikeHW(ismember([allclusters.location],'tail')),20,'filled');
 xlabel('Frequency (Hz)');
 ylabel('Spike Half Width (ms)');
 set(gca,'xscale','log')
 title('Frequency vs Spike Half Width')
 
subplot(3,3,9)
scatter(allpeakISI(ismember([allclusters.location],'cortex')),allspikeHW(ismember([allclusters.location],'cortex')),20,'filled');
hold on
scatter(allpeakISI(ismember([allclusters.location],'tuber')),allspikeHW(ismember([allclusters.location],'tuber')),20,'filled');
scatter(allpeakISI(ismember([allclusters.location],'tail')),allspikeHW(ismember([allclusters.location],'tail')),20,'filled');

 xlabel('PeakISI (ms)');
 ylabel('Spike Half Width (ms)');
 set(gca,'xscale','log')
 title('Peak ISI vs Spike Half Width')
 
  set(gcf,'position',[20 20 1000 800]);
  saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/ClusterClusters05_bylocation.png');
  
  %% Plot filters by location
  
  allPSF = [allclusters.PSF];
  
  subplot(3,1,1)
  plot(allPSF(:,ismember([allclusters.location],'cortex')'));
  title('Cortical PSFs')
  ylim([-0.5 1.5])
  xlabel('time (ms)')
  ylabel('gain')
  
  subplot(3,1,2)
  plot(allPSF(:,ismember([allclusters.location],'tuber')'));
  title('Tuber PSFs')
  ylim([-0.5 1.5])
  xlabel('time (ms)')
  ylabel('gain')
  
  subplot(3,1,3)
  plot(allPSF(:,ismember([allclusters.location],'tail')'));
  title('Tail PSFs')
  ylim([-0.5 1.5]) 
  xlabel('time (ms)')
  ylabel('gain')
  
  set(gcf,'position',[20 20 1000 800]);
  saveas(gcf,'/home/aswinchari/Documents/NNAIS_0506/PrelimFigs/PSFs_bylocation.png');