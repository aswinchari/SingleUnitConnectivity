% %% Get meta-data for batch processing
% xls_file = '~/Desktop/Listallcells - Copy.xls';
% [~, ~, xls_data] = xlsread(xls_file);
% 
% %% Tetrode files
% % NTT files
% filepath = xls_data(2:end, strcmp(xls_data(1, :), 'Path'));
% ntt_file = xls_data(2:end, strcmp(xls_data(1, :), 'NTT'));
% 
% data_files = cellfun(@(x, y) [x, '\', y, '.ntt'], filepath, ntt_file,...
%     'UniformOutput', false);
% 
% % Clusters
% clust = cell2mat(xls_data(2:end, strcmp(xls_data(1, :), ' Cluster')));
% 
% %%

%load ~/Desktop/test_model.mat

fm = fitted_models;
dt = fm.dt;
rt = fm.rt;
samp_ts = fm.samp_ts;

event_data_file = '~/Desktop/newevent.csv';
event_data = dlmread(event_data_file);
event_data(:, 1) = event_data(:, 1) * 1e-6;

% 4: Nose Poke
ev = bin_raster({event_data(event_data(:, 2) == 4, 1)'}, dt, 0,...
    minmax(samp_ts));

%% Get event response basis

evprs.ncols = 5;
evprs.hpeaks = [0.05, 2];
evprs.b = 0.5;
[evt, evbasis] = makeBasis_PostSpike(evprs, dt);

%% Plot peri-nosepoke time histogram

max_lag = ceil(max(evt) * 1e3); % in milliseconds
%max_lag = 5e3; % in milliseconds
sm_width = 1e2; % in milliseconds

% Peri-stimulus time histogram
[psth, lags] = xcorr(rt, ev, max_lag);
psth = gaussian_smooth(psth, sm_width);

% Presses-only model prediction
intns = fm.poiss_intensity * dt;
psth_poiss = xcorr(intns, ev, max_lag);

% Full model prediction
intns = fm.full_intensity * dt;
psth_full = xcorr(intns, ev, max_lag)';
psth_full = gaussian_smooth(psth_full, sm_width);

% Plot
figure;
hold on;
h = zeros(3, 1);
h(1) = plot(lags * dt, psth, 'b', 'linewidth', 2);
h(2) = plot(lags * dt, psth_poiss, 'm', 'linewidth', 2);
h(3) = plot(lags * dt, psth_full, 'g', 'linewidth', 2);

title('Peri-nosepoke histogram', 'fontsize', 16)
xlabel('Time (s)', 'fontsize', 16)
ylabel('Smoothed counts', 'fontsize', 16)
legend(h, {'Smoothed PSTH', 'Poisson', 'GLM'})

set(gca, 'linewidth', 2, 'fontsize', 16,...
    'ylim', [0, max(max(psth), max(psth_full))])

%%

np_coefs = fm.full.x(strcmp(fm.basis_labels, 'Nosepoke'));
np_coefs_for = np_coefs(1:evprs.ncols);
np_coefs_back = np_coefs((evprs.ncols + 1):length(np_coefs));

np_filter = [flipud(evbasis * np_coefs_back); evbasis * np_coefs_for];
npt = [-flipud(evt); evt];

figure;
plot(npt, exp(np_filter))

%%

[cp, lags] = xcorr(fm.correct_press, ev, max_lag);
cp_times = lags(cp == 1) * dt;

[ip, lags] = xcorr(fm.incorrect_press, ev, max_lag);
ip_times = lags(ip == 1) * dt;

[lp, lags] = xcorr(fm.left_press, ev, max_lag);
lp_times = lags(lp == 1) * dt;

[rp, lags] = xcorr(fm.right_press, ev, max_lag);
rp_times = lags(rp == 1) * dt;

%%

% Plot
figure;
hold on;
h = zeros(3, 1);
h(1) = plot(lags * dt, psth, 'b', 'linewidth', 2);
h(2) = plot(lags * dt, psth_poiss, 'm', 'linewidth', 2);
h(3) = plot(lags * dt, psth_full, 'g', 'linewidth', 2);
h(4) = plot(npt, psth_poiss(1) * exp(np_filter), 'c', 'linewidth', 2);
h(5) = plot(cp_times, zeros(size(cp_times)), 'g^', 'linewidth', 2);
h(6) = plot(ip_times, zeros(size(ip_times)), 'r^', 'linewidth', 2);
h(7) = plot(lp_times, zeros(size(lp_times)), 'c.', 'linewidth', 2);
h(8) = plot(rp_times, zeros(size(rp_times)), 'm.', 'linewidth', 2);

title('Peri-nosepoke histogram', 'fontsize', 16)
xlabel('Time (s)', 'fontsize', 16)
ylabel('Smoothed counts', 'fontsize', 16)
legend(h, {'Smoothed PSTH', 'Poisson', 'GLM', 'NP filter',...
    'Correct press', 'Incorrect press', 'Left press', 'Right press'})

set(gca, 'linewidth', 2, 'fontsize', 16,...
    'ylim', [0, max(max(psth), max(psth_full))])








