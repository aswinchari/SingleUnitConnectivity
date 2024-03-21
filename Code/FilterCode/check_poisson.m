function [p, firing_rate] = check_poisson(raster, dt, TS)
% p = check_poisson(raster, dt, TS)
% Compute the p-value of the time-rescaled ISI probablities under the null
% (Poisson) model.

% Bin raster
rt = bin_raster(raster, dt, false, minmax(TS));

% Compute firing rate
num_spikes = sum(rt);
firing_rate = num_spikes / diff(minmax(TS));

% Intensity of Poisson
intensity = dt * firing_rate * ones(size(rt))';

% Kolmogorov-Smirnov stats for poisson models
[fig, z] = time_rescale_qqplot(intensity, rt,...
    {'Poisson'});
close(fig);

[~, p] = kstest(norminv(z));