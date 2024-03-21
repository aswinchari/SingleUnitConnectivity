function [data, min_time, max_time] = bin_raster(raster, time_scale,...
    sparse_tf, time_bnds)
% [binned_raster, min_time, max_time] = bin_raster(raster, time_scale, sparse_tf, time_bnds)
% Bin raster data over given time scale.

%% Parse input
if nargin < 3
    sparse_tf = 0;
end

%% Get range of time stamps
if nargin < 4
    max_time = 0;
    min_time = 1e10;
    for i = 1:length(raster)
        if ~isempty(raster{i})
            max_time = max(max_time, max(raster{i}));
            min_time = min(min_time, min(raster{i}));
        end
    end
else
    min_time = time_bnds(1);
    max_time = time_bnds(2);
end

%% Bin over time scale

bins = (min_time - time_scale):time_scale:(max_time + time_scale);

if sparse_tf
    data = sparse(length(raster), length(bins));
else
    data = zeros(length(raster), length(bins));
end

for i = 1:length(raster)
    if ~isempty(raster{i})
        data(i,:) = histc(raster{i}, bins);
    end
end