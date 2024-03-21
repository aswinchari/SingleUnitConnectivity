function [fig, z] = time_rescale_qqplot(intensity, rt, leg_names, fig_title)
% [fig, z] = time_rescale_qqplot(intensity, rt,  leg_names, fig_title))
% Make time rescaling QQ-plots for point process models

%%

if nargin < 3
    leg_names = {};
end

if nargin < 4
    fig_title = '';
end

%%

aug_rt = [1, rt];

n = sum(rt  > 0);
k = 1:n;
bk = (k - 0.5) / n;

conf = 1.63  /sqrt(n);

fig = figure;
hold on;

plot(bk, bk, 'k', 'linewidth', 2)
plot(bk, bk + conf, 'k--', 'linewidth', 2)
plot(bk, bk - conf, 'k--', 'linewidth', 2)
set(gca, 'xlim', [0, 1], 'ylim', [0, 1], 'linewidth', 2, 'fontsize', 16)

% Compute emperical quantiles
z = [];
for i = 1:size(intensity, 2)
    cum_intensity = cumsum(intensity(:, i));
    lam = cum_intensity(aug_rt > 0);
    z = [z, sort(1 - exp(-diff(lam)))];
end

% Plot lines
h = plot(bk, z, 'linewidth', 2);

legend(h, leg_names, 'Location', 'NorthWest', 'linewidth', 2)
xlabel('Theoretical quantiles')
ylabel('Model quantiles')
title(fig_title);

