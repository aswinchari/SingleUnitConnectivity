function sm_data = gaussian_smooth(data, sig)
% sm_data = gaussian_smooth(data, sig)
% Smooth the counts in DATA. SIG is the input argument for 'gwin.m'. It
% corresponds to the number of points that constitute one standard deviation
% of the window.

% function sm_data = gaussian_smooth(data,width,alpha)
% % sm_data = gaussian_smooth(data,window,alpha)
% % Smooth the counts in DATA. WIDTH is the number of points in the smoothing
% % window. ALPHA (optional) is the same argument for 'gausswin.m'.
% 
% % Parse inputs
% if nargin < 3
%     alpha = 2.5;
% end
% 
% % % Construct smoothing window
% half_width = ceil(width/2);
% gauss_filter = gausswin(width,alpha)';
% gauss_filter = gauss_filter/sum(gauss_filter); % Normalize.

%% Construct smoothing window
half_width = ceil(2*sig);
gauss_filter = gwin(sig);
gauss_filter = gauss_filter/sum(gauss_filter); % Normalize.

%% Loop over rows and smooth with gauss filter
sm_data = data;

for i = 1:size(data,1)
    
    % Smooth the current row
    sm_row = conv(data(i,:),gauss_filter);
    
    % Trim edge effect (ceiling of half the window width on either side)
    sm_data(i,:) = sm_row((half_width + 1):(end - half_width));
    
end