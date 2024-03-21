function win = gwin(sig)
% win = gwin(sig)
% Returns a Gaussian window of length 4*N+1, where N = ceil(2*sig). The window spans two standard
% deviations on either side of the center. The input SIG is the number of
% timepoints that comprise a standard deviation of the window.

% Normal distribution
f = @(x) (1/(sig*sqrt(2*pi)))*exp(-0.5*(x/sig).^2);

% Compute window
N = ceil(2*sig);
eval_pts = -N:N;

win = f(eval_pts);

% Normalize window
win = win/sum(win);