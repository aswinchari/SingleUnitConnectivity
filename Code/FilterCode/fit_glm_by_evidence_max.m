function [x, eta_opt, eta, evd_val, fitting_output] =...
    fit_glm_by_evidence_max(dt, rt, B, shrink_idx, eta, x0)
% [x, eta_opt, eta, evd_val, fitting_output] =...
%     fit_glm_by_evidence_max(dt, rt, B, shrink_idx, eta, x0)
% Fit spike train model with by maximum likelihod with ridge prior. Choose
% ridge penalty by evidence maximization.

%% Loop over grid and update evidence value
% Grid to search
if (nargin < 5) || isempty(eta)
    %eta = logspace(-4, 4, 5);
    eta = logspace(-10, 0, 6);
end

% Initial point for fitter
if nargin < 6 || isempty(x0)
    
    % Initialize at zeros
    x0 = zeros(size(B, 2), 1);
    
    % Initialize by linearization
%     A = B' * B * dt;
%     b = (sum(B(rt == 1, :)) - sum(B * dt))';
%     x0 = linsolve(A, b);
end

% Initialize optimal model
evd_opt = -Inf;
eta_opt = -Inf;

%% Compute optimal model for path of 
fitting_output = struct('x', [], 'fval', [], 'exitflag', [], 'output', [],...
    'grad', [], 'hessian', []);
fitting_output(length(eta)).x = [];

evd_val = zeros(size(eta));
h = waitbar(0, 'Fitting GLM by evidence maximization');
for i = 1:length(eta)
    % Fit model at current shrinkage intensity
    [x, fval, exitflag, output, grad, hessian] =...
        fit_glm_nested(dt, rt, B, eta(i), shrink_idx, x0);
    
    fitting_output(i).x = x;
    fitting_output(i).fval = fval;
    fitting_output(i).exitflag = exitflag;
    fitting_output(i).output = output;
    fitting_output(i).grad = grad;
    fitting_output(i).hessian = hessian;
    
    % Compute evidence at current shrinkage intensity
    evd_val(i) =...
        log(fval) + log(4 * pi * eta(i)) / 2 - log(det(hessian)) / 2;
    
    
    % Update maximum
    if evd_val(i) > evd_opt
        evd_opt = evd_val(i);
        eta_opt = eta(i);
    end
    
    % Update initialization point for solver. (We know that the solution
    % for the next shrinkage intensity will be close to this one.)
    x0 = x;
    
    % Update waitbar
    waitbar(i / length(eta));
end
try
close(h);
end

% Report coefficient values for optimal model
if length(eta) > 1
    
    eta,
    eta_opt,
    
    try
        x = fitting_output(eta == eta_opt).x;
    catch
        x = fitting_output(end).x;
    end
else
    x = fitting_output(1).x;
end

%% % Use derivative-free optimization
%     % Function handle
%     fun = @(eta) -post_evidence(dt, rt, B, 10^(eta), shrink_idx);
% 
%     % Optimize
%     l = -5;
%     u = 5;
% 
%     eta_opt = fminbnd(fun, l, u);
%     eta_opt = 10 ^ eta_opt;
%     
%     % Get optimal model
%     [x, fval, exitflag, output, grad, hessian] =...
%         fit_glm_nested(dt, rt, B, eta_opt, shrink_idx);
    
end