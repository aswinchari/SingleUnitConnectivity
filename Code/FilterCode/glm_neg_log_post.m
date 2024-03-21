function [f, g, H] = glm_neg_log_post(x, dt, rt, B, eta, shrink_idx)
% [f, g, H] = glm_neg_log_post(x, dt, rt, B, eta, shrink_idx)
% Compute the value and derivatives of the negative log-posterior.

%% Parse inputs
if isempty(eta)
    eta = 1;
end

if nargin < 6
    shrink_idx = 1:length(x);
end

num_bas = size(B, 2);

%% Compute log-likelihood
log_lam = B * x;
lam = exp(log_lam);

log_lik = sum(rt .* (log_lam + log(dt)) - lam * dt);

% Compute log-prior
log_prior = -eta * sum(x(shrink_idx) .^ 2);

% Compute log-posterior
log_post = log_lik + log_prior;

%% Compute gradient
% Gradient of log-likelihood
g = sum(repmat(rt - dt * lam, [1, num_bas]) .* B);

% Shrinkage term
shrink_g = zeros(size(g));
shrink_g(shrink_idx) = x(shrink_idx);

% Combine
g = g - 2 * eta * shrink_g;

%% Compute Hessian
% Hessian of log-likelihood
H = -(repmat(lam * dt, [1, num_bas]) .* B)' * B;

% Shrinkage term
shrink_H = zeros(size(H));
shrink_H(shrink_idx, shrink_idx) = eye(length(shrink_idx));

% Combine
H = H - 2 * eta * shrink_H;

%% Negate output

f = -log_post;
g = -g;
H = -H;




