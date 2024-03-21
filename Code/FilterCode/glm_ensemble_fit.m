function fitted_models = glm_ensemble_fit(raster, TS, eta, target_neuron)
% fitted_models = dnms_ensemble_fit(raster, TS, eta, target_neuron)
% Fits coupling filters for all cells in an ensemble onto a target neuron.
% raster is cell array with each cell being spikes (in seconds)
% TS is timestamps
% eta can set to 0 to begin with


%% Initializations
num_neuron = length(raster);

%% Bin raster into 1 ms time bins
dt = 1e-3;
rt = bin_raster(raster, dt, 0, minmax(TS));
samp_ts = (0:(size(rt, 2) - 1)) * dt - dt / 2 + min(TS);

%% Basis convolutions with spike trains
% Get post spike basis
ihprs.ncols = 10;           % number of basis functions
ihprs.hpeaks = [.01, 0.5];  % where they are
ihprs.b = 0.5;
[iht, ihbasis] = makeBasis_PostSpike(ihprs, dt);

nrefrac = 1;                % basis that enforces refractory period
refrac_basis = [eye(nrefrac); zeros(length(iht) - nrefrac, nrefrac)];

%% Construct post-spike filter convolutions for neuron 1
spike_basis = cell(1, num_neuron);
spike_basis{target_neuron} = spike_train_conv(rt(target_neuron, :), ihbasis);

% Refractory basis for spike
spike_basis{target_neuron} =...
    [spike_train_conv(rt(target_neuron, :), refrac_basis),...
    spike_basis{target_neuron}];

pspk_labels = cell(1, num_neuron);
pspk_labels{target_neuron} = [{'Refractory'},...
    repmat({'PSF'}, [1, ihprs.ncols])];

%% Construct post-spike filter convolutions for source neurons

for nrn = setdiff(1:num_neuron, target_neuron)
    spike_basis{nrn} = spike_train_conv(rt(nrn, :), ihbasis);
    pspk_labels{nrn} = repmat({['CPF', num2str(nrn)]}, [1, ihprs.ncols]);
end

%% Build basis
back = ones(size(rt, 2), 1);
covar_basis = back;
covar_labels = {'Background'};

%% Assemble full basis
B = [covar_basis, cell2mat(spike_basis)];
basis_labels = [covar_labels, cat(2, pspk_labels{:})];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimize by evidence maximization

shrink_idx = 2:size(B, 2);

fitted_models = [];
fitted_models.dt = dt;
fitted_models.rt = rt;
fitted_models.samp_ts = samp_ts;
fitted_models.basis_labels = basis_labels;
fitted_models.shrink_idx = shrink_idx;

%% Full model for target neuron

% Regularization grid
if (nargin < 6) || (~isempty(eta))
    eta = logspace(-10, -2, 3);
end

%% Starting point
x0 = zeros(size(B, 2), 1);

% Fit model for neuron 1
[x, eta_opt, eta, evd_val, fitting_output] =...
    fit_glm_by_evidence_max(dt, rt(target_neuron, :)', B,...
    2:size(B, 2), eta, x0);

%%
fitted_models.x= x;
fitted_models.eta_opt = eta_opt; 
fitted_models.eta = eta;
fitted_models.evd_val = evd_val;
fitted_models.fitting_output = fitting_output;
fitted_models.full_intensity = exp(B * fitted_models.x);

%% Poisson firing
firing_rate = sum(rt, 2) / diff(minmax(TS));
poiss_intensity = ones(size(rt')) * diag(firing_rate);

fitted_models.poiss_intensity = poiss_intensity;

