function fitted_models = glm_single_neuron_fit(raster, TS)
% fitted_models = dnms_single_neuron_fit(raster, TS, event_data)
% Fit GLMs for multiple sources of spike dependence including post-spike
% filters 

%% Bin raster into 1 ms time bins
dt = 1e-3;
rt = bin_raster(raster, dt, 0, minmax(TS));
samp_ts = (0:(size(rt, 2) - 1)) * dt - dt / 2 + min(TS);

%% Basis convolutions with spike trains
% Get post spike basis
ihprs.ncols = 10;
ihprs.hpeaks = [.01, 0.5];
ihprs.b = 0.5;
[iht, ihbasis] = makeBasis_PostSpike(ihprs, dt);

%% Construct post-spike filter convolutions
bas_spike = cell(1, length(raster));
bas_spike{1} = spike_train_conv(rt, ihbasis);

% Refractory basis for spike
nrefrac = 1;
refrac_basis = [eye(nrefrac); zeros(length(iht) - nrefrac, nrefrac)];
refrac_spike = spike_train_conv(rt, refrac_basis);

back = ones(size(rt, 2), 1);

covar_basis = [back];

covar_labels = ['Background'];

%% Post-spike bases
pspk_basis = [refrac_spike, cell2mat(bas_spike)];
pspk_labels = [repmat({'Refractory'}, [1, size(refrac_basis, 2)]),...
    repmat({'Postspike'}, [1, ihprs.ncols])];

%% Assemble full basis
B = [covar_basis, pspk_basis];
basis_labels = [covar_labels, pspk_labels];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimize by evidence maximization

shrink_idx = 2:size(B, 2);

fitted_models = [];
fitted_models.dt = dt;
fitted_models.rt = rt;
fitted_models.samp_ts = samp_ts;
%fitted_models.B = B;
fitted_models.basis_labels = basis_labels;
fitted_models.shrink_idx = shrink_idx;

%% Full model
[x, eta_opt, eta, evd_val, fitting_output] =...
    fit_glm_by_evidence_max(dt, rt', B,...
    2:size(B, 2));

fitted_models.full = [];
fitted_models.full.x = x;
fitted_models.full.eta_opt = eta_opt; 
fitted_models.full.eta = eta;
fitted_models.full.evd_val = evd_val;
fitted_models.full.fitting_output = fitting_output;

fitted_models.full_intensity = exp(B * x);

%% Poisson firing
firing_rate = sum(rt) / diff(minmax(TS));
poiss_intensity = ones(size(rt')) * firing_rate;

fitted_models.poiss_intensity = poiss_intensity;


