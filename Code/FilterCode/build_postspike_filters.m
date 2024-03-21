function [psf, iht] = build_postspike_filters(fm)
% [psf, iht] = build_postspike_filters(fm)
% Build post-spike filter for a fitted model.

% Post-spike basis
ihprs.ncols = 10;
ihprs.hpeaks = [.01, 0.5];
ihprs.b = 0.5;
[iht, ihbasis] = makeBasis_PostSpike(ihprs, fm.dt);

% Make PSF    
psf_coefs = fm.full.x(strcmp(fm.basis_labels, 'Postspike'));
psf = ihbasis * psf_coefs;
    