function [lp_filter, rp_filter, cp_filter, ip_filter, np_filter,...
    lp_coefs_for, lp_coefs_back,...
    rp_coefs_for, rp_coefs_back,...
    cp_coefs_for, cp_coefs_back,...
    ip_coefs_for, ip_coefs_back,...
    np_coefs_for, np_coefs_back, evt] = build_event_filters(fm)
% [lp_filter, rp_filter, cp_filter, ip_filter, np_filter,...
%     lp_coefs_for, lp_coefs_back,...
%     rp_coefs_for, rp_coefs_back,...
%     cp_coefs_for, cp_coefs_back,...
%     ip_coefs_for, ip_coefs_back,...
%     np_coefs_for, np_coefs_back, evt] = build_event_filters(fm) 
% Build event response filters for fitted GLM.

evprs.ncols = 5;
evprs.hpeaks = [0.05, 2];
evprs.b = 0.5;
[for_evt, evbasis] = makeBasis_PostSpike(evprs, fm.dt);

ev_list = {'Left', 'Right', 'Correct', 'Incorrect', 'Nosepoke'};

for i = 1:length(ev_list)
    
    ev_coefs = fm.full.x(strcmp(fm.basis_labels, ev_list{i}));
    ev_coefs_for = ev_coefs(1:evprs.ncols);
    ev_coefs_back = ev_coefs((evprs.ncols + 1):length(ev_coefs));
    
    curr_filter = [flipud(evbasis * ev_coefs_back);...
        evbasis(2:end, :) * ev_coefs_for];
    evt = [-flipud(for_evt); for_evt(2:end)];
    
    if i == 1
        lp_filter = curr_filter;
        lp_coefs_for = ev_coefs_for;
        lp_coefs_back = ev_coefs_back;
    elseif i == 2
        rp_filter = curr_filter;
        rp_coefs_for = ev_coefs_for;
        rp_coefs_back = ev_coefs_back;
    elseif i == 3
        cp_filter = curr_filter;
        cp_coefs_for = ev_coefs_for;
        cp_coefs_back = ev_coefs_back;
    elseif i == 4
        ip_filter = curr_filter;
        ip_coefs_for = ev_coefs_for;
        ip_coefs_back = ev_coefs_back;
    else
        np_filter = curr_filter;
        np_coefs_for = ev_coefs_for;
        np_coefs_back = ev_coefs_back;
    end
end