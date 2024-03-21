function bas_spike = spike_train_conv(rt, ihbas, for_back)
% ihbas_spike = spike_train_conv(rt, ihbas, for_or_back)
% Compute convolution of spike trains with post-spike basis. 

%% Parse inputs
if nargin < 3
    for_back = [true, false];
end

%% Post-spike filter basis convolution
bas_spike_for = [];
bas_spike_back = [];
if for_back(1)
    
    num_ihbas = size(ihbas, 2);
    bas_spike_for = zeros(length(rt), num_ihbas);
%     for i = 1:num_ihbas
%         basis_conv = conv([0, rt], ihbas(:, i));
%         bas_spike_for(:, i) = basis_conv(1:length(rt));
%     end

    for i = 1:num_ihbas

        filt = ihbas(:, i);
        filt = [zeros(size(filt)); filt];

        basis_conv = conv([0, rt], filt, 'same');
        bas_spike_for(:, i) = basis_conv(1:length(rt));
        
    end
end

if for_back(2)
    
    num_ihbas = size(ihbas, 2);
    bas_spike_back = zeros(length(rt), num_ihbas);
    for i = 1:num_ihbas
        
        filt = ihbas(:, i);
        filt = [flipud(filt); zeros(size(filt))];
        
        basis_conv = conv([0, rt], filt, 'same');
        bas_spike_back(:, i) = basis_conv(1:length(rt));
        
    end
end

bas_spike = [bas_spike_for, bas_spike_back];
