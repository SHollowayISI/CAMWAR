function [SNR] = CalculateSNR(scenario,RCS,Range)
%CALCULATESNR Calculates SNR of target for CAMWAR scenario
%   Takes radar scenario object, RCS (absolute) value, and Range 
%   as input, provides SNR value as output.

%% Unpack Variables
radarsetup = scenario.radarsetup;

%% Calculate SNR

lambda = physconst('LightSpeed')/radarsetup.f_c;
t_cpi = radarsetup.t_ch*radarsetup.n_p*radarsetup.n_tx_sub;
Lr = db2pow(radarsetup.rf_sys_loss);
NF = db2pow(radarsetup.rx_nf);
c = (4*pi)^3;
n = physconst('Boltzmann')*290;

Gt = radarsetup.tx_gain*radarsetup.n_tx_ant*radarsetup.n_tx_sub;
Gr = radarsetup.rx_gain*radarsetup.n_rx_ant;


SNR_abs = (radarsetup.tx_pow * Gt * Gr * t_cpi * lambda * lambda * RCS) ...
    ./ (c * (Range.^4) * n * NF * Lr);

SNR = pow2db(SNR_abs);

end

