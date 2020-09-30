function [detection] = GenerateRCSCube(scenario)
%GENERATERCSCUBE Calculates RCS of targets for CAMWAR project
%   Takes scenario object as input, provides scenario.detection object as
%   output, containing RCS estimation of detected targets.

%% Unpack Variables

detection = scenario.detection;
radarsetup = scenario.radarsetup;
cube = scenario.cube;
simsetup = scenario.simsetup;
sim = scenario.sim;

%% Calculate RCS Cube

% Set up constants for calculation
lambda = physconst('LightSpeed')/radarsetup.f_c;
t_cpi = radarsetup.t_ch*radarsetup.n_p*radarsetup.n_tx_sub;
Lr = db2pow(radarsetup.rf_sys_loss);
NF = db2pow(radarsetup.rx_nf);
const = (4*pi)^3;
np = physconst('Boltzmann')*290;

% Set up angular axes for antenna gain calculation
az_idx = round(cube.azimuth_axis*10) + 1801;
el_idx = round(simsetup.elevation_slices*10) + 901;

% Calculate antenna gain
tx_gain = ones(length(el_idx), length(az_idx));
rx_gain = ones(length(el_idx), length(az_idx));

for n = 1:length(el_idx)
    tx_gain(n,:) = sim.tx_pattern(el_idx(n), az_idx);
    rx_gain(n,:) = sim.rx_pattern(el_idx(n), az_idx);
end

tx_gain = db2pow(tx_gain - max(sim.tx_pattern, [], 'all'));
rx_gain = db2pow(rx_gain - max(sim.rx_pattern, [], 'all'));

Gt = radarsetup.tx_gain*radarsetup.n_tx_ant*radarsetup.n_tx_sub .* ...
    permute(tx_gain, [3 2 1]);
Gr = radarsetup.rx_gain*radarsetup.n_rx_ant .* ...
    permute(rx_gain, [3 2 1]);

% Calculate scale factor between power and RCS
detection.scale_factor = (const * (cube.range_axis'.^4) * np * NF * Lr) ./ ...
    (radarsetup.tx_pow * Gt .* Gr * t_cpi * lambda * lambda);

% Calculate full RCS cube
detection.RCS_cube = detection.SNR_cube .* detection.scale_factor;

%% Calculate Surface Reflectivity Cube

% Apply scaling factor to approximate surface reflectivity
range_diff = diff(scenario.cube.range_axis'.^2);
range_diff = [scenario.cube.range_axis(1)^2; range_diff];

detection.surface_RCS = squeeze(sum(detection.RCS_cube./range_diff, 1));
detection.surface_RCS = pow2db(detection.surface_RCS .* sind(-scenario.simsetup.elevation_slices));
detection.surface_RCS = detection.surface_RCS - pow2db(2*deg2rad(scenario.cube.azimuth_res)) - 3;

end

