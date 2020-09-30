function [scenario_out] = RadarSimulation_CAMWAR(scenario)
%RADARSIMULATION_CAMWAR Generates simulated radar response for CAMWAR
%   Takes scenario.sim, .simsetup, .radarsetup, .environment as inputs and
%   outputs scenario.rx_sig struct containing the received signal.

%% Unpack Variables

scenario_out = scenario;
target_list = scenario.target_list;
sim = scenario.sim;
simsetup = scenario.simsetup;
radarsetup = scenario.radarsetup;
flags = scenario.flags;

%% Setup

% Physical constants
c = physconst('LightSpeed');

% Derived variables
lambda = c/radarsetup.f_c;


%% Phase Shift Setup

% Set up elevation beam steering
sv = steervec(getElementPosition(sim.tx_subarray)/lambda, ...
    [0; -1*simsetup.elevation_slices(flags.slice)]);

% Set up MIMO antenna weighting
code = hadamard(radarsetup.n_tx_sub);
weights = sv .* permute(code, [3 1 2]);

%% Allocate Targets in Slice

if simsetup.allocate_slice || simsetup.allocate_az
    
    % Determine target angles
    [tgt_rng, tgt_ang] = rangeangle(target_list.pos, simsetup.radar_pos);
    
    % Generate indices of targets in desired area
    tgt_ind = true(size(tgt_ang(1,:)));
    
    % Remove targets outside of elevation slice
    if simsetup.allocate_slice
        tgt_ind = (abs(tgt_ang(2,:) - simsetup.elevation_slices(flags.slice)) ...
            < simsetup.slice_bw/2);
    end
    
    % Remove targets outside of azimuth FOV
    if simsetup.allocate_az
        tgt_ind = tgt_ind & (abs(tgt_ang(1,:)) < ...
        (simsetup.fov/2 + simsetup.fov_bw));
    end
    
    % Remove targets outside of maximum range
    r_max = radarsetup.t_ch*c/2;
    tgt_ind = tgt_ind & (tgt_rng < r_max);
    
    % Generate list of targets in elevation slice
    target_list.pos = target_list.pos(:,tgt_ind');
    target_list.vel = target_list.vel(:,tgt_ind');
    target_list.rcs = target_list.rcs(tgt_ind');
    
    % Reset radiator phased array object
    targetListReset(scenario);
    
    if numel(target_list.pos) == 0
        
        % If no target, create dummy target with 0 RCS
        target_list.pos = [1000; 0; 0];
        target_list.vel = [0; 0; 0];
        target_list.rcs = 0;
    end 
    
end
    

% Create platform object with target list
sim.target_plat = phased.Platform( ...
    ...
    'MotionModel',              'Velocity', ...
    'InitialPosition',          target_list.pos, ...
    'Velocity',                 target_list.vel);

%% Simulation

% Allocate size for signals
rx_sig = zeros(radarsetup.n_s, radarsetup.n_p, radarsetup.n_rx_ant*radarsetup.n_tx_sub);

% Generate the pulse
tx_sig = sim.waveform();

% Transmit the pulse
tx_sig = sim.transmitter(tx_sig);

% Get range and angle to all targets
[~, tgt_ang] = rangeangle(target_list.pos, simsetup.radar_pos);

for block = 1:radarsetup.n_p
    
    for chirp = 1:radarsetup.n_tx_sub
        
        % Update target position and velocities
        [target_list.pos, target_list.vel] = sim.target_plat(radarsetup.t_ch);
        
        % Radiate signal towards all targets
        sig = sim.radiator(tx_sig, tgt_ang, weights(:,:,chirp));
        
        % Propogate signal to the target through two-way environmental channel
        sig = sim.channel(sig, ...
            simsetup.radar_pos, target_list.pos, ...
            simsetup.radar_vel, target_list.vel);
        
        % Reflect the pulse off of the target
        sig = sim.target(sig, target_list.rcs, true);
        
        % Collect the reflected target at the antenna array
        sig = sim.collector(sig, tgt_ang);
        
        % Apply receiver noise and gains
        sig = sim.receiver(sig);
        
        % Dechirp signal
        sig = dechirp(sig,tx_sig);
        
        % Apply phase coding to result signal
        mimo_sig = reshape(sig.*permute(code(:,chirp), [2 3 1]), length(sig), 1, []);
        
        % Save Rx signal by fast time x slow time x virtual antenna
        rx_sig(:,block,:) = rx_sig(:,block,:) + mimo_sig;
        
    end
    
end


%% Re-pack Variables

scenario_out.rx_sig = rx_sig;

% Note: Do not export target list, since it is modified during


end