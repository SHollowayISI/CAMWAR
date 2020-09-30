%% Radar System Modeling

%{
    Sean Holloway
    2/6/2019
    Version 3
    Simulating radar system response for CAMWAR project.

    Simulator should be able to take RCS environmental scene as
    input and export signal to Radar Processing program.

    Version 3 restructures primary loop, speeding up process by factor of ~35 

    TODO: Implement correct chirp parameters
    TODO: Implement parameter import
    
    NOTE: Follows "End-to-End Radar System" guide and "Automotive Adaptive
    Cruise Control Using FMCW Technology" example on Mathworks site

    Working as of 2/6/2020

%}

%% Housekeeping
clear variables;
close all;
tic;

addpath(genpath('Functions'));
addpath(genpath('Scene Assets'));
addpath(genpath('MAT Files'));

c = physconst('LightSpeed');

%% Signal & Radar Variables

% Assigned Variables
fc = 77e9;                      % Operating frequency in Hz
fs = 37.5e6;                    % Sample frequency in Hz
tm = 62.5e-6;                   % Sweep time in seconds
bw = 0.088e9;                   % Sweep bandwidth in Hz
num_blocks = 8;                 % Number of chirps per frame
chirps_block = 2;               % Number of chirps per MIMO block
% samples_per_chirp = 2048;       % Number of samples to keep

% multiplexingMethod = 'TDM';
multiplexingMethod = 'CDM';
% multiplexingMethod = 'none';

% Correction
fs = floor(tm*fs)/tm;
samples_per_chirp = 2^floor(log2(tm*fs));

% Derived Variables
lambda = c/fc;
ts = 1/fs;
num_samples = floor(fs*tm);
drop_samples = num_samples - samples_per_chirp;
sweep_slope = bw/tm;
chirp_rate = 1/tm;
block_time = tm*chirps_block;
frame_time = block_time*num_blocks;
t_p = 0:ts:(tm-ts);
code = hadamard(chirps_block);

% Save Radar Variables
filename = 'MAT Files/SignalParams.mat';
save(filename);

%% Target Variables

% Manually insert targets
%{
% Static target placeholder position and velocity
% dist = 1000;
% dist_ang = [30, -45, 0];
% 
% target_pos{1} = [dist*cosd(dist_ang(1)); dist*sind(dist_ang(1)); 0];
% target_vel{1} = [0;0;0];
% target_rcs{1} = 1;
% 
% target_pos{2} = [dist*cosd(dist_ang(2)); dist*sind(dist_ang(2)); 0];
% target_vel{2} = [4; 0; 0];
% target_rcs{2} = 1;
% 
% target_pos{3} = [dist*cosd(dist_ang(3)); dist*sind(dist_ang(3)); 0];
% target_vel{3} = [-4; 0; 0];
% target_rcs{3} = 1;
%}

% Load target list from file
filename = 'ImageTargets3D.mat';
load(filename);

%% Transceiver Variables

% Transceiver centroid position and velocity
radar_pos = [0; 0; 0];
radar_vel = [0; 0; 0];
tx_pos{1} = [0; -8*lambda; 0];
tx_pos{2} = [0;  8*lambda; 0];

% Elevation sweep angle setup
if isempty(elevation_axis)
    N_e = 4;
    angle_res = 10;
    maximum_angle = 0;
    elevation_slices = maximum_angle:-angle_res:(maximum_angle - angle_res*(N_e-1));
else
    elevation_slices = elevation_axis;
end
        
% Antenna element parameters
tx_gain = db2pow(8.86);         % Transmit antenna gain in absolute
tx_az_bw = 74.61;               % Tx azimuth beamwidth in degrees
tx_el_bw = 57.86;               % Tx elevation beamwidth in degrees
tx_az_pow = 4;                  % Tx azimuth sinc power
tx_el_pow = 4;                  % Tx elevation sinc power

rx_gain = db2pow(12.3);         % Receive antenna gain in absolute
rx_az_bw = 69.17;               % Rx azimuth beamwidth in degrees
rx_el_bw = 26.09;               % Rx elevation beamwidth in degrees
rx_az_pow = 2;                  % Rx azimuth sinc power
rx_el_pow = 4;                  % Rx elevation sinc power
rx_el_tilt = -10;                 % Rx elevation tilt in degrees

% Transmitter parameters
tx_power = db2pow(12-30);       % Total transmit power in Watts
tx_ant_gain = 0;                % NOTE: Custom antenna gain not normalized

% Receiver parameters
rx_dev_gain = 21;               % Receiver chip gain (N/A?)
rx_nf = 4;                      % Receiver noise figure in dB
rx_ant_gain = 0;                % NOTE: Custom antenna gain not normalized

%% Signal Setup

% Generate FMCW waveform using Phased Array Toolbox
waveform = phased.FMCWWaveform('SampleRate', fs, 'SweepTime', tm, ...
    'SweepBandwidth', bw, 'NumSweeps', 1);
%TODO: Correct number of sweeps to provide

%% Target & Channel Setup

% Create objects for targets

% Model target position and motion
target_plat = phased.Platform('MotionModel', 'Velocity', 'InitialPosition', ...
    target_pos, 'Velocity', target_vel);


% Set up target to take RCS value at each object call
target = phased.RadarTarget('EnablePolarization', false, ...
    'MeanRCSSource', 'Input port', 'Model', 'Nonfluctuating', ...
    'OperatingFrequency', fc, 'PropagationSpeed', c);

% TODO: Implement non-free space environments

% Set up Tx and Rx free space path loss
tx_chan = phased.LOSChannel('PropagationSpeed', c, 'OperatingFrequency', fc, ...
    'TwoWayPropagation', false, 'SampleRate', fs);
rx_chan = phased.LOSChannel('PropagationSpeed', c, 'OperatingFrequency', fc, ...
    'TwoWayPropagation', false, 'SampleRate', fs);
% Set property 'SpecifyAtmosphere' to true to enable temperature, pressure,
% environment, etc.

%% Antenna Setup

% TODO: Update tilt shifting

az_axis = -180:0.1:180;
el_axis = -90:0.1:90;
phase_pattern = zeros(length(el_axis),length(az_axis));

% Create Tx vertical antenna pattern
tx_pattern = pow2db(sincAntennaPattern(az_axis, el_axis, 0, 0, ...
    tx_az_bw, tx_el_bw, tx_az_pow, tx_el_pow, tx_gain));

% Create single Tx antenna element
tx_antenna = phased.CustomAntennaElement('AzimuthAngles', az_axis, ...
    'ElevationAngles', el_axis, 'MagnitudePattern', tx_pattern, ...
    'PhasePattern', phase_pattern);

% Create 24-element vertical array of Tx antennas
tx_subarray = phased.ULA('Element', tx_antenna, 'NumElements', 24, ...
    'ArrayAxis', 'z', 'ElementSpacing', lambda/2);

% Set up beamformer weights
tx_angle = 0;
sv = steervec(getElementPosition(tx_subarray)/lambda, [0; -tx_angle]);
tx_subarray.Taper = sv;

% Create 2-element horizontal array of Tx subarrays
tx_array = phased.ReplicatedSubarray('Subarray', tx_subarray, ...
    'Layout', 'Rectangular', 'GridSize', [1,2], 'GridSpacing', ...
    32*lambda/2);%, 'SubarraySteering', 'Custom');
% tx_array = phased.ConformalArray('Element', tx_subarray, 'ElementPosition', ...
%     [tx_pos{1}, tx_pos{2}]);

% Create Rx horizontal antenna pattern
rx_pattern = pow2db(sincAntennaPattern(az_axis, el_axis, 0, rx_el_tilt, ...
    rx_az_bw, rx_el_bw, rx_az_pow, rx_el_pow, rx_gain));

% Create single Rx antenna element
rx_antenna = phased.CustomAntennaElement('AzimuthAngles', az_axis, ...
    'ElevationAngles', el_axis, 'MagnitudePattern', rx_pattern, ...
    'PhasePattern', phase_pattern);

% Create 32-element horizontal array of Rx antennas
rx_array = phased.ULA('Element', rx_antenna, 'NumElements', 32, ...
    'ElementSpacing', lambda/2);


%% Transceiver Setup

% Transceiver position and motion
radar_plat = phased.Platform('MotionModel', 'Velocity', 'InitialPosition', ...
    radar_pos, 'Velocity', radar_vel);

% Set up transmitter parameters
% Parameters from TI AWR1243P T/R Module
transmitter = phased.Transmitter('PeakPower', tx_power, 'Gain', tx_ant_gain);
radiator = phased.Radiator('Sensor', tx_subarray, 'PropagationSpeed', c, ...
    'OperatingFrequency', fc, 'CombineRadiatedSignals', true);

% Set up receiver parameters
% Parameters from TI AWR1243P T/R Module
collector = phased.Collector('Sensor', rx_array, 'PropagationSpeed', c, ...
    'OperatingFrequency', fc, 'Wavefront', 'Plane');
receiver = phased.ReceiverPreamp('Gain', rx_dev_gain + rx_ant_gain,  ...
    'NoiseFigure', rx_nf, 'SampleRate', fs);

%% Simulation

% Allocate array for received signal, 1 frame long
% Dimension 1 is fast time, dimension 2 is slow time, dimension 3 is
% receive element, dimension 4 is transmit element

% When using subarrays, each subarray is single channel.
rx_sig = zeros(num_samples, num_blocks, prod(rx_array.NumElements), ...
    chirps_block, length(elevation_slices));

% Allocate data for signal variables
rad_sig = zeros(num_samples, 1, chirps_block);
ref_sig = zeros(num_samples, length(target_pos), chirps_block);
col_sig = zeros(num_samples, prod(rx_array.NumElements), chirps_block);

% Set up time reporting variables
startTime = tic;
startTimeDate = now;
timeGate = 0;

% Override elevation data
elevation_slices = -30:2:0;

% Loop through MIMO blocks in full frame
for j = 1:length(elevation_slices)
    
    % Set up beamformer weights
    tx_angle = elevation_slices(j);
    sv = steervec(getElementPosition(tx_subarray)/lambda, [0; -tx_angle]);
    tx_subarray.Taper = sv;
    
    for n = 1:num_blocks
        
        % Loop through chirps in each MIMO block
        for m = 1:chirps_block
            
            % Update the target position
            [tgt_pos,tgt_vel] = target_plat(tm);
      
            % Loop through transmitters in each chirp
            for k = 1:chirps_block
                
                % Get the range and angle to the target
                [tgt_rng,tgt_ang] = rangeangle(tgt_pos,tx_pos{k});
                
                % Generate the pulse
                rad_sig(:,:,k) = code(m,k)*waveform();
                
                % Transmit the pulse. Output transmitter status
                rad_sig(:,:,k) = transmitter(rad_sig(:,:,k));
                
                % Radiate the pulse toward the target
                ref_sig(:,:,k) = radiator(rad_sig(:,:,k),tgt_ang);
                
                % Propagate the pulse to the target in free space
                ref_sig(:,:,k) = tx_chan(ref_sig(:,:,k),tx_pos{k},tgt_pos,radar_vel,tgt_vel);
                
                % Reflect the pulse off the target
                ref_sig(:,:,k) = target(ref_sig(:,:,k), target_rcs);
                
                % Propagate the echo to the antenna in free space
                ref_sig(:,:,k) = rx_chan(ref_sig(:,:,k),tgt_pos,tx_pos{k},tgt_vel,radar_vel);
                
                % Collect the echo from the incident angle at the antenna
                col_sig(:,:,k) = collector(ref_sig(:,:,k),tgt_ang);
                
            end
            
            % Add together two code-division signals
            col_sig = sum(col_sig,3);
            
            % Receive the echo at the antenna
            col_sig = receiver(col_sig);
            
            % Dechirp the received signal
            rx_sig(:,n,:,m,j) = dechirp(col_sig, waveform());
            
            % Estimate time of completion
            loops = (num_blocks*chirps_block*(j-1) + chirps_block*(n-1) + m);
            completePercent = 100*loops/(chirps_block*num_blocks*length(elevation_slices));
            nowTimeDate = now;
            elapsedTime = nowTimeDate - startTimeDate;
            estComplete = nowTimeDate + ((100/completePercent)-1)*elapsedTime;
            
            % Display current progress
            minInterval = 2;
            if (timeGate == 0) || (toc > timeGate + minInterval)
                toc(startTime);
                
                fprintf('Number of iterations complete: %d of %d\n', ...
                    loops, chirps_block*num_blocks*length(elevation_slices));
                timeGate = toc;
            
                dispMessage = ['Current time: ' datestr(now)];
                disp(dispMessage);
                
                dispMessage = ['Estimated time of completion: ' datestr(estComplete)];
                disp(dispMessage);
            end
            
        end
    end
end




disp('Simulation complete, saving file...')

%% Visualization

% Plot I/Q return signal of first chirp
%{
figure('Name', 'Response of first chirp')
plot(t_p, sum(real(rx_sig(:,:,1)),2), t_p, sum(imag(rx_sig(:,:,1)),2));
xlabel = 'Time [s]';
ylabel = 'Amplitude';
%}

%% Export Receive Signal
%
filename = 'MAT Files/SimRxSignal.mat';
save(filename, 'rx_sig', 'elevation_slices', '-v7.3');
toc;
disp(['Signal saved as "', filename, '"'])
%}