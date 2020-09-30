%% Radar Signal Processing for CAMWAR Project
%{
    
    Sean Holloway
    2/7/2020
    Version 3
    Radar signal processing for CAMWAR project.

    Range-Doppler-Angle processing for CAMWAR signal
    Uses result of RadarTransceiverModel_v3.m

    Version 3 separates signal processing and data processing. Output is
    full Range-Doppler-Azimuth-Elevation data cube, with detection 
    performed in RadarDataProcessing_v1.m

    Working as of 2/7/2020

%}

%% Housekeeping
% clear variables
clear xlabel;
clear ylabel;
close all;
tic

addpath(genpath('Functions'));
addpath(genpath('Scene Assets'));
addpath(genpath('MAT Files'));

c = physconst('LightSpeed');

%% Variables

% Assigned Variables
fc = 77e9;                      % Operating frequency in Hz
fs = 37.5e6;                    % Sample frequency in Hz
tm = 62.5e-6;                   % Sweep time in seconds
bw = 0.088e9;                   % Sweep bandwidth in Hz
num_blocks = 8;                 % Number of chirps per frame
chirps_block = 2;               % Number of chirps per MIMO block
% samples_per_chirp = 2048;       % Number of samples to keep

% MIMO Settings
num_rx_channels_x = 32;
num_rx_channels_y = 1;

% multiplexingMethod = 'TDM';
multiplexingMethod = 'CDM';
% multiplexingMethod = 'none';


% Overwrite variables from file
filename = 'SignalParams.mat';
load(filename);

% Derived Variables
lambda = c/fc;
ts = 1/fs;
num_samples = floor(fs*tm);
drop_samples = num_samples - samples_per_chirp;
sweep_slope = bw/tm;
chirp_rate = 1/tm;
block_time = tm*chirps_block;
frame_time = block_time*num_blocks;

range_res = (num_samples/samples_per_chirp)*c/(2*bw);
vel_res = lambda/(2*frame_time);

switch multiplexingMethod
    case {'TDM', 'CDM'}
        num_tx_channels_x = 2;
        num_tx_channels_y = 1;
    case 'none'
        num_tx_channels_x = 1;
        num_tx_channels_y = 1;
end

%% Load signal from file

% Import signal from RadarTransceiverModel
cube = load('SimRxSignal.mat','rx_sig','elevation_slices');
elevation_slices = cube.elevation_slices;
elev_axis = elevation_slices;
cube = cube.rx_sig;


%% Run signal processing loop for each elevation slice

for slice = 1:length(elevation_slices)

    % Ignore additional Tx channels if not using MIMO
    switch multiplexingMethod
        case {'TDM', 'CDM'}
            chirps=cube(end-samples_per_chirp+1:end,1:num_blocks,:,:,slice);
        case 'none'
            chirps=cube(end-samples_per_chirp+1:end,1:num_blocks,:,1,slice);
    end

    
    %% Range FFT & Windowing
    
    % Set Range FFT Size
    %N_r = 1024;
    N_r = samples_per_chirp;
    
    % Perform Range FFT in Range_Calc
    [range_bins, k_r] = Range_Calc(chirps, N_r);
    range_axis = (k_r-0.5) * range_res;
    % toc;
    
    %% Process MIMO Blocks
    
    % Break blocks into virtual channels
    
    switch multiplexingMethod
        case 'TDM'
            mimo_range_bins = reshape(range_bins, size(range_bins,1), ...
                size(range_bins,2), []);
            
        case 'CDM'
            mimo_range_bins(:,:,1:size(range_bins,3)) = ...
                range_bins(:,:,:,1) + range_bins(:,:,:,2);
            mimo_range_bins(:,:,(size(range_bins,3)+1):(2*size(range_bins,3))) = ...
                range_bins(:,:,:,1) - range_bins(:,:,:,2);
            
        case 'none'
            mimo_range_bins = range_bins;
    end
    
    % toc;
    
    %% Doppler FFT & Windowing
    
    % Set Doppler FFT Size
    N_d = num_blocks;
    % N_d = num_blocks * 4;
    % N_d = 128;
    
    % Perform Doppler FFT in Doppler_Calc
    [doppler_bins, k_d] = Doppler_Calc(mimo_range_bins, N_d);
    doppler_axis = k_d * vel_res * num_blocks/N_d;
    % toc;
    
    %% Angle FFT & Windowing
    
    % Set Angle FFT Size
    % N_a = 4;
    N_a = 64;
    % N_a = num_rx_channels_x * num_tx_channels_x;
    % N_a = num_rx_channels_x * num_tx_channels_x * 4;
    
    % Perform Angle FFTs in Angle_Calc
    [angle_bins(:,:,:,slice), k_a] = Angle_Calc(doppler_bins, N_a);
    azimuth_axis = -1*asin((k_a)*(2/N_a))*(180/pi);
    % toc;
    
end


%% Plot Range-Doppler Intensity Heatmap
%{
figure('Name', 'Range-Doppler Intensity Heatmap')
imagesc(doppler_axis, range_axis, abs(angle_bins(:,:,ceil(end/2))).^2)
set(gca,'YDir','normal')
xlabel('Velocity [m/s]')
ylabel('Range [m]')
%}

%% Plot Range-Azimuth Intensity Heatmap
%{
view_slice = min(1, length(elevation_slices));
x_axis = range_axis'.*sind(azimuth_axis);
y_axis = range_axis'.*cosd(azimuth_axis);

figure('Name', 'Range-Azimuth Intensity Heatmap')
surf(x_axis, y_axis, permute(sum(abs(angle_bins(:,:,:,view_slice)).^2, 2), [1 3 2]), ...
    'EdgeColor', 'none')
% surf(x_axis, y_axis, 10*log10(permute(sum(abs(angle_bins).^2, 2), [1 3 2])), ...
%     'EdgeColor', 'none')
set(gca,'YDir','normal')
xlabel('Cross Range [m]')
ylabel('Range [m]')
xlim([-2000 2000])
ylim([0 2000])
view(2)
%}

%% Save full cube to file
%
disp('Processing complete, saving file...')
filename = 'MAT Files/processedCube.mat';
save(filename,'angle_bins', 'range_axis', 'doppler_axis', 'azimuth_axis', 'elev_axis', '-v7.3');
toc;
disp(['Signal saved as "', filename, '"'])
%}



