%% Radar Detection for CAMWAR Project
%{
    
    Sean Holloway
    2/7/2020
    Version 1
    Radar detection system for CAMWAR project.

    Target detection using Range-Doppler-Azimuth-Elevation cube provided by
    RadarSignalProcessing_v3.m

    TODO: CFAR detection
    TODO: Binary integration
    TODO: Doppler beam sharpening?

    Working as of 2/18/2020

%}

%% Housekeeping
clear variables
clear xlabel;
clear ylabel;
close all;
tic

addpath(genpath('Functions'));
addpath(genpath('Scene Assets'));
addpath(genpath('MAT Files'));

c = physconst('LightSpeed');

%% Variables

threshold = 14.4;               % Detection threshold in dB

%% Load signal from file

% Import signal from RadarSignalProcessing
load('processedCube.mat');

% Optional: Import original target list
tgt_imp = load('ImageTargets3D.mat', 'tgt_exp');
tgt_imp = tgt_imp.tgt_exp;

%% Signal Pre-Processing

% Determine power in each angle bin
power_bins = abs(angle_bins).^2;

% Incoherent sum over Doppler bins
dopsum_bins = squeeze(sum(power_bins, 2));

% Zero-Doppler radar cube
% zerodop_bins = squeeze(power_bins(:,ceil(end/2),:,:));

%% Threshold Target Detection

% Determine noise power
noise_power = -80;

% Set target detection threshold
thresh_abs = db2pow(noise_power + threshold);

% Detect bins above threshold
targets_threshold = double(power_bins >= thresh_abs);
zerodop_targets = squeeze(targets_threshold(:,ceil(end/2),:,:));


%% CFAR Target Detection

% TODO: Implement

targets_CFAR = zeros(size(power_bins));

%% Data Visualization

% Plot Range-Azimuth Intensity Heatmap
%{
view_slice = ceil(length(elevation_slices)/2);
x_axis = range_axis'.*sind(azimuth_axis);
y_axis = range_axis'.*cosd(azimuth_axis);

figure('Name', 'Range-Azimuth Intensity Heatmap')
subplot(121)
surf(x_axis, y_axis, dopsum_bins(:,:,view_slice), ...
    'EdgeColor', 'none')
% surf(x_axis, y_axis, 10*log10(permute(sum(abs(angle_bins).^2, 2), [1 3 2])), ...
%     'EdgeColor', 'none')
set(gca,'YDir','normal')
xlabel('Cross Range [m]')
ylabel('Range [m]')
xlim([-2000 2000])
ylim([0 2000])
view(2)


% Plot Detected Targets in Range-Azimuth Intensity Heatmap

subplot(122)
surf(x_axis, y_axis, zerodop_targets(:,:,view_slice), 'EdgeColor', 'none')
set(gca,'YDir','normal')
xlabel('Cross Range [m]')
ylabel('Range [m]')
xlim([-2000 2000])
ylim([0 2000])
view(2)
%}

%% Display 3D Scatter Plot of Detections
%{

x3_axis = range_axis'.*sind(azimuth_axis).*cosd(permute(elev_axis, [1 3 2]));
y3_axis = range_axis'.*cosd(azimuth_axis).*cosd(permute(elev_axis, [1 3 2]));
z3_axis = range_axis'.*ones(1,length(azimuth_axis)).*sind(permute(elev_axis, [1 3 2]));

target_x = x3_axis(zerodop_targets == 1);
target_y = y3_axis(zerodop_targets == 1);
target_z = z3_axis(zerodop_targets == 1);

target_list = [target_x, target_y, target_z];


figure('Name', '3D Detection Scatter Plot')
subplot(122)
scatter3(target_x, target_y, target_z)
xlim([-500 500])
ylim([0 1000])
zlim([-500 500])
% xlim([-200 200])
% ylim([300 700])
% zlim([-200 200])
title('Target Detections')

subplot(121)
scatter3(tgt_imp.x, tgt_imp.y, tgt_imp.z);
xlim([-500 500])
ylim([0 1000])
zlim([-500 500])
% xlim([-200 200])
% ylim([300 700])
% zlim([-200 200])
title('Distributed Target Locations')
%}

%% Save target list to file
%
disp('Processing complete, saving file...')
filename = 'MAT Files/detectionCube.mat';
save(filename,'targets_CFAR', 'targets_threshold', 'range_axis', 'doppler_axis', 'azimuth_axis', 'elev_axis', '-v7.3');
toc;
disp(['Signal saved as "', filename, '"'])
%}



