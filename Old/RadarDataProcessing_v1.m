%% Radar Data Processing for CAMWAR Project
%{
    
    Sean Holloway
    2/7/2020
    Version 1
    Radar data processing for CAMWAR project.

    Data processing using target list from RadarDetection_v1 and full radar
    cube from RadarSignalProcessing_v3

    TODO: Terrain map algorithm
    TODO: Monopulse angle processing

    Not working as of 2/18/2020

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


%% Load signal from file

% Import signal from RadarSignalProcessing
load('detectionCube.mat');

% Optional: Import original target list
tgt_imp = load('ImageTargets3D.mat', 'tgt_exp');
tgt_imp = tgt_imp.tgt_exp;


%% Height above terrain (HAT) algorithm



%% Data Visualization


%% Save data to file
%{
disp('Processing complete, saving file...')
filename = 'MAT Files/processedCube.mat';
save(filename,'angle_bins', 'range_axis', 'doppler_axis', 'azimuth_axis', 'elev_axis', '-v7.3');
toc;
disp(['Signal saved as "', filename, '"'])
%}



