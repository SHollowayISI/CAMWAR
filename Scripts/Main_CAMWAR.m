%% CAMWAR Radar System - Main Simulation Loop
%{

    Sean Holloway
    CAMWAR Main Simulation Loop
    
    This file specifies performs simulation, signal processing, detection,
    data processing, and results collection for CAMWAR system.

    Use script 'FullSystem_CAMWAR.m' to run scenarios.
    
%}

%% Main Loop

% Start timing for estimation
timeStart(scenario);

% Loop through elevation slices
for elevation_slice = 1:length(scenario.simsetup.elevation_slices)
    
    % Set current slice flag
    scenario.flags.slice = elevation_slice;
    
    %% Radar Simulation
    
    % Run simulation to retrieve fast time x slow time x antenna Rxsignal
    scenario = RadarSimulation_CAMWAR(scenario);
    
    %% Signal Processing
    
    % Perform signal processing on received signal
    scenario.cube = SignalProcessing_CAMWAR(scenario);
    
    % Generate 2-D Cartesian Mesh Grid
    generateCoordinates2D(scenario);
    
    %% Single Slice Data Processing
    
    % Perform single-frame radar detection
    scenario.detection = DetectionSingle_CAMWAR(scenario);
    
    % Read out peak SNR
    readOut(scenario, 'max');
    
    %% Single Slice Visualization
    
    % View Range-Doppler heat map of center azimuth direction
    % viewRDCube(scenario, 'heatmap')
    
    % View Range-Doppler surface of center azimuth direction
    % viewRDCube(scenario, 'surface')
    
    % View Range-Angle heat map of zero-doppler slice
    % viewRACube(scenario, 'heatmap')
    
    % View Range-Angle surface of zero-doppler slice
    % viewRACube(scenario, 'surface')
    
    % View Range-Angle PPI of zero-doppler slice
    % viewRACube(scenario, 'PPI')
    
    % View single-frame detection heatmap
    % viewDetections(scenario, 'heatmap')

    % View single-frame detection PPI
    % viewDetections(scenario, 'PPI')

    %% Loop Update Procedures
    
    % Read out estimated time of completion
    timeUpdate(scenario, 1, 'loops')
    
end

%% Full Cube Data Processing

% Generate 3-D Cartesian Mesh Grid
generateCoordinates3D(scenario);

% View 3-D Scatter Plot of Detections
% viewDetectionsFull(scenario);

% View 3-D Scatter Plot of CFAR Detections
% viewDetectionsCFARFull(scenario);

% Perform HAT algorithm
scenario.terrain = HAT_CAMWAR(scenario);

% View Estimated Terrain
% viewTerrain(scenario);

% Separate Non-Ground Detections
scenario.detection = Declutter_CAMWAR(scenario);

% View 3-D Scatter Plot of Over-Terrain Detections
% viewDetectionsDiscreteFull(scenario);

% View 3-D Scatter Plot of Centroids of Detection Regions
% viewDetectionsDiscreteCentroid(scenario);

% Generate RCS Cube
scenario.detection = GenerateRCSCube(scenario);

%% Data Visualization

% Display result visualization plots
Visualization_CAMWAR(scenario);








