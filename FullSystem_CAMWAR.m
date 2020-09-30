%% CAMWAR Radar System
%{

    Sean Holloway
    CAMWAR (Collision Avoidance Millimeter Wave Radar) System
    MATLAB Simulation & Processing

    This shell file runs successive scripts and gauges progress.

    TODO: DetectionSingle
            - Range integration

    TODO: SetupEnvironment:
            - Object generation
            - Channel effects

    TODO: Detection:
            - m-of-n Processing

    TODO: DataProcessing:
            - Monopulse Processing
            - DBS

%}

%% Housekeeping
clear variables
close all
addpath(genpath(pwd));
tic

%% Initialize Scenario Object

% Initialization
scenario = RadarScenario_CAMWAR;

%% Setup Structures for Simulation

% Set up simulation parameters
SetupSimulation_CAMWAR

% Set up radar environment, returning target list
SetupEnvironment_CAMWAR

% Set up transceiver parameters
SetupRadarScenario_CAMWAR

%% Run Simulation & Signal Processing

% Perform main processes of simulation, signal and data processing
Main_CAMWAR

%% Save and Package Resultant Data

% Run all end-of-simulation tasks
EndSimulationSingle_CAMWAR

















