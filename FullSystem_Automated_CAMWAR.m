%% CAMWAR Radar System
%{

    Sean Holloway
    CAMWAR (Collision Avoidance Millimeter Wave Radar) System
    MATLAB Simulation & Processing

    This shell file runs successive scripts and gauges progress.

    TODO's logged in FullSystem_CAMWAR.m

%}

%% Housekeeping
clear variables
close all
addpath(genpath(pwd));
tic

%% Loop Through Simulation Setup Files

% Pull list of files from directory
dir_list = dir('Automated Testing/To Run/*.m');

for file_index = 1:length(dir_list)
    
    % Clear scenario object
    clear scenario
    
    % Initialization of scenario object
    scenario = RadarScenario_CAMWAR;
    
    % Pull filename from directory
    scenario_filename = dir_list(file_index).name;

    % Display current test
    disp(['Beginning Scenario: ', scenario_filename]);
    
    % Run setup files
    run(['Automated Testing/To Run/', scenario_filename]);
    
    % Perform main processes of simulation, signal and data processing
    Main_CAMWAR

    % Run all end-of-simulation tasks
    EndSimulationAuto_CAMWAR

end























