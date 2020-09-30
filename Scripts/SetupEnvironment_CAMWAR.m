%% CAMWAR Radar System - Environment Initialization File
%{

    Sean Holloway
    CAMWAR Environment Init File
    
    This file specifies the environmental scene model used in CAMWAR
    simulation, including path losses due to low visibility conditions, and
    creates a list of targets based on terrain, target, and object
    locations and types.

    Use script 'FullSystem_CAMWAR.m' to run scenarios.
    
%}

%% Environment Parameter Setup

% Environment parameter setup
scenario.environment = struct( ...
    ...
    ... % Channel Properties
    'channel_loss_on',      false, ...      % Whether to include environment
    ...                                     % modeling in LOS channel loss
    ... % Terrain Properties
    'gnd_on',               true, ...       % Simulate terrain T/F
    'gnd_res',              2, ...         % Distance between targets
    'gnd_xlim',             [0 1500], ...  % Down-range min-max limits
    'gnd_ylim',             [-1000 1000], ...   % Cross-range min-max limits
    'gnd_zavg',             -100, ...       % Average ground elevation
    'gnd_var',              2, ...         % Ground height variance
    'gnd_cor_l',            500, ...         % Ground correlation length
    'gnd_gamma',            db2pow(-10), ...  % Surface gamma
    ...
    ... % Swath Properties
    'swath_on',             false, ...           % Simulate different gamma swath T/F
    'swath_diff',           db2pow(-20), ...    % Difference between swath RCS and ground RCS in dB
    'swath_width',          30, ...             % Width of swath
    ...
    ... % Discrete Target Properties (DEBUG)
    'pt_tgt',               false, ...          % Single point target T/F
    'tgt_pos',              [42.6434; 7.5192; -15], ...    % Discrete target location
    'tgt_vel',              [10; 0; 0], ...      % Discrete target velocity
    'tgt_rcs',              db2pow(0));         % Discrete target RCS


%% Generate Target List from Terrain Parameters

if scenario.environment.gnd_on
    % Generate elevation data
    scenario = TerrainToTarget(scenario);
end

if scenario.environment.swath_on && scenario.environment.gnd_on
    % Determine targets in swath
    swath_index = (abs(scenario.target_list.pos(2,:)) < ...
        scenario.environment.swath_width/2);
    
    % Modify RCS of targets in swath
    scenario.target_list.rcs(swath_index) = ...
        scenario.target_list.rcs(swath_index) ...
        * scenario.environment.swath_diff;
end

%% Generate Target List for Discrete Targets

if scenario.environment.pt_tgt
    % Generate single point target
    scenario.target_list.pos = [scenario.target_list.pos, scenario.environment.tgt_pos];
    scenario.target_list.vel = [scenario.target_list.vel, scenario.environment.tgt_vel];
    scenario.target_list.rcs = [scenario.target_list.rcs, scenario.environment.tgt_rcs];
end

%% Display scatter plot of targets

% Visualize target locations
% viewTargets(scenario)

% Visualize terrain RCS
% viewTerrainRCS(scenario)

% Visualize generated terrain
% viewRealTerrain(scenario)




