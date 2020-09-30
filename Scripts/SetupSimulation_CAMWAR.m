%% CAMWAR Radar System - Example Simulation Initialization File
%{

    Sean Holloway
    CAMWAR Simulation Init File
    
    This file specifies simulation parameters for CAMWAR simulation.

    Use script 'FullSystem_CAMWAR.m' to run scenarios.
    
%}

%% Simulation Parameter Setup

save_format.list = {'.png','.fig'};

% Radar simulation and processing setup
scenario.simsetup = struct( ...
    ...
    ... % Transceiver Trajectory Properties
    'radar_pos',    [0; 0; 0], ...              % Position of radar unit
    'radar_vel',    [0; 0; 0], ...              % Velocity of radar unit
    ...
    ... % Radar Scan Settings
    'elevation_slices', -30:2:-10, ...          % Elevation angles to simulate
    'allocate_slice',   true, ...               % Only simulate targets in slice Y/N
    'slice_bw',         8, ...                  % Beamwidth for slice allocation
    'fov',              90, ...                 % Field of view in Azimuth direction
    'allocate_az',      true, ...               % Simulate targets in FOV
    'fov_bw',           5, ...                  % Extra beamwidth for azimuth allocation
    'HAT_max',          -10, ...                % Highest elevation for HAT algorithm
    ...
    ... % Visualization Settings
    'max_range_method', 'fit', ...              % 'max', 'set', 'fit', or 'none'
    'set_max',          500, ...                % Max range setting for if 'set' is used
    'min_range_method', 'fit', ...              % 'min', 'set', 'fit', or 'none'
    'set_min',          100, ...                % Min range setting for if 'set' is used
    'contour_scale',    5, ...                  % Scale for contour lines
    'interpolate_nan',  true, ...               % Interpolate through NaN values
    ...
    ... % Simulation Properties
    'sim_rate',     2^0, ...                    % Rate to divide up fast x 
    ...                                           slow time simulation
    'send_alert',   false, ...                  % Send email alert T/F
    'attach_zip',   false, ...
    'alert_address', 'sholloway@intellisenseinc.com', ...
    ...                                         % Email address for status
    ...                                            updates
    'filename',     'RCSTest_CAMWAR', ...     % Filename to save data as
    'save_format',  save_format, ...            % File types to save figures
    'save_figs',    false, ...                   % Save figures T/F
    'save_mat',     false, ...                   % Save mat file T/F
    'reduce_mat',   false);                      % Reduce mat file for saving




