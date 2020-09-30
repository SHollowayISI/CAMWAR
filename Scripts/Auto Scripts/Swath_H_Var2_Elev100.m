%% CAMWAR Radar System - Example Automated Testing Initialization File
%{

    Sean Holloway
    CAMWAR Automated Testing Init File
    
    This file specifies simulation, environment, and radar scenario
    parameters for automated script testing.

    Use script 'FullSystem_Automated_CAMWAR.m' to run scenarios.
    
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
    'elevation_slices', -30:1:0, ...          % Elevation angles to simulate
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
    'send_alert',   true, ...                  % Send email alert T/F
    'attach_zip',   false, ...
    'alert_address', 'sholloway@intellisenseinc.com', ...
    ...                                         % Email address for status
    ...                                            updates
    'filename',     'RCSTest_CAMWAR', ...     % Filename to save data as
    'save_format',  save_format, ...            % File types to save figures
    'save_figs',    true, ...                   % Save figures T/F
    'save_mat',     true, ...                   % Save mat file T/F
    'reduce_mat',   true);                      % Reduce mat file for saving

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
    swath_index = (abs(scenario.target_list.pos(1,:) - 400) < ...
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

%% Radar Parameter Setup

% Radar simulation and processing setup
scenario.radarsetup = struct( ...
    ...
    ... % Waveform Properties
    'f_c',          77e9, ...               % Operating frequency in Hz
    'f_s',          37.5e6, ...             % ADC sample frequency in Hz
    't_ch',         31.25e-6, ...              % Chirp duration in seconds
    'bw',           88e6, ...               % Chirp bandwidth in Hz
    'n_p',          16, ...                  % Number of (MIMO) chirps per frame
    ...
    ... % Transceiver Properties
    'n_tx_ant',     24, ...             % Number of vertical elements in Tx subarray
    'n_tx_sub',     2, ...              % Number of horizontal Tx subarrays
    'n_rx_ant',     32, ...             % Number of horizontal elements in Rx array
    'tx_pow',       0.76, ...           % Transmit power per antenna in Watts
    'rf_sys_loss',  4, ...              % RF system loss in dB
    'rx_nf',        4, ...              % Rx noise figure in dB
    ...
    ... % Antenna Properties
    'tx_gain',      100/48, ...         % Transmit antenna gain in absolute
    'tx_az_bw',     74.61, ...          % Tx azimuth beamwidth in degrees
    'tx_el_bw',     57.86, ...          % Tx elevation beamwidth in degrees
    'tx_az_pow',    4, ...              % Tx azimuth sinc power
    'tx_el_pow',    4, ...              % Tx elevation sinc power
    ...
    'rx_gain',      640/32, ...         % Receive antenna gain in absolute
    'rx_az_bw',     69.17, ...          % Rx azimuth beamwidth in degrees
    'rx_el_bw',     26.09, ...          % Rx elevation beamwidth in degrees
    'rx_az_pow',    2, ...              % Rx azimuth sinc power
    'rx_el_pow',    4, ...              % Rx elevation sinc power
    'rx_el_tilt',   0, ...-10, ...            % Rx elevation tilt in degrees
    ...
    ... % Processing Properties
    'r_win',        'hanning', ...          % Window for range processing
    'd_win',        'blackmanharris', ...   % Window for range processing
    'a_win',        'hanning', ...          % Window for range processing
    ...
    ... % Detection Properties
    'thresh',       14.4, ...           % Detection threshold above noise power in dB
    'CFAR_Pfa',     1e-4, ...           % False alarm probability for CFAR detection
    'int_start',    460, ...            % Bin to start performing range bin combination
    'range_int',    false, ...           % Range combining Y/N
    'num_guard',    [10 2], ...         % Number of R-D guard cells for CFAR detection
    'num_train',    [10 2]);            % Number of R-D training cells for CFAR detection


%% Run Setup Scripts

% Set up Phased Array Toolbox system objects
scenario = PhasedSetup_CAMWAR(scenario);



