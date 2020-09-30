%% CAMWAR Radar System - Example Radar Initialization File
%{

    Sean Holloway
    CAMWAR Init File
    
    This file specifies radar parameters for CAMWAR simulation.

    Use script 'FullSystem_CAMWAR.m' to run scenarios.
    
%}

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





