function [cube] = SignalProcessing_CAMWAR(scenario)
%SIGNALPROCESSING_CAMWAR Performs signal processing for CAMWAR
%   Takes scenario struct as input, retuns scenario.cube struct containing
%   processed Range-Doppler cube

%% Unpack Variables

radarsetup = scenario.radarsetup;
simsetup = scenario.simsetup;

%% Define Constants

c = physconst('LightSpeed');
lambda = c/radarsetup.f_c;

%% Perform Range FFT

% Calculate FFT Size
N_r = 2^ceil(log2(size(scenario.rx_sig,1)));

% Apply windowing
expression = '(size(scenario.rx_sig,1)).*scenario.rx_sig;';
expression = [ radarsetup.r_win, expression];
cube.range_cube = eval(expression);

% FFT across fast time dimension
cube.range_cube = fft(scenario.rx_sig, N_r, 1);

% Remove negative complex frequencies
cube.range_cube = cube.range_cube(1:ceil(end/2),:,:);


%% Perform Doppler FFT

% Calculate FFT size
N_d = 2^ceil(log2(size(cube.range_cube,2)));

% Apply windowing
expression = '(size(cube.range_cube,2))).*cube.range_cube;';
expression = ['transpose(', radarsetup.d_win, expression];
cube.rd_cube = eval(expression);

% FFT across slow time dimension
cube.rd_cube = fftshift(fft(cube.rd_cube, N_d, 2), 2);

% Wrap max negative frequency and positive frequency
cube.rd_cube(:,(end+1),:) = cube.rd_cube(:,1,:);

%% Perform Angle FFT

% Calculate FFT size 
%TODO: CHANGE TO SETTING OR RESOLUTION SETTING
N_a = 2^ceil(log2(size(cube.rd_cube, 3)));

% Apply windowing
expression = '(size(cube.rd_cube,3)), [2 3 1]).*cube.rd_cube;';
expression = ['permute(', radarsetup.a_win, expression];
cube.angle_cube = eval(expression);

% FFT across angle dimension
cube.angle_cube = fftshift(fft(cube.angle_cube, N_a, 3), 3);

% Wrap max negative frequency and positive frequency
cube.angle_cube(:,:,(end+1)) = cube.angle_cube(:,:,1);

%% Calculate Power Cube

% Take square magnitude of radar cube
cube.pow_cube = abs(cube.angle_cube).^2;

% Sum over doppler axis
cube.pow_cube_nodop = squeeze(cube.pow_cube(:,ceil(end/2),:));

%% Derive Axes

% Derive Range axis
cube.range_res = (size(scenario.rx_sig,1)/N_r)*(c/(2*radarsetup.bw));
cube.range_axis = ((1:(N_r/2))-0.5)*cube.range_res;

% Derive Doppler axis
cube.vel_res = lambda/(2*radarsetup.t_ch*radarsetup.n_p*radarsetup.n_tx_sub);
cube.vel_axis = ((-N_d/2):(N_d/2))*cube.vel_res;

% Derive Azimuth axis
cube.azimuth_axis = -asind(((-N_a/2):(N_a/2))*(2/N_a));
cube.azimuth_res = min(abs(diff(cube.azimuth_axis)));

%% Trim Cube to Field of View

% Determine Azimuth bins in FOV
[~, az_idx] = find(cube.azimuth_axis < simsetup.fov/2, 1);

% Reduce azimuth axis
cube.azimuth_axis = cube.azimuth_axis(az_idx:(end-az_idx+1));

% Reduce power cubes
cube.pow_cube = cube.pow_cube(:,:,az_idx:(end-az_idx+1));
cube.pow_cube_nodop = cube.pow_cube_nodop(:,az_idx:(end-az_idx+1));


end