%% Terrain to Target System
%{

    Sean Holloway
    2/3/2020
    Version 1
    Separation of a terrain surface into targets for CAMWAR project.

    Uses surface roughness/feature/rcs description and produces a list of
    targets for use in RadarTransceiverModel_v2.m

    TODO: Implement terrain generation
    TODO: Implement correct gamma model

    Working as of 2/3/2020
%}

%% Housekeeping

tic;
clear variables;
close all;

addpath(genpath('Functions'));
addpath(genpath('Scene Assets'));
addpath(genpath('MAT Files'));


%% Set up terrain parameters

% Assigned variables
mean_elev = -100;           % Average height above radar POV


%% Set up axes in each dimension

% Set up range bin limits
range_res = 1.948724168518067;
N_r = 1024;
range_axis = ((1:N_r)-0.5)*range_res;

% Set up azimuth angle bin limits
N_a = 64;
k_a = -N_a/2:N_a/2;
azimuth_axis = asind((k_a)*(2/N_a));

% Set up elevation angle bin limits
N_e = 12;
angle_res = 2;
maximum_angle = 6;
elevation_axis = maximum_angle:-angle_res:(maximum_angle - angle_res*(N_e-1));

% Overwrite with linspace elevation axis
% elevation_axis = 90:-2:-90;
elevation_axis = -5:-10:-15;
N_e = length(elevation_axis);

% Set up 3D axis grid
elevation_axis = permute(elevation_axis, [1 3 2]);
x3_axis = range_axis'.*sind(azimuth_axis).*cosd(elevation_axis);
y3_axis = range_axis'.*cosd(azimuth_axis).*cosd(elevation_axis);
z3_axis = range_axis'.*ones(1,length(azimuth_axis)).*sind(elevation_axis);


%% Set up 3D model parameters

% Set up patch area
% Y is down range, X is cross range
terrain_res = 10;
x_min = -400;
x_max = 400;
y_min = 0;
y_max = 800;

% Create axes
image_x = x_min:terrain_res:x_max;
image_y = y_min:terrain_res:y_max;
image_z = mean_elev;

% Create mesh grid of targets
[xcoord, ycoord, zcoord] = meshgrid(image_x, image_y, image_z);

%%  Distribute RCS and make target list

% Set target position
target_pos = [ycoord(:)'; xcoord(:)'; zcoord(:)'];

% Set target velocity
target_vel = repmat([0;0;0], 1, numel(xcoord));

% Set target RCS
rcs_per_bin = terrain_res^2;
[~, grazing_angle] = rangeangle(target_pos);

target_rcs = rcs_per_bin.*abs(sind(grazing_angle(2,:)));


%% Display scatter plot of targets

% Show 3-D scatter plot of target locations
%{
figure('Name', 'Target 3D Scatter Plot')
scatter3(target_pos(2,:), target_pos(1,:), target_pos(3,:))
xlim([-500 500])
ylim([0 1000])
zlim([-500 500])
%}


%% Export Target List to File

% Export real target list for use in Data Processing comparison
tgt_exp.x = target_pos(2,:);
tgt_exp.y = target_pos(1,:);
tgt_exp.z = target_pos(3,:);

message = sprintf('%d targets saved', size(target_pos, 2));
disp(message)

filename = 'MAT Files/ImageTargets3D.mat';
save(filename, 'target_pos', 'target_vel', 'target_rcs', ...
    'elevation_axis', 'tgt_exp', '-v7.3');
toc;
disp(['Target list saved as "', filename, '"'])





