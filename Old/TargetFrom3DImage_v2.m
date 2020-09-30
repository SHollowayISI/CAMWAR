%% Image to Target System
%{

    Sean Holloway
    Version 2
    Separation of a discrete objects into targets for CAMWAR project.

    Uses STL 3D model file and produces a list of targets for use in 
    RadarTransceiverModel_v2.m

    Version 2 changes output format to work with new RadarTransceiverModel
    simulation loop.

    TODO: Correct RCS models

    Working as of 1/16/2020
%}

%% Housekeeping

tic;
clear variables;
close all;

addpath(genpath('Functions'));
addpath(genpath('Scene Assets'));
addpath(genpath('MAT Files'));

%% Read STL File

% Import 3D image from file
fv = stlread('sphere.stl');
% fv = stlread('3D Images/femur.stl');

% Scale if needed
% fv.vertices = fv.vertices*0.25;


%% Convert STL file to 3D binary image

% Creates 3D binary image
% imageSize does not scale automatically, must be set
imageSize = [50, 50, 50];
image_raw = createBinImageFrom3DpointMesh(fv.vertices', imageSize, fv.faces);

image_raw(end) = false;

%% Plot 3D binary image

% [x y z] = ind2sub(size(bin_image), find(bin_image));
% plot3(x, y, z, 'k.')

% volumeViewer(image_raw)


%% Set up 3D model parameters

% Set image size and position parameters
image_bin_size = 2;         % Size in meters of pixel
image_pos_x = 0;          % Location of (0,0) center of image
image_pos_y = 500;          % In this case, y is downrange, x crossrange
image_pos_z = 0;

% Translate image to real location
image_x = (size(image_raw,1)/2)*linspace(-1,1,size(image_raw,1))*image_bin_size + image_pos_x;
image_y = (size(image_raw,2)/2)*linspace(-1,1,size(image_raw,2))*image_bin_size + image_pos_y;
image_z = (size(image_raw,3)/2)*linspace(-1,1,size(image_raw,3))*image_bin_size + image_pos_z;

% Create mesh grid of image coordinates
[xcoord, ycoord, zcoord] = meshgrid(image_x, image_y, image_z);


%% Slice 3D binary image and collect data

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

elevation_axis = 90:-2:-90;
N_e = length(elevation_axis);


%% Translate image pixels to radar bins

radar_bins = false(N_r, N_a+1, N_e);

% Loop through all pixels in image
for pixel = 1:numel(image_raw)
    if image_raw(pixel)
        
        % Determine range and angle to pixel
        [rng, ang] = rangeangle([ycoord(pixel); xcoord(pixel); zcoord(pixel)]);
        
        % If the pixel is in the field of view
        if((rng>0) && (rng<range_axis(end)) && ...
                (ang(1) > -90) && (ang(1) < 90) && ...
                (ang(2) < elevation_axis(1)) && (ang(2) > elevation_axis(end)))
            
            % Determine minimum distance radar bin
            [~, r_idx] = min(abs(range_axis - rng));
            [~, az_idx] = min(abs(azimuth_axis - ang(1)));
            [~, el_idx] = min(abs(elevation_axis - ang(2)));
            
            % Mark location as "true"
            radar_bins(r_idx, az_idx, el_idx) = true;
            
        end
        
    end
end


%% Display slices

% Create 2-D axis
x_axis = range_axis'.*sind(azimuth_axis);
y_axis = range_axis'.*cosd(azimuth_axis);

% Display target locations in Range-Azimuth domain
figure('Name', 'Range-Azimuth Target Locations')
subplot(121)
surf(x_axis, y_axis, double(radar_bins(:,:,ceil(end/2))), 'EdgeColor', 'none')
% surf(x_axis, y_axis, double(radar_bins(:,:,13)), 'EdgeColor', 'none')
set(gca,'YDir','normal')
xlabel('Cross Range [m]')
ylabel('Range [m]')
xlim([-2000 2000])
ylim([0 2000])
view(2)


%% Implement shadowing

% Remove targets which are obscured by closer targets in angle-angle slice

% Set shadowing depth, in # of bins
% Shadowing depth of 0 indicates no shadowing
shadow_depth = 1;

if shadow_depth > 0
    for az_slice = 1:length(azimuth_axis)
        for el_slice = 1:length(elevation_axis)
            
            % Find index of closest target
            sh_idx = find(radar_bins(:,az_slice,el_slice),1);
            
            % Set all further bins to false
            radar_bins((sh_idx + shadow_depth):end, az_slice, el_slice) = 0;
        end
    end
end


%% Plot Range-Azimuth Target Locations

% Plot location of targets in Range-Azimuth domain, after shadowing
%
subplot(122)
surf(x_axis, y_axis, double(radar_bins(:,:,ceil(end/2))), 'EdgeColor', 'none')
% surf(x_axis, y_axis, double(radar_bins(:,:,8)), 'EdgeColor', 'none')
set(gca,'YDir','normal')
xlabel('Cross Range [m]')
ylabel('Range [m]')
xlim([-2000 2000])
ylim([0 2000])
view(2)
%}


%% Display 3D Scatter Plot of Detections

% Create 3-D non-Cartesian axes
elevation_axis = permute(elevation_axis, [1 3 2]);
x3_axis = range_axis'.*sind(azimuth_axis).*cosd(elevation_axis);
y3_axis = range_axis'.*cosd(azimuth_axis).*cosd(elevation_axis);
z3_axis = range_axis'.*ones(1,length(azimuth_axis)).*sind(elevation_axis);

% Create list of target coordinates 
target_x = x3_axis(radar_bins == 1);
target_y = y3_axis(radar_bins == 1);
target_z = z3_axis(radar_bins == 1);

% Set list of original targets for export to detection and processing
tgt_exp.x = target_x;
tgt_exp.y = target_y;
tgt_exp.z = target_z;

% Show 3-D scatter plot of target locations
%
figure('Name', '3D Detection Scatter Plot')
scatter3(target_x, target_y, target_z)
xlim([-500 500])
ylim([0 1000])
zlim([-500 500])
%}

%%  Distribute RCS and make target list

% Set total RCS of target, distribute to bins
% TODO: Correct methodology
total_rcs = 10000;
rcs_per_bin = total_rcs/nnz(radar_bins);

% Remove elevation slices which contain no targets
elevation_axis = elevation_axis(any(radar_bins, [1 2]));
radar_bins = radar_bins(:,:,any(radar_bins,[1 2]));

% Collect indices
[r_idx, az_idx, el_idx] = ind2sub(size(radar_bins), find(radar_bins));

% Loop through targets
for target = 1:nnz(radar_bins)
    
    % Set position of individual target
    target_pos(:,target) = ...
        [range_axis(r_idx(target))*...
        cosd(azimuth_axis(az_idx(target)))*...
        cosd(elevation_axis(el_idx(target))); ...
        range_axis(r_idx(target))*...
        sind(azimuth_axis(az_idx(target)))*...
        cosd(elevation_axis(el_idx(target))); ...
        range_axis(r_idx(target))*...
        sind(elevation_axis(el_idx(target)))];
    
    % Set velocity of individual target
    target_vel(:,target) = [0;0;0];
    
    % Set RCS of individual target
    target_rcs(:,target) = rcs_per_bin;
    
end


%% Export Target List

message = sprintf('%d targets saved', nnz(radar_bins));
disp(message)

filename = 'MAT Files/ImageTargets3D.mat';
save(filename, 'target_pos', 'target_vel', 'target_rcs', ...
    'elevation_axis', 'tgt_exp', '-v7.3');
toc;
disp(['Target list saved as "', filename, '"'])





