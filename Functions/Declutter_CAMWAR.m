function [detection] = Declutter_CAMWAR(scenario)
%DECLUTTER_CAMWAR Performs clutter rejection for CAMWAR project
%   Takes scenario object as input, provides scenario.detection object as
%   output, containing target list of non-terrain detections

%% Unpack Variables

detection = scenario.detection;
radarsetup = scenario.radarsetup;
simsetup = scenario.simsetup;
terrain = scenario.terrain;
res = scenario.cube.range_res;
raxis = scenario.cube.range_axis;
results = scenario.results;

ff = 5;

%% Calculate Non-Terrain Detections

% Calculate array beamwidth
el_bw = rad2deg(1/radarsetup.n_tx_ant);

% Calculate height-to-range coefficient
trig_factor = 1 - ff*sind(el_bw)./tand(-simsetup.elevation_slices);
         
% Calculate maximum non-terrain range
max_range = permute(terrain.range_cube.*trig_factor, [3 1 2]);
% max_bin = round(max_range/res + 0.5);
ranges = repmat(raxis', size(max_range));
keep_logical = (ranges < max_range);

% Remove terrain detections from detection cube
detection.discrete_cube_nodop = detection.detect_cube_nodop;
detection.discrete_cube_nodop(~keep_logical) = 0;

detection.discrete_cube = detection.detect_cube;
detection.discrete_cube( ...
    repmat(permute(~keep_logical, [1 4 2 3]), 1, size(detection.CFAR_cube, 2), 1, 1)) ...
    = 0;

detection.discrete_cube_any = squeeze(any(detection.discrete_cube, 2));


%% Find Centroids of Detection Regions

% Find connected components and pull region data
components = bwconncomp(detection.discrete_cube_nodop);
centroid = regionprops(components, 'Centroid');



% Create list of Cartesian coordinates of detections
detection.detect_list = [];
detection.centroid_angles = [];
for n = 1:length(centroid)
    
    cent_loc = centroid(n).Centroid;
    cent_x = interp3(results.x_grid, cent_loc(1), cent_loc(2), cent_loc(3)); 
    cent_y = interp3(results.y_grid, cent_loc(1), cent_loc(2), cent_loc(3)); 
    cent_z = interp3(results.z_grid, cent_loc(1), cent_loc(2), cent_loc(3)); 
    
    cent_az = interp1(scenario.cube.azimuth_axis, cent_loc(1));
    cent_el = interp1(scenario.simsetup.elevation_slices, cent_loc(3));

    detection.centroid_angles(:,(end+1)) = [cent_az; cent_el];
    detection.detect_list(:,(end+1)) = [cent_x; cent_y; cent_z];
end
   
end

