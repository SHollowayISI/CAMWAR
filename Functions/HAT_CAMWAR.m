function [terrain] = HAT_CAMWAR(scenario)
%HAT_CAMWAR Perform Height-Above-Terrain (HAT) algorithm for CAMWAR
%   Takes radar scenario object as input, returns 

%% HAT Algorithm

% Convert detection cube to cube of ranges
ranges = scenario.detection.detect_cube_nodop.*scenario.cube.range_axis';

% Remove outliers from data set
ranges(ranges == 0) = nan;
keep_logical = ~isoutlier(ranges, 1);
range_sum = sum(ranges .* keep_logical .* ~isnan(ranges), 1, 'omitnan');
num_keep = sum(keep_logical .* ~isnan(ranges), 1);

% Estimate average ranges for each angle bin
terrain.range_cube = squeeze(range_sum./num_keep);

% Reshape if single range slice
if length(scenario.simsetup.elevation_slices) == 1
    terrain.range_cube = terrain.range_cube';
end

% Remove returns over the maximum angle
% terrain.range_cube(:,(max_idx+1):end) = nan;

% Estimate terrain location from range cube
terrain.locations(1,:,:) = terrain.range_cube .* ...
    cosd(scenario.cube.azimuth_axis') .* ...
    cosd(scenario.simsetup.elevation_slices);
terrain.locations(2,:,:) = terrain.range_cube .* ...
    sind(scenario.cube.azimuth_axis') .* ...
    cosd(scenario.simsetup.elevation_slices);
terrain.locations(3,:,:) = terrain.range_cube .* ...
    ones(size(scenario.cube.azimuth_axis')) .* ...
    sind(scenario.simsetup.elevation_slices);

% Save height cube of Z-values
terrain.height_cube = squeeze(terrain.locations(3,:,:));


end

