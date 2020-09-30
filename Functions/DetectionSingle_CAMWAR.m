function [detection] = DetectionSingle_CAMWAR(scenario)
%DETECTION_CAMWAR Performs single-frame target detection for CAMWAR project
%   Takes scenario object as input, provides scenario.detection object as
%   output, containing information about detected targets.

%% Unpack Variables

detection = scenario.detection;
radarsetup = scenario.radarsetup;
cube = scenario.cube;

%% Perform Threshold Detection

% Estimate noise power
detection.noise_pow = pow2db(median(median(mean(scenario.cube.pow_cube))));

% Calculate threshold in absolute
abs_thresh = db2pow(radarsetup.thresh + detection.noise_pow);

% Calculate start bin
[~, start_bin] = find(scenario.cube.range_axis > scenario.radarsetup.int_start, 1);

% Combine range bins in cube
if radarsetup.range_int
    comb_cube = cube.pow_cube;
    comb_cube(start_bin:end) = ...
        comb_cube(start_bin:end) + ...
        comb_cube((start_bin-1):(end-1));
else
    comb_cube = cube.pow_cube;
end

% Perform detection
detection.detect_cube(:,:,:,scenario.flags.slice) = ...
    (comb_cube > abs_thresh);

% % Generate zero doppler detection cube
% detection.detect_cube_nodop(:,:,scenario.flags.slice) = ...
%     squeeze(any(detection.detect_cube(:,:,:,scenario.flags.slice), 2));

detection.detect_cube_nodop(:,:,scenario.flags.slice) = ...
    squeeze(detection.detect_cube(:,ceil(end/2),:,scenario.flags.slice));

% Generate cube of SNR
detection.SNR_cube(:,:,scenario.flags.slice) = detection.detect_cube_nodop(:,:,scenario.flags.slice) .* ...
    cube.pow_cube_nodop / db2pow(detection.noise_pow);

%% Read SNR and Power Values

% Determine peak power in cube
detection.max_SNR = 10*log10(max(cube.pow_cube, [], 'all', 'linear')) ...
    - detection.noise_pow;

% Determine total power in cube
detection.total_SNR = 10*log10(sum(cube.pow_cube, 'all')) ...
    - detection.noise_pow;

% Determine total power in detection bins
detection.detect_SNR = 10*log10(sum(cube.pow_cube(detection.detect_cube(:,:,:,scenario.flags.slice)))) ...
    - detection.noise_pow;

%% Perform CFAR Detection

% Set up detectionwindow sizes
[r_sz, d_sz, a_sz] = size(scenario.cube.pow_cube);

tr_x = radarsetup.num_guard(1) + radarsetup.num_train(1);
tr_y = radarsetup.num_guard(2) + radarsetup.num_train(2);

cut_range = ((1+tr_x):(r_sz-tr_x))';
cut_dop = ((1+tr_y):(d_sz-tr_y))';

cut_range2D = reshape(repmat(cut_range, [1 length(cut_dop)])', ...
    length(cut_range)*length(cut_dop), []);
cut_dop2D = reshape(repmat(cut_dop, [1 length(cut_range)]), ...
    length(cut_range)*length(cut_dop), []);

cut2d_RD = [cut_range2D, cut_dop2D]';

% Loop through azimuth angles
for az = 1:a_sz
    % Run 2-D CFAR Detection
    result_2D = scenario.sim.CFAR(scenario.cube.pow_cube(:,:,az), cut2d_RD);
    
    % Reshape to original data cube size
    result_2D_shaped = reshape(result_2D, length(cut_range), []);
    result_2D_shaped = [zeros(length(result_2D_shaped), tr_y), ...
        result_2D_shaped, zeros(length(result_2D_shaped), tr_y)];
    result_2D_shaped = [zeros(tr_x, size(result_2D_shaped,2)); ...
        result_2D_shaped; zeros(tr_x, size(result_2D_shaped,2))];
    
    % Save CFAR detections per slice
    detection.CFAR_cube(:,:,az, scenario.flags.slice) = result_2D_shaped;

end

% % Generate zero doppler CFAR detection cube
% detection.CFAR_cube_nodop(:,:,scenario.flags.slice) = ...
%     squeeze(any(detection.detect_cube(:,:,:,scenario.flags.slice), 2));

detection.CFAR_cube_nodop(:,:,scenario.flags.slice) = ...
    squeeze(detection.detect_cube(:,ceil(end/2),:,scenario.flags.slice));

end

