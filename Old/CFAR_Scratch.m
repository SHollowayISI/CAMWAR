%% Bookkeeping

addpath(genpath('MAT Files'));
load('processedCube.mat');


%% Set up CFAR

power_bins = abs(angle_bins).^2;
sum_dop_bins = squeeze(sum(power_bins,2));

center1D = power_bins(:,ceil(end/2),ceil(end/2),ceil(end/2));
center2D_RD = power_bins(:,:,ceil(end/2),ceil(end/2));
center2D_RA = squeeze(power_bins(:,ceil(end/2),:,ceil(end/2)));

centersum1D = sum_dop_bins(:,ceil(end/2),ceil(end/2));
centersum2D = sum_dop_bins(:,:,ceil(end/2));

center3D_RD_A = power_bins(:,:,:,ceil(end/2));


Pfa = 1e-4;

numGuard1D = 4;
numTraining1D = 6;

detector1D = phased.CFARDetector('Method', 'CA', 'ProbabilityFalseAlarm', Pfa, ...
    'NumGuardCells', numGuard1D, 'NumTrainingCells', numTraining1D);


numGuard2D = [1 2];
numTraining2D = [1 2];

detector2D = phased.CFARDetector2D('Method', 'CA', 'ProbabilityFalseAlarm', Pfa, ...
    'GuardBandSize', numGuard2D, 'TrainingBandSize', numTraining2D);


%% Run CFAR Tests

% close all;
% 1 Dimensional CFAR test
%{

cutCenter = (1:length(center1D))';
cutAzimuth = repmat(cutCenter,[1 size(power_bins,3)]);
cutDoppler = repmat(cutCenter,[1 size(power_bins,2)]);


subplot(211)
plot(range_axis, detector1D(center1D, cutCenter), range_axis, detector1D(centersum1D, cutCenter))
xlim([400 500])

subplot(212)
plot(range_axis, center1D, range_axis, centersum1D)
xlim([400 500])
%}

% 2 Dimensional Range-Doppler CFAR test
%

tr_x = numGuard2D(1) + numTraining2D(1);
tr_y = numGuard2D(2) + numTraining2D(2);

cut_range = ((1+tr_x):(length(center1D)-tr_x))';
cut_dop = ((1+tr_y):(size(center2D_RD,2)-tr_y))';

cut_range2D = reshape(repmat(cut_range, [1 length(cut_dop)])', ...
    length(cut_range)*length(cut_dop), []);
cut_dop2D = reshape(repmat(cut_dop, [1 length(cut_range)]), ...
    length(cut_range)*length(cut_dop), []);

cut2d_RD = [cut_range2D, cut_dop2D]';


result_2D = detector2D(center2D_RD, cut2d_RD);
result_2D_shaped = reshape(result_2D, length(cut_range), []);
result_2D_shaped = [zeros(length(result_2D_shaped), tr_y), ...
    result_2D_shaped, zeros(length(result_2D_shaped), tr_y)];
result_2D_shaped = [zeros(tr_x, size(result_2D_shaped,2)); ...
    result_2D_shaped; zeros(tr_x, size(result_2D_shaped,2))];  

subplot(121)
imagesc(doppler_axis, range_axis, result_2D_shaped)
set(gca,'YDir','normal')

subplot(122)
imagesc(doppler_axis, range_axis, center2D_RD)
set(gca, 'YDir','normal')
%}

% 2 Dimensional Azimuth-Doppler CFAR test
%{
tr_x = numGuard2D(1) + numTraining2D(1);
tr_y = numGuard2D(2) + numTraining2D(2);

cut_range = ((1+tr_x):(length(center1D)-tr_x))';
cut_az = ((1+tr_y):(size(centersum2D,2)-tr_y))';

cut_range2D = reshape(repmat(cut_range, [1 length(cut_az)])', ...
    length(cut_range)*length(cut_az), []);
cut_az2D = reshape(repmat(cut_az, [1 length(cut_range)]), ...
    length(cut_range)*length(cut_az), []);

cut2d_RA = [cut_range2D, cut_az2D]';


result_2D = detector2D(centersum2D, cut2d_RA);
result_2D_shaped = reshape(result_2D, length(cut_range), []);
result_2D_shaped = [zeros(length(result_2D_shaped), tr_y), ...
    result_2D_shaped, zeros(length(result_2D_shaped), tr_y)];
result_2D_shaped = [zeros(tr_x, size(result_2D_shaped,2)); ...
    result_2D_shaped; zeros(tr_x, size(result_2D_shaped,2))];  

subplot(121)
imagesc(azimuth_axis, range_axis, result_2D_shaped)
set(gca,'YDir','normal')

subplot(122)
imagesc(azimuth_axis, range_axis, centersum2D)
set(gca, 'YDir','normal')
%}

% 2 Dimensional Range-Doppler CFAR across Azimuth pages
%
tr_x = numGuard2D(1) + numTraining2D(1);
tr_y = numGuard2D(2) + numTraining2D(2);

cut_range = ((1+tr_x):(length(center1D)-tr_x))';
cut_dop = ((1+tr_y):(size(center2D_RD,2)-tr_y))';

cut_range2D = reshape(repmat(cut_range, [1 length(cut_dop)])', ...
    length(cut_range)*length(cut_dop), []);
cut_dop2D = reshape(repmat(cut_dop, [1 length(cut_range)]), ...
    length(cut_range)*length(cut_dop), []);

cut2d_RD = [cut_range2D, cut_dop2D]';

az_num = size(center3D_RD_A,3);

result_3D = detector2D(center3D_RD_A, cut2d_RD);
result_3D_shaped = reshape(result_3D, length(cut_range), [], az_num);
result_3D_shaped = [zeros(length(result_3D_shaped), tr_y, az_num), ...
    result_3D_shaped, zeros(length(result_3D_shaped), tr_y, az_num)];
result_3D_shaped = [zeros(tr_x, size(result_3D_shaped,2), az_num); ...
    result_3D_shaped; zeros(tr_x, size(result_3D_shaped,2), az_num)];  

viewangle = 35;

subplot(121)
imagesc(azimuth_axis, range_axis, result_3D_shaped(:,:,viewangle))
set(gca,'YDir','normal')

subplot(122)
imagesc(azimuth_axis, range_axis, center3D_RD_A(:,:,viewangle))
set(gca, 'YDir','normal')
%}

