function [angle_bins, k_a] = Angle_Calc(doppler_bins,N)
%ANGLE_CALC Applies Hanning Window and performs N-point Angle 2D-FFT 
%along dimensions 3 and 4 of "doppler_bins" input
%   "angle_bins" is FFT output
%   "k_a" is bin index

% Apply windowing
h = hanning(size(doppler_bins,3));
doppler_bins = permute(h, [3 2 1]).*doppler_bins;

% Set up axis
k_a = -N/2:N/2;

% Perform FFT and shift
angle_bins = fft(doppler_bins, N, 3);
angle_bins = fftshift(angle_bins, 3);

% Wrap around FFT shift
angle_bins(:,:,(end+1)) = angle_bins(:,:,1);