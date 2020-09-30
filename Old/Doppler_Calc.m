function [doppler_bins, k_d] = Doppler_Calc(slow_bins,N)
%DOPPLER_CALC Applies Blackman Harris Window and performs N-point Range FFT 
%along dimension 2 of "slow_bins" input.
%Note that this is Doppler_Calc for TRASAT/SEMTA projects.
%   "doppler_bins" is FFT output
%   "k_d" is bin index

% Apply windowing
h = blackmanharris(size(slow_bins,2));
slow_bins = h'.*slow_bins;

% Create axis
k_d = (-N/2:N/2);

% Perform FFT and shift
doppler_bins = fft(slow_bins, N, 2);
doppler_bins = fftshift(doppler_bins, 2);

% Wrap end of FFT shfit
doppler_bins(:,(end+1),:) = doppler_bins(:,1,:);
