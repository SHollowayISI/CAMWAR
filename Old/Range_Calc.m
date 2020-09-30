function [range_bins, k_r] = Range_Calc(chirps, N)
%RANGE_CALC Applies Hanning Window and performs N-point Range FFT along
%dimension 1 of "chirps" input
%   "range_bins" is FFT output
%   "k_r" is bin index

h = hanning(size(chirps,1));
chirps_hanning = h.*chirps;


k_r = 1:N/2;
% k_r = 0:N-1;
range_bins_complex = fft(chirps_hanning,N,1);
% range_bins = range_bins_complex;
range_bins = range_bins_complex(1:end/2,:,:,:);
%range_bins = flip(range_bins_complex(N/2+1:end,:,:),1);
%range_bins = fftshift(range_bins,1);