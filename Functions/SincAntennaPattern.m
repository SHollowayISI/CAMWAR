function [pattern] = SincAntennaPattern(az, el, az_tilt, el_tilt, ...
    az_bw, el_bw, az_pow, el_pow, gain)
%SINCANTENNAPATTERN Returns antenna pattern that is sinc function in each
%direction.
%   az:         Azimuth angle(s) in degrees
%   el:         Elevation angle(s) in degrees
%   az_tilt:    Azimuth angle of main beam in degrees
%   el_tilt:    Elevation angle of main beam in degrees
%   az_bw:      Azimuthal beamwidth in degrees
%   el_bw:      Elevation beamwidth in degrees
%   az_pow:     Power to raise azimuth sinc function to. Must be even number.
%               Each power of 2 decreases first sidelobe by 13dB.
%   el_pow:     Power to raise elevation sinc function to. See az_pow.
%   gain:       Maximum gain at main beam direction in absolute

phi = (az-az_tilt)/(0.5*pi*az_bw);
theta = ((el-el_tilt)/(0.5*pi*el_bw))';

pattern = gain*(sinc(phi).^az_pow).*(sinc(theta).^el_pow);

end

