function [az_out, el_out, pattern] = rxPatternCAMWAR(tilt)
%RXPATTERNCAMWAR Returns Rx antenna pattern for CAMWAR
%   Input tilt is the angle above boresight at which the antenna beam is
%   directed.
%   
%   Outputs az_out, el_out are axes.
%   
%   Pattern can be used in phased.CustomAntennaElement object in the Phased
%   Array Toolbox

az_out = -180:0.1:180;
el_out = -90:0.1:90;

Gr = 640;
del_theta = 30*pi/180;
del_phi = 2*pi/180;

tilt_rad = tilt*pi/180;

phi = az_out*pi/180;
theta = (el_out*pi/180)';

pattern = Gr*(sinc(phi/(0.5*pi*del_phi)).^4).*(sinc((theta-tilt_rad)/(0.5*pi*del_theta)).^4);

end