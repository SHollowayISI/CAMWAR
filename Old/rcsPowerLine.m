function [rcs] = rcsPowerLine(inc_angle)
%RCSPOWERLINE Returns RCS of power line
%   Input is incident angle in degrees, output is RCS in meter^2

K = 25;
N = 10;
M = 200;

phi = inc_angle*pi/180;

rcs = (sin(K*N*phi).^2).*...
    ((sinc(N*(phi+0.12)/pi).^4)+(sinc(N*(phi-0.12)/pi).^4))./...
    (M*sin(K*phi).^2);

rcs(mod(K*inc_angle, 180) == 0) = ((sinc(N*(phi(mod(K*inc_angle, 180) == 0)+0.12)/pi).^4)...
    +(sinc(N*(phi(mod(K*inc_angle, 180) == 0)-0.12)/pi).^4))*(N^2)/M;

end

