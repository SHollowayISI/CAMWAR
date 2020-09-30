
res = 2;

theta = -30:1:0;
theta_step = diff(theta);
theta_step = theta_step(1);

sweeps = 1:length(theta);

red_fac = sweeps;
% red_fac = 2*ones(size(sweeps));

x = -1000:res:1000;
y = 0:res:1500;

[y_grid, x_grid] = meshgrid(y,x);

z_grid = -100*ones(size(x_grid));
r_grid = res*res*ones(size(x_grid));

xl = x_grid(:)';
yl = y_grid(:)';
zl = z_grid(:)';
rl = r_grid(:)';

rem = zeros(size(xl));

% Calculate closest elevation bin of each target
[rng, ang] = rangeangle([xl; yl; zl]);
ang_bin = round((ang(2,:) - min(theta))/theta_step) + 1;

% Set targets outside of elevation range to be removed
rem((ang_bin <= 0) | (ang_bin > length(theta))) = 1;
ang_bin(rem == 1) = 1;

rl_red = rl.*red_fac(ang_bin).^2;

for n = 1:length(rl_red)
    
    [r, c] = ind2sub(size(x_grid), n);
    
    rem(n) = (rem(n) || ...
             (mod(r, red_fac(ang_bin(n))) > 0) || ...
             (mod(c, red_fac(ang_bin(n))) > 0));
    
end    

xl = xl(~rem);
yl = yl(~rem);
zl = zl(~rem);
rl_red = rl_red(~rem);
    
scatter3(xl, yl, zl)
figure;
scatter3(xl, yl, rl_red)


