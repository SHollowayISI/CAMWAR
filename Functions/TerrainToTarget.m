function [scenario] = TerrainToTarget(scenario_in)
%TERRAINTOTARGET Generates target list from terrain parameters
%   Takes radar scenario object as input, outputs target list with terrain
%   targets appended

%% Unpack Variables

scenario = scenario_in;
target_list = scenario.target_list;
environment = scenario.environment;
simsetup = scenario.simsetup;

%% Generate Mesh Grid for Targets

% Create axes
x_lin = environment.gnd_xlim(1):environment.gnd_res:environment.gnd_xlim(2);
y_lin = environment.gnd_ylim(1):environment.gnd_res:environment.gnd_ylim(2);

% Create mesh grid of targets
[xcoord, ycoord] = meshgrid(x_lin, y_lin);

%% Generate Elevation Data

% Generate initial variable surface, with padding at edges
pad_size = 10;
elevs_seed = randn(size(xcoord) + 2*pad_size);

% Calculate normalized frequency cutoff
f_cut = environment.gnd_res/environment.gnd_cor_l;

% If cutoff below Nyquist frequency, perform low pass filtering
if f_cut >= 1
    elevs_smooth = elevs_seed;
else
    elevs_smooth = lowpass(elevs_seed, f_cut);
    elevs_smooth = lowpass(elevs_smooth', f_cut)';
end

% Apply smoothing window
n = 3;
elevs_smooth = movmean(movmean(elevs_smooth, ...
    environment.gnd_cor_l/(environment.gnd_res*n))', environment.gnd_cor_l/(environment.gnd_res*n))';

% Remove padding
elevs_smooth = elevs_smooth((pad_size+1):(end-pad_size), (pad_size+1):(end-pad_size));

% Correct variance and mean, and export z-positions
z_pos = elevs_smooth * environment.gnd_var/std(elevs_smooth, [], 'all', 'omitnan') ...
    + environment.gnd_zavg;

% Save surface to environment object
environment.surface.x = xcoord;
environment.surface.y = ycoord;
environment.surface.z = z_pos;

% Generate RCS for broadside terrain
rcs_broadside = environment.gnd_gamma * environment.gnd_res^2;

%% Reduce Number of Targets By Ground Range

% DEBUG: WRITE IN REDUCTION FACTOR
scenario.simsetup.reduce_factor = 1:length(simsetup.elevation_slices);

% DEBUG: CLEAN THIS ALL UP
red_fac = scenario.simsetup.reduce_factor;

r_grid = rcs_broadside * ones(size(xcoord));

theta = simsetup.elevation_slices;
theta_step = diff(theta);
theta_step = theta_step(1);

xl = xcoord(:)';
yl = ycoord(:)';
zl = z_pos(:)';
rl = r_grid(:)';

rem = zeros(size(xl));

% Calculate closest elevation bin of each target
[~, ang] = rangeangle([xl; yl; zl]);
ang_bin = round((ang(2,:) - min(theta))/theta_step) + 1;

% Set targets outside of elevation range to be removed
rem((ang_bin <= 0) | (ang_bin > length(theta))) = 1;
ang_bin(rem == 1) = 1;

rl = rl.*red_fac(ang_bin).^2;

for n = 1:length(rl)
    
    [r, c] = ind2sub(size(xcoord), n);
    
    rem(n) = (rem(n) || ...
             (mod(r, red_fac(ang_bin(n))) > 0) || ...
             (mod(c, red_fac(ang_bin(n))) > 0));
    
end    


%% Generate Location and Velocity Target List

% Reshape target grid into list
tgt_pos = [xl; yl; zl];

% Set all velocity to zero
tgt_vel = zeros(size(tgt_pos));

%% Generate RCS Data

% Generate surface normals for landscape
[Nx, Ny, Nz] = surfnorm(xcoord, ycoord, z_pos);
norm_list = [Nx(:)'; Ny(:)'; Nz(:)'];

% Generate normalized direction vectors from surface to radar system
dir_list = simsetup.radar_pos - tgt_pos;
dir_list = dir_list ./ sqrt(sum(dir_list.^2,1));

% Generate fading factor for each target
fading_list = dot(norm_list, dir_list);

% Generate RCS using fading factor
tgt_rcs = rl .* fading_list;

% Remove occluded targets
rem = (rem | (tgt_rcs <= 0));

% Remove all designated targets
tgt_pos = tgt_pos(:, ~rem);
tgt_vel = tgt_vel(:, ~rem);
tgt_rcs = tgt_rcs(~rem);

%% Re-pack Variables

% Append target lists
target_list.pos = [target_list.pos, tgt_pos];
target_list.vel = [target_list.vel, tgt_vel];
target_list.rcs = [target_list.rcs, tgt_rcs];

% Return variables
scenario.target_list = target_list;
scenario.environment = environment;

end

