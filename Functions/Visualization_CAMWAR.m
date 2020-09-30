function [] = Visualization_CAMWAR(scenario)
%VISUALIZATION_CAMWAR Displays "Plan Position Indicator" and POV view for CAMWAR project
%   Take radar scenario data object as input, returns graph of PPI view of
%   HAT results.

%% Unpack Variables

surface_RCS = scenario.detection.surface_RCS;
detection = scenario.detection;

%% User Options

% Discrete Point Settings
mkr = '+';
sz = 200;

% Surface Visualization Settings
interp_method = 'nearest';
extrap_method = 'nearest';
show_real = true;

% Color Map Settings
color_min = [0, 0.5, 0];
color_max = [1, 1, 0];
color_lim = [-30 0];

%% PPI Calculation

% Determine maximum elevation bin to use
[~, max_idx] = find(scenario.simsetup.elevation_slices <= ...
    scenario.simsetup.HAT_max, 1, 'last');

% Create condensed terrain map
terr_x = squeeze(scenario.terrain.locations(1,:,:));
terr_y = squeeze(scenario.terrain.locations(2,:,:));
terr_z = squeeze(scenario.terrain.locations(3,:,:));

% Determine ground distances from radar system
terr_ground_distance = sqrt(terr_x.^2 + terr_y.^2);

% Set maximum visualization distance based on method selection
skip_max = false;
switch scenario.simsetup.max_range_method
    case 'max'
        r_max = scenario.cube.range_axis(end);
    case 'fit'
        r_max = max(terr_ground_distance(:,1:max_idx), [], 'all');
    case 'set'
        r_max = scenario.simsetup.set_max;
    case 'none'
        r_max = max(terr_ground_distance(:,1:max_idx), [], 'all');
        skip_max = true;
end
r_max_bin = round(r_max/scenario.cube.range_res + 0.5);


% Set minimum visualization distance based on method selection
skip_min = false;
switch scenario.simsetup.min_range_method
    case 'min'
        r_min = scenario.cube.range_axis(1);
    case 'fit'
        r_min = min(terr_ground_distance, [], 'all');
    case 'set'
        r_min = scenario.simsetup.set_min;
    case 'none'
        r_min = min(terr_ground_distance, [], 'all');
        skip_min = true;
end
r_min_bin = round(r_min/scenario.cube.range_res + 0.5);

% Remove HAT results over the maximum angle
terr_x(:, (max_idx+1):end) = nan;
terr_y(:, (max_idx+1):end) = nan;
terr_z(:, (max_idx+1):end) = nan;


% Interpolate throught NaN values
if scenario.simsetup.interpolate_nan
    
    az = scenario.cube.azimuth_axis;
    el = scenario.simsetup.elevation_slices;
    [el_grid, az_grid] = meshgrid(el, az);
    
    % Determine indices of NaNs
    [nan_r, nan_c] = find(isnan(terr_x));
    nan_az = az(nan_r)';
    nan_el = el(nan_c)';
    
    % Determine indices of non-NaNs
    [good_r, good_c] = find(~isnan(terr_x));
    good_az = az(good_r)';
    good_el = el(good_c)';
    
    % Interpolate elevation data through NaN values
    interp_x = scatteredInterpolant(good_az, good_el, terr_x(~isnan(terr_x)), 'linear');
    interp_y = scatteredInterpolant(good_az, good_el, terr_y(~isnan(terr_y)), 'linear');
    interp_z = scatteredInterpolant(good_az, good_el, terr_z(~isnan(terr_z)), interp_method);
    
    terr_x = interp_x(az_grid, el_grid);
    terr_y = interp_y(az_grid, el_grid);
    terr_z = interp_z(az_grid, el_grid);
    
end

% % Perform smoothing on surface (TEST)
% el_smooth = 4;
% az_smooth = 8;
% terr_z = movmean(movmean(terr_z, el_smooth, 2), az_smooth, 1); 

% Pre-allocation
extrap_x_max = [];
extrap_y_max = [];
extrap_z_max = [];
extrap_rcs_max = [];
extrap_x_min = [];
extrap_y_min = [];
extrap_z_min = [];
extrap_rcs_min = [];

if ~skip_max
    
    % Pare down results if exceed maximum distance
    [~, col] = find(terr_ground_distance > r_max, 1, 'first');
    if ~isempty(col)
        
        terr_x = terr_x(:,1:(col-1));
        terr_y = terr_y(:,1:(col-1));
        terr_z = terr_z(:,1:(col-1));
%         surface_RCS = surface_RCS(:,1:(col-1));
    end
    
    % Extrapolate out to maximum range
    az = scenario.cube.azimuth_axis;
    rax = scenario.cube.range_axis(r_max_bin:end);
    
    extrap_x_max = rax .* cosd(az');
    extrap_y_max = rax .* sind(az');
    
    extrap_z_obj = scatteredInterpolant(terr_x(:), terr_y(:), terr_z(:), interp_method, extrap_method);
%     extrap_rcs_obj = scatteredInterpolant(terr_x(:), terr_y(:), surface_RCS(:), interp_method, extrap_method);
    
    extrap_z_max = extrap_z_obj(extrap_x_max, extrap_y_max);
%     extrap_rcs_max = extrap_rcs_obj(extrap_x_max, extrap_y_max);
end

% Pare down results if exceed minimum distance
if ~skip_min
    [~, col] = find(terr_ground_distance < r_min, 1, 'last');
    if ~isempty(col)
        
        terr_x = terr_x(:,(col+1):end);
        terr_y = terr_y(:,(col+1):end);
        terr_z = terr_z(:,(col+1):end);
%         surface_RCS = surface_RCS(:,(col+1):end);
    end
    
    % Extrapolate out to minimum range
    az = scenario.cube.azimuth_axis;
    rax = scenario.cube.range_axis(1:r_min_bin);
    
    extrap_x_min = rax .* cosd(az');
    extrap_y_min = rax .* sind(az');
    
    extrap_z_obj = scatteredInterpolant(terr_x(:), terr_y(:), terr_z(:), interp_method, extrap_method);
%     extrap_rcs_obj = scatteredInterpolant(terr_x(:), terr_y(:), surface_RCS(:), interp_method, extrap_method);
    
    extrap_z_min = extrap_z_obj(extrap_x_min, extrap_y_min);
%     extrap_rcs_min = extrap_rcs_obj(extrap_x_min, extrap_y_min);
end

% Combine extrapolated surfaces
terr_x = [extrap_x_min, terr_x, extrap_x_max];
terr_y = [extrap_y_min, terr_y, extrap_y_max];
terr_z = [extrap_z_min, terr_z, extrap_z_max];
surface_RCS = [ones(size(extrap_x_min))*(-Inf), surface_RCS];
surface_RCS = [surface_RCS, ones(size(surface_RCS,1), ...
    size(terr_x,2) - size(surface_RCS,2))*-Inf];

%% Set up color map

map = repmat(color_min, 100, 1);

map = [map; ...
       linspace(color_min(1),color_max(1))', ...
       linspace(color_min(2),color_max(2))', ...
       linspace(color_min(3),color_max(3))'];
map = [[0, 0, 0]; map];

%% Create surface plot
figure('Name', 'Top Down PPI Display');
if show_real
    subplot(121)
end
surf(terr_x, terr_y, terr_z, ...
    surface_RCS, ...
    'EdgeColor', 'none', 'FaceColor', 'interp');
hold on

% Set scale for contour lines
contour_scale = scenario.simsetup.contour_scale;
max_z = contour_scale*ceil(max(terr_z, [], 'all')/contour_scale);
min_z = contour_scale*floor(min(terr_z, [], 'all')/contour_scale);
scale = [min_z:contour_scale:max_z];

% Plot contour lines over surface
contour(terr_x, terr_y, terr_z, ...
    scale, ...
    'k', ...
    'ShowText', 'on');
view(-90,90)

% Set correct color bar scale
colormap(map)
set(gca, 'CLim', color_lim)

% Append discrete target list
hold on
if ~isempty(detection.detect_list)
    scatter3( ...
        detection.detect_list(1,:), ...
        detection.detect_list(2,:), ...
        detection.detect_list(3,:), ...
        sz, mkr, 'r');
end

% Set up labels
ylabel('Cross-Range Distance [m]', 'FontWeight', 'bold')
xlabel('Down-Range Distance [m]', 'FontWeight', 'bold')
title('Top Down PPI Display: Radar Response')

% Set up limits
y_max = 2*r_max*sind(scenario.simsetup.fov/2);
lim_max = max([y_max, r_max]);
lim_max = 100*ceil(lim_max/100);
rlim_max = 100*ceil(r_max/100);
xlim([0 lim_max])
ylim([-lim_max/2 lim_max/2])

% PPI for real terrain
if show_real
    subplot(122)
    real_t = scenario.environment.surface;
    [rng, ang] = rangeangle([real_t.x(:)'; real_t.y(:)'; real_t.z(:)']);
    real_t.x(abs(ang(1,:)) > scenario.simsetup.fov/2) = nan;
    real_t.y(abs(ang(1,:)) > scenario.simsetup.fov/2) = nan;
    real_t.z(abs(ang(1,:)) > scenario.simsetup.fov/2) = nan;
    real_t.x(rng > r_max) = nan;
    real_t.y(rng > r_max) = nan;
    real_t.z(rng > r_max) = nan;
    
    surf(real_t.x, real_t.y, real_t.z, ...
        -10*ones(size(real_t.x)), ...
        'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(map)
    hold on
    
    % Plot contour lines over surface
    contour(real_t.x, real_t.y, real_t.z, ...
        scale, ...
        'k', ...
        'ShowText', 'on');
    view(-90,90)
    
    % Set correct color bar scale
    colormap(map)
    set(gca, 'CLim', color_lim)
    
    % Set up labels
    ylabel('Cross-Range Distance [m]', 'FontWeight', 'bold')
    xlabel('Down-Range Distance [m]', 'FontWeight', 'bold')
    title('Top Down PPI Display: Real Terrain')
    
    % Set up limits
    xlim([0 lim_max])
    ylim([-lim_max/2 lim_max/2])
end

%% Create Pilot Perspective Plot

% Create same PPI plot as before
figure('Name', 'Perspective PPI Display Radar Response')
title('Perspective PPI Display: Radar Response')
surf(terr_x, terr_y, terr_z, ...
    surface_RCS, ...
    'EdgeColor', 'none', 'FaceColor', 'interp');
hold on

% Plot contour lines at exact z-levels
for lev = 1:length(scale)
    % Apply single contour
    c = contour(terr_x, terr_y, terr_z, ...
        [scale(lev) scale(lev)], ...
        'k', ...
        'LineWidth', 1e-9);
    hold on
    
    % Set plot Z-Level
    ax = gca;
    ax.Children(1).ContourZLevel = scale(lev);
    
    % Apply text to contour
    t = clabel(c,scale(lev));
    hold on
    
    % Loop through text objects
    for tn = 2:2:length(t)
        t(tn-1).Visible = 'off';
        t(tn).Position(3) = scale(lev) + 2;
    end
    
end

% Set correct color bar scale
colormap(map)
set(gca, 'CLim', color_lim)

% Append discrete target list
hold on
if ~isempty(detection.detect_list)
scatter3( ...
    detection.detect_list(1,:), ...
    detection.detect_list(2,:), ...
    detection.detect_list(3,:), ...
    sz, mkr, 'r');
end

% Set limits to very large numbers
s_n = r_max;
xlim([-s_n s_n])
ylim([-s_n s_n])
zlim([-s_n s_n])

% Set camera properties for POV perspective
l_n = 1e20;
% set(gca, 'Visible', 'off')
camproj('perspective')
campos(scenario.simsetup.radar_pos)
camtarget([l_n 0 0])
camva(60)

% Display real terrain
if show_real
    figure('Name', 'Perspective PPI Display Real Terrain')
    title('Perspective PPI Display: Real Terrain')
    surf(real_t.x, real_t.y, real_t.z, ...
        -10*ones(size(real_t.x)), ...
        'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on
    
    % Plot contour lines at exact z-levels
    for lev = 1:length(scale)
        % Apply single contour
        c = contour(real_t.x, real_t.y, real_t.z, ...
            [scale(lev) scale(lev)], ...
            'k', ...
            'LineWidth', 1e-9);
        hold on
        
        % Set plot Z-Level
        ax = gca;
        ax.Children(1).ContourZLevel = scale(lev);
        
        % Apply text to contour
        t = clabel(c,scale(lev));
        hold on
        
        % Loop through text objects
        for tn = 2:2:length(t)
            t(tn-1).Visible = 'off';
            t(tn).Position(3) = scale(lev) + 2;
        end        
    end
    
    % Set limits to very large numbers
    s_n = r_max;
    xlim([-s_n s_n])
    ylim([-s_n s_n])
    zlim([-s_n s_n])
    
    % Set camera properties for POV perspective
    l_n = 1e20;
    % set(gca, 'Visible', 'off')
    camproj('perspective')
    campos(scenario.simsetup.radar_pos)
    camtarget([l_n 0 0])
    camva(60)
    
    % Set correct color bar scale
    colormap(map)
    set(gca, 'CLim', color_lim)
end

%% Create C-Scope Plot

% Set up axes for C-Scope
el_axis = scenario.simsetup.elevation_slices;
surf_RCS = [scenario.detection.surface_RCS, ...
    nan(size(scenario.detection.surface_RCS, 1), ...
    length(el_axis) - size(scenario.detection.surface_RCS, 2))];

% Plot C-Scope
figure('Name', 'C Scope')
if show_real
    subplot(121)
end
imagesc(scenario.cube.azimuth_axis, ...
    el_axis, ...
    surf_RCS')
set(gca, 'YDir', 'normal')
xlabel('Azimuth Angle [degree]', 'FontWeight', 'bold')
ylabel('Elevation Angle [degree]', 'FontWeight', 'bold')
title('C-Scope Display: Radar Response')

% Append discrete target list
hold on
if ~isempty(detection.detect_list)
    scatter(detection.centroid_angles(1), detection.centroid_angles(2), ...
        sz, mkr, 'r');
end

% Set correct color bar scale
colormap(map)
set(gca, 'CLim', color_lim)

% Plot C-Scope for real terrain
if show_real
    subplot(122)
    imagesc(scenario.cube.azimuth_axis, ...
        el_axis, ...
        -10*ones(size(surf_RCS)))
    colormap(map)
    set(gca, 'YDir', 'normal')
    xlabel('Azimuth Angle [degree]', 'FontWeight', 'bold')
    ylabel('Elevation Angle [degree]', 'FontWeight', 'bold')
    title('C-Scope Display: Real Terrain')
    
    % Set correct color bar scale
colormap(map)
set(gca, 'CLim', color_lim)
end

end









