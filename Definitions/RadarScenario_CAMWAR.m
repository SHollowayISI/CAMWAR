% ClassDef File for CAMWAR Radar Scenario

classdef RadarScenario_CAMWAR < handle
    properties
        environment
        target_list
        simsetup
        radarsetup
        sim
        rx_sig
        cube
        detection
        terrain
        flags
        timing
        results
        calc
    end
    
    methods
        
        function RadarScenario = RadarScenario_CAMWAR()
            % Initialize structure of target list
            RadarScenario.target_list = struct( ...
                'pos',      [], ...
                'vel',      [], ...
                'rcs',      []);
        end
        
        function targetListReset(RadarScenario)
            release(RadarScenario.sim.radiator);
            RadarScenario.sim.radiator = phased.Radiator( ...
                'Sensor',                   RadarScenario.sim.tx_array, ...
                'PropagationSpeed',         physconst('LightSpeed'), ...
                'OperatingFrequency',       RadarScenario.radarsetup.f_c, ...
                'CombineRadiatedSignals',   true);
            
            release(RadarScenario.sim.channel);
            RadarScenario.sim.channel = phased.LOSChannel( ...
                'PropagationSpeed',         physconst('LightSpeed'), ...
                'OperatingFrequency',       RadarScenario.radarsetup.f_c, ...
                'SpecifyAtmosphere',        RadarScenario.environment.channel_loss_on, ...
                'TwoWayPropagation',        true, ...
                'SampleRate',               RadarScenario.radarsetup.f_s);
            
            release(RadarScenario.sim.target);
            RadarScenario.sim.target = phased.RadarTarget( ...
                'EnablePolarization',       false, ...
                'MeanRCSSource',            'Input Port', ...
                'Model',                    'Swerling1', ...
                'OperatingFrequency',       RadarScenario.radarsetup.f_c, ...
                'PropagationSpeed',         physconst('LightSpeed'));
            
            release(RadarScenario.sim.collector);
            RadarScenario.sim.collector = phased.Collector( ...
                'Sensor',                   RadarScenario.sim.rx_array, ...
                'PropagationSpeed',         physconst('LightSpeed'), ...
                'OperatingFrequency',       RadarScenario.radarsetup.f_c, ...
                'Wavefront',                'Plane');
        end
        
        function generateCoordinates2D(RadarScenario)
            %Generate Input Coordinate Grid
            [range_grid, azimuth_grid] = meshgrid( ...
                RadarScenario.cube.range_axis, ...
                RadarScenario.cube.azimuth_axis);
            %Generate Output Coordinate Grid
            RadarScenario.cube.x_grid = transpose(range_grid .* cosd(azimuth_grid));
            RadarScenario.cube.y_grid = transpose(range_grid .* sind(azimuth_grid));
        end
        
        function generateCoordinates3D(RadarScenario)
            %Generate Input Coordinate Grid
            [range_grid, azimuth_grid, elevation_grid] = meshgrid( ...
                RadarScenario.cube.range_axis, ...
                RadarScenario.cube.azimuth_axis, ...
                RadarScenario.simsetup.elevation_slices);
            %Generate Output Coordinate Grid
            RadarScenario.results.x_grid = ...
                permute( ...
                range_grid .* cosd(azimuth_grid) .* cosd(elevation_grid), ...
                [2 1 3]);
            RadarScenario.results.y_grid = ...
                permute( ...
                range_grid .* sind(azimuth_grid) .* cosd(elevation_grid), ...
                [2 1 3]);
            RadarScenario.results.z_grid = ...
                permute( ...
                range_grid .* sind(elevation_grid), ...
                [2 1 3]);
        end
        
        function timeStart(RadarScenario)
            % Begin timing for progress readout
            RadarScenario.timing.timing_logical = true;
            RadarScenario.timing.startTime = tic;
            RadarScenario.timing.TimeDate = now;
            RadarScenario.timing.numLoops = ...
                length(RadarScenario.simsetup.elevation_slices);
            RadarScenario.timing.timeGate = 0;
            
            % OBSELETE: Time estimation from before resolution correction
            %{
            % Estimate progress for each elevation slice
            el = RadarScenario.simsetup.elevation_slices;
            bw = RadarScenario.simsetup.slice_bw;
            FOV = RadarScenario.simsetup.fov;
            res = RadarScenario.environment.gnd_res;
            h = RadarScenario.environment.gnd_zavg;

            total_area = min([2*pi*(FOV/360)*(h^2)*...
                (cotd(el+bw/2).^2 - cotd(el-bw/2).^2); ...
                2*pi*(FOV/360)*((1000^2) - (h*cotd(el-bw/2).^2))], [], 1);
            total_targets = (total_area/(res^2))./(RadarScenario.simsetup.reduce_factor.^2);
            runtime = (20.124*(total_targets/1000).^2 + 2.8422*(total_targets/1000)+ 5.0627)/60;
            total = sum(runtime);
            RadarScenario.timing.progress = cumsum(runtime)/total;
            %}
            
        end
        
        function timeUpdate(RadarScenario, repetition, rep_method)
            
            if ~RadarScenario.timing.timing_logical
                error('Must use method timeStart() before timeUpdate()');
            end
            
            % Calculate progress through simulation
            loops_complete = RadarScenario.flags.slice;
            percent_complete = 100*loops_complete/RadarScenario.timing.numLoops;
            
            % Calculate remaining time in simulation
            nowTimeDate = now;
            elapsedTime = nowTimeDate - RadarScenario.timing.TimeDate;
            estComplete = nowTimeDate + ((100/percent_complete)-1)*elapsedTime;
            
            % Form message to display in command window
            message_l = sprintf('%d Loops complete out of %d', loops_complete, RadarScenario.timing.numLoops);
            message_p = [sprintf('Percent complete: %0.0f', percent_complete), '%'];
            message_t = ['Estimated time of completion: ', datestr(estComplete)];
            
            % Display current progress
            switch rep_method
                case 'loops'
                    if (mod(loops_complete, repetition) == 1) || (repetition == 1)
                        disp(message_l);
                        disp(message_p);
                        disp(message_t);
                        disp('');
                    end
                    
                case 'time'
                    if ((RadarScenario.timing.timeGate == 0) || ...
                            (toc > repetition + RadarScenario.timing.timeGate))
                        disp(message_p);
                        disp(message_t);
                        disp('');
                        RadarScenario.timing.timeGate = toc;
                    end
                    
            end
        end
        
        function readOut(RadarScenario, result)
            switch result
                case 'max'
                    fprintf('\nMaximum SNR: %0.1f [dB]\n', ...
                        RadarScenario.detection.max_SNR);
                case 'total'
                    fprintf('\nTotal Power SNR: %0.1f [dB]\n', ...
                        RadarScenario.detection.total_SNR);
                case 'detect'
                    fprintf('\nDetection Power SNR: %0.1f [dB]\n', ...
                        RadarScenario.detection.detect_SNR);
                case 'noise'
                    fprintf('\nNoise Power: %0.1f [dBW]\n', ...
                        RadarScenario.detection.noise_pow);
                case 'all'
                    fprintf('\nMaximum SNR: %0.1f [dB]\n', ...
                        RadarScenario.detection.max_SNR);
                    fprintf('Total Power SNR: %0.1f [dB]\n', ...
                        RadarScenario.detection.total_SNR);
                    fprintf('Detection Power SNR: %0.1f [dB]\n', ...
                        RadarScenario.detection.detect_SNR);
                    fprintf('Noise Power: %0.1f [dBW]\n', ...
                        RadarScenario.detection.noise_pow);
            end     
        end
        
        function viewTargets(RadarScenario)
            % Show 3-D scatter plot of target locations
            figure('Name', 'Target 3D Scatter Plot')
            scatter3(RadarScenario.target_list.pos(1,:), ...
                RadarScenario.target_list.pos(2,:), ...
                RadarScenario.target_list.pos(3,:))
        end
        
        function viewTerrainRCS(RadarScenario)
            % Show 2-D spatial scatter plot with RCS as 3rd dimension
            figure('Name', 'Target RCS Scatter Plot')
            scatter3(RadarScenario.target_list.pos(1,:), ...
                RadarScenario.target_list.pos(2,:), ...
                RadarScenario.target_list.rcs)
        end
        
        function viewRDCube(RadarScenario, graphType)
            if strcmp(graphType, 'heatmap')
                figure('Name', 'Range-Doppler Heat Map');
                imagesc(RadarScenario.cube.vel_axis, ...
                    RadarScenario.cube.range_axis, ...
                    10*log10(RadarScenario.cube.pow_cube(:,:,ceil(end/2))))
                title('Range-Doppler Heat Map')
                set(gca,'YDir','normal')
                xlabel('Velocity [m/s]','FontWeight','bold')
                ylabel('Range [m]','FontWeight','bold')
            else
                figure('Name', 'Range-Doppler Surface');
                surf(RadarScenario.cube.vel_axis, ...
                    RadarScenario.cube.range_axis, ...
                    10*log10(RadarScenario.cube.pow_cube(:,:,ceil(end/2))), ...
                    'EdgeColor', 'none')
                title('Range-Doppler Surface')
                set(gca,'YDir','normal')
                xlabel('Velocity [m/s]','FontWeight','bold')
                ylabel('Range [m]','FontWeight','bold')
                zlabel('FFT Log Intensity [dB]','FontWeight','bold')
            end
            
        end
        
        function viewRACube(RadarScenario, graphType)
            switch graphType
                case 'heatmap'
                    figure('Name', 'Range-Azimuth Heat Map');
                    imagesc(RadarScenario.cube.azimuth_axis, ...
                        RadarScenario.cube.range_axis, ...
                        10*log10(squeeze(RadarScenario.cube.pow_cube(:,ceil(end/2),:))))
                    title('Range-Azimuth Heat Map')
                    set(gca,'YDir','normal')
                    xlabel('Azimuth Angle [degree]','FontWeight','bold')
                    ylabel('Range [m]','FontWeight','bold')
                case 'surface'
                    figure('Name', 'Range-Azimuth Surface');
                    surf(RadarScenario.cube.azimuth_axis, ...
                        RadarScenario.cube.range_axis, ...
                        10*log10(squeeze(RadarScenario.cube.pow_cube(:,ceil(end/2),:))), ...
                        'EdgeColor', 'none')
                    title('Range-Azimuth Surface')
                    xlabel('Azimuth Angle [degree]','FontWeight','bold')
                    ylabel('Range [m]','FontWeight','bold')
                    zlabel('FFT Log Intensity [dB]','FontWeight','bold')
                case 'PPI'
                    figure('Name', 'Range-Azimuth PPI');
                    surf(RadarScenario.cube.x_grid, ...
                        RadarScenario.cube.y_grid, ...
                        10*log10(squeeze(RadarScenario.cube.pow_cube(:,ceil(end/2),:))), ...
                        'EdgeColor', 'none')
                    title('Range-Azimuth PPI')
                    xlabel('Cross-Range Distance [m]','FontWeight','bold')
                    ylabel('Down-Range Distance [m]','FontWeight','bold')
                    zlabel('FFT Log Intensity [dB]','FontWeight','bold')
            end
            
        end
        
        function viewDetections(RadarScenario, graphType)
            switch graphType
                case 'heatmap'
                    figure('Name', 'Detection Heatmap')
                    imagesc(RadarScenario.cube.azimuth_axis, ...
                        RadarScenario.cube.range_axis, ...
                        RadarScenario.detection.detect_cube_nodop( ...
                        :, :, RadarScenario.flags.slice))
                    set(gca, 'YDir', 'Normal')
                    xlabel('Azimuth Angle [degree]', 'FontWeight', 'bold')
                    ylabel('Ragne [m]', 'FontWeight', 'bold')
                case 'PPI'
                    figure('Name', 'Detection PPI');
                    surf(RadarScenario.cube.x_grid, ...
                        RadarScenario.cube.y_grid, ...
                        RadarScenario.detection.detect_cube_nodop( ...
                        :, :, RadarScenario.flags.slice), ...
                        'EdgeColor', 'none')
                    view(270,90)
                    title('Range-Azimuth PPI')
                    xlabel('Cross-Range Distance [m]','FontWeight','bold')
                    ylabel('Down-Range Distance [m]','FontWeight','bold')
                    zlabel('FFT Log Intensity [dB]','FontWeight','bold')
            end
        end
        
        function viewDetectionsCFAR(RadarScenario, graphType)
            switch graphType
                case 'heatmap'
                    figure('Name', 'Detection Heatmap')
                    imagesc(RadarScenario.cube.azimuth_axis, ...
                        RadarScenario.cube.range_axis, ...
                        RadarScenario.detection.CFAR_cube_nodop( ...
                        :, :, RadarScenario.flags.slice))
                    set(gca, 'YDir', 'Normal')
                    xlabel('Azimuth Angle [degree]', 'FontWeight', 'bold')
                    ylabel('Ragne [m]', 'FontWeight', 'bold')
                case 'PPI'
                    figure('Name', 'Detection PPI');
                    surf(RadarScenario.cube.x_grid, ...
                        RadarScenario.cube.y_grid, ...
                        RadarScenario.detection.CFAR_cube_nodop( ...
                        :, :, RadarScenario.flags.slice), ...
                        'EdgeColor', 'none')
                    view(270,90)
                    title('Range-Azimuth PPI')
                    xlabel('Cross-Range Distance [m]','FontWeight','bold')
                    ylabel('Down-Range Distance [m]','FontWeight','bold')
                    zlabel('FFT Log Intensity [dB]','FontWeight','bold')
            end
        end
        
        function viewDetectionsFull(RadarScenario)
            x_list = RadarScenario.results.x_grid( ...
                RadarScenario.detection.detect_cube_nodop);
            y_list = RadarScenario.results.y_grid( ...
                RadarScenario.detection.detect_cube_nodop);
            z_list = RadarScenario.results.z_grid( ...
                RadarScenario.detection.detect_cube_nodop);
            
            figure('Name', 'Detection 3-D Scatter Plot')
            scatter3(x_list, y_list, z_list)
            title('Detection 3-D Scatter Plot')
            xlabel('Cross-Range Distance [m]','FontWeight','bold')
            ylabel('Down-Range Distance [m]','FontWeight','bold')
            zlabel('Elevation [m]','FontWeight','bold')
        end
        
        function viewDetectionsCFARFull(RadarScenario)
            x_list = RadarScenario.results.x_grid( ...
                RadarScenario.detection.CFAR_cube_nodop);
            y_list = RadarScenario.results.y_grid( ...
                RadarScenario.detection.CFAR_cube_nodop);
            z_list = RadarScenario.results.z_grid( ...
                RadarScenario.detection.CFAR_cube_nodop);
            
            figure('Name', 'Detection 3-D Scatter Plot')
            scatter3(x_list, y_list, z_list)
            title('Detection 3-D Scatter Plot')
            xlabel('Cross-Range Distance [m]','FontWeight','bold')
            ylabel('Down-Range Distance [m]','FontWeight','bold')
            zlabel('Elevation [m]','FontWeight','bold')
        end
        
        function viewDetectionsDiscreteFull(RadarScenario)
            x_list = RadarScenario.results.x_grid( ...
                RadarScenario.detection.discrete_cube_any);
            y_list = RadarScenario.results.y_grid( ...
                RadarScenario.detection.discrete_cube_any);
            z_list = RadarScenario.results.z_grid( ...
                RadarScenario.detection.discrete_cube_any);
            
            figure('Name', 'Detection 3-D Scatter Plot')
            scatter3(x_list, y_list, z_list)
            title('Detection 3-D Scatter Plot')
            xlabel('Cross-Range Distance [m]','FontWeight','bold')
            ylabel('Down-Range Distance [m]','FontWeight','bold')
            zlabel('Elevation [m]','FontWeight','bold')
            
        end
        
        function viewDetectionsDiscreteCentroid(RadarScenario)
            x_list = RadarScenario.detection.detect_list(1,:);
            y_list = RadarScenario.detection.detect_list(2,:);
            z_list = RadarScenario.detection.detect_list(3,:);
            
            figure('Name', 'Detection 3-D Scatter Plot')
            scatter3(x_list, y_list, z_list)
            title('Detection 3-D Scatter Plot')
            xlabel('Cross-Range Distance [m]','FontWeight','bold')
            ylabel('Down-Range Distance [m]','FontWeight','bold')
            zlabel('Elevation [m]','FontWeight','bold')
            
        end
        
        function viewRealTerrain(RadarScenario)
            figure('Name', 'Generated Ground Surface');
            surf(RadarScenario.environment.surface.x, ...
                RadarScenario.environment.surface.y, ...
                RadarScenario.environment.surface.z, ...
                'EdgeColor', 'none')
            title('Generated Ground Surface')
            xlabel('Cross-Range Distance [m]','FontWeight','bold')
            ylabel('Down-Range Distance [m]','FontWeight','bold')
            zlabel('Elevation [m]','FontWeight','bold')
            zlim([2 * RadarScenario.environment.gnd_zavg 0])
        end
        
        function viewTerrain(RadarScenario)
            if length(RadarScenario.simsetup.elevation_slices) > 1
                figure('Name', 'Estimated Ground Surface');
                surf(squeeze(RadarScenario.terrain.locations(1,:,:)), ...
                    squeeze(RadarScenario.terrain.locations(2,:,:)), ...
                    squeeze(RadarScenario.terrain.locations(3,:,:)), ...
                    'EdgeColor', 'none')
                title('Estimated Ground Surface')
                xlabel('Cross-Range Distance [m]','FontWeight','bold')
                ylabel('Down-Range Distance [m]','FontWeight','bold')
                zlabel('Elevation [m]','FontWeight','bold')
                zlim([2 * RadarScenario.environment.gnd_zavg 0])
            end
        end
        
    end
end






