%% CAMWAR Radar System - End of Simulation Tasks
%{

    Sean Holloway
    CAMWAR End-of-Simulation Tasks

    Container script which saves output files from CAMWAR simulation,and
    sends alert email.

    For use in FullSystem_CAMWAR without automated simulation.
    
%}

%% Announce Elapsed Time

toc

%% Save Files and Figures

%TEMP: OVERWRITE

scenario.simsetup.set_filename = false;

% Establish file name
if scenario.simsetup.set_filename
    save_name = [scenario.simsetup.filename, '_', datestr(now, 'mmddyy_HHMM')];
else
    save_name = [scenario_filename(1:end-2), '_', datestr(now, 'mmddyy_HHMM')];    
end

% Establish filepaths for saving
mat_path = 'MAT Files\Scenario Objects\';
fig_path = ['Figures\', save_name, '\'];

% Save scenario object if chosen
if scenario.simsetup.save_mat
    SaveScenario(scenario, save_name, mat_path);
end

% Save open figures if chosen
if scenario.simsetup.save_figs
    for ftype = 1:length(scenario.simsetup.save_format.list)
        SaveFigures(save_name, fig_path, scenario.simsetup.save_format.list{ftype});
    end
end

% Close figures after saving
close all


%% Send Email Alert

% Send email alert with attachment if chosen
if scenario.simsetup.send_alert
    
    % Set up email process
    EmailSetup();
    
    % Send email
    EmailAlert( ...
        scenario.simsetup.alert_address, ...
        save_name, ...
        scenario.simsetup.attach_zip);
end


%% File Management 

% Move file to "Complete" folder
movefile(['Automated Testing/To Run/', scenario_filename], ...
    'Automated Testing/Complete/');

% Display message to command window
disp(['Simulation scenario complete: ', scenario_filename]);










