%% Module for fitting isotherm parameters to batch isotherm data

% Run by pressing start in the editor tab or by pressing F5.
% Individual functions run from this script are not required to be
% modified.

% Start timer and program starting notification
tic; fprintf("\n---------------------------PROGRAM STARTING" + ...
    "---------------------------\n");

% Clear variables and plots
clear; close all;

% Change directory to location of this script
cd(fileparts(which('isofit_module')));


%% Import and process batch isotherm data

% Import and bin isotherm data and parameters into structures
[iso_data,iso_pars] = gen_struct();

% Generate bound concentration data (q) from batch data measurements
[iso_data,iso_pars] = isodata_processing(iso_data,iso_pars);


%% Import parameter bounds

% Parse parameter search range from input sheet
iso_pars = parse_pars(iso_pars,iso_data);


%% Begin optimization to fit parameters

% Run fitting routine using genetic algorithm
[iso_data,iso_pars] = run_opt(iso_data,iso_pars);


%% Set directory for exporting results

% Set directory for output
iso_pars = set_dir(iso_pars);


%% Calculate fitting statistics

% Calculate statistics for isotherm fit results
[iso_data,iso_pars] = calc_stats(iso_data,iso_pars);


%% Generate plots

% Generate 2D and 3D plots of fits to batch data and export results
plot_iso(iso_data,iso_pars);


%% Export fitting results and data

% Write parameters and fit statistics to excel file
export_results(iso_data,iso_pars);

% End timer and program ending notification
fprintf("\n---------------------------PROGRAM ENDING" + ...
    "---------------------------\n"); toc;



