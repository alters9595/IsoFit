%% Generates structures containing isotherm data and parameters
function [iso_data,iso_pars] = gen_struct()

% Load adsorption isotherm data from excel sheet
[data,~,~] = xlsread('Sheet for adsorption isotherm data entry',...
    'Adsorption Isotherm Data Sheet');

% Load fitting options from excel sheet
[~,~,opts] = xlsread('Sheet for adsorption isotherm data entry',...
    'Fitting Options');

% Bin data into isotherm data structure
iso_data.pH = data(:,1);             % pH
iso_data.salt = data(:,2);           % Salt concentration (mM)
iso_data.load_conc = data(:,3);      % Load concentration (mg/mL)
iso_data.flowthrough = data(:,4);    % Flowthrough concentration (mg/mL)
iso_data.resin_vol = data(:,5);      % Resin volume (uL)
iso_data.load_vol = data(:,6);       % Load volume (uL)
iso_data.ft_frac = data(:,7:11);     % Flowthrough mass fractions (%)
iso_data.feed_frac = data(:,14);     % Feed mass fractions (%)
iso_data.mol_weight = data(:,15);    % Molecular weight (kDa)

% Bin options into parameters structure
iso_pars.num_comp = opts{2,2};       % Number of components (#)
iso_pars.isotherm = opts{3,2};       % Isotherm model (name)
iso_pars.pH_dep = opts{4,2};         % pH dependency (Y/N)
iso_pars.ref_pH = opts{5,2};         % Reference pH (pH)
iso_pars.holdup_frac = opts{6,2};    % Hold-up fraction (#) 
iso_pars.meas_error = opts{7,2};     % Batch data measurement error (mg/mL) 
iso_pars.use_parallel = opts{8,2};   % Option to use parallel processing
iso_pars.resin_name = opts{9,2};     % Name of resin sample
iso_pars.protein_names = opts{10,2}; % Names of protein components
iso_pars.mode_oper = opts{11,2};     % Resin mode of operation
iso_pars.ionic_cap = opts{12,2};     % Ionic capacity (M)
iso_pars.porosity_Ec = opts{13,2};   % Interstitial porosity Ee
iso_pars.porosity_Ep = opts{14,2};   % Intraparticle porosity Ep
iso_pars.pH_asmn = opts{15,2};       % pH assumption for multicomp.
iso_pars.size_asmn = opts{16,2};     % Size assumption for multicomp.
iso_pars.ks_asmn = opts{17,2};       % Ks assumption for multicomp.
iso_pars.kp_asmn = opts{18,2};       % Kp assumption for multicomp.
iso_pars.beta_asmn = opts{19,2};     % Beta assumption for multicomp.
iso_pars.plot_lines = opts{20,2};    % Setting for plot lines (Y/N)
iso_pars.export_folder = opts{21,2}; % Name of export folder for results
iso_pars.opt_tolerance = opts{22,2}; % Optimization tolerance cutoff
iso_pars.ga_max_gens = opts{23,2};   % Max number of generations for GA

end

