%% Processing of adsorption isotherm data
function [iso_data,iso_pars] = isodata_processing(iso_data,iso_pars)

% Load variables from input data structure
c_load = iso_data.load_conc;                 % Load concentration (mg/ml)
c_ft = iso_data.flowthrough;                 % FT concentration (mg/ml)
v_resin = iso_data.resin_vol;                % Resin volume (ul)
v_load = iso_data.load_vol;                  % Load volume (ul)
x_ft = iso_data.ft_frac;                     % Mass fraction in FT (%)
x_load = iso_data.feed_frac;                 % Mass fraction in load (%) 
salt = iso_data.salt;                        % Salt concentration (mM)
pH = iso_data.pH;                            % pH
mol_weight = iso_data.mol_weight;            % Molecular weight (kDa)

% Load variables from input parameters structure
f_holdup = iso_pars.holdup_frac;             % Hold-up volume fraction (#)
n_comp = iso_pars.num_comp;                  % Number of components
p_names = iso_pars.protein_names;            % Protein names

% Check for missing data, throw error if any are missing
if any(isnan(salt))
    error('ERROR: Missing salt concentrations in data sheet!')
elseif any(isnan(pH))
    error('ERROR: Missing pH values in data sheet!');
elseif any(isnan(c_load))
    error('ERROR: Missing load concentrations in data sheet!');
elseif any(isnan(c_ft))
    error('ERROR: Missing flowthrough concentrations in data sheet!');
elseif any(isnan(v_resin))
    error('ERROR: Missing resin volumes in data sheet!');
elseif any(isnan(v_load))
    error('ERROR: Missing load volumes in data sheet!');
end

% Calculate hold-up volume and provide volumetric correction
v_holdup = v_resin*f_holdup;
corr_holdup = v_load./(v_load + v_holdup);

% Adjust load concentration and phase ratio
c_load = c_load.*corr_holdup;
phase_ratio = (v_load + v_holdup)./v_resin;

% Follow multicomponent case if more than one component
if n_comp > 1                           
    % Match mass fractions to number of components, convert to decimal
    x_ft = x_ft(:,1:n_comp)/100;
    x_load = x_load(1:n_comp)'/100;
    
    % Remove data points if a component is missing mass fraction data
    i = 1;
    while i < size(x_ft,1)                    % For each data point
        if any(isnan(x_ft(i,:)))
            c_ft(i,:) = []; c_load(i) = []; 
            v_resin(i) = []; v_load(i) = [];
            phase_ratio(i) = [];
            salt(i) = []; pH(i) = [];
            x_ft(i,:) = [];
        else
            i = i + 1;
        end
    end
    
    % Calculate fluid phase concentrations for each component
    c_load = c_load.*x_load;
    c_ft = c_ft.*x_ft;  
elseif n_comp == 1
    % Set homogeneous composition for single component
    x_ft = ones(size(c_load));
    x_load = 1;
elseif n_comp < 1
    % Throw error if number of components is not greater than zero
    error('ERROR: Number of components must be greater than zero!')
end

% Calculate bound concentration
q = (c_load - c_ft).*phase_ratio;

% Set negative bound conc. to zero
q(q < 0) = 0;

% Remove properties for unused component
mol_weight = mol_weight(1:n_comp)';
p_names = strsplit(p_names,', ');
p_names = string(p_names(1:n_comp));

% Repacking variables into structure for output
iso_data.load_conc = c_load;                % Load concentration (mg/ml)
iso_data.flowthrough = c_ft;                % FT concentration (mg/ml)
iso_data.bound_conc = q;                    % Bound concentration (mg/ml)
iso_data.ft_frac = x_ft;                    % Mass fraction in flowthrough
iso_data.feed_frac = x_load;                % Mass fraction in load
iso_data.resin_vol = v_resin;               % Resin volume (ul)
iso_data.load_vol = v_load;                 % Load volume (ul)
iso_data.perc_comps = x_ft;                 % Mass fraction in FT (%)
iso_data.salt = salt;                       % Salt concentration (mM)
iso_data.pH = pH;                           % pH
iso_pars.mol_weight = mol_weight;           % Molecular weight (kDa)
iso_pars.protein_names = p_names;           % Protein names

end




