%% Parsing of parameter ranges 
function iso_pars = parse_pars(iso_pars,iso_data)

    % Load parameter ranges from excel sheet (depending on resin mode)
    switch iso_pars.mode_oper
        % For cation exchange (CEX)
        case 'CEX'
            [~,~,par_sheet] = ...
                xlsread('Sheet for adsorption isotherm data entry',...
                'Parameter Range (CEX)'); % #ok

        % For anion exchange (AEX)
        case 'AEX'
            [~,~,par_sheet] = ...
                xlsread('Sheet for adsorption isotherm data entry',...
                'Parameter Range (AEX)'); % #ok

        % For multimodal cation exchange (MMCEX)
        case 'MMCEX'
            [~,~,par_sheet] = ...
                xlsread('Sheet for adsorption isotherm data entry',...
                'Parameter Range (MMCEX)'); % #ok

        % For multimodal anion exchange (MMAEX)
        case 'MMAEX'
            [~,~,par_sheet] = ...
                xlsread('Sheet for adsorption isotherm data entry',...
                'Parameter Range (MMAEX)'); % #ok

        % For invalid selection
        otherwise
            % Throw error if mode of operation does not match any option
            error("ERROR: Invalid mode of operation. " + ...
                "Please choose an option from the excel sheet.");
    end
    
    % Enter parameter sheet into structure 
    iso_pars.par_sheet = par_sheet;

    % Load variables from input structure
    pars = iso_pars.par_sheet;                   % Parameter sheet
    isotherm = iso_pars.isotherm;                % Name of isotherm
    n_comp = iso_pars.num_comp;                  % Number of components
    pH_dep = iso_pars.pH_dep;                    % pH dependent parameters
    n_cmp = iso_pars.num_comp;                   % Number of components
    
    % Mark if multicomponent and if parameters have been transformed
    if n_cmp > 1
        transform = 0;
        iso_pars.multicomp_transform = transform;
    end

    % Search for mode of operation specification
    i = 1;
    while ~strcmp(pars{i,1},isotherm)
        i = i + 1;
    end
    
    % Throw error if isotherm is not found in sheet
    if i > size(pars,1)
        error("ERROR: Invalid isotherm selection! " + ...
            "Please choose an option from the excel sheet.");
    end
    
    % Parse parameter names
    par_names = string(pars(i-1,3:end));
    
    % Parse parameter upper bound
    ub = cell2mat(pars(i,3:end));
    
    % Parse parameter lower bound
    lb = cell2mat(pars(i+1,3:end));
    
    % Average to obtain initial parameters
    pg = (ub+lb)/2;
    
    % Remove nans
    ub = ub(~isnan(ub));
    lb = lb(~isnan(lb));
    pg = pg(~isnan(pg));
    par_names = par_names(~ismissing(par_names));
    
    % Remove pH dependent parameters (if applicable)
    if strcmpi(pH_dep,'no')
        pH_idx = contains(par_names,'1');
        ub = ub(~pH_idx);
        lb = lb(~pH_idx);
        pg = pg(~pH_idx);
        par_names = par_names(~pH_idx);
    elseif strcmpi(pH_dep,'yes')
    else
        error("ERROR: Invalid selection of pH dependency. " + ...
            "Please select Yes or No for Add pH dependency.");
    end

    % Transform parameter set for multi-component
    ub = repmat(ub,[n_comp 1]);
    lb = repmat(lb,[n_comp 1]);
    pg = repmat(pg,[n_comp 1]);
    par_names = repmat(par_names,[n_comp 1]);
    
    % Store parameter info into parameters structure
    iso_pars.par_upper_bound = ub;             % Parameter set upper bound
    iso_pars.par_lower_bound = lb;             % Parameter set lower bound
    iso_pars.par_guess = pg;                   % Parameter intitial guess 
    iso_pars.par_names = par_names;            % Names of parameters

    % Perform transformation on parameters for multicomponent
    if n_cmp > 1
        iso_pars = multicomp_par_transform(iso_pars,iso_data);
    end
    
    % Remove parameter sheet from structure
    iso_pars = rmfield(iso_pars,'par_sheet');
    
end


