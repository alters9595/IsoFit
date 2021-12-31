%% Transformation of parameter set for multicomponent formalisms
function iso_pars = multicomp_par_transform(iso_pars,iso_data)

% Load from input parameters structure
pH_asmn = iso_pars.pH_asmn;               % pH assumption for multicomp.
size_asmn = iso_pars.size_asmn;           % Size assumption for multicomp.
ks_asmn = iso_pars.ks_asmn;               % Ks assumption for multicomp.
kp_asmn = iso_pars.kp_asmn;               % Kp assumption for multicomp.
beta_asmn = iso_pars.beta_asmn;           % Beta assumption for multicomp.
pH_dep = iso_pars.pH_dep;                 % pH dependent parameters
transform = iso_pars.multicomp_transform; % Flag for multicomp. transform
mw = iso_pars.mol_weight;                 % Molecular weight
ub = iso_pars.par_upper_bound;            % Parameter set upper bound
lb = iso_pars.par_lower_bound;            % Parameter set lower bound
pg = iso_pars.par_guess;                  % Parameter intitial guess 
p_names = iso_pars.par_names;             % Names of parameters
n_cmp = iso_pars.num_comp;                % Number of components

% Load from input data structure
feed_comp = iso_data.feed_frac;           % Composition of feedstock

% Load current isotherm parameters for solver
if isfield(iso_pars,'current_par')
    par_0 = iso_pars.current_par;
end

% Load deletion matrix post-transformations
if isfield(iso_pars,'deletion_mat')
    ptfm = iso_pars.deletion_mat;
end

% Generate scaling factor for size assumption, calculated from relative 
% surface area using hydrodynamic radii (Rh) obtained from molecular weight
% (MW) correlation. Surface area here is relative to the component in 
% highest concentration in the load material. Correlation was obtained
% from a power law fit using Rh and MW data for 25 proteins between
% 5 and 1000 kDa.
r_h = 0.7429.*mw.^0.3599;
[~,loc] = max(feed_comp);
size_scale = r_h.^2./r_h(loc)^2;
    
% Handling error for invalid pH dependence options selections
if strcmpi(pH_dep,'no') && strcmpi(pH_asmn,'yes')
    error("ERROR: Setting for pH dependent parameters must be enabled" + ...
        " if multicomponent pH assumption is enabled! Please disable both" + ...
        " settings if pH dependence is not desired.");
end

if transform == 1
    % Load original parameter names
    p_names0 = iso_pars.par_names_original;
    
    % Initialize parameter matrix
    par_i = zeros(size(ptfm));
    p_names2 = strings(size(ptfm));
    
    % Find where parameter names match between transformed and original set
    % Start with parameters with no applied assumptions
    p_asmn = zeros(size(ptfm));
    for i=1:size(p_names0,1)
        for j=1:size(p_names0,2)
            for k=1:numel(p_names)
                % Check if there is a match
                if strcmpi(p_names0(i,j),p_names(k))
                    % Check if pH assumption was selected
                    if strcmpi(pH_asmn,'yes')
                        % Check if parameter is pH dependent/
                        if strcmpi(p_names(k),'Keq1') || ...
                            strcmpi(p_names(k),'Nu1') || ...
                            strcmpi(p_names(k),'N1')
                            % Assign pH dependent parameters
                            if all(par_i(:,j) == 0)
                                par_i(:,j) = par_0(k);
                                p_names2(:,j) = p_names(k);
                                p_asmn(:,j) = 1;
                            end
                        end
                    end
                    % Check if Ks assumption was selected
                    if strcmpi(ks_asmn,'yes') && contains(p_names(k),'Ks')
                        % Assign Ks parameters
                        if all(par_i(:,j) == 0)
                            par_i(:,j) = par_0(k);
                            p_names2(:,j) = p_names(k);
                            p_asmn(:,j) = 1;
                        end
                    end
                    % Check if Kp assumption was selected
                    if strcmpi(kp_asmn,'yes') && contains(p_names(k),'Kp')
                        % Assign Kp parameters
                        if all(par_i(:,j) == 0)
                            par_i(:,j) = par_0(k);
                            p_names2(:,j) = p_names(k);
                            p_asmn(:,j) = 1;
                        end
                    end
                    % Check if Beta assumption was selected
                    if strcmpi(beta_asmn,'yes')
                        if strcmpi(p_names(k),'Beta0') || ...
                            contains(p_names(k),'Beta1')
                            % Assign beta parameters
                            if all(par_i(:,j) == 0)
                                par_i(:,j) = par_0(k);
                                p_names2(:,j) = p_names(k);
                                p_asmn(:,j) = 1;
                            end
                        end
                    end
                    % Check if size assumption was selected
                    if strcmpi(size_asmn,'yes')
                        if strcmpi(p_names(k),'Sigma') || ...
                            strcmpi(p_names(k),'S')   
                            % Assign steric parameters
                            if all(par_i(:,j) == 0)
                                par_i(:,j) = size_scale.*par_0(k);
                                p_names2(:,j) = p_names(k);
                                p_asmn(:,j) = 1;
                            end
                        end
                    end  
                end
            end
        end
    end
         
    % Find remaining matches for parameters without assumptions
    for i=1:size(p_names0,1)
        for j=1:size(p_names0,2)
            % Find indices of parameters to match
            n_asmn = length(p_asmn(p_asmn(i,1:j) == 1));
            if p_asmn(i,j) == 0
                p_idx = i + (j-1)*n_cmp + n_asmn*(1 - n_cmp);
            end
            % Proceed to assignment if index is matching
            for k=1:numel(p_names)  
                if strcmpi(p_names0(i,j),p_names(k)) && par_i(i,j) == 0
                    if p_idx == k
                        % Assign remaining parameters
                        par_i(i,j) = par_0(k);
                        p_names2(i,j) = p_names(k);
                    end
                end
            end
        end
    end    
    
    % Save transformed parameter set
    iso_pars.current_par = par_i;
    iso_pars.par_names = p_names2;
    
elseif transform == 0
    % Save original parameter names
    iso_pars.par_names_original = p_names;
    
    % Initialize transformation matrix
    ptfm = ones(size(p_names));
    
    % Handle selection for pH assumption for multicomponent
    if ~strcmpi(pH_asmn,'yes') && ~strcmpi(pH_asmn,'no') 
        error("ERROR: Please select yes or no for multicomponent" + ...
            " assumption on pH dependent isotherm parameters.")
    end

    % Handle selection for size assumption for multicomponent
    if ~strcmpi(size_asmn,'yes') && ~strcmpi(size_asmn,'no')
        error("ERROR: Please select yes or no for multicomponent" + ...
            " assumption on size (steric) isotherm parameters.")
    end

    % Handle selection for Ks assumption for multicomponent
    if ~strcmpi(ks_asmn,'yes') && ~strcmpi(ks_asmn,'no')
        error("ERROR: Please select yes or no for multicomponent" + ...
            " assumption on Mollerup Ks isotherm parameters.")
    end

    % Handle selection for Kp assumption for multicomponent
    if ~strcmpi(kp_asmn,'yes') && ~strcmpi(kp_asmn,'no')
        error("ERROR: Please select yes or no for multicomponent" + ...
            " assumption on Mollerup Kp isotherm parameters.")
    end

    % Handle selection for Beta assumption for multicomponent
    if ~strcmpi(beta_asmn,'yes') && ~strcmpi(beta_asmn,'no')
        error("ERROR: Please select yes or no for multicomponent" + ...
            " assumption on Beta isotherm parameters.")
    end
    
    % Mark parameters for deletion depending on selected assumption
    for i=1:size(ptfm,2)
        if strcmpi(pH_asmn,'yes')
        % Mark Keq1, Nu1, N1 for deletions
            if strcmpi(p_names(1,i),'Keq1') || ...
                    strcmpi(p_names(1,i),'Nu1') || ...
                    strcmpi(p_names(1,i),'N1')
                ptfm(2:end,i) = 0;
            end
        end
        % Mark Ks for deletions
        if strcmpi(p_names(1,i),'Ks') && strcmpi(ks_asmn,'yes')
            ptfm(2:end,i) = 0;
        end
        % Mark Kp for deletions
        if strcmpi(p_names(1,i),'Kp') && strcmpi(kp_asmn,'yes')
            ptfm(2:end,i) = 0;
        end
        % Mark Beta0 and Beta1 for deletions
        if strcmpi(beta_asmn,'yes')
            if strcmpi(p_names(1,i),'Beta0') || ...
                    strcmpi(p_names(1,i),'Beta1')
                ptfm(2:end,i) = 0;
            end
        end
        % Mark Sigma and S for deletions
        if strcmpi(size_asmn,'yes')
            if strcmpi(p_names(1,i),'Sigma') || ...
                    strcmpi(p_names(1,i),'S')
                ptfm(2:end,i) = 0;
            end
        end
    end
    
    % Perform deletions
    ub = ub(ptfm == 1);
    lb = lb(ptfm == 1);
    pg = pg(ptfm == 1);
    p_names = p_names(ptfm == 1);
    
    % Store parameter info into parameters structure
    iso_pars.par_upper_bound = ub;             % Parameter set upper bound
    iso_pars.par_lower_bound = lb;             % Parameter set lower bound
    iso_pars.par_guess = pg;                   % Parameter intitial guess 
    iso_pars.par_names = p_names;              % Names of parameters
    iso_pars.deletion_mat = ptfm;              % Deletion matrix
    
    % Mark that transformation has been performed
    transform = 1;
end

% Save transformation state
iso_pars.multicomp_transform = transform; 

end

