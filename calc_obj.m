%% Function for evaluation of objective value used in parameter estimation
function obj = calc_obj(par_i,iso_data,iso_pars)
   
    % Save current search parameters into data structure
    iso_pars.current_par = par_i;
    
    % Load parameters from input structures
    n_comp = iso_pars.num_comp;
    pH = iso_data.pH;
    salt = iso_data.salt;
       
    % Check if initial objective has been determined
    if isfield(iso_pars,'initial_obj')
        % Declare initial objective value for iterations beyond start
        obj0 = iso_pars.initial_obj;
    else
        % Declare initial objective value for first iteration
        obj0 = 1;
    end
    
    % Perform transformation on parameters for multicomponent scenario 
    if n_comp > 1
        iso_pars = multicomp_par_transform(iso_pars,iso_data);
    end
    
    % Solve for fitted q values
    [iso_data,iso_pars] = solve_iso(iso_data,iso_pars);

    % Load variables from data structure
    q = iso_data.bound_conc;                       % Experimental q values
    q_fit = iso_data.bound_conc_eval;              % Fitted q values
    
    % Load variables from parameters structure
    n_pts = iso_pars.num_points;                   % Number of data points
    stoich_flag = iso_pars.stoich_flag;            % Stoichiometric flag
    solver_flag = iso_pars.solver_flag;            % Failed solver flag
    
    % Count number of data points per each pH and salt
    pH_unq = unique(pH);
    pts_ct = 0*salt;
    for i=1:numel(pH_unq)
        pH_idx = pH == pH_unq(i);
        salt_tmp = salt(pH_idx);
        salt_unq = unique(salt_tmp);       
        for j=1:numel(salt_unq)
            ct = numel(find(salt_tmp == salt_unq(j)));
            for k=1:numel(salt)
                if pH(k) == pH_unq(i) && salt(k) == salt_unq(j) 
                    pts_ct(k) = ct;
                end
            end
        end
    end
    
    % If all stoichiometric parameters are positive
    if stoich_flag == 0       
        % Calculate optimization objective value
        obj = sqrt(sum(abs((q_fit - q).^2)./pts_ct).*...
            (var(q)./sum(q.^2))/n_pts);
              
    % If any stoichiometric parameters are negative or solver failed
    elseif stoich_flag == 1 || solver_flag == 1
        % Set arbitarily high objective value
        obj = 2;
    end
    
    % Normalize objective values (to start from obj = 1)
    if obj0 ~= 1
         obj = mean(obj./obj0);
    end  
end   