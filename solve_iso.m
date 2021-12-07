%% Solve for q with isotherm function based on selected model
function [iso_data,iso_pars] = solve_iso(iso_data,iso_pars)

    % Load variables from input structure (data)
    c_p = iso_data.flowthrough;       % Liquid phase protein conc. (mg/ml)
    c_s = iso_data.salt;              % Salt concentrations (mM)
    q = iso_data.bound_conc;          % Solid phase protein conc. (mg/ml) 
    pH = iso_data.pH;                 % pH values

    % Load variables from input structure (options)
    par_names = iso_pars.par_names;   % Parameter names
    par = iso_pars.current_par;       % Parameter values
    pH_dep = iso_pars.pH_dep;         % pH dependency check
    pH_ref = iso_pars.ref_pH;         % Reference pH
    lm = iso_pars.ionic_cap;          % Ionic capacity (M)
    Ec = iso_pars.porosity_Ec;        % Interstitial porosity
    Ep = iso_pars.porosity_Ep;        % Intraparticle porosity
    mw = iso_pars.mol_weight;         % Molecular weight (kDa)
    prl_p = iso_pars.use_parallel;    % Parallel processing check
    n_cmp = iso_pars.num_comp;        % Number of components
    
    % Unit conversions and calculations
    c_p = c_p./mw;                    % Convert to mM basis
    lm = lm*1000;                     % Convert to mM basis
    pH = pH - pH_ref;                 % Reference pH adjustment
    Et = Ec + Ep*(1 - Ec);            % Total porosity calculation
    q = q./(mw.*(1-Et));              % Convert to mM basis & solid volume
    n_pts = size(c_p,1);              % Number of data points
    
    % Transform pH and salt arrays for multi-component
    pH = repmat(pH,[1 n_cmp]);
    c_s = repmat(c_s,[1 n_cmp]);
    
    % Remove units for later parameter checks
    p_names = strrep(par_names,' (1/M)','');
    
    % Set parameter values if present in input
    for i=1:size(par,2)
        switch p_names(1,i)
            case 'Keq0';  keq0 = par(:,i); % Base equilibrium constant
            case 'Keq1';  keq1 = par(:,i); % pH-dep equilibrium constant
            case 'Nu0';   nu0 = par(:,i);  % Base characteristic charge
            case 'Nu1';   nu1 = par(:,i);  % pH-dep characteristic charge
            case 'N0';    n0 = par(:,i);   % Base char. hydrophobicity
            case 'N1';    n1 = par(:,i);   % pH-dep char. hydrophobicity
            case 'Sigma'; sh0 = par(:,i);  % Charged sites steric factor
            case 'S';     s0 = par(:,i);   % Hydrophobic steric factor
            case 'Ks';    ks = par(:,i);   % Mollerup salt-prot activity
            case 'Kp';    kp = par(:,i);   % Mollerup prot-prot activity
            case 'Beta0'; b0 = par(:,i);   % Base number of released water
            case 'Beta1'; b1 = par(:,i);   % Salt modification on beta
        end
    end
    
    % Set pH dependent parameters (keq, nu, n)
    if strcmpi(pH_dep,'yes')
        if exist('keq1','var') 
            % Initialize and set Keq[pH]
            keq = zeros(n_pts,n_cmp);
            for i=1:n_cmp
                keq(:,i) = keq0(i)*exp(pH(:,i).*keq1(i)); 
            end
        end
        if exist('nu1','var')
            % Initialize and set Nu[pH]
            nu = zeros(n_pts,n_cmp);
            for i=1:n_cmp
                nu(:,i) = nu0(i) + pH(:,i).*nu1(i); 
            end
        end
        if exist('n1','var')
            % Initialize and set N[pH]
            n = zeros(n_pts,n_cmp);
            for i=1:n_cmp
                n(:,i) = n0(i) + pH(:,i).*n1(i); 
            end
        end
        
    % Set non pH dependent parameters (keq, nu, n)
    elseif strcmpi(pH_dep,'no')
        % Initialize and set Keq
        if exist('keq0','var') 
            keq = ones(n_pts,n_cmp);
            for i=1:n_cmp
                keq(:,i) = keq(:,i).*keq0(i); 
            end
        end
        % Initialize and set Nu
        if exist('nu0','var')
            nu = ones(n_pts,n_cmp);
            for i=1:n_cmp
                nu(:,i) = nu(:,i).*nu0(i); 
            end
        end
        % Initialize and set N
        if exist('n0','var')
            n = ones(n_pts,n_cmp);
            for i=1:n_cmp
                n(:,i) = n(:,i).*n0(i); 
            end
        end
    end
    
    % Initialize and set Sigma
    if exist('sh0','var') 
        sh = ones(n_pts,n_cmp);
        for i=1:n_cmp
            sh(:,i) = sh(:,i).*sh0(i); 
        end
    end
    
    % Initialize and set S
    if exist('s0','var') 
        s = ones(n_pts,n_cmp);
        for i=1:n_cmp
            s(:,i) = s(:,i).*s0(i); 
        end
    end
    
    % Initialize and set B[salt]
    if exist('b0','var') 
        b = zeros(n_pts,n_cmp);
        for i=1:n_cmp
            b(:,i) = b0(i).*exp(b1(i).*c_s(:,i)/1000); 
        end
    end
    
    % Initialize and set Gamma[salt,protein]
    if exist('ks','var') && exist('kp','var')
        g = zeros(n_pts,n_cmp);
        for i=1:n_cmp
            g(:,i) = exp(ks(i).*c_s(:,i)/1000 + kp(i).*c_p(:,i)/1000);
        end
    end
    
    % Perform checks on stoichiometric parameters (must be non-negative)
    stoich_flag = 0;
    if exist('nu','var') && any(nu < 0,'all'); stoich_flag = 1; end
    if exist('n','var') && any(n < 0,'all');   stoich_flag = 1; end   
    
    % Define initial guess for function solver
    q0 = ceil(max(q,[],1))/2;
    
    % Initialize function handles
    isofxn = cell(n_pts,1);
    
    % Select function depending on isotherm model
    for i=1:n_pts
        switch iso_pars.isotherm
            case 'SMA'      
                isofxn{i,1} = @(q0)sma_solve(q0,lm,keq(i,:),nu(i,:),...
                    sh(i,:),c_p(i,:),c_s(i,:));
            case 'SMA Extended'     
                isofxn{i,1} = @(q0)sma_ext_solve(q0,lm,keq(i,:),nu(i,:),...
                    sh(i,:),g(i,:),c_p(i,:),c_s(i,:));
            case 'Ottens'           
                isofxn{i,1} = @(q0)ottens_solve(q0,lm,keq(i,:),nu(i,:),...
                    sh(i,:),n(i,:),s(i,:),c_p(i,:),c_s(i,:));
            case 'Ottens Extended'  
                isofxn{i,1} = @(q0)ottens_ext_solve(q0,lm,keq(i,:),nu(i,:),...
                    sh(i,:),n(i,:),s(i,:),g(i,:),c_p(i,:),c_s(i,:));
            case 'SMAHIC'        
                isofxn{i,1} = @(q0)smahic_solve(q0,lm,keq(i,:),nu(i,:),...
                    sh(i,:),n(i,:),s(i,:),b(i,:),c_p(i,:),c_s(i,:));
            case 'SMAHIC Extended' 
                isofxn{i,1} = @(q0)smahic_ext_solve(q0,lm,keq(i,:),nu(i,:),...
                    sh(i,:),n(i,:),s(i,:),g(i,:),b(i,:),c_p(i,:),c_s(i,:));
        end
    end
    
    % Initialize q for solver
    q_fit = zeros(n_pts,n_cmp);
    solver_flag = zeros(n_pts,1);
    
    % Function solver options
    options = optimoptions('fsolve','Display','off');
    
    % Solve for q values if stoichiometric parameters are non-negative
    if stoich_flag == 0
        if strcmpi(prl_p,'yes')
            % Solve for q values in parallel (faster, but higher CPU load)
            parfor i=1:n_pts
                % Try to solve for q, failure if solution is undefined
                try
                    q_fit(i,:) = fsolve(isofxn{i,1},q0,options);
                catch ME
                    % Set flag if solver failed
                    if contains(ME.message,"undefined values")
                        solver_flag(i) = 1;
                    else
                        throw(ME)
                    end
                end
            end
        elseif strcmpi(prl_p,'no')
            % Solve for q values in sequence
            for i=1:n_pts
                % Try to solve for q, failure if solution is undefined
                try
                    q_fit(i,:) = fsolve(isofxn{i,1},q0,options);
                catch ME
                    % Set flag if solver failed
                    if contains(ME.message,"undefined values")
                        solver_flag(i) = 1;
                    else
                        throw(ME)
                    end
                end
            end
        else 
            % Throw error if invalid parallel processing option
            error("ERROR: Invalid option selected for use parallel " + ...
                "processing. " + ...
                "Please select an option from the excel sheet."); 
        end
    end
    
    % Set flag if solver has failed at any data point
    if any(solver_flag == 1)
        solver_flag = 1;
    end
    
    % Fix pH dependent parameter names    
    if strcmp(iso_pars.pH_dep,'No') || strcmp(iso_pars.pH_dep,'no')
        iso_pars.par_names = erase(par_names,'0');
    end
    
    % Unit conversions
    q_fit = q_fit.*(mw.*(1-Et));                  % Convert to mg/ml basis
        
    % Store data into output structure
    iso_data.bound_conc_eval = q_fit;             % Fitted q values
    iso_pars.num_points = n_pts;                  % Number of data points
    iso_pars.stoich_flag = stoich_flag;           % Stoich parameter flag
    iso_pars.solver_flag = solver_flag;           % Flag for failed solver
end

%% Isotherm functions
% SMA isotherm
function q_solve = sma_solve(q0,lm,keq,nu,sh,c_p,c_s)
    q_ch = (lm - sum((nu+sh).*q0)).^nu;
    q_solve = c_p.*keq.*q_ch - c_s.^nu.*q0;
end

% Extended SMA isotherm
function q_solve = sma_ext_solve(q0,lm,keq,nu,sh,g,c_p,c_s)
    q_ch = (lm - sum((nu+sh).*q0)).^nu;
    q_solve = c_p.*g.*keq.*q_ch - c_s.^nu.*q0;
end

% Ottens isotherm
function q_solve = ottens_solve(q0,lm,keq,nu,sh,n,s,c_p,c_s)
    q_ch = (lm - sum((nu+sh).*q0)).^nu;
    q_hp = (lm - sum((n+s).*q0)).^n;
    q_solve = c_p.*keq.*q_ch.*q_hp - c_s.^nu.*q0;
end

% Extended ottens isotherm
function q_solve = ottens_ext_solve(q0,lm,keq,nu,sh,n,s,g,c_p,c_s)
    q_ch = (lm - sum((nu+sh)*q0)).^nu;
    q_hp = (lm - sum((n+s).*q0)).^n;
    q_solve = c_p.*g.*keq.*q_ch.*q_hp - c_s.^nu.*q0;
end

% SMA/HIC isotherm 
function q_solve = smahic_solve(q0,lm,keq,nu,sh,n,s,b,c_p,c_s)
    q_ch = (lm - sum((nu+sh).*q0)).^nu;
    q_hp = (lm - sum((n+s).*q0)).^n;
    q_solve = c_p.*keq.*q_ch.*q_hp - c_s.^nu.*q0.^(1+n.*b);
end

% Extended SMA/HIC isotherm 
function q_solve = smahic_ext_solve(q0,lm,keq,nu,sh,n,s,g,b,c_p,c_s)
    q_ch = (lm - sum((nu+sh).*q0)).^nu;
    q_hp = (lm - sum((n+s).*q0)).^n;
    q_solve = c_p.*g.*keq.*q_ch.*q_hp - c_s.^nu.*q0.^(1+n.*b);
end
