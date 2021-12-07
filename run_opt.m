%% Runs optimization routine to fit isotherm parameters
function [iso_data,iso_pars] = run_opt(iso_data,iso_pars)

% Load variables from parameters structure
ub = iso_pars.par_upper_bound;        % Upper bound on isotherm parameters
lb = iso_pars.par_lower_bound;        % Lower bound on isotherm parameters
pg = iso_pars.par_guess;              % Parameter guess

% Determine initial objective value
obj0 = calc_obj(pg,iso_data,iso_pars); iso_pars.initial_obj = obj0;

% Set options for parameter estimation 
opt_tol = iso_pars.opt_tolerance; max_gens = iso_pars.ga_max_gens;
options = optimoptions('ga','FunctionTolerance',opt_tol,'Display','iter',...
    'PlotFcn',@gaplotbestfun,'MaxGenerations',max_gens);

% Declare beginning of parameter fitting routine
fprintf("\n---------------------------FIT STARTING" + ...
    "---------------------------\n");

% Run parameter estimation routine using genetic algorithm
fit_pars = ga(@(par_i)calc_obj(par_i,iso_data,iso_pars),...
    numel(pg),[],[],[],[],lb,ub,[],options);
    
% Evaluate q values using fitted parameters
iso_pars.current_par = fit_pars;
if iso_pars.num_comp > 1
    iso_pars = multicomp_par_transform(iso_pars,iso_data);
end
[iso_data,iso_pars] = solve_iso(iso_data,iso_pars);
iso_data.bound_conc_eval = real(iso_data.bound_conc_eval);

% Declare end of parameter fitting routine
fprintf("\n---------------------------FIT COMPLETED" + ...
    "---------------------------\n");

end

