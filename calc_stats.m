%% Calculate statistics for fit accuracy and confidence intervals 
function [iso_data,iso_pars] = calc_stats(iso_data,iso_pars)

% Load variables from input data structure
q = iso_data.bound_conc;
q_fit = iso_data.bound_conc_eval;
x_load = iso_data.feed_frac;

% Load variables from input parameters structure
protein_names = iso_pars.protein_names;
pars_fit = iso_pars.current_par;
meas_err = iso_pars.meas_error;

% Declare beginning of statistics calculation
fprintf("\n------------------------CALCULATING STATISTICS" + ...
    "------------------------\n\n");

% Number of data points, components, and parameters
n_pts = size(q,1);
n_cmp = size(q,2);
n_par = numel(iso_pars.par_guess);

% Statistical parameters
sig_level = 0.05;                          % Significance level (95%)
deg_f = n_pts*n_cmp - n_par;               % Degrees of freedom
if deg_f <= 0
    deg_f = 1;
end
crit = tinv(1 - sig_level/2,deg_f);        % Critical level
meas_err = meas_err*x_load;

% Mean values for each component
q_mean = mean(q,1);                    

% Squared statistics for each component
sst = sum((q - q_mean).^2,1);             % Total sum of squares (SST)
sse = sum((q_fit - q).^2,1);              % Sum squared error (SSE)
mse = sse/n_pts;                          % Mean squared error (MSE)
rmse = sqrt(mse);                         % Root mean squared error (RMSE)
nrmse = rmse./q_mean;                     % Normalized RMSE (NRMSE)
r2 = 1 - (sse./sst);                      % R squared

% Absolute statistics for each component
sae = sum(abs(q_fit - q),1);              % Sum absolute error (SAE)
mae = sae/n_pts;                          % Mean absolute error (MAE)
nmae = mae./q_mean;                       % Normalized MAE

% Calculate pearson coefficients for each component
pearson = zeros(size(1,n_cmp));
for i=1:n_cmp
    pson_tmp = corrcoef(q(:,i),q_fit(:,i));
    pearson(i) = pson_tmp(2);
end

% Initialize AIC and BIC arrays
aic_r = zeros(1,n_cmp); aic_small = aic_r;
bic_r = aic_r; bic_small = bic_r;

for i=1:n_cmp
    % Prediction error
    pred_err = q(:,i) - q_fit(:,i);
    
    % Residual sum of squares
    rss = sum(pred_err.^2,'all')/n_pts;
    
    % Calculate Akaike information criteria (AIC) from RSS
    aic_r(i) = ...       % AIC from residual sum of squares (RSS)
        n_pts*log(rss) + 2*n_par;
    aic_small(i) = ... % AIC RSS for small sample size (n_pts/n_par < 40)
        aic_r(i) + 2*n_par*(n_par + 1)/(n_pts - n_par - 1);
    
    % Calculate Bayesian information criteria (BIC) from RSS
    bic_r(i) = ...       % BIC from residual sum of squares (RSS)
        n_pts*log(rss) + n_par*log(n_pts);
    bic_small(i) = ... % BIC RSS for small sample size (n_pts/n_par < 40)
        bic_r(i) + n_par*log(n_pts)*(n_par + 1)/(n_pts - n_par - 1);
end

% Collate statistics
stats_names = ["SSE" "MSE" "RMSE" "NRMSE" "R2" "SAE" "MAE" "NMAE" "Pearson" ...
    "AIC" "AIC small" "BIC" "BIC small"]';
stats_all = [sse; mse; rmse; nrmse; r2; sae; mae; nmae; pearson; ...
    aic_r; aic_small; bic_r; bic_small];

% Overall scores for multicomponent model
if n_cmp > 1
    % Mean values
    q_mean = mean(q,'all');                    

    % Squared statistics
    sst = sum((q - q_mean).^2,'all');     % Total sum of squares (SST)
    sse = sum((q_fit - q).^2,'all');      % Sum squared error (SSE)
    mse = sse/(n_pts*n_cmp);              % Mean squared error (MSE)
    rmse = sqrt(mse);                     % Root mean squared error (RMSE)
    nrmse = rmse./q_mean;                 % Normalized RMSE (NRMSE)
    r2 = 1 - (sse./sst);                  % R squared

    % Absolute statistics
    sae = sum(abs(q_fit - q),'all');      % Sum absolute error (SAE)
    mae = sae/(n_pts*n_cmp);              % Mean absolute error (MAE)
    nmae = mae./q_mean;                   % Normalized MAE

    % Calculate pearson coefficients
    pson_tmp = corrcoef(q,q_fit);
    pearson = pson_tmp(2);
    
    % Prediction error
    pred_err = q - q_fit;
    
    % Residual sum of squares
    rss = sum(pred_err.^2,'all')/(n_pts*n_cmp);
    
    % Calculate Akaike information criteria (AIC) from RSS
    aic_r = ...           % AIC from residual sum of squares (RSS)
        n_pts*n_cmp*log(rss) + 2*n_par;
    aic_small = ...       % AIC for small sample size (n_pts/n_par < 40)
        aic_r + 2*n_par*(n_par + 1)/(n_pts*n_cmp - n_par - 1);
    
    % Calculate Bayesian information criteria (BIC) from RSS
    bic_r = ...           % BIC from residual sum of squares (RSS)
        n_pts*n_cmp*log(rss) + n_par*log(n_pts*n_cmp);
    bic_small = ...       % BIC for small sample size (n_pts/n_par < 40)
        bic_r + n_par*log(n_pts*n_cmp)*(n_par + 1)/...
        (n_pts*n_cmp - n_par - 1);
    
    % Collate statistics
    stats_tmp = [sse; mse; rmse; nrmse; r2; sae; mae; nmae; pearson; ...
        aic_r; aic_small; bic_r; bic_small];
    stats_all = [stats_all stats_tmp];
    
    protein_names = [protein_names "All Components"];
end

% Final table of statistics
stats_table = ["Statistic Name" protein_names; stats_names string(stats_all)];

% Initialize parameters for CI calculation
pars_ptb = pars_fit;
ci = zeros(size(pars_fit));

% Calculate confidence intervals on model parameters via uncertainty method
for i=1:size(pars_fit,1)
    for j=1:size(pars_fit,2)
        % Apply perturbation to parameter (1%)
        pars_ptb(i,j) = pars_fit(i,j)*1.01;      
        iso_pars.current_par = pars_ptb;
        
        % Solve for q using perturbed parameter set
        [iso_data,iso_pars] = solve_iso(iso_data,iso_pars);        
        q_ptb = iso_data.bound_conc_eval;
        
        % Calculate jacobian matrix
        jac = (q_fit - q_ptb)/abs(pars_ptb(i,j) - pars_fit(i,j));
        
        % Calculate hessian matrix
        hessian = jac'*jac;
        
        % Clear last warning
        lastwarn('');
        
        % Calculate covariance matrix
        covm = (hessian\eye(size(hessian)))*meas_err(i)^2; 
       
        % Detect if inv threw 'nearly singular' warning.
        [~, warnId] = lastwarn;  

        % Compute covm using pseudoinverse if singular
        if strcmp(warnId,'MATLAB:nearlySingularMatrix')
            covm = pinv(hessian)*meas_err(i)^2;
        end
        
        % Parse diagonal values
        dgl = abs(diag(covm));
        
        % Calculate confidence interval and normalize value to percentage
        ci(i,j) = 100*crit*sqrt(dgl(i))/abs(pars_fit(i,j)); 
        
        % Reset parameter value
        pars_ptb(i,j) = pars_fit(i,j);
    end
end

% Save statistics and confidence intervals
iso_pars.current_par = pars_fit;
iso_pars.stats = stats_table;
iso_pars.confidence_int = ci;

end









