%% Function for plotting of experimental and fitted isotherm data
function plot_iso(iso_data,iso_pars)

    % Load data from input data structure
    pH = iso_data.pH;                 % pH values
    salt = iso_data.salt;             % Salt concentrations (mM)
    c = iso_data.flowthrough;         % Mobile phase concentration (mg/ml)
    q = iso_data.bound_conc;          % Bound concentrations (mg/ml)
    
    % Load data from input options structure
    isotherm = iso_pars.isotherm;     % Name of isotherm
    resin = iso_pars.resin_name;      % Name of resin
    comp = iso_pars.protein_names;    % Name of protein component
    n_comp = iso_pars.num_comp;       % Number of components
    lines_on = iso_pars.plot_lines;   % Setting for plot lines
    directory = iso_pars.export_dir;  % Export directory
    
    % Number of data points for plot
    n_pts = 500;
    
    % Sort based on pH
    [pH,pH_idx] = sort(pH);
    salt = salt(pH_idx);
    c = c(pH_idx,:);
    q = q(pH_idx,:);
    
    % Determine unique pH values 
    [pH_unq,pH_idx] = unique(pH);
    pH_idx = [pH_idx; numel(pH)+1];
    n_pH = numel(pH_unq);
    
    % Initialize temporary variables
    salt_tmp = 0*salt; salt_idx_tmp = 0*salt;
    
    % Sort salt values for each unique pH value
    for i=1:n_pH
        % Set temporary index for pH changes
        pH_idx_tmp = pH_idx(i):pH_idx(i+1)-1;
        
        % Sort salt concentrations
        [salt_tmp1,salt_idx_tmp1] = sort(salt(pH_idx_tmp));
        
        % Generate temporary salt array and index
        salt_idx_tmp1 = salt_idx_tmp1 + pH_idx(i) - 1;
        salt_tmp(pH_idx_tmp) = salt_tmp1;
        salt_idx_tmp(pH_idx_tmp) = salt_idx_tmp1;
    end

    % Adjust c and q values to match sorted salt
    salt = salt_tmp;
    c = c(salt_idx_tmp,:);
    q = q(salt_idx_tmp,:);
 
    % Create markers for scatter plot
    mkr = {'o','s','d','^','v','h','p','*','x',...
        'o','s','d','^','v','h','p','*','x'};
    
    % Declare beginning of plot generation
    fprintf("\n---------------------------GENERATING PLOTS" + ...
    "---------------------------\n");

    % Go through plot scheme for each component at each pH
    for n=1:n_comp
        for i=1:n_pH
            % Initialize temporary variables
            c_tmp = c(pH_idx(i):pH_idx(i+1)-1,:);
            salt_tmp = salt(pH_idx(i):pH_idx(i+1)-1);
            q_tmp = q(pH_idx(i):pH_idx(i+1)-1,:);

            % Initialize figure for 3D plot
            figure(i+(n-1)*n_pH+1);

            % Initialize figure for 2D plot
            figure(i+(n-1)*n_pH+100);

            % Determine unique salt values at each pH
            [salt_unq,salt_idx] = unique(salt_tmp);
            salt_idx = [salt_idx; numel(salt_tmp)+1];

            % Develop c array for model-generated curves with lines
            if n_comp > 1 && strcmpi(lines_on,'yes')
                c_plot_tmp = ...
                    linspaceNDim(zeros(1,n_comp),max(c_tmp,[],1),n_pts)';
            elseif n_comp == 1 && strcmpi(lines_on,'yes')
                c_plot_tmp = linspace(0,max(c_tmp),n_pts)';
            end
      
            for j=1:numel(salt_unq)
                % Develop arrays for experimental data
                c_exp = c_tmp(salt_idx(j):salt_idx(j+1)-1,n);
                q_exp = real(q_tmp(salt_idx(j):salt_idx(j+1)-1,n));
                salt_exp = salt_unq(j).*ones(size(c_exp,1),1);
                
                % Develop c array for model-generated curves without lines
                if strcmpi(lines_on,'no')
                    c_plot_tmp = c_tmp(salt_idx(j):salt_idx(j+1)-1,:);
                    n_pts = size(c_plot_tmp,1);
                end
                
                % Develop pH array for model-generated curves
                pH_plot = pH_unq(i).*ones(n_pts,1);

                % Develop salt array for model-generated curves
                salt_plot = salt_unq(j).*ones(n_pts,1);        

                % Prepare plot inputs for q evaluation
                iso_data.flowthrough = c_plot_tmp;
                iso_data.salt = salt_plot; 
                iso_data.pH = pH_plot;

                % Determine q values for plot
                [iso_data,iso_pars] = solve_iso(iso_data,iso_pars);
                q_plot_tmp = iso_data.bound_conc_eval;
                
                % Set final plot inputs
                c_plot = c_plot_tmp(:,n); q_plot = real(q_plot_tmp(:,n));
                
                % Add line to 3D plot for simulated curve
                figure(i+(n-1)*n_pH+1);
                if strcmpi(lines_on,'yes')
                    plot3(c_plot,salt_plot,q_plot,'LineStyle','-',...
                        'LineWidth',2); 
                elseif strcmpi(lines_on,'no')
                    scatter3(c_plot,salt_plot,q_plot,'Marker',mkr{j},...
                        'LineWidth',1);
                end
                hold on;

                % Get info on individual curves
                h = get(gca,'Children');
                if strcmpi(lines_on,'yes')
                    clr = h(1).Color;
                elseif strcmpi(lines_on,'no')
                    clr = h(1).CData;
                end

                % Add to 3D scatter for experimental data points
                scatter3(c_exp,salt_exp,q_exp,'Marker',mkr{j},'LineWidth',1,...
                    'MarkerEdgeColor','k','MarkerFaceColor',clr); 
                hold on;

                % Add line to 2D plot for simulated curves
                figure(i+(n-1)*n_pH+100);
                if strcmpi(lines_on,'yes')
                    plot(c_plot,q_plot,'LineStyle','-','LineWidth',1.5);
                elseif strcmpi(lines_on,'no')
                    scatter(c_plot,q_plot,'Marker',mkr{j},'LineWidth',1);
                end 
                hold on;

                % Add to 2D scatter for experimental data points
                scatter(c_exp,q_exp,'Marker',mkr{j},'LineWidth',1,...
                    'MarkerEdgeColor','k','MarkerFaceColor',clr); 
                hold on;

                % Generate legend for plots (referring to salt conc.)
                lgd(2*j-1:2*j) = string([salt_unq(j); salt_unq(j)]) + " mM";
            end

            % Set name for pH condition
            pH = num2str(pH_unq(i));

            % Call 3D plot
            figure(i+(n-1)*n_pH+1);

            % Reverse direction for 3D
            set(gca,'ydir','reverse');

            % Set title
            fig_title = strcat(comp(n),"  ",resin,"  ",isotherm,...
                " Isotherm  ","pH ",pH);
            title(fig_title);

            % Set axes
            xlabel('c_p (mg/ml)');      % Protein concentration
            ylabel('c_s (mM)');         % Salt concentration
            zlabel('q_p (mg/ml)');      % Bound concentration

            % Set axes limits
            xlim([0 inf])
            ylim([min(salt_tmp) max(salt_tmp)])
            zlim([0 inf])

            % Set legend
            legend(lgd)

            % Grid on
            grid on

            % Remove legend entries for lines
            p = get(gca,'Children');
            for k=1:numel(p)
                if ~mod(k,2)
                    set(get(get(p(k),'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                end
            end
            
            % Export 3D plot
            fig_name = strcat(fig_title,' 3D.png');
            saveas(gcf,fullfile(directory,fig_name));
            delete(gcf);
            
            % Call 2D plot
            figure(i+(n-1)*n_pH+100);

            % Set title
            title(fig_title);

            % Set axes
            xlabel('c_p (mg/ml)');      % Protein concentration
            ylabel('q_p (mg/ml)');      % Bound concentration

            % Set axes limits
            xlim([0 inf])
            ylim([0 inf])

            % Set legend
            legend(lgd,'location','bestoutside')

            % Grid on
            grid on

            % Remove legend entries for lines
            s = get(gca,'Children');
            for k=1:numel(s)
                if ~mod(k,2)
                    set(get(get(s(k),'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                end
            end
            
            % Export 2D plot
            fig_name = strcat(fig_title,' 2D.png');
            saveas(gcf,fullfile(directory,fig_name));
            delete(gcf);
                       
            clear('lgd')
        end
    end
    
end
