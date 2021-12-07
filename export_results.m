%% Exports isotherm parameters and fit statistics to excel
function export_results(iso_data,iso_pars)

% Initialize variables from data structure
c_p = iso_data.flowthrough;                  % Flowthrough concentrations
q = iso_data.bound_conc;                     % Bound concentrations
q_fit = iso_data.bound_conc_eval;            % Fitted bound concentrations
c_s = iso_data.salt;                         % Salt concentrations
pH = iso_data.pH;                            % pH values

% Initialize variables from parameters structure
pars = iso_pars.current_par;                 % Fitted isotherm parameters
ci = iso_pars.confidence_int;                % Confidence intervals
par_names = iso_pars.par_names;              % Parameter names
stats = iso_pars.stats;                      % Fit statistics
directory = iso_pars.export_dir;             % Export directory
protein_names = iso_pars.protein_names';     % Protein names
n_cmp = iso_pars.num_comp;                   % Number of components

% Initialize parameter table for output
par_table = strings(3*n_cmp,size(par_names,2)+1);

% Declare beginning of results export
fprintf("\n---------------------------EXPORTING RESULTS" + ...
    "--------------------------\n");

% Collate parameter names, values, and confidence intervals (CI) per comp
for i=1:n_cmp
    % Collect names
    name_out = [protein_names(i) par_names(i,:)];
    
    % Collect parameters
    par_out = ["Parameter Values" string(round(pars(i,:),5))];
    
    % Collect CIs
    ci_out = ["Confidence Intervals (%)" string(round(ci(i,:),5))];
    
    % Collate into table
    par_table(3*i-2:3*i,:) = [name_out; par_out; ci_out];
    
    % Clear temporary variables
    clear('name_out','par_out','ci_out');
end

% Generate names for output data
out_data_str = ["Cp (mg/ml) " + protein_names', ...
    "q Actual (mg/ml) " + protein_names', ...
    "q Fitted (mg/ml) " + protein_names', "Salt (mM)", "pH"];

% Collate isotherm data (Cp, q actual, q fit, salt, pH)
out_data = [out_data_str; ...
    round(c_p,5) round(q,5) round(real(q_fit),5) round(c_s,5) round(pH,5)];    

% Set file location for excel export
filename = 'Isotherm fit results.xlsx';
fileloc = fullfile(directory,filename);

% Suppress added new sheet warning
warning('off','MATLAB:xlswrite:AddSheet');

% Record instances of excel files open before export process
[~,tasks] = system('tasklist/fi "imagename eq Excel.exe"');
tasks = sscanf(tasks,'%s');

% Establish process IDs (PIDs) of open excel processes
exe = strfind(tasks,'.EXE');
console = strfind(tasks,'Console');
pid_before = zeros(1,size(exe,2));

% Populate list of excel PIDs
for i=1:size(exe,2)
   pid_before(1,i) = str2double(tasks((exe(i)+4):(console(i)-1)));    
end

% Write parameters table to excel sheet
sheetname = 'Parameters and CIs';
xlswrite(fileloc,par_table,sheetname);

% Perform autofit for parameters table
for i=1:size(par_table,2)
    % Find cell index of the longest string in parameters table
    [~,idx] = ...
        max(strlength(par_table(:,i)).*isnan(str2double(par_table(:,i))));
    cell_idx = strcat(idx2A1(i),num2str(idx)); 
    
    % Autofit column with respect to cell with longest string
    xlsAutoFitCol(fileloc,sheetname,cell_idx);
end

% Write statistics table to excel sheet
sheetname = 'Fit Statistics';
xlswrite(fileloc,stats,sheetname);

% Perform autofit for statistics table
for i=1:size(stats,2)
    % Find cell index of the longest string in parameters table
    [~,idx] = ...
        max(strlength(stats(:,i)).*isnan(str2double(stats(:,i))));
    cell_idx = strcat(idx2A1(i),num2str(idx)); 
    
    % Autofit column with respect to cell with longest string
    xlsAutoFitCol(fileloc,sheetname,cell_idx);
end

% Write output data to excel sheet
sheetname = 'Isotherm Data';
xlswrite(fileloc,out_data,sheetname);

% Perform autofit for output data table
for i=1:size(out_data,2)
    % Find cell index of the longest string in output data table
    [~,idx] = ...
        max(strlength(out_data(:,i)).*isnan(str2double(out_data(:,i))));
    cell_idx = strcat(idx2A1(i),num2str(idx)); 
    
    % Autofit column with respect to cell with longest string
    xlsAutoFitCol(fileloc,sheetname,cell_idx);
end

% Open Excel file
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(fileloc);

% Try to delete unnecessary "Sheet 1" in the excel file
try
  % Avoid error if sheet does not exist
  objExcel.ActiveWorkbook.Worksheets.Item('Sheet1').Delete;
catch
  % Do nothing
end

% Save, close, and clean up
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;

% Parse PIDs of newly generated processes
[~,tasks] = system('tasklist/fi "imagename eq Excel.exe"');
tasks = sscanf(tasks,'%s');
exe = strfind(tasks,'.EXE');
console = strfind(tasks,'Console');
pid_after = zeros(1,size(exe,2));

% Populate with these processes
for i=1:size(exe,2)
    pid_after(1,i) = str2double(tasks((exe(i)+4):(console(i)-1)));    
end

% We track excel file PIDs so that they can be closed without closing other
% excel files that are open. This has to be done because when excel is
% closed through the ActiveX server in MATLAB, the background, "ghost",
% processes are not closed. These will propagate and consume system memory,
% so they must be tracked and closed specifically.

% Close processes
if size(pid_after,2) > size(pid_before,2)
    command = ['taskkill /f /PID ',mat2str(pid_after(1,end))];
    system(command);
end

end




