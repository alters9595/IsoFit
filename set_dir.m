%% Set directory for output of isotherm fitting results
function iso_pars = set_dir(iso_pars)

% Make directory for output folder
if ~exist('Output','dir')
    mkdir('Output')
end

% Make directory for isotherm fit figures
filename = iso_pars.export_folder; 
directory = fullfile(pwd,'Output',filename);

% Append (1) to directory if it already exists (to prevent overwrite)
i = 1;
while exist(directory,'dir') == 7
    filename = strcat(filename,' (1)');
    directory = fullfile(pwd,'Output',filename);
    i = i + 1;
end

% Set final output directory
mkdir(directory); iso_pars.export_dir = directory;

end

