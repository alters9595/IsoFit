IsoFit Software 
Copyright (c) 2023, Scott Altern
All rights reserved.

Contact: alters9595@gmail.com for any questions, suggestions, or request for additional features and/or bug fixes.

----------------------------------------------------INSTRUCTIONS-------------------------------------------------------------

1) Open the excel file titled "Sheet for adsorption isotherm data entry". Enter the batch data into the "Adsorption Isotherm 
Data Sheet" tab and specify the pH, salt (mM), protein load concentration (mg/ml), protein flowthrough (mobile phase) 
concentration (mg/ml), resin volume (uL), load volume (uL) and mass fractions (%) of each component in the system. 
Enter the load material composition and the molecular weights of each component. Component mass fractions should add up to 
100%. It is not required to use all five components here. Do not leave any other fields blank. A example data set is 
provided in the sheet for the user. Be sure to save any changes.

2) Switch to the sheet titled "Fitting Options" in the same excel file. Enter all the fields in the sheet according 
to the possible selections and the purpose of each field. Do not leave any fields blank. Be sure to save any changes.

3) Switch to the sheet titled "Parameter Range (...)" in the same excel file. Enter the desired upper and lower bounds 
of each parameter corresponding to the desired isotherm model and resin mode (CEX, AEX, MMCEX, MMAEX). Parameter bounds
are provided to start but can be modified if desired. Be sure to save any changes.

4) Open up the isofit_module.m file in the IsoFit folder. Ensure that the modified excel sheet is present in the folder
along with all functions listed in the below section. Start the program by hitting start in the editor tab or by
pressing F5. It is not required to edit any code written in this script or any other functions. Look out for any 
errors that appear in the command window and correct any errors in the excel sheet as suggested by the error.
Please report unforseen errors to alters9595@gmail.com along with the excel sheet.

5) Monitor the parameter fitting process in the window showing the genetic algorithm progress will appear. The optimization 
can be prematurely stopped by pressing the STOP button, however it is recommended to allow the fitting to run to completion.
The expected time for fitting may be between 10-60 minutes depending on the number of isotherm parameters, number of data
points, optimization settings, and computational power (and whether parallel processing is enabled in the fitting options).

6) Once the program has finished running (PROGRAM ENDING in command window), open the Output folder in the IsoFit program
folder. Open the Isotherm Fits folder, will have (1) if this folder was already existing (to prevent overwrite). The plots
should be located in this folder in individual .png files with the names corresponding to the protein name, resin name,
and pH. Another export containing the parameters, confidence intervals, statistics, and isotherm data will be present
here in the same folder titled "Isotherm fit results".  

------------------------------------------------SUMMARY OF FUNCTIONS---------------------------------------------------------

isofit_module

Main script from which the program is run by hitting start in the editor tab or pressing F5. 
It is not required to edit the other functions.

-----------------------------------------------------------------------------------------------------------------------------

[iso_data,iso_pars] = gen_struct(data,opts)

This function generates two data structures (structs): iso_data and iso_pars. These structs are used throughout the 
remainder of functions to pass around isotherm data (in iso_data), fitting options (in iso_pars), 
and isotherm parameters (in iso_pars).

-----------------------------------------------------------------------------------------------------------------------------

[iso_data,iso_pars] = isodata_processing(iso_data,iso_pars)

This function processes the raw batch data provided in the Excel sheet "Sheet for adsorption isotherm data entry" under
the "Adsorption Isotherm Data Sheet" by calculating bound concentrations (q) of each component from the flowthrough (c_ft)
and load concentrations (c_load). Phase ratio is employed in this calculation and is corrected for hold-up volume using
the resin hold-up fraction (recommended value of 0.6).

-----------------------------------------------------------------------------------------------------------------------------

iso_pars = parse_pars(iso_pars,iso_data)

This function parses the parameter ranges specified in the Excel sheet for the desired isotherm formalism and
resin "Parameter Ranges (insert ligand mode, e.g. MMCEX)".

-----------------------------------------------------------------------------------------------------------------------------

iso_pars = multicomp_par_transform(iso_pars,iso_data)

This function transforms the provided isotherm parameter set into a form suitable for multicomponent isotherm. 
This may simply be reshaping the parameter array for all components or applying the multicomponent assumptions
described in the Excel sheet "Fitting Options".

-----------------------------------------------------------------------------------------------------------------------------

obj = calc_obj(pg,iso_data,iso_pars)

This function calculates the fitting objective value which is based on normalized root-mean squared error (NRMSE)
between the actual q values and the fitted q values. The objective value is normalized with respect to the starting point
value and should max at around unity. A negative weight is given to parameter sets that have negative stoichiometric
terms (nu and n) at any pH in the provided batch data.

-----------------------------------------------------------------------------------------------------------------------------

[iso_data,iso_pars] = run_opt(iso_data,iso_pars);

This function runs the parameter fitting routine (optimization) using a genetic algorithm with the specified 
optimization tolerance and maximum number of generations.

-----------------------------------------------------------------------------------------------------------------------------

[iso_data,iso_pars] = solve_iso(iso_data,iso_pars)

This function solves for q using the salt, pH, and c_ft for the selected isotherm formalism for single component 
and multicomponent forms.

-----------------------------------------------------------------------------------------------------------------------------

iso_pars = set_dir(iso_pars);

This function sets the final directory that the outputs will be saved to.

-----------------------------------------------------------------------------------------------------------------------------

[iso_data,iso_pars] = calc_stats(iso_data,iso_pars)

This function calculates the fitting statistics using the fitted q values calculated from the final isotherm parameter 
set obtained from the optimization routine. Confidence intervals are also calculated here using the Fischer information
criteria method. A perturbation of 1% is applied to each parameter value one at a time and a jacobian matrix is developed.
The covariance matrix is derived from here and used to calculate the confidence intervals (expressed in % of the
corresponding parameter value).

-----------------------------------------------------------------------------------------------------------------------------

plot_iso(iso_data,iso_pars)

This function generates plots of the resulting isotherm data using the parameters obtained from fitting. Fitted parameters
can be used to generate model curves (arbitrary input c_ft) or data points evaluated at the same c_ft from the experimental
data. Plots are generated at each pH with a q vs c curve per each salt concentration in both 2D and 3D with labels
according to the protein and resin used. Plots are exported in .png format to the specified export folder.

-----------------------------------------------------------------------------------------------------------------------------

export_results(iso_data,iso_pars)

This function exports the parameter values, confidence intervals, statistics, and isotherm data to an excel sheet
"Isotherm fit results" located in the output folder. The cells in each sheet are formatted automatically to select
appropriate column sizes.

-----------------------------------------------------------------------------------------------------------------------------

FROM TAL SHIR ON MATHWORKS FORUM, COPYRIGHT (c) 2008
xlsAutoFitCol(filename,sheetname,varargin)

This function changes the specified cell on an excel sheet to be autofit to its appropriate column width.

-----------------------------------------------------------------------------------------------------------------------------

FROM STEVE AMBROISE ON MATHWORKS FORUM, COPYRIGHT (c) 2009
y = linspaceNDim(d1,d2,n)

This function creates a linearly spaced matrix in N dimensions. It is analogous to linspace, except scaled to higher
dimensional space.

-----------------------------------------------------------------------------------------------------------------------------

FROM MATT BRUNNER ON MATHWORKS FORUM, COPYRIGHT (c) 2010
a1String = idx2A1(idx)

This function converts an array index to a excel cell string.

-----------------------------------------------------------------------------------------------------------------------------
