% calculate quasiphasematching period as function of wavelength for several temperatures.

% This example demonstrates how to set up the problem parameters for the quasiphasematching model (qpm),
% call the SNLO qpm function with the parameters, run the model, load and process the model
% outputs. We demonstrate several sets of similar input parameters, with a single input parameter varying
% from run to run. In this case, we'll vary the temperature, and examine its effect on the qpm period
% calculated as a function of wavelength.

calculate_curves = true;
if calculate_curves
    temperature_list = 250:25:400; % generate vector of equally spaced temperature from 250 K to 400 K in increments of 25 K
    
    % then define the full set of input parameters for qpm; these correspond to the input boxes
    % visible when you run qpm. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. calling the qpm function
    % with a set of input parameters requires specifying a single input argument to the function
    % snlo_qpm_func, where the input parameter is a data structure with the fieldnames specified
    % below, and values in the form specified below. 
    
    % In this case, the crystal is passed to the snlo_qpm_func as its number in the list of crystals
    % (which changes if we add crystals). We use the function 'find_selection_number_from_string' to
    % determine which number corresponds to the string currentCrystalName.
    qpm_crystal_list = importdata('crystal_list_qpm.txt'); % load list of crystals
    inputs.qpm_temperature = 250;               % temperature (in K)
    inputs.qpm_wavelengths = [752,900,3000];    % pump, bluest, reddest wavelengths [nm]
    inputs.currentCrystalPopupName = 'KNBO3';   % specify one of the strings which are contained in qpm crystal list (crystal_list_qpm.txt)
    inputs.qpm_currentCrystalPopupValue = find_selection_number_from_string(qpm_crystal_list,inputs.currentCrystalPopupName); % find the position in the list of crystals which corresponds to the name of the crystal you specified
    inputs.qpm_currentPolValue = 4;             % position in list of polarizations which corresponds to the one you want to use (run qpm and select your crystal to find the list of polarizations, and which set this value to the position in the list of polarizations)
    inputs.qpm_temptune_temp_range = [273,450]; % temperature range for temperature tuning [K] - not used in this example
    inputs.qpm_temptune_period = 10;            % period for temperature tuning [um] - not used in this example
    inputs.qpm_pumptune_wavelength_range = [796,810]; % wavelength range for pump tuning [nm] - not used in this example
    inputs.qpm_pumptune_period = 28;            % period for pump tuning [um] - not used in this example
    
    
    %% call qpm, run model, load output
    clear problem;
    outputs = cell(size(temperature_list)); % pre-allocate a cell array to stick contents of model output file into
    for K = 1:length(temperature_list) % loop over vector of temperatures 
        problem(K) = inputs;  % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).qpm_temperature = temperature_list(K); % modify this problem's temperature to be an element of the vector of temperatures generated earlier
        fcn_handles = snlo_qpm_func(problem(K)); % call the SNLO qpm function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_fcn = fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        run_fcn(); % call the 'Run' button callback - equivalent to clicking the run button
        outputs{K} = importdata('qpm.dat'); % load the model output from the file qpm.dat, stick it in an element in the cell array 'output'
    end
    
end

%% plot output
fh = figure; % make new output figure
ax = axes('Parent',fh); % make new axes in the output figure
phs = zeros(length(temperature_list),1); % pre-allocate array to store plot handles in
legend_strings = cell(size(temperature_list)); % pre-allocate cell array to store legend label strings in
for K = 1:length(temperature_list) % loop over vector of temperatures 
    % plot wavelengths as a function of period: join the first and second columns (with second column
    % order reversed), use it as y data; join third and a replica of reverse order third column, used it
    % as x data
    phs(K) = plot(abs([outputs{K}(:,3);outputs{K}(end:-1:1,3)]),... 
        [outputs{K}(:,1);outputs{K}(end:-1:1,2)]);
    hold(ax,'on'); % tell axis the next curves plotted should be added to axis instead of replace existing curves
    legend_strings{K} = sprintf('T=%i',temperature_list(K)); % generate string for the legend indicating the temperature for this curve
end
ylh = ylabel(ax,'Wavelength [nm]'); % label y axis
xlh = xlabel(ax,'Period [{\mu}m]'); % label x axis
lh  = legend(ax,legend_strings,'Location','Best'); % add legend describing the curves to this axis
set(ax,'FontName','Times New Roman'); % set the font used in this axis to time new roman