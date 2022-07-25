% calculate signal wavelength vs internal signal angle for several values of pump tilt angle.

% This example demonstrates how to set up the problem parameters for the opo angles model, call the SNLO
% opoangles function with the parameters, run the model, load and process the model outputs. We
% demonstrate several sets of similar input parameters, with a single input parameter varying from run to
% run. In this case, we'll vary the pump tilt angle, and examine its effect on the signal wavelength
% calculated as a function of internal signal angle.

calculate_curves = true;
if calculate_curves
    pump_tilt_list = 0:0.5:3.5; % generate vector of equally spaced pump tilts from 0 deg to 3.5 deg in increments of 0.5 deg
    
    % then define the full set of input parameters for opoangles; these correspond to the input boxes
    % visible when you run opoangles. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. calling the opoangles function
    % with a set of input parameters requires specifying a single input argument to the function
    % snlo_opoangles_func, where the input parameter is a data structure with the fieldnames specified
    % below, and values in the form specified below. 
    
    % In this case, the crystal is passed to the snlo_opoangles_func as its number in the list of crystals
    % (which changes if we add crystals). We use the function 'find_selection_number_from_string' to
    % determine which number corresponds to the string currentCrystalName.
    crystal_list = importdata('crystal_list.txt'); % load list of crystals
    % example 1 for opoangles (#3 from snlo examples button)
    inputs.opoangles_currentCrystalPopupName = 'BBO'; % string of crystal name
    inputs.opoangles_currentCrystalPopupValue = find_selection_number_from_string(crystal_list,inputs.opoangles_currentCrystalPopupName); % which position in the list of crystals this string appears at
    inputs.opoangles_temperature = 293; % temperature (in K)
    inputs.opoangles_lambda = [650,1200,532]; % pump, reddest, bluest wavelengths (in nm)
    inputs.opoangles_delta = 0;         % pump tilt angle (in degrees)
    inputs.opoangles_type = 'Type 1';   % mixing type ('Type 1' or 'Type 2')
    inputs.opoangles_plane = 'XZ';      % plane (ignored for uniaxial crystals)
    inputs.opoangles_signal_pol = 'o';  % signal polarization (ignored for type 1 mixing)
    
    
    %% call opoangles, run model, load output
    clear problem;
    outputs = cell(size(pump_tilt_list)); % pre-allocate a cell array to stick contents of model output file into
    for K = 1:length(pump_tilt_list) % loop over vector of temperatures 
        problem(K) = inputs;  % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).opoangles_delta = pump_tilt_list(K); % modify this problem's temperature to be an element of the vector of temperatures generated earlier
        fcn_handles = snlo_opoangles_func(problem(K)); % call the SNLO opoangles function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_fcn = fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        run_fcn(); % call the 'Run' button callback - equivalent to clicking the run button
        outputs{K} = importdata('opoangle.dat'); % load the model output from the file opoangle.dat, stick it in an element in the cell array 'output'
    end
    
end

%% plot output
fh = figure; % make new output figure
ax = axes('Parent',fh); % make new axes in the output figure
phs = zeros(length(pump_tilt_list),1); % pre-allocate array to store plot handles in
legend_strings = cell(size(pump_tilt_list)); % pre-allocate cell array to store legend label strings in
% define an array of line colors to use in plotting several curves on a single axis object
co = [0   0.447   0.741;
   0.850   0.325   0.098;
   0.929   0.694   0.125;
   0.494   0.184   0.556;
   0.466   0.674   0.188;
   0.301   0.745   0.933;
   0.635   0.078   0.184;
   0.0     0.0     0.0;];
for K = 1:length(pump_tilt_list) % loop over vector of temperatures 
    % plot wavelengths as a function of internal signal angle
    phs(K) = plot(ax,outputs{K}(:,1),outputs{K}(:,3));
    set(ax,'LineStyleOrder',{'-','--'},'ColorOrder',co);
    hold(ax,'on'); % tell axis the next curves plotted should be added to axis instead of replace existing curves
    legend_strings{K} = sprintf('pump tilt=%.1f deg.',pump_tilt_list(K)); % generate string for the legend indicating the temperature for this curve
end
ylh = ylabel(ax,'Wavelength [nm]');                 % label y axis
xlh = xlabel(ax,'Internal signal angle [degrees]'); % label x axis
th = title(ax,'BBO Type 1 - 532 nm pump','FontWeight','Normal'); % title the plot
set(ax,'FontName','Times New Roman');               % set the font used in this axis to time new roman

ylim(ax,[650,1200]);    % set y-axis limits
xlim(ax,[20.5,23.5]);   % set x-axis limits
drawnow;
lh  = legend(ax,legend_strings,'Location','northeast');  % add legend describing the curves to this axis