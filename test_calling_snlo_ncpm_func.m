
% Example: using a batch script to call NCPM function

% This performs a simple optimization to find the temperature required to produce a pair
% of specified noncritical phase match wavelengths for a pair of target blue and red1
% waves (the wavelength of the red2 wave is determined by the other 2). Loop through
% several temperatures, keeping the contents of the ncpm.dat file which contains an N x 3
% array (1st column is blue wavelength, 2nd & 3rd columns are red wavelengths). Use these
% in a 1-dimensional interpolation to find the red1 wavelength that phase matches with the
% target blue wavelength.

target_blue = 1600; % desired blue wavelength in nm
target_red1 = 5700; % desired red1 wavelength in nm
temperatures = 270:10:370; % vector of temperatures we'll loop through to find the temperature that produces a phasematch closest to our target wavelengths

% First, make a datastructure with fieldnames as shown
%  For more info, see the ncpm tab in Help (launched from main menu panel)

% Note: NCPM requires an integer to be passed in the ncpm_currentCrystalPopupValue
% field; this integer is the line where the crystal specified resides in 
% crystals.txt. The MATLAB file find_selection_number_from_string.m will find this for
% you automatically if you give it two arguments: a cell array of the strings from
% each line of crystal_list.txt and the crystal in question.

clear inputs
inputs.ncpm_temperature = 300;
inputs.ncpm_currentCrystal = 'LNB_S';
inputs.ncpm_type = 'Type 1';
inputs.ncpm_propagation_axis = 'X';
inputs.ncpm_wavelengths = [1200, 200];

% Find index for ncpm_currentCrystalPopupValue
% crystal_list = importdata('crystal_list.txt');
% inputs.ncpm_currentCrystalPopupValue = find_selection_number_from_string(crystal_list, ncpm_crystal);
inputset = cell(size(temperatures));
outputs = cell(size(temperatures));

for K = 1:length(temperatures)
    inputset{K} = inputs;
    inputset{K}.ncpm_temperature = temperatures(K);
    fcns = snlo_ncpm_func(inputset{K});
    run_fcn = fcns{1}; % run button callback is the first element of cell array fcns
    run_fcn(); % call the run button function
    outputs{K} = importdata('ncpm.dat');
end

% Now to deal with the outputs cell array: 
%  Unfortunately, one complication here is the function ncpm combines two sets of phase
%  matching curves for each run in the output file ncpm.dat. Luckily, in our simple case
%  it is not strictly necessary to separate the two sets in order to just use a
%  1-dimensional interpolation away from the discontinuity.

red1_matches = zeros(size(temperatures));
for K = 1: length(temperatures)
    % interpolate to find what red1 wavelength phase matches with blue wave at target wavelength at each temperature
    red1_matches(K) = interp1(outputs{K}(:,1), outputs{K}(:,2), target_blue);
end

% plot the output of this:
fh = figure;
ax = axes(fh);
plot(ax,temperatures, red1_matches);
xlabel(ax,'Temperatures [K]');
ylabel(ax,'Red1 phase match wavelength');

% Find optimal temperature: interpolate to find temperature
T_target = interp1(red1_matches, temperatures, target_red1); % this time, use red1 phase match wavelength as x and temperature as y

% Indicate this optimal temperature on the figure
hold(ax,'on');                                      % add to current contents of figure with these next calls
plot(T_target, target_red1, '*', 'markersize', 10); % Plot the optimal temperature at the target red1 wavelength

% add title to indicate what crystal this optimization is for, and for what target wavelengths
title(ax,sprintf('Temperature to phase match %s with\n %g nm (blue) and %g nm (red) wavelengths',inputs.ncpm_currentCrystal, ...
    target_blue, target_red1),'fontweight','normal'); 

% also add a text label to the plot near the optimal temperature and target red1 wavelength
th = text(T_target-10, target_red1-50, sprintf('T = %.0f K', T_target)); 
