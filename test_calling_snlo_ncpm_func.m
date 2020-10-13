
% Example: using a batch script to call NCPM function

% First, make a datastructure with fieldnames as shown
%  For more info, see the ncpm tab in Help (launched from main
%  menu panel)

% Note: NCPM requires an integer to be passed in the ncpm_currentCrystalPopupValue
% field; this integer is the line where the crystal specified resides in 
% crystals.txt. The MATLAB file find_selection_number_from_string.m will find this for
% you automatically if you give it two arguments: a cell array of the strings from
% each line of crystal_list.txt and the crystal in question.
ncpm_crystal = 'LNB_S'; 
clear inputs
inputs.ncpm_temperature = 300;
inputs.ncpm_currentCrystal = 'BBO';
inputs.ncpm_type = 'Type 1';
inputs.ncpm_propagation_axis = 'X';
inputs.ncpm_wavelengths = [1200, 200];

temperatures = 270:10:370;


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
    keyboard
end




