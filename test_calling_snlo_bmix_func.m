% This is an example of how to call the bmix function in mlSNLO from the
% command-line or through a script.

% First we need to generate a set of inputs to be passed to snlo_bmix_func.
% The inputs must be in a data structure containing these fields (with
% values as described)
% bmix_selected_crystal     % an index of the cell array with all crystals listed. 
%   Needs to be an integer, chooses the crystal out of a list in crystal_listb.txt. Use 
%   the function find_selection_number_from_string to find the correct value, or
%   specify it manually.
% bmix_temperature          % scalar input with 
% bmix_wavelengths          % 3 element vector with wavelengths in nm (one must be 0)
% bmix_rotation_axis        % choose between 'X', 'Y', 'Z'
% bmix_type                 % choose between 'Type 1','Type 2','Type 3'

crystal_name = 'LBO';
crystal_list = importdata('crystal_listb.txt');
ind = find_selection_number_from_string(crystal_list, crystal_name);
inputs.bmix_selected_crystal = ind;     % Index of crystal from list in crystal_listb.txt
inputs.bmix_temperature = 300;          % temperature in K
inputs.bmix_wavelengths = [1500, 1500, 0];
inputs.bmix_rotation_axis = 'Y';
inputs.bmix_type = 'Type 2';

% Call bmix function, get a cell array containing all the internal functions it uses
fcns = snlo_bmix_func(inputs);

% The first element in the fcns cell array is the callback function to the
% 'Run' button; calling it lets you simulate pressing 'Run'.
run_fcn = fcns{1};

% You can simulate pressing the 'Deff' callback using the function which is
% sixth in the list
deff_fcn = fcns{6};

run_fcn(); % Equivalent to clicking the 'Run' button
deff_fcn([],[]); % Equivalent to clicking the 'deff' button. **Calling this function needs two input variables due to a bug; if none are passed to it then it crashes, though these variables aren't actually used


% The code below basically replicates what happens when you press
% deff button or run deff_fcn.
bmixdata = load('bmix.dat'); % load the data Bmix saves to file

% bmixdata will be a 2d array that has size n_phi x 12. 
% bmix(:,1) are phi angles, bmix(:,2) theta angles (in degrees)

% bmix(:,[3,4,5,6]) are deffs for four quadrants (phi1, phi4, phi2, phi3)
% vectors of phi angles (with 0.1 degree steps), where
%   0 deg < phi1 <  90 deg,
%  90 deg < phi2 < 180 deg,
%-180 deg < phi3 < -90 deg,
% -90 deg < phi4 <   0 deg,

% bmix(:,[7,8]) are walkoffs for lo, hi refractive indices
% bmix(:,[9,10,11]) are red, middle, and blue group velocity indices
% bmix(:,[12]) is temperature range

% In order to use these 

% concatenate four vectors of angles to make an array with size 4 in 2nd dimension;
% phi3, phi4, phi1, phi2 in in ascending order
phi_array = cat(2, -180+bmixdata(:,1), -bmixdata(:,1), ...
    bmixdata(:,1), 180-bmixdata(:,1));

% concatenate four vectors of deff values to make an array with size 4 in
% 2nd dimension; these are paired with the phi_array above
deff_array = cat(2, abs(bmixdata(:,6)), abs(bmixdata(:,4)), abs(bmixdata(:,3)), abs(bmixdata(:,5)));

fh = figure;
ax = axes('parent',fh);
ph = plot(ax, phi_array(:,1), deff_array(:,1), '-k', ...  % phi3 (-180 to -90 deg)
    phi_array(:,2), deff_array(:,2), '--k', ...           % phi4 (-90 to 0 deg)
    phi_array(:,3), deff_array(:,3), ':k', ...            % phi1 (0 to 90 deg)
    phi_array(:,4), deff_array(:,4), '-.k');              % phi2 (90 to 180 deg)
xlabel(ax,'\phi [deg]');
ylabel(ax,'d_{eff} [pm/V]');
legend(ax,'{\phi}_3', '{\phi}_4', '{\phi}_1', '{\phi}_2', 'location', 'best');
