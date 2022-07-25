% this is an example of how to call the function snlo_ref_ind_func from a
% MATLAB script. first, we'll define a set of inputs the function needs to
% populate the input edit boxes. then, we'll call the snlo_ref_ind_func
% with those inputs and get back a list of functions internal to the GUI
% which we can use to simulate pressing the 'Run' button. We'll
% simulate pressing the run button and read the function's resulting
% output. However, since snlo_ref_ind_func doesn't save the results of its
% calculations to a file, but it does display those results in a text box
% in the GUI, we'll try to read the contents of the output text box and
% save the string as the variable output_string which will be a cell array
% of differently sized strings.

% first define some values we'll use in specifying set of inputs for
% calling snlo_ref_ind_func
crystal_name = 'BBO';           % name of crystal. string. must match an entry in the file crystal_list.txt.
temperature = 350;              % temperature in kelvin for use in temperature dependent Sellmeier calculations.
propagation_angle_theta = 30;   % (in degrees) direction you want to find the refractive index for.
propagation_angle_phi = 30;     % (in degrees) direction you want to find the refractive index for; irrelevant for uniaxial crystals, but must be specified anyway
wavelength = 632.8;             % wavelength in nanometers

% find which element in the list of crystals matches the crystal named
% above, which is what snlo_ref_ind_func needs to properly select the
% crystal in question
ref_ind_crystal_list = importdata('crystal_list.txt');
crystal_list_number = find_selection_number_from_string(ref_ind_crystal_list,crystal_name); % find the position in the list of crystals which corresponds to the name of the crystal you specified

% construct data structure to use as input when calling snlo_ref_ind_func
input_set.ref_ind_selected_crystal = crystal_list_number;
input_set.ref_ind_temperature = temperature;
input_set.ref_ind_theta = propagation_angle_theta;
input_set.ref_ind_phi = propagation_angle_phi;
input_set.ref_ind_wavelength = wavelength;

% call snlo_ref_ind_func with set of inputs, get list of function handles
% for run, print buttons and other functions we don't care about here
function_handles = snlo_ref_ind_func(input_set);
run_function = function_handles{1};

% simulate pressing the 'Run' button
run_function();

% find snlo_ref_ind_func's window which we'll try to read the output text
% from
ref_ind_fig = findobj('Tag','RefInd');
% get the list of elements in this figure
ref_ind_fig_children = get(ref_ind_fig,'children');

% try to find which item in the list of elements is the output text box

% for every element in the list of the children of the figure, do this
for K = 1:length(ref_ind_fig_children)
    try
        % if we can, get the 'String' property of the object we're looking
        % at. this will work for the figure's child objects which are the
        % pushbuttons, and the output text edit box.
        temp_output_string = get(ref_ind_fig_children(K),'String');
        % let's try to isolate only the output text edit box and assign its
        % contents to the variable output_string. we know that the output
        % text box from a calculation with valid inputs will always be a
        % cell array of differently sized strings, but that sometimes the
        % inputs contain a wavelength out of the range of Sellmeier
        % validity, in which case the output box will contain only a string
        % reporting the wavelength is out of range. however, the pushbutton
        % objects also contain only one string in their 'String' property,
        % but its always set to 'Run' or 'Print', and we don't want those
        % strings saved.
        
        if iscell(temp_output_string)
%             if ~(strcmp(temp_output_string{1},'Run') || strcmp(temp_output_string{1},'Print'))
            output_string = temp_output_string; % if we've found the output text box, save its contents as output_string
%             end
        else
            if ~(strcmp(temp_output_string,'Run') || strcmp(temp_output_string,'Print'))
                output_string = temp_output_string;
            end
        end
    catch
        % do nothing if error encountered, which will happen if you can't
        % get the 'String' property of the object we're looking at. this
        % will happen for the children of the GUI's figure which are
        % uipanel objects used to organize input boxes.
    end
end