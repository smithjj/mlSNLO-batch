% In this example batch file, we'll use the Focus function in mlSNLO to
% generate a curve of beam size vs distance from the beam waist. This isn't
% the most intelligent way to calculate this information, but it serves as
% an example of how to construct a set of inputs for Focus and how to read
% the function's output.

clear inputs

inputs.focus_wavelength = 1000;         % wavelength (in nm)
inputs.focus_ref_ind = 1.6;             % refractive index
inputs.focus_waist = 10e-3;             % waist (in mm)
inputs.focus_dist = 0;                  % distance to crystal input face (in mm)
inputs.focus_face_to_focus = 500e-3;    % distance from input face to focus waist in crystal (in mm)
% initially the value of focus_face_to_focus is unimportant since we'll set it in the for
% loop below

face_to_focus_list = linspace(0.05,1,20); % in mm
waist_list = 0.05:-0.01:0.01;

% Pre allocate some arrays to store results of Focus calculations
rayleigh_range_mm   = zeros(length(face_to_focus_list),length(waist_list));
size_mm             = zeros(length(face_to_focus_list),length(waist_list));
curvature_mm        = zeros(length(face_to_focus_list),length(waist_list));
farfield_angle_mrad = zeros(length(face_to_focus_list),length(waist_list));

clear problem
for J = 1:length(waist_list)
    for K = 1:length(face_to_focus_list)
        problem(K,J) = inputs;
        problem(K,J).focus_dist = face_to_focus_list(K);
        problem(K,J).focus_waist = waist_list(J);
        funcs = snlo_focus_func(problem(K,J));

        % funcs is a cell array containing function handles in the snlo_focus_func
        % function you can call from the commandline or this script. The first element of
        % funcs is the run button callback.

        focus_run_func = funcs{1};  % Choose the run button callback
        focus_run_func();           % Call the run button callback

        % We'll read the text written in the four output edit boxes in the figure. First,
        % get a handle to the figure using findobj looking for any items with a 'tag'
        % property equal to 'focus'
        focus_fig = findobj('tag','focus');

        % Then get the guidata stored in the focus figure
        D = guidata(focus_fig);
        
        % The guidata D contains 1 data structure with fieldname focus. That structure
        % contains fields corresponding to each uicomponent created by the Focus function.
        % Choose the four output edit boxes:
        rayleigh_range_box  = D.focus.edit_rayleigh;
        size_box            = D.focus.edit_size;
        curvature_mm_box    = D.focus.edit_rad_curv;
        farfield_angle_box  = D.focus.edit_farfield_ang;

        % Get the string contents of each of the four output edit boxes
        rayleigh_range_mm_string    = get(rayleigh_range_box,   'string');
        size_mm_string              = get(size_box,             'string');
        curvature_mm_string         = get(curvature_mm_box,     'string');
        farfield_angle_mrad_string  = get(farfield_angle_box,   'string');

        % Convert strings into numerical data, store in arrays
        farfield_angle_mrad(K,J)    = str2num(farfield_angle_mrad_string);
        curvature_mm(K,J)           = str2num(curvature_mm_string);
        size_mm(K,J)                = str2num(size_mm_string);
        rayleigh_range_mm(K,J)      = str2num(rayleigh_range_mm_string);
        fprintf(1,'Run %i of %i\n', (J-1)*length(face_to_focus_list)+K, length(face_to_focus_list)*length(waist_list));
    end
end

% Plot resulting data: beam width vs position, one curve per waist
fh = figure;
ax = axes('parent',fh);
ph = zeros(1, length(waist_list));
legendstrings = cell(length(waist_list),1);

for J = 1:length(waist_list)
    ph(J) = plot(ax,[-face_to_focus_list(end:-1:1),face_to_focus_list],...
        [size_mm(end:-1:1,J);size_mm(:,J)]);
    legendstrings{J} = sprintf('Waist = %.3G mm',waist_list(J));
    hold(ax,'on');
end
xlabel(ax,'Face to focus [mm]');
ylabel(ax,'Beam size [mm]');
legend(ax, legendstrings, 'location','SE');
