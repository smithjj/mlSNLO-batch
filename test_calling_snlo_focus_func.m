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

face_to_focus_list = linspace(0.05,1,20); % in mm
waist_list = 0.05:-0.01:0.01;
% dist_to_focus_list = dist_to_focus_list(2:end);

farfield_angle_mrad = zeros(length(face_to_focus_list),length(waist_list));
curvature_mm = zeros(length(face_to_focus_list),length(waist_list));
size_mm = zeros(length(face_to_focus_list),length(waist_list));
rayleigh_range_mm = zeros(length(face_to_focus_list),length(waist_list));

clear problem
for J = 1:length(waist_list)
    for K = 1:length(face_to_focus_list)
        problem(K,J) = inputs;
        problem(K,J).focus_dist = face_to_focus_list(K);
        problem(K,J).focus_waist = waist_list(J);
        funcs = snlo_focus_func(problem(K,J));
        focus_run_func = funcs{1};
        focus_run_func();

        % The edit boxes in the Focus function aren't normal MATLAB
        % uicomponents; instead they are Java Swing components, which I chose
        % to use because doing so lets you define callbacks to fire on
        % keypresses and also when the cursor selects an edit box. We use this
        % to draw the diagrams you see when selcting 'Waist,' 'Face to focus,'
        % or 'Dist. to focus'.
        focus_figure = findobj('tag','focus');
        children = get(focus_figure,'children');

        farfield_angle_box = get(children(1),'JavaPeer');
        curvature_mm_box = get(children(2),'JavaPeer');
        size_box = get(children(3),'JavaPeer');
        rayleigh_range_box = get(children(4),'JavaPeer');

        curvature_mm_string = curvature_mm_box.getText;
        size_mm_string = size_box.getText;
        farfield_angle_mrad_string = farfield_angle_box.getText;
        rayleigh_range_mm_string = rayleigh_range_box.getText;

        farfield_angle_mrad(K,J) = str2num(farfield_angle_mrad_string);
        curvature_mm(K,J) = str2num(curvature_mm_string);
        size_mm(K,J) = str2num(size_mm_string);
        rayleigh_range_mm(K,J) = str2num(rayleigh_range_mm_string);

    end
end

% keyboard;
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
