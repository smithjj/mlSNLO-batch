% This example demonstrates how to set up the problem parameters for the 2d model for
% optical parametric generation. The inputs are collected in a data structure which is fed
% to the SNLO OPG function. The script then run the model, loads and processes the model
% outputs. We demonstrate several sets of similar input parameters, with a single input
% parameter varying from run to run.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    input_pulse_energies = linspace(1e-6,60e-6,15); % make vector of equally spaced pulse energies between 100 nanojoules and 5 microjoules
    
    % then define the full set of input parameters for opg; these correspond to the input boxes
    % visible when you run opg. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. Calling the opg function
    % with a set of input parameters requires selecting a single set of inputs and sending it to
    % the snlo_opg_func as an intput argument. 
    
    inputs.opg_2d_wavelengths = [1550,1550,775];     % wavelengths (in nm)
    inputs.opg_2d_ref_inds = [2.142,2.142,2.184];    % refractive indices
    inputs.opg_2d_gvi = [2.1875,2.1875,2.2811];      % group velocity indices
    inputs.opg_2d_gdd = [1.04e2,1.04e2,3.99e2];      % group delay dispersions (in fs^2/mm)
    inputs.opg_2d_phase = [0,0,0];                   % input phases (in radians)
    inputs.opg_2d_input_refl = [0,0,0];              % crystal input face power reflectivity coefficients (0-1)
    inputs.opg_2d_output_refl = [0,0,0];             % crystal output face power reflectivity coefficients (0-1)
    inputs.opg_2d_crystal_losses = [0,0,0];          % linear crysta losses (in 1/mm)
    inputs.opg_2d_n2_red1 = [0,0,0];                 % red1 nonlinear refractive indices (in sq cm/W)
    inputs.opg_2d_n2_red2 = [0,0,0];                 % red2 nonlinear refractive indices (in sq cm/W)
    inputs.opg_2d_n2_blue = [0,0,0];                 % blue nonlinear refractive indices (in sq cm/W)
    inputs.opg_2d_beta_red1 = [0,0,0];               % red1 two photon absorption coefficients (in cm/W)
    inputs.opg_2d_beta_red2 = [0,0,0];               % red2 two photon absorption coefficients (in cm/W)
    inputs.opg_2d_beta_blue = [0,0,0];               % blue two photon absorption coefficients (in cm/W)
    inputs.opg_2d_pulseenergy = [0,0,1e-5];          % input pulse energies (in J)
    inputs.opg_2d_pulse_durations = [0.1,0.1,0.1;1,1,1].'; % input pulse durations (in fwhm ps); with temporal supergaussian coefficients (positive integers)
    inputs.opg_2d_pulse_delays = [0,0];              % temporal delays relative to blue pulse (in ps)
    inputs.opg_2d_pulse_chirps = [0,0,0];            % linear chrips (in THz/ps)
    inputs.opg_2d_beam_diameters = [0.5,0.5,0.5;0.5,0.5,0.5].'; % input beam diameters (in fwhm mm), in walkoff direction; and perpendicular to walkoff direction
    inputs.opg_2d_supergaussian_coeff = [1,1,1];     % spatial supergaussian coefficients (integers 1-10)
    inputs.opg_2d_wo_angles = [0,0,0];               % spatial walkoff angles (in milliradians)
    inputs.opg_2d_offset_wodir = [0,0];              % spatial offsets in walkoff directions (in mm)
    inputs.opg_2d_rad_curv = [1e7,1e7,1e7;1e7,1e7,1e7].';  % input beam radii of curvature (in mm in air) for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction]
    inputs.opg_2d_nt = 512;                           % number of time points to consider
    inputs.opg_2d_nxny = [32,32];                    % number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
    inputs.opg_2d_grid_duration = 0.3;               % total time to model (in ps)
    inputs.opg_2d_crystal_length = 0.5;              % length of nonlinear crystal (in mm)
    inputs.opg_2d_lx_ly = [2,2];                     % size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
    inputs.opg_2d_deff = 30;                         % crystal nonlinear coefficient (in pm/V)
    inputs.opg_2d_deltak = 0;                        % phase mismatch (in rad/mm)
    inputs.opg_2d_nz = 30;                           % number of (z) integration steps through the crystal
    inputs.opg_2d_dist_to_image = 0;                 % Distance from the crystal exit face to image plane (in mm) of the output beam profile displayed during the run and by the post-run fluence and movie buttons.
    
    %% call 2d mix sp, run model, load output
    clear problem

    powers = cell(size(input_pulse_energies));
    output_energies = zeros(length(input_pulse_energies),3);

    for K = 1:length(input_pulse_energies)
        problem(K) = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).opg_2d_pulseenergy(3) = input_pulse_energies(K); % modify this problem's phase mismatch to be an element of the vector of phase mismatches generated earlier
        fcn_handles = snlo_opg_func(problem(K)); % call the SNLO 2D-MIX-SP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1};    % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        accept_handle = fcn_handles{2}; % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        close_handle = fcn_handles{end};% the function handle of the 'Close figure' callback is last returned
        accept_handle();    % call the accept button callback - equivalent to clicking the accept button
        run_handle();       % call the 'Run' button callback - equivalent to clicking the run button
        data = load('opg_2d_output.mat');
        powers{K} = data.power;
        output_energies(K,1) = trapz(powers{K}(:,1),powers{K}(:,2));
        output_energies(K,2) = trapz(powers{K}(:,1),powers{K}(:,3));
        output_energies(K,3) = trapz(powers{K}(:,1),powers{K}(:,4));
%         close_handle(); % call the 'Close figure' callback. it's not necessary to include this but might help limit the clutter of figures we open
    end
end

%% plot output energies as a function of input energies
fh = figure; % make new output figure
ax = axes('Parent',fh); % make new axes in the new figure
ph = plot(ax, input_pulse_energies, output_energies(:,1), 'r-', ...
    input_pulse_energies, output_energies(:,2), 'g-', ...
    input_pulse_energies, output_energies(:,3), 'b-');
xlabel(ax, 'Input pulse energies [J]');
ylabel(ax, 'Output pulse energies [J]');
legend(ax, 'Red1', 'Red2', 'Blue');
