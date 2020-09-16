% This example of batch processing in the PW-CAV-LP function in SNLO is based on Exercise 1 of Chapter 11
% of "Crystal Nonlinear Optics: with SNLO examples," by Arlee Smith.

% This example demonstrates how to set up the problem parameters for the plane-wave model for long pulses
% in intracavity mixing pulses (PW-CAV-LP), call the SNLO PW-CAV-LP function with the parameters, run the
% model, load and process the model outputs. We demonstrate several sets of similar input parameters,
% with a single input parameter varying from run to run. In this case, we'll vary the blue wave input
% pulse energy and examine its affect on the output fluences.

calculate_curves = true;

if calculate_curves
    %% generate inputs
    blue_pulse_energy_list = linspace(1e-3,20e-3,10); % vary input pulse energy of blue wave: make equally spaced vector of energies from 1 mJ to 20 mJ with 10 elements
    % based on pw opo lp example 10 (#69 from snlo examples button)
    
    % define the full set of input parameters for PW-CAV-LP; these parameters correspond to the input
    % boxes visible when you run the PW-CAV-LP function. for most parameters, a row vector with three
    % elements is defined and corresponds to the values for the red1, red2, and blue waves. calling the
    % PW-CAV-LP function with a set of input parameters requires specifying a single input argument to
    % the function snlo_pw_cav_lp_func, where the input parameter is a data structure with the fieldnames
    % specified below, and values in the form specified below.
    
    inputs.pw_cav_lp_wavelengths = [800,1588.1,532];    % wavelengths (in nm)
    inputs.pw_cav_lp_ref_inds = [1.817,1.736,1.79];     % refractive inidices
    inputs.pw_cav_lp_left_refl_c = [0.013,0.013,0];     % left crystal face power reflectivities (0-1)
    inputs.pw_cav_lp_right_refl_c = [0.013,0.013,0];    % right crystal face power reflectivities (0-1)
    inputs.pw_cav_lp_crystal_loss = [0,0,0];            % linear crystal losses (in 1/mm)
    inputs.pw_cav_lp_left_energy_pwr = [1e-9,0,4.236e-3]; % left input pulse energies (in J) or powers (in W)
    inputs.pw_cav_lp_right_energy_pwr = [1e-9,0];       % right input pulse energies (in J) or powers (in W)
    inputs.pw_cav_lp_pulse_duration = [0,0,5];          % pulse durations (in FWHM ns)
    inputs.pw_cav_lp_beam_diameters = [1,1,1];          % beam diameters (in FWHM mm)
    inputs.pw_cav_lp_left_refl_m = [0.5,0.0,0.0];       % left mirror power reflectivities 
    inputs.pw_cav_lp_right_refl_m = [1,1,1];            % right mirror power reflectivities
    inputs.pw_cav_lp_phase_left = [0,0,0];              % left mirror to crystal phases (in radians)
    inputs.pw_cav_lp_phase_right = [0,0,0];             % right mirror to crystal phases (in radians)
    inputs.pw_cav_lp_third_leg_phase = [0,0,0];         % third leg phases (in radians) if ring cavity
    inputs.pw_cav_lp_cavity_length = 25.00;             % total length of cavity (in mm)
    inputs.pw_cav_lp_crystal_length = 5.00;             % length of nonlinear crystal (in mm)
    inputs.pw_cav_lp_deff = 3.21;                       % crystal nonlinear coefficient (in pm/V)
    inputs.pw_cav_lp_cavity_type = 1;                   % integer boolean for whether cavity is linear or ring type (1=lin, 0=ring)
    inputs.pw_cav_lp_deltak = 0;                        % phase mismatch (in radians/mm)
    inputs.pw_cav_lp_nz = 50;                           % number of z integration points through the nonlinear crystal
    inputs.pw_cav_lp_n_starts = 3;                      % time resolution variable (1 implies resolution 1 round trip time, 2 is half round trip time, etc)

    inputs.max_model_time = 50e-9;                      % stops model after this time, if all duration inputs are 0
    
    %% call 2d cav lp, run model, load output
    clear problem % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    outputs = cell(size(blue_pulse_energy_list)); % pre-allocate a cell array to store the contents of the output files in
    for K = 1:length(blue_pulse_energy_list) % loop over list of pulse energies for blue wave
        problem(K) = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).pw_cav_lp_left_energy_pwr(3) = blue_pulse_energy_list(K); % modify this problem's blue input pulse energy to be those defined in the vector of blue pulse energies generated earlier
        fcn_handles = snlo_pw_cav_lp_func(problem(K)); % call the SNLO pw-cav-lp function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned
        close_handle = fcn_handles{end}; % the function handle of the 'Close figure' callback is last returned
        run_handle(); % call the 'Run' button callback - equivalent to clicking the run button
        outputs{K} = load('pw_cav_lp.mat'); % load the pw cav lp output file into memory, stick into the output cell array
%         close_handle(); % call the 'Close figure' callback
    end
end

%% process output
left_red1_fluences = zeros(size(blue_pulse_energy_list));
for K = 1:length(blue_pulse_energy_list)
    % calculate the fluences for red1 wave by integrating in time the irradiance vs time (using
    % trapezoidal method). the first row of the left_output field of the outputs data structure is
    % time (in seconds), second row is red1 irradiance.
    left_red1_fluences(K) = trapz(outputs{K}.left_output(1,:),outputs{K}.left_output(2,:)); 
end

%% plot processed model outputs
fh = figure;            % make new figure
ax = axes('Parent',fh); % make new axes in that figure
ph = plot(ax,blue_pulse_energy_list,left_red1_fluences*1e-4,'ro-'); % plot left red1 fluences as function of input blue pulse energy & convert from J/m^2 to J/cm^2 with 1e-4 factor
xlabel(ax,'Input blue pulse energy [J]');   % label x axis
ylabel(ax,'Left red1 fluence [J/cm^2]');    % label y axis
lh = legend(ax,'Red1','Location','Best');   % use legend to label curve
