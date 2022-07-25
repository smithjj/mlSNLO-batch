% This example of batch processing in the PW-OPO-BB function in SNLO is based on Exercise 2 of Chapter 9
% of "Crystal Nonlinear Optics: with SNLO examples," by Arlee Smith.

% This example demonstrates how to set up the problem parameters for the 2d model for long pulses in
% intracavity mixing pulses (2d-cav-lp), call the SNLO 2d-cav-lp function with the parameters, run the
% model, load and process the model outputs. We demonstrate several sets of similar input parameters,
% with a single input parameter varying from run to run. In this case, we'll vary the phase between the
% crystal and right mirror, and examine its affect on the output powers when the crystal nonlinearity is
% set to zero. This is supposed to show that resonating light in the cavity requires the correct phase to
% counteract the Gouy phase. It also shows that higher reflectivity mirrors means narrower resonance or
% more sensitivity to the phase - compare the sensitivity of the transmitted red1 light to that of the
% red2.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    phase_cr_list = linspace(1.546,1.746,21); % make vector of equally spaced phases between crystal and right mirror, between 1.546 and 1.746 radians with 21 elements
    
    % then define the full set of input parameters for 2d-cav-lp; these correspond to the input boxes
    % visible when you run 2d-cav-lp. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. calling the 2d-cav-lp function
    % with a set of input parameters requires specifying a single input argument to the function
    % snlo_2d_cav_lp_func, where the input parameter is a data structure with the fieldnames specified
    % below, and values in the form specified below.
    
    % based on 2d opo lp example 2 (#34 from snlo examples button)
    inputs.cav_2d_lp_wavelengths = [800,800,400];   % wavelengths (in nm)
    inputs.cav_2d_lp_ref_inds = [1.6,1.6,1.625];    % refractive indices
    inputs.cav_2d_lp_left_refl_c = [0,0,0];         % left crystal power reflectivity coefficient (0-1)
    inputs.cav_2d_lp_right_refl_c = [0,0,0];        % right crystal power reflectivity coefficient (0-1)
    inputs.cav_2d_lp_crystal_loss = [0,0,0];        % linear crystal losses (in 1/mm)
    inputs.cav_2d_lp_left_energy_pwr = [0.001,0.001,0.001]; % left input pulse energy (in J) for pulses or power (in W) for cw
    inputs.cav_2d_lp_right_energy_pwr = [0,0]; % red1 and red2 right input pulse energy (in J) for pulses or power (in W) for cw
    inputs.cav_2d_lp_pulse_duration = [0,0,0;1,1,1].'; % input pulse durations (in fwhm ns), with temporal super Gaussian coeffcient (integer 1-10)
    inputs.cav_2d_lp_pulse_delay = [0,0]; % temporal delays relative to blue input pulse for red1,red2 waves (in ns)
    inputs.cav_2d_lp_beam_diameters  = [0.13787,0.13787,0.15];  % input beam diameters (in FWHM mm)
    inputs.cav_2d_lp_spatial_supergaussian = [1,1,1];           % spatial input beam super gaussian coefficient (integer 1-10)
    inputs.cav_2d_lp_walkoff_angles = [0,0,0];                  % spatial walkoff angles (in milliradians)
    inputs.cav_2d_lp_beam_offset = [0,0];                       % spatial input beam offset in the walkoff direction relative to blue beam center(in mm)
    inputs.cav_2d_lp_beam_roc = [50,50,50];                     % input beam radii of curvature (in mm in air) for [red1,red2,blue]
    inputs.cav_2d_lp_left_refl_m = [0.95,0.9,0.1];              % left mirror power reflectivity coefficient (0-1)
    inputs.cav_2d_lp_right_refl_m = [0.95,0.9,0.9];             % right mirror power reflectivity coefficient (0-1)
    inputs.cav_2d_lp_phase_lc = [0,0,0];                        % single pass phase for left mirror to left crystal face (in radians)
    inputs.cav_2d_lp_phase_cr = [1.646,1.646,0];                % single pass phase for right crystal face to right mirror (in radians)
    inputs.cav_2d_lp_phase_rl = [0,0,0];                        % single pass phase for right mirror to left mirror if ring cavity (in radians)
    inputs.cav_2d_lp_mirror_rocs = [50,50];                     % left and right mirrors radii of curvature (in mm)
    inputs.cav_2d_lp_leg_lengths = [17.5,17.5,30];              % cavity leg lengths for left mirror to left crystal face, right crystal face to right mirror, and right mirror to left mirror for ring cavity (in mm)
    inputs.cav_2d_lp_nz = 30;                                   % number of (z) integration steps through nonlinear crystal
    inputs.cav_2d_lp_nxny = [32,32];                            % number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
    inputs.cav_2d_lp_crystal_length = 30;                       % length of nonlinear crystal (in mm)
    inputs.cav_2d_lp_lxly = [0.75,0.75];                        % size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
    inputs.cav_2d_lp_cavity_type = 1;                           % 0 for ring cavity;1 for linear cavity.
    inputs.cav_2d_lp_cavity_inversion = 0;                      % 1 if ring inverts beam in walk off direction on each cavity pass. 0 otherwise.
    inputs.cav_2d_lp_deff = 0;                                  % crystal nonlinear coefficient (in pm/V)
    inputs.cav_2d_lp_deltak = 0;                                % phase mismatch (in radians/mm)
    inputs.cav_2d_lp_nstarts = 2;                               % number of starts-inverse of time resolution; number of time points per cavity round trip time; increase to improve time resolution
    
    inputs.max_model_time = 50e-9;                              % stop model execution after modeling up to this point in time
    
    %% call 2d cav lp, run model, load output
    clear problem % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    output = cell(size(phase_cr_list)); % pre-allocate a cell array to store the contents of the output files in
    for K = 1:length(phase_cr_list)     % for each element in the vector of phases between right crystal face and right mirror
        problem(K) = inputs;            % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).cav_2d_lp_phase_cr(1:2) = phase_cr_list(K); % modify this problem's C-R phase to be an element of the vector of C-R phases generated earlier
        fcn_handles = snlo_2d_cav_lp_func(problem(K)); % call the SNLO 2d-cav-lp function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1};    % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        accept_handle = fcn_handles{2}; % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        close_handle = fcn_handles{end};% the function handle of the 'Close figure' callback is last returned
        accept_handle();                % call the 'accept' button's callback
        run_handle();                   % call the 'Run' button callback - equivalent to clicking the run button
        output{K} = importdata('PWR_R.DAT'); % load the 2d-cav-lp right output powers from file, stick them in an element in the output cell array
%         close_handle();                 % call the 'Close figure' callback
    end
end

%% process output
red1_output_power = zeros(size(phase_cr_list)); % pre-allocate vector of red1 output powers 
red2_output_power = zeros(size(phase_cr_list)); % pre-allocate vector of red1 output powers 
for K = 1:length(phase_cr_list)
    red1_output_power(K) = output{K}(end,2); % select the last element in the time dimension ('end') for the red1 wave (2nd column)
    red2_output_power(K) = output{K}(end,3); % select the last element in the time dimension ('end') for the red2 wave (3rd column)
end

%% plot processed output
fh = figure; % make new figure
ax = axes('Parent',fh); % make new axes in that figure
ph = plot(ax,phase_cr_list,red1_output_power,'r-',...
    phase_cr_list,red2_output_power,'b-');  % plot red1 and red2 output powers as a function of right crystal to right mirror phases
xlabel(ax,'Phase (C-R) [rad]'); % label x axis
ylabel(ax,'Transmitted power [W]'); % label y axis
xlim(ax,[min(phase_cr_list),max(phase_cr_list)]); % set x axis limits to be lowest and highest value in vector of C-R phases
legend(ax,'Red1','Red2','Location','Best'); % make legend in the axis with appropriate strings