% This example demonstrates how to set up the problem parameters for the 2d short-pulse mixing function
% 2D-MIX-SP, call the SNLO 2D-MIX-SP function with the parameters, run the model, load the model outputs,
% process the model outputs and plot them. We demonstrate several sets of similar input parameters, with
% the phase mismatch parameter varying from run to run. The output energies of the red1, red2, and blue
% waves are then plotted as a function of phase mismatch.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    delta_k_list = linspace(-15,15,31); % make a vector of phase mismatches, spanning -15 to 15 radians with 31 elements.
    % based on 2d mix sp example 2 (#66 from snlo examples button)
    
    % define the full set of input parameters for 2D-MIX-SP; these parameters correspond to the input
    % boxes visible when you run the 2D-MIX-SP function. for most parameters, a row vector with three
    % elements is defined and corresponds to the values for the red1, red2, and blue waves. calling the
    % 2D-MIX-SP function with a set of input parameters requires specifying a single input argument to
    % the function snlo_2d_mix_sp_func, where the input parameter is a data structure with the fieldnames
    % specified below, and values in the form specified below.
    
    inputs.mix_2d_sp_wavelengths = [1550,1550,775];     % wavelengths (in nm)
    inputs.mix_2d_sp_ref_inds = [2.142,2.142,2.184];    % refractive indices
    inputs.mix_2d_sp_gvi = [2.1875,2.1875,2.2811];      % group velocity indices
    inputs.mix_2d_sp_gdd = [1.04e2,1.04e2,3.99e2];      % group delay dispersion (in fs^2/mm)
    inputs.mix_2d_sp_phase = [0,0,0];                   % phases at input (in rad)
    inputs.mix_2d_sp_input_refl = [0,0,0];              % input face power reflectivity (0-1)
    inputs.mix_2d_sp_output_refl = [0,0,0];             % output face power reflectivity (0-1)
    inputs.mix_2d_sp_crystal_losses = [0,0,0];          % linear power absorption coefficient (in 1/mm)
    inputs.mix_2d_sp_n2_red1 = [0,0,0];                 % red1 nonlinear refractive index (in sq cm/W)
    inputs.mix_2d_sp_n2_red2 = [0,0,0];                 % red2 nonlinear refractive index (in sq cm/W)
    inputs.mix_2d_sp_n2_blue = [0,0,0];                 % blue nonlinear refractive index (in sq cm/W)
    inputs.mix_2d_sp_beta_red1 = [0,0,0];               % red1 two-photon absorption coefficient (in cm/W)
    inputs.mix_2d_sp_beta_red2 = [0,0,0];               % red2 two-photon absorption coefficient (in cm/W)
    inputs.mix_2d_sp_beta_blue = [0,0,0];               % blue two-photon absorption coefficient (in cm/W)
    inputs.mix_2d_sp_pulseenergy = [5e-8,5e-8,0];       % input pulse energies (in J)
    inputs.mix_2d_sp_pulse_durations = [0.1,0.1,0.1;1,1,1].'; % input pulse durations  (in fwhm ps), with optional (super) gaussian coefficient (
    inputs.mix_2d_sp_pulse_delays = [0,0];              % pulse delay relative to the blue pulse (in ps)
    inputs.mix_2d_sp_pulse_chirps = [0,0,0];            % Linear chirp of input pulse (in THz/ps)
    inputs.mix_2d_sp_beam_diameters = [0.5,0.5,0.5;0.5,0.5,0.5].'; % input beam diameter (in fwhm mm), specified in walkoff direction, and perpendicular to walkoff direction
    inputs.mix_2d_sp_supergaussian_coeff = [1,1,1];     % spatial super gaussian coefficients
    inputs.mix_2d_sp_wo_angles = [0,0,0];               % spatial walkoff angles (in mrad)
    inputs.mix_2d_sp_offset_wodir = [0,0];              % spatial offset relative to blue beam center in walk off direction (in mm)
    inputs.mix_2d_sp_rad_curv = [1e7,1e7,1e7;1e7,1e7,1e7].';  % input beam radii of curvature (in mm in air) for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction]
    inputs.mix_2d_sp_nt= 64;                            % number of time points to consider
    inputs.mix_2d_sp_nxny = [64,64];                    % number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
    inputs.mix_2d_sp_crystal_length = 0.5;              % length of nonlinear crystal (in mm)
    inputs.mix_2d_sp_lx_ly = [2,2];                     % size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
    inputs.mix_2d_sp_deff = 30;                         % crystal nonlinear coefficient (in pm/V)
    inputs.mix_2d_sp_deltak = 0;                        % phase mismatch (in rad/mm)
    inputs.mix_2d_sp_nz = 30;                           % number of (z) integration steps through the crystal
    inputs.mix_2d_sp_dist_to_image = 0;                 % Distance from the crystal exit face to image plane (in mm) of the output beam profile displayed during the run and by the post-run fluence and movie buttons.
    
    %% call 2d mix sp, run model, load output
    clear problem   % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    output = cell(size(delta_k_list)); % pre-allocate a cell array to store the contents of the output files in
    for K = 1:length(delta_k_list) % loop over each element in the vector of phase mismatches
        problem(K) = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).mix_2d_sp_deltak = delta_k_list(K);  % modify this problem's phase mismatch to be an element of the vector of phase mismatches generated earlier
        fcn_handles = snlo_2d_mix_sp_func(problem(K));  % call the SNLO 2D-MIX-SP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1};        % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        accept_handle = fcn_handles{2};     % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        close_handle = fcn_handles{end};    % the function handle of the 'Close figure' callback is last returned
        accept_handle();                    % call the 'accept' button's callback - equivalent to clicking the accept button
        run_handle();                       % call the 'Run' button callback - equivalent to clicking the run button
        output{K} = load('mix_2d_sp_output.mat'); % load the 2d mix sp output file into memory, stick into the output cell array
%         close_handle(); % call the 'Close figure' callback
    end
end

%% process output
red1_energies = zeros(size(delta_k_list)); % pre-allocate a vector for red1 output energies as fcn of phase mismatch
red2_energies = zeros(size(delta_k_list)); % pre-allocate a vector for red2 output energies as fcn of phase mismatch
blue_energies = zeros(size(delta_k_list)); % pre-allocate a vector for blue output energies as fcn of phase mismatch
for K = 1:length(delta_k_list) % loop over the phase mismatch vector
    red1_energies(K) = trapz(output{K}.power(:,1),output{K}.power(:,2))*1e3; % calculate red1 energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
    red2_energies(K) = trapz(output{K}.power(:,1),output{K}.power(:,3))*1e3; % calculate red2 energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
    blue_energies(K) = trapz(output{K}.power(:,1),output{K}.power(:,4))*1e3; % calculate blue energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
end

%% plot processed output
fh = figure; % make new outut figure
ax = axes('Parent',fh); % make new axes in that output figure
% plot the red1, red2, and blue energies as a function of phase mismatch
ph = plot(ax,delta_k_list,red1_energies,'ro-',...
    delta_k_list,red2_energies,'bo-',...
    delta_k_list,blue_energies,'go-');
xlabel(ax,'\Delta k');                                  % label x-axis 
ylabel(ax,'Output energy [mJ]');                        % label y-axis 
lh = legend(ax,'Red1','Red2','Blue','Location','Best'); % make legend to describe each curve on the plot