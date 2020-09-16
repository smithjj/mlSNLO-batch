% This example demonstrates how to set up the problem parameters for the 2d long-pulse mixing function
% 2D-MIX-LP, call the SNLO 2D-MIX-LP function with the parameters, run the model, load the model outputs,
% process the model outputs and plot them. We demonstrate several sets of similar input parameters, with
% the input pulse energies for the red1 and red2 waves varying from run to run. The output energies of
% the red1, red2, and blue waves are then plotted as a function of the sum of red1 and red2 input pulse
% energy.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    red1red2_pulse_energy_list = linspace(1e-3,10e-2,10); % make a vector of input pulse energies between 0.001 J and 0.1 J with 10 elements
    % based on 2d mix lp example 4 (#26 from snlo examples button)
    
    % define the full set of input parameters for 2D-MIX-LP; these parameters correspond to the input
    % boxes visible when you run the 2D-MIX-LP function. for most parameters, a row vector with three
    % elements is defined and corresponds to the values for the red1, red2, and blue waves. calling the
    % 2D-MIX-LP function with a set of input parameters requires specifying a single input argument to
    % the function snlo_2d_mix_lp_func, where the input parameter is a data structure with the fieldnames
    % specified below, and values in the form specified below.
    inputs.mix_2d_lp_wavelengths = [1000,1000,500];         % wavelengths (in nm)
    inputs.mix_2d_lp_ref_inds = [1.655,1.655,1.655];        % refractive indices
    inputs.mix_2d_lp_phase = [0,0,0];                       % input phases (in radians)
    inputs.mix_2d_lp_input_refl = [0,0,0];                  % power reflectivity at crystal input face
    inputs.mix_2d_lp_output_refl = [0,0,0];                 % power reflectivity at crystal output face
    inputs.mix_2d_lp_crystal_losses = [0,0,0];              % linear absorption (in 1/mm)
    inputs.mix_2d_lp_pulseenergy = [8e-2,8e-2,0];           % input pulse energies (in joules)
    inputs.mix_2d_lp_pulse_durations = [1,1,1];             % input pulse durations (in fwhm ns)
    inputs.mix_2d_lp_beam_diameters = [10,10,10;10,10,10].';% input beam diameters (in fwhm mm) [red1 in walkoff dir, red2 in walkoff dir, blue in walkoff dir;red1 perp to walkoff, red2 perp to walkoff, blue perp to walkoff]
    inputs.mix_2d_lp_supergaussian_coeff = [10,10,10];      % supergaussian coefficient for spatial input beam profile
    inputs.mix_2d_lp_n2_red1 = [0,0,0];                     % red1 nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_n2_red2 = [0,0,0];                     % red2 nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_n2_blue = [0,0,0];                     % blue nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_beta_red1 = [0,0,0];                   % red1 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_beta_red2 = [0,0,0];                   % red2 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_beta_blue = [0,0,0];                   % blue 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_wo_angles = [0,0,0];                   % spatial walkoff angles (in milliradians)
    inputs.mix_2d_lp_offset_wodir = [0,0,0];                % spatial offset in walkoff direction (in mm)
    inputs.mix_2d_lp_rad_curv = [1e9,1e9,1e9;1e9,1e9,1e9].';% input beam radii of curvature for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction] (in mm in air)
    inputs.mix_2d_lp_nz = 100;                              % number of z integration steps
    inputs.mix_2d_lp_nxny = [64,64];                        % number of grid points in the transverse (x,y) directions
    inputs.mix_2d_lp_crystal_length = [10,1];               % crystal length (in mm); optionally, the second parameter is number of walkoff compensating sections
    inputs.mix_2d_lp_lx_ly = [30,30];                       % size of grid in walkoff direction, perpendicular to walkoff direction (in mm)
    inputs.mix_2d_lp_deff = 2.01;                           % crystal nonlinear coefficient (in pm/V)
    inputs.mix_2d_lp_deltak = 0;                            % phase mismatch (in rad/mm)
    inputs.mix_2d_lp_dist_to_image = 0;                     % Distance from the crystal exit face to image plane  (in mm) of the output beam profile displayed during the run and by the post-run fluence and movie buttons.
    inputs.mix_2d_lp_nt = 30;                               % number of time steps to include (Recommend starting with 32 and varying to check convergence.)
    
    %% call 2d mix lp, run model, load output
    output_red1 = cell(size(red1red2_pulse_energy_list)); % pre-allocate a cell array to store the contents of the red1 output files in
    output_red2 = cell(size(red1red2_pulse_energy_list)); % pre-allocate a cell array to store the contents of the red2 output files in
    output_blue = cell(size(red1red2_pulse_energy_list)); % pre-allocate a cell array to store the contents of the blue output files in
    clear problem % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    for K = 1:length(red1red2_pulse_energy_list) % loop over each element in the vector of input pulse energies
        problem(K) = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).mix_2d_lp_pulseenergy = [red1red2_pulse_energy_list(K),red1red2_pulse_energy_list(K),0]; % modify this problem's red1 and red2 input pulse energies to be those defined in the vector of pulse energies generated earlier
        fcn_handles = snlo_2d_mix_lp_func(problem(K)); % call the SNLO 2D-MIX-LP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        accept_handle = fcn_handles{2}; % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        close_handle = fcn_handles{end}; % the function handle of the 'Close figure' callback is last returned
        accept_handle();    % call the 'accept' button's callback - equivalent to clicking the accept button
        run_handle();       % call the 'Run' button callback - equivalent to clicking the run button
        output_red1{K} = importdata('mix_2d_lp_red1_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        output_red2{K} = importdata('mix_2d_lp_red2_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        output_blue{K} = importdata('mix_2d_lp_blue_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
%         close_handle(); % call the 'Close figure' callback
    end
end
%% process output
red1_energies = zeros(size(red1red2_pulse_energy_list)); % pre-allocate a vector for red1 output energies as fcn of input pulse energy
red2_energies = zeros(size(red1red2_pulse_energy_list)); % pre-allocate a vector for red2 output energies as fcn of input pulse energy
blue_energies = zeros(size(red1red2_pulse_energy_list)); % pre-allocate a vector for blue output energies as fcn of input pulse energy
for K = 1:length(red1red2_pulse_energy_list) % loop over the input pulse energy vector
    red1_energies(K) = trapz(output_red1{K}.data(:,1),output_red1{K}.data(:,2))*1e3; % calculate red1 energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
    red2_energies(K) = trapz(output_red2{K}.data(:,1),output_red2{K}.data(:,2))*1e3; % calculate red2 energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
    blue_energies(K) = trapz(output_blue{K}.data(:,1),output_blue{K}.data(:,2))*1e3; % calculate blue energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
end

%% plot processed output
fh = figure;  % make new outut figure
ax = axes('Parent',fh); % make new axes in that output figure

% plot the red1, red2, and blue output energies as a function of total input pulse energy (red1+red2)
ph = plot(ax,2*red1red2_pulse_energy_list,red1_energies,'ro-',...
    2*red1red2_pulse_energy_list,red2_energies,'bo-',...
    2*red1red2_pulse_energy_list,blue_energies,'go-');
xlabel(ax,'Input red1+red2 pulse energy [J]');  % label x-axis 
ylabel(ax,'Output energy [mJ]');                % label y-axis 
lh = legend(ax,'Red1','Red2','Blue','Location','Best'); % make legend to describe each curve on the plot
