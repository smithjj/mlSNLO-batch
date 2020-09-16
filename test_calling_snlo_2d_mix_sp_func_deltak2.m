% Phase mismatch tuning of the second harmonic light in short pulse frequency doubling

% This example demonstrates how to set up the problem parameters for the 2d model for short pulses mixing
% (2d-mix-sp), call the SNLO 2d-mix-sp function with the parameters, run the model, load and process the
% model outputs. We demonstrate several sets of similar input parameters, with a single input parameter
% varying from run to run. In this case, we'll vary the phase mismatch, and examine its effect on the
% output spectra for the blue wave.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    delta_k_list = linspace(-15,15,5); % make vector of equally spaced phase mismatches between -15 radians and 15 radians with 5 elements
    
    % then define the full set of input parameters for 2d-mix-sp; these correspond to the input boxes
    % visible when you run 2d-mix-sp. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. calling the 2d-mix-sp function
    % with a set of input parameters requires specifying a single input argument to the function
    % snlo_2d_mix_sp_func, where the input parameter is a data structure with the fieldnames specified
    % below, and values in the form specified below.
    
    % based on 2d mix sp example 2 (#66 from snlo examples button)
    inputs.mix_2d_sp_wavelengths = [1550,1550,775];     % wavelengths (in nm)
    inputs.mix_2d_sp_ref_inds = [2.142,2.142,2.184];    % refractive indices
    inputs.mix_2d_sp_gvi = [2.1875,2.1875,2.2811];      % group velocity indices
    inputs.mix_2d_sp_gdd = [1.04e2,1.04e2,3.99e2];      % group delay dispersions (in fs^2/mm)
    inputs.mix_2d_sp_phase = [0,0,0];                   % input phases (in radians)
    inputs.mix_2d_sp_input_refl = [0,0,0];              % crystal input face power reflectivity coefficients (0-1)
    inputs.mix_2d_sp_output_refl = [0,0,0];             % crystal output face power reflectivity coefficients (0-1)
    inputs.mix_2d_sp_crystal_losses = [0,0,0];          % linear crysta losses (in 1/mm)
    inputs.mix_2d_sp_n2_red1 = [0,0,0];                 % red1 nonlinear refractive indices (in sq cm/W)
    inputs.mix_2d_sp_n2_red2 = [0,0,0];                 % red2 nonlinear refractive indices (in sq cm/W)
    inputs.mix_2d_sp_n2_blue = [0,0,0];                 % blue nonlinear refractive indices (in sq cm/W)
    inputs.mix_2d_sp_beta_red1 = [0,0,0];               % red1 two photon absorption coefficients (in cm/W)
    inputs.mix_2d_sp_beta_red2 = [0,0,0];               % red2 two photon absorption coefficients (in cm/W)
    inputs.mix_2d_sp_beta_blue = [0,0,0];               % blue two photon absorption coefficients (in cm/W)
    inputs.mix_2d_sp_pulseenergy = [5e-8,5e-8,0];       % input pulse energies (in J)
    inputs.mix_2d_sp_pulse_durations = [0.1,0.1,0.1;1,1,1].'; % input pulse durations (in fwhm ps); with temporal supergaussian coefficients (positive integers)
    inputs.mix_2d_sp_pulse_delays = [0,0];              % temporal delays relative to blue pulse (in ps)
    inputs.mix_2d_sp_pulse_chirps = [0,0,0];            % linear chrips (in THz/ps)
    inputs.mix_2d_sp_beam_diameters = [0.5,0.5,0.5;0.5,0.5,0.5].'; % input beam diameters (in fwhm mm), in walkoff direction; and perpendicular to walkoff direction
    inputs.mix_2d_sp_supergaussian_coeff = [1,1,1];     % spatial supergaussian coefficients (integers 1-10)
    inputs.mix_2d_sp_wo_angles = [0,0,0];               % spatial walkoff angles (in milliradians)
    inputs.mix_2d_sp_offset_wodir = [0,0];              % spatial offsets in walkoff directions (in mm)
    inputs.mix_2d_sp_rad_curv = [1e7,1e7,1e7;1e7,1e7,1e7].';  % input beam radii of curvature (in mm in air) for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction]
    inputs.mix_2d_sp_nt= 64;                            % number of time points to consider
    inputs.mix_2d_sp_nxny = [32,32];                    % number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
    inputs.mix_2d_sp_crystal_length = 0.5;              % length of nonlinear crystal (in mm)
    inputs.mix_2d_sp_lx_ly = [2,2];                     % size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
    inputs.mix_2d_sp_deff = 30;                         % crystal nonlinear coefficient (in pm/V)
    inputs.mix_2d_sp_deltak = 0;                        % phase mismatch (in rad/mm)
    inputs.mix_2d_sp_nz = 30;                           % number of (z) integration steps through the crystal
    inputs.mix_2d_sp_dist_to_image = 0;                 % Distance from the crystal exit face to image plane (in mm) of the output beam profile displayed during the run and by the post-run fluence and movie buttons.
    
    %% call 2d mix sp, run model, load output
    clear problem
    spectra = cell(size(delta_k_list)); % pre-allocate cell array to load spectra from file into
    for K = 1:length(delta_k_list)
        problem(K) = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).mix_2d_sp_deltak = delta_k_list(K); % modify this problem's phase mismatch to be an element of the vector of phase mismatches generated earlier
        fcn_handles = snlo_2d_mix_sp_func(problem(K)); % call the SNLO 2D-MIX-SP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1};    % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        accept_handle = fcn_handles{2}; % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        close_handle = fcn_handles{end};% the function handle of the 'Close figure' callback is last returned
        accept_handle();    % call the accept button callback - equivalent to clicking the accept button
        run_handle();       % call the 'Run' button callback - equivalent to clicking the run button
        spectra{K} = load('BEAM_3WP.DAT'); % load the blue spectra file
        spectra{K}(:,1) = spectra{K}(:,1)*1E-12; % convert frequencies from Hz to THz
        spectra{K}(:,2) = spectra{K}(:,2)./max(spectra{K}(:,2)); % normalize power spectral densities to be maximum of 1
%         close_handle(); % call the 'Close figure' callback
    end
end

%% plot output spectra
fh = figure; % make new output figure
ax = axes('Parent',fh); % make new axes in the new figure
ph = zeros(size(delta_k_list)); % pre-allocate vector for plot handles
lstrings = cell(size(delta_k_list)); % pre-allocate cell vector to stick strings for legend in
bx1 = zeros(size(delta_k_list)); % pre-allocate vector to store lowest frequencies substantial power in
bx2 = zeros(size(delta_k_list)); % pre-allocate vector to store highest frequencies substantial power in
for K = 1:length(delta_k_list) % loop over vector of phase mismatches 
    ph(K) = plot(ax,spectra{K}(:,1),spectra{K}(:,2)); % plot output power spectral density as function of frequency in THz
    hold(ax,'on'); % tell axis the next plots will add additional curves instead of replacing them
    lstrings{K} = ['{\Delta}k = ',sprintf('%.2f',delta_k_list(K))]; % generate string to describe this iteration's phase mismatch value, stick in cell vector element
    bx1(K) = spectra{K}(find(spectra{K}(:,2)>1E-4,1,'first'),1); % calculate lowest frequency with substantial power (>0.01% maximum power)
    bx2(K) = spectra{K}(find(spectra{K}(:,2)>1E-4,1,'last'),1);  % calculate highest frequency with substantial power (>0.01% maximum power)
end
lh = legend(ax,lstrings,'Location','Best'); % make legend in this axis using the strings we generated earlier
xlabel(ax,'Frequency [THz]'); % label x axis
ylabel(ax,'Blue wave spectral power density [arb]'); % label y axis
xlim(ax,[min(bx1),max(bx2)]); % set x axis limits to be lowest and highest frequencies with substantial powers