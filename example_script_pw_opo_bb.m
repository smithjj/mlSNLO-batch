% This example of batch processing in the PW-OPO-BB function in SNLO is based on Exercise
% 12 of Chapter 11 of "Crystal Nonlinear Optics: with SNLO examples," by Arlee Smith.
% Here, we examine how the output spectra change as Delta k changes.

% This example demonstrates how to set up the problem parameters for the plane-wave model for long pulses
% in intracavity mixing pulses (PW-OPO-BB), call the SNLO PW-OPO-BB function with the parameters, run the
% model, load and process the model outputs. We demonstrate several sets of similar input parameters,
% with a single input parameter varying from run to run. In this case, we'll vary the phase mismatch and
% examine its affect on the output spectra. This is supposed to show that seeding fails if the phase
% mismatch gets too large.

calculate_curves = true; % set to true to run the model several times and process and plot the output, to false if you just want to tweak the processing and plotting without rerunning the model
if calculate_curves
    %% generate inputs
    delta_k_list = linspace(0.25,0.30,6); % make equally spaced vector of phase mismatches from 0.25 to 0.30 with 6 elements
    % the values used here are based on pw opo bb example 1 (#48 from snlo examples button)
    
    % then define the full set of input parameters for PW-OPO-BB; these correspond to the input boxes
    % visible when you run PW-OPO-BB. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. calling the PW-OPO-BB function
    % with a set of input parameters requires specifying a single input argument to the function
    % snlo_pw_opo_bb_func, where the input parameter is a data structure with the fieldnames specified
    % below, and values in the form specified below.
    
    inputs.pw_opo_bb_wavelengths = [800,1588.1,532];    % wavelengths (in nm)
    inputs.pw_opo_bb_ref_inds = [1.817,1.736,1.79];     % refractive indices
    inputs.pw_opo_bb_gvi = [1.95,2.05,2.0];             % group velocity indices
    inputs.pw_opo_bb_left_refl_c = [0,0,0];             % left crystal face power reflectivity coefficient (0-1)
    inputs.pw_opo_bb_right_refl_c = [0,0,0];            % right crystal face power reflectivity coefficient (0-1)
    inputs.pw_opo_bb_crystal_loss = [0,0,0];            % linear crystal loss coefficient (in 1/mm)
    inputs.pw_opo_bb_left_energy_pwr = [0,0,0.005];     % left input pulse energy (in J) for pulses or input power (in W) for cw
    inputs.pw_opo_bb_right_energy_pwr = [1e-4,0];       % right input pulse energy (in J) for pulses or input power (in W) for cw; only red1 & red2
    inputs.pw_opo_bb_pulse_duration = [0,0,3];          % pulse durations (in FWHM ps) (0 for cw)
    inputs.pw_opo_bb_beam_diameters= [1,1,1];           % beam diameters (in FWHM mm)
    inputs.pw_opo_bb_left_refl_m = [0.99,0,0];          % left mirror power reflectivity coefficient (0-1)
    inputs.pw_opo_bb_right_refl_m = [0.7,0,0];          % right mirror power reflectivity coefficient (0-1)
    inputs.pw_opo_bb_phase_left = [0,0,0];              % left input phases (in radians)
    inputs.pw_opo_bb_phase_right = [0,0,0];             % right input phases (in radians)
    inputs.pw_opo_bb_cavity_length  = 30.00;            % total length of cavity (in mm)
    inputs.pw_opo_bb_crystal_length = 10.00;            % length of nonlinear crystal (in mm)
    inputs.pw_opo_bb_pump_bandwidth = 0.00;             % Bandwidth (in fwhm MHZ).
    inputs.pw_opo_bb_pump_mode_spacing = 250;           % Frequency spacing of longitudinal modes (in MHz)
    inputs.pw_opo_bb_pump_fm = 0;                       % integer boolean for consider frequency modulated light (1) or chaotic light (0)
    inputs.pw_opo_bb_deff = 3.2;                        % crystal nonlinear coefficient (in pm/V)
    inputs.pw_opo_bb_cavity_type  = 0;                  % integer boolean for whether cavity is linear or ring type (1=lin, 0=ring)
    inputs.pw_opo_bb_deltak = 0;                        % phase mismatch (in radians/mm)
    inputs.pw_opo_bb_nz = 30;                           % number of z integration points through the nonlinear crystal
    inputs.pw_opo_bb_new_noise = 1;                     % integer boolean for calculating new noise (1) or re-using noise calculated in previous run (0)
    inputs.r1 = 1; % show red1 curves in plots (1=true or 0=false)
    inputs.r2 = 1; % show red2 curves in plots (1=true or 0=false)
    inputs.bl = 1; % show blue curves in plots (1=true or 0=false)
    inputs.leftright = 'Right';
    
    %% call pw opo bb, run model, load output
    clear problem % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    outputs = cell(size(delta_k_list));         % pre-allocate a cell array to store the contents of the output files in
    output_spectra = cell(size(delta_k_list));  % pre-allocate a cell array to store the contents of the output spectra files in
    for K = 1:length(delta_k_list)              % loop over each element in the vector of phase mismatches
        problem(K) = inputs;                    % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).pw_opo_bb_deltak = delta_k_list(K); % modify this problem's phase mismatch to be those defined in the vector of phase mismatches generated earlier
        fcn_handles = snlo_pw_opo_bb_func(problem(K)); % call the SNLO pw-opo-bb function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1};            % the function handle of the 'Run' button callback is the first returned
        close_handle = fcn_handles{end};        % the function handle of the 'Close figure' callback is last returned
        run_handle();                           % call the 'Run' button callback - equivalent to clicking the run button
        outputs{K} = load('pw_opo_bb_output.mat'); % load the pw opo bb output file into memory, stick into the output cell array
        output_spectra{K} = importdata('BBSPECR.DAT'); % load the pw-opo-bb right spectra output file
%         close_handle();                         % call the 'Close figure' callback
    end
end

%% plot the model outputs
fh = figure;                % make new figure
ax = axes('Parent',fh);     % make new axis in that figure
lstrings = cell(size(delta_k_list)); % pre-allocate cell array to stick strings for use in figure legend's labels
r1x1 = zeros(size(delta_k_list));    % pre-allocate vector to stick lowest frequencies with significant spectral power 
r1x2 = zeros(size(delta_k_list));    % pre-allocate vector to stick highest frequencies with significant spectral power
for K = 1:length(delta_k_list) 
    % make new plot curve in axis for each phase mismatch
    ph = plot(ax,output_spectra{K}(1,:)*1e-9,(K-1)+output_spectra{K}(2,:)./max(output_spectra{K}(2,:))); % plot normalized output power spectral densities as a function of frequency (plot freq in GHz, and we have it in Hz, so factor of 1E-9)
    hold(ax,'on'); % make next plot in loop add another curve to the axis instead of replace curve we just plotted
    lstrings{K} = ['{\Delta}k=',sprintf('%.3f rad/mm',delta_k_list(K))]; % generate string for this phase mismatch to use in figure legend
    r1x1(K) = output_spectra{K}(1,find(output_spectra{K}(2,:)./max(output_spectra{K}(2,:))>=1E-2,1,'first')); % find lowest frequency for this phase mismatch which has >=1% peak power spectral density
    r1x2(K) = output_spectra{K}(1,find(output_spectra{K}(2,:)./max(output_spectra{K}(2,:))>=1E-2,1,'last'));  % find highest frequency for this phase mismatch which has >=1% peak power spectral density
end
xlabel(ax,'Frequency [GHz]'); % label x axis
ylabel(ax,'Red1 power spectral density [arb]'); % label y axis
lh = legend(ax,lstrings,'Location','Best'); % add legend to the axis using the cell array of strings we used earlier
xlim(ax,[min(r1x1),max(r1x2)]*1e-9); % set x axis limits to the lowest and highest frequencies which has >=1% peak power spectral density