% This example of batch processing in the PW-OPO-SP function in SNLO is based on Exercise
% 6 of Chapter 10 of "Crystal Nonlinear Optics: with SNLO examples," by Arlee Smith. Here,
% we look at how changing cavity delay alters the output spectra.

% This example demonstrates how to set up the problem parameters for the plane-wave OPO model for short
% pulses (PW-OPO-SP), call the SNLO PW-OPO-SP function with the parameters, run the model, load and
% process the model outputs. We demonstrate several sets of similar input parameters, with a single input
% parameter varying from run to run. In this case, we'll vary the cavity delay and examine its  affect on
% the output spectrum.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    
    cavity_delay_list = linspace(0.067,0.060,5); % first make a vector of the values to use for cavity delay
    
    % the values used here are based on pw opo sp example 3 (#77 from the snlo examples button)
    
    % then define the full set of input parameters for PW-OPO-SP; these correspond to the input boxes
    % visible when you run PW-OPO-SP. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. calling the PW-OPO-SP function
    % with a set of input parameters requires specifying a single input argument to the function
    % snlo_pw_opo_sp_func, where the input parameter is a data structure with the fieldnames specified
    % below, and values in the form specified below.
    inputs.pw_opo_sp_wavelengths = [1200,2400,800]; % wavelengths (in nm)
    inputs.pw_opo_sp_ref_inds = [2.15,2.115,2.176]; % refractive indices
    inputs.pw_opo_sp_gvi = [2.199,2.183,2.265];     % group velocity index
    inputs.pw_opo_sp_gdd = [1.95e2,-1.7e2,3.71e2];  % group delay dispersion (in fs^2/mm)
    inputs.pw_opo_sp_crystal_loss = [0,0,0];        % linear absorption (in 1/mm)
    inputs.pw_opo_sp_n2_red1 = [0,0,0];     % Red1 n2 (in sq. cm/W)
    inputs.pw_opo_sp_n2_red2 = [0,0,0];     % Red2 n2 (in sq. cm/W)
    inputs.pw_opo_sp_n2_blue = [0,0,0];     % Blue n2 (in sq. cm/W)
    inputs.pw_opo_sp_beta_red1 =  [0,0,0];  % Red1 beta (in cm/W)
    inputs.pw_opo_sp_beta_red2 = [0,0,0];   % Red2 beta (in cm/W)
    inputs.pw_opo_sp_beta_blue = [0,0,0];   % Blue beta (in cm/W)
    inputs.pw_opo_sp_pulse_energy = 0.03;   % Blue pulse energy (in microjoules)
    inputs.pw_opo_sp_pulse_duration = 0.15; % Blue pulse duration (in FWHM picoseconds)
    inputs.pw_opo_sp_chirp = 0;             % Blue pulse chirp (in THz/ps)
    inputs.pw_opo_sp_beam_diameter = 0.1;   % Blue beam diameter (in FWHM mm)
    inputs.pw_opo_sp_output_refl = 0.9;     % Output coupler power reflectivity (for red1 wave only)
    inputs.pw_opo_sp_loss = 0.0;            % other cavity loss
    inputs.pw_opo_sp_cavity_delay = 0.067;  % Cavity delay (in ps)
    inputs.pw_opo_sp_cavity_gdd = 0;        % Cavity group delay dispersion (in fs^2)
    inputs.pw_opo_sp_lcrystal = 0.3;        % crystal length (in mm)
    inputs.pw_opo_sp_deff = 15;             % Nonlinear coefficient (in pm/V)
    inputs.pw_opo_sp_deltak = 0;            % Phase mismatch (in rad/mm)
    inputs.pw_opo_sp_nz = 30;               % number of z steps
    inputs.pw_opo_sp_nt = 1024;             % number of time steps
    inputs.display = 'Irrad';               % 'Irrad' radio button selected in the 'Display' panel
        
    inputs.max_num_passes = 500;            % Number of cavity passes to model (this doesn't correspond 
                                            % to an input box, and is used only for batching)

	% pre-allocate a cell array to store the contents of the output files in
    outputs = cell(size(cavity_delay_list));
    output_spectra = cell(size(cavity_delay_list));
    
    clear problem % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    %% loop over the cavity delay list:
    % call the function with a set of inputs, call the run-button callback, load the function's output, and close the gui figure
    for K = 1:length(cavity_delay_list)
        problem(K) = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).pw_opo_sp_cavity_delay = cavity_delay_list(K); % modify set the cavity delay to be an element of the vector of cavity delays generated earlier
        fcn_handles = snlo_pw_opo_sp_func(problem(K)); % call the SNLO PW-OPO-SP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles        
        close_handle = fcn_handles{end}; % the function handle of the 'Close figure' callback is the last element in the cell array of function handles
        run_handle(); % call the 'Run' button callback - equivalent to clicking the run button
        outputs{K} = importdata('PWOPO_SP.DAT'); % load the pw cav lp output file into memory, stick into the output cell array 'outputs'
        output_spectra{K} = importdata('PWOPOSPS.DAT');  % load the pw cav lp spectral output file into memory, stick into the output cell array 'output_spectra'
%         close_handle(); % call the 'Close figure' callback
    end
end

%% plot the output of the batch run
fh = figure; % make new figure
ax = axes('Parent',fh); % make new axes in the new figure
lstrings = cell(size(cavity_delay_list)); % pre-allocate a cell-array for sticking legend strings into
r1x1 = zeros(size(cavity_delay_list)); % pre-allocate a vector for lower x-axis limits for the plot
r1x2 = zeros(size(cavity_delay_list)); % pre-allocate a vector for upper x-axis limits for the plot
for K = 1:length(cavity_delay_list) % loop over the elements in the vector of cavity delays
    ph = plot(ax,output_spectra{K}(:,1),fftshift(output_spectra{K}(:,2)./max(output_spectra{K}(:,2)))); % add a plot curve to the axes we created earlier
    hold(ax,'on'); % make sure the next plot in the axes gets added to the axes instead of replacing curves already on that axes
    lstrings{K} = sprintf('%.3f fs',cavity_delay_list(K)*1e3); % generate a string to use in the legend for this newly plotted curve
    r1x1(K) = output_spectra{K}(find(fftshift(output_spectra{K}(:,2)./max(output_spectra{K}(:,2)))>=1E-4,1,'first'),1); % calculate the lowest frequency which has 0.01% of the peak power spectral density or larger
    r1x2(K) = output_spectra{K}(find(fftshift(output_spectra{K}(:,2)./max(output_spectra{K}(:,2)))>=1E-4,1,'last'),1);  % calculate the highest frequency which has 0.01% of the peak power spectral density or larger
end
xlabel(ax,'Frequency [THz]');                       % label the x-axis 
ylabel(ax,'Red1 power spectral density [arb]');     % label the y-axis
lh = legend(ax,lstrings,'Location','Best');         % make a legend to identify the plot curves in the axes, using the cell array of strings we generated earlier
xlim(ax,[min(r1x1),max(r1x2)]);                     % set the x-axis limits to be the lowest and highest frequencies that any curve has >= 0.01% of peak power spectral density at
title(ax,'Signal spectra vs. cavity delay','FontWeight','normal');