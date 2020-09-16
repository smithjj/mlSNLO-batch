% This example demonstrates how to set up the problem parameters for the 2d model for long pulses in
% intracavity mixing pulses (pw-mix-sp), call the SNLO pw-mix-sp function with the parameters, run the
% model, load and process the model outputs. We demonstrate several sets of similar input parameters,
% with a single input parameter varying from run to run. In this case, we'll model a short-pulse
% frequency-doubling device and determine the output pulse energy with varying input pulse energies.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    % make a vector of equally spaced input pulse energies between 0 and 250 nanojoules with 15 elements.
    % each element in this vector will be used in an iteration of the model as the value for each of the
    % red input pulse energies.
    red1red2_input_energy_list = linspace(0,2.5e-7,15); 
    
    % then define the full set of input parameters for pw-mix-sp; these correspond to the input boxes
    % visible when you run pw-mix-sp. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. calling the pw-mix-sp function
    % with a set of input parameters requires specifying a single input argument to the function
    % snlo_pw_mix_sp_func, where the input parameter is a data structure with the fieldnames specified
    % below, and values in the form specified below.
    
    % based on pw mix sp example 1 (#1 from snlo examples button)
    input.pw_mix_sp_wavelengths = [1064,1064,532];  % wavelengths (in nm)
    input.pw_mix_sp_ref_inds = [1.654,1.654,1.654]; % refractive indices
    input.pw_mix_sp_gvi = [1.60,1.60,1.61];         % group velocity indices
    input.pw_mix_sp_gdd = [0,0,0];                  % group delay dispersion (in fs^2/mm)
    input.pw_mix_sp_phases = [0,0,0];               % input phases (in radians)
    input.pw_mix_sp_n2_red1 = [0,0,0];              % nonlinear refractive indices from red1 wave (in cm^2/W)
    input.pw_mix_sp_n2_red2 = [0,0,0];              % nonlinear refractive indices from red2 wave (in cm^2/W)
    input.pw_mix_sp_n2_blue = [0,0,0];              % nonlinear refractive indices from blue wave (in cm^2/W)
    input.pw_mix_sp_beta_red1 = [0,0,0];            % two-photon absorptions from red1 wave (in cm/W)
    input.pw_mix_sp_beta_red2 = [0,0,0];            % two-photon absorptions from red2 wave (in cm/W)
    input.pw_mix_sp_beta_blue = [0,0,0];            % two-photon absorptions from blue wave (in cm/W)
    input.pw_mix_sp_input_refl = [0,0,0];           % input face power reflectivity coefficient (0-1)
    input.pw_mix_sp_output_refl = [0,0,0];          % output face power reflectivity coefficient (0-1)
    input.pw_mix_sp_crystal_losses = [0,0,0];       % linear crystal absorption (in 1/mm)
    input.pw_mix_sp_pulseenergy = [1e-8,1e-8,0];    % input pulse energies (in J)
    input.pw_mix_sp_beam_diameters = [4,4,4];       % beam diameters (in fwhm mm)
    input.pw_mix_sp_pulse_durations = [0.05,0.05,0.05;1,1,1].'; % input pulse durations (in fwhm ps); and temporal supergaussian coefficient (1-10) or 0 for hyperbolic secant shape
    input.pw_mix_sp_pulse_chirp = [0,0,0];          % linear chirps of the input pulses (in THz/ps)
    input.pw_mix_sp_pulse_delay = [0,0];            % pulse delays relative to pump (in ps)
    input.pw_mix_sp_crystal_length = 10;            % nonlinear crystal length (in mm)
    input.pw_mix_sp_deff = 5;                       % crystal nonlinear coefficient (in pm/V)
    input.pw_mix_sp_deltak = 0;                     % phase mismatch (in rad/mm)
    input.pw_mix_sp_nz = 100;                       % number of z points
    input.pw_mix_sp_nt = 512;                       % number of time points
    
    %% call pw mix sp, run model, load output
    
    % loop over blue input energies: launch pw_mix_lp, run the model, load the output into memory, and close
    % the figure
    clear problem; % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    output = cell(size(red1red2_input_energy_list)); % pre-allocate a cell array to store the contents of the output files in
    for K = 1:length(red1red2_input_energy_list)
        problem(K) = input;             % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).pw_mix_sp_pulseenergy(1:2) = red1red2_input_energy_list(K); % modify this problem's red1 and red2 input pulse energies to be an element of the vector of input pulse energies generated earlier
        fcn_handles = snlo_pw_mix_sp_func(problem(K)); % call the SNLO pw-mix-sp function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1};    % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        close_handle = fcn_handles{end};% the function handle of the 'Close figure' callback is last returned
        run_handle(); % call the 'Run' button callback - equivalent to clicking the run button
        output{K} = importdata('pw_mix_sp.dat'); % load the pw_mix_sp output file into memory, stick into the output cell array
%         close_handle(); % call the 'Close figure' callback
    end
end
%% process output
output_fluences = zeros(3,length(red1red2_input_energy_list)); % pre-allocate vector of red1 output powers 
for K = 1:length(red1red2_input_energy_list)
    % calculate the on-axis output fluences for the red1, red2, and blue waves
    times = output{K}(:,1); % the first column of the output file is the time in seconds
    output_red1_irrads = output{K}(:,2); % the second column is the red1 on-axis irradiance at that time
    output_red2_irrads = output{K}(:,3); % the third column is the red2 on-axis irradiance at that time
    output_blue_irrads = output{K}(:,4); % the fourth column is the blue on-axis irradiance at that time
    output_fluences(1,K) = trapz(times,output_red1_irrads)*1e-4; % use the trapezoidal method to integrate the red1 irradiance in time, multiply by 10^-4 to get fluence in J/cm^2
    output_fluences(2,K) = trapz(times,output_red2_irrads)*1e-4; % use the trapezoidal method to integrate the red2 irradiance in time, multiply by 10^-4 to get fluence in J/cm^2
    output_fluences(3,K) = trapz(times,output_blue_irrads)*1e-4; % use the trapezoidal method to integrate the blue irradiance in time, multiply by 10^-4 to get fluence in J/cm^2
end


%% plot the processed model outputs
fh = figure;            % make new figure
ax = axes('Parent',fh); % make new axis in that figure
% plot the output on-axis fluence as a function of the total input pulse energy. since half the input
% pulse energy goes in red1 and half in red2, multiply red1red2_input_energy_list by 2
ph = plot(ax,2*red1red2_input_energy_list,output_fluences(1,:),'ro-',...     % plot the red1 fluence as fcn of blue input energy
    2*red1red2_input_energy_list,output_fluences(2,:),'bo-',...              % plot the red2 fluence as fcn of blue input energy
    2*red1red2_input_energy_list,output_fluences(3,:),'go-');                % plot the blue fluence as fcn of blue input energy
xlh = xlabel(ax,'Input red1+red2 energy [J]');              % label the x-axis
ylh = ylabel(ax,'Output on-axis fluences [J/cm^2]');        % label the y-axis
lh  = legend(ax,'Red1','Red2','Blue','Location','Best');    % make legend for this axes
