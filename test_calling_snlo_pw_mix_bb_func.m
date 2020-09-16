% This example demonstrates how to set up the problem parameters for the plane-wave model a broad
% bandwidth mixer including group velocities but no spatial walkoff or diffraction, PW-MIX-BB, call the
% SNLO PW-MIX-BB function with the parameters, run the model, load the model outputs, process the model
% outputs and plot them. We demonstrate several sets of similar input parameters, with the input pulse
% energy parameter varying from run to run. The output energies of the red1, red2, and blue waves are
% then plotted as a function of input pulse energy. In this example, since the broad bandwidth model
% includes a noise term, each set of parameters is run several times and we plot the average output
% on-axis fluence with error bars set by the standard deviation of the on-axis fluence.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    red1red2_energy_list = linspace(1e-10,100E-8,11); % make an evenly spaced vector of red1 and red2 input pulse energies, from 0.1 nanojoule to 1 microjoule with 11 elements
    N_run_per_energy = 15; % run each set of input parameters this many times and use the results to calculate mean output fluences and their spread, which we'll plot
    
    % based on pw mix bb example 1 (#12 from snlo examples button)
    
    % define the full set of input parameters for PW-MIX-BB; these parameters correspond to the input
    % boxes visible when you run the PW-MIX-BB function. for most parameters, a row vector with three
    % elements is defined and corresponds to the values for the red1, red2, and blue waves. calling the
    % PW-MIX-BB function with a set of input parameters requires specifying a single input argument to
    % the function snlo_pw_mix_bb_func, where the input parameter is a data structure with the fieldnames
    % specified below, and values in the form specified below.
    
    inputs.pw_mix_bb_wavelengths = [800,800,400];       % wavelengths (in nm) [of red1,red2,blue]
    inputs.pw_mix_bb_ref_inds = [1.8,1.8,1.8];          % refractive indices
    inputs.pw_mix_bb_gvi = [2,2,2];                     % group velocity indices
    inputs.pw_mix_bb_gdd  = [0,0,0];                    % group delay dispersion (in fs^2/mm)
    inputs.pw_mix_bb_input_refl = [0,0,0];              % input face power reflectivies (0-1)
    inputs.pw_mix_bb_output_refl = [0,0,0];             % output face power reflectivies (0-1)
    inputs.pw_mix_bb_crystal_loss = [0,0,0];            % linear crystal losses (in 1/mm)
    inputs.pw_mix_bb_pulse_energy_power = [5e-9,5e-9,0];% input pulse energies (in J)
    inputs.pw_mix_bb_pulse_duration = [3,3,3];          % pulse durations (in fwhm ns)
    inputs.pw_mix_bb_beam_diameter = [1,1,1];           % input beam diameters (in fwhm mm)
    inputs.pw_mix_bb_bandwidth = [100,100,0];           % Bandwidth for multi-longitudinal mode pulses (in fwhm GHZ). Set to 0 for single longitudinal mode pulses.
    inputs.pw_mix_bb_mode_spacing = [0.5,0.5,0.25];     % Frequency spacing of longitudinal modes (in GHz)
    inputs.pw_mix_bb_freq_mod  = [0,0,0];               % integer boolean for consider frequency modulated light (1) or chaotic light (0)
    inputs.pw_mix_bb_crystal_length = 20;               % length of crystal (in mm)
    inputs.pw_mix_bb_deff = 10;                         % crystal nonlinear coefficient (in pm/V)
    inputs.pw_mix_bb_deltak = 0;                        % phase mismatch (in rad/mm)
    inputs.pw_mix_bb_nz = 100;                          % Number of integration steps through the crystal (in z direction)

    
    %% call pw mix bb, run model, load output
    clear problem % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    output = cell(length(red1red2_energy_list),N_run_per_energy);
    for K = 1:length(red1red2_energy_list)  % loop over each element in the vector of input pulse energies
        for J = 1:N_run_per_energy          % run each input pulse energy N_run_per_energy times times 
            problem(K,J) = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
            problem(K,J).pw_mix_bb_pulse_energy_power = [red1red2_energy_list(K),red1red2_energy_list(K),0]; % modify this problem's input pulse energies to be an element of the vector of input pulse energies generated earlier
            fcn_handles = snlo_pw_mix_bb_func(problem(K));  % call the SNLO 2D-MIX-SP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
            run_handle = fcn_handles{1};                    % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
            close_handle = fcn_handles{end};                % the function handle of the 'Close figure' callback is last returned
            run_handle();   % call the 'Run' button callback - equivalent to clicking the run button
            output{K,J} = importdata('pw_mix_bb.dat'); % load the pw_mix_bb output file into memory, stick into the output cell array
%             close_handle(); % call the 'Close figure' callback
        end
    end
    
end

%% process output files
red1_fluences = zeros(length(red1red2_energy_list),N_run_per_energy); % pre-allocate a array of red1 output energies as fcn of input pulse energies, and with several runs per input energy (N_run_per_energy)
red2_fluences = zeros(length(red1red2_energy_list),N_run_per_energy); % pre-allocate a array of red2 output energies as fcn of input pulse energies, and with several runs per input energy (N_run_per_energy)
blue_fluences = zeros(length(red1red2_energy_list),N_run_per_energy); % pre-allocate a array of blue output energies as fcn of input pulse energies, and with several runs per input energy (N_run_per_energy)
red1_fluence_mean = zeros(length(red1red2_energy_list),1); % pre-allocate a vector of mean output pulse energies for red1 wave
red2_fluence_mean = zeros(length(red1red2_energy_list),1); % pre-allocate a vector of mean output pulse energies for red2 wave
blue_fluence_mean = zeros(length(red1red2_energy_list),1); % pre-allocate a vector of mean output pulse energies for blue wave
red1_fluence_std = zeros(length(red1red2_energy_list),1);  % pre-allocate a vector of standard deviation of output pulse energies for red1 wave
red2_fluence_std = zeros(length(red1red2_energy_list),1);  % pre-allocate a vector of standard deviation of output pulse energies for red2 wave
blue_fluence_std = zeros(length(red1red2_energy_list),1);  % pre-allocate a vector of standard deviation of output pulse energies for blue wave
for K = 1:length(red1red2_energy_list)  % loop over the input pulse energies vector
    for J = 1:N_run_per_energy          % run each input pulse energy N_run_per_energy times times 
        t = output{K,J}(:,1)*1e-9;      % time vector is saved as the first column in the output file, specified in nanoseconds, so convert to seconds
        I_red1 = output{K,J}(:,2);      % red1 irradiance vs time is saved as the second column in the output file (in W/m^2)
        I_red2 = output{K,J}(:,3);      % red2 irradiance vs time is saved as the second column in the output file (in W/m^2) 
        I_blue = output{K,J}(:,4);      % blue irradiance vs time is saved as the second column in the output file (in W/m^2)
        red1_fluences(K,J) = trapz(t,I_red1)*1e-4; % find red1 fluence by integrating in time the irradiance using trapezoidal method, convert to J/cm^2 
        red2_fluences(K,J) = trapz(t,I_red2)*1e-4; % find red2 fluence by integrating in time the irradiance using trapezoidal method, convert to J/cm^2 
        blue_fluences(K,J) = trapz(t,I_blue)*1e-4; % find blue fluence by integrating in time the irradiance using trapezoidal method, convert to J/cm^2 
    end
    red1_fluence_mean(K) = mean(red1_fluences(K,:)); % find the mean output fluence for red1 at each input pulse energy
    red2_fluence_mean(K) = mean(red2_fluences(K,:)); % find the mean output fluence for red2 at each input pulse energy
    blue_fluence_mean(K) = mean(blue_fluences(K,:)); % find the mean output fluence for blue at each input pulse energy
    red1_fluence_std(K) = std(red1_fluences(K,:));   % find the standard deviation of the output fluence for red1 at each input pulse energy
    red2_fluence_std(K) = std(red2_fluences(K,:));   % find the standard deviation of the output fluence for red2 at each input pulse energy
    blue_fluence_std(K) = std(blue_fluences(K,:));   % find the standard deviation of the output fluence for blue at each input pulse energy
end

%% plot processed output
fh = figure; % make new output figure
ax = axes('Parent',fh); % make new axes in that output figure
phs = zeros(3,1); % pre-allocate a vector with 3 elements for the 3 plot curves we'll stick in the new axes (red1, red2, and blue)
phs(1) = errorbar(ax,2*red1red2_energy_list,red1_fluence_mean,red1_fluence_std,'r-'); % in the new axes, plot red1 mean fluence as function of total input pulse energy
hold(ax,'on'); % tell the axes the next plot will be added to the existing plot, instead of replacing it
phs(2) = errorbar(ax,2*red1red2_energy_list,red2_fluence_mean,red2_fluence_std,'b-'); % in the same axes, plot red2 mean fluence as function of total input pulse energy
phs(3) = errorbar(ax,2*red1red2_energy_list,blue_fluence_mean,blue_fluence_std,'g-'); % in the same axes, plot blue mean fluence as function of total input pulse energy
xlabel(ax,'Input red1+red2 energy [J]'); % label the x-axis
ylabel(ax,'Output fluence [J/cm^2]');    % label the y-axis
legend(ax,'Red1 fluence','Red2 fluence','Blue fluence','Location','Best'); % make legend to describe each curve on the plot, locate the legend in the 'best' place