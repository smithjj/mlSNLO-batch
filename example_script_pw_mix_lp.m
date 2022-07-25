% This example of batch processing in the PW-MIX-LP function in SNLO is based on Exercise
% 5 of Chapter 3 of "Crystal Nonlinear Optics: with SNLO examples," by Arlee Smith. This
% example demonstrates an example where a small absorption on the red 1
% wave increases the output energy of the red 2 wave because it acts to limit back conversion.

% We show here how to set up the problem parameters for the plane-wave model for short
% pulses (PW-MIX-LP), call the SNLO PW-MIX-LP function with the parameters, run the model, load and
% process the model outputs. We demonstrate several sets of similar input parameters, with a single input
% parameter varying from run to run. In this case, we'll vary the red1 crystal loss and examine its
% affect on the output fluences.

calculate_curves = true;
if calculate_curves
    %% generate inputs
    red1_crystal_loss_list = linspace(0,1.5,10); % vary red1 crystal losses input energies
    % based on pw mix lp example 5 (#21 from snlo examples button)
     
    % define the full set of input parameters for PW-MIX-LP; these parameters correspond to the input
    % boxes visible when you run the PW-MIX-LP function. for most parameters, a row vector with three
    % elements is defined and corresponds to the values for the red1, red2, and blue waves. calling the
    % PW-MIX-LP function with a set of input parameters requires specifying a single input argument to
    % the function snlo_pw_mix_lp_func, where the input parameter is a data structure with the fieldnames
    % specified below, and values in the form specified below.
    
    input.pw_mix_lp_wavelengths = [1500,750,500];   % wavelengths (in nm)
    input.pw_mix_lp_ref_inds = [1.5,1.5,1.5];       % refractive indicies
    input.pw_mix_lp_phases = [0,0,0];               % phases at input (in radians)
    input.pw_mix_lp_input_reflectivities = [0,0,0]; % power reflectivities at crystal input face
    input.pw_mix_lp_output_reflectivities = [0,0,0]; % power reflectivities at crystal output face
    input.pw_mix_lp_crystal_losses = [0,0,0];       % linear crystal loss (in 1/mm)
    input.pw_mix_lp_energy_power = [0,1e-7,0.1];    % input energy (in J) if duration is nonzero or power (in W) if duration is zero
    input.pw_mix_lp_pulse_durations = [1,1,1];      % pulse durations (in ns fwhm)
    input.pw_mix_lp_beam_diameters = [1,1,1];       % beam diameters (in fwhm mm)
    input.pw_mix_lp_crystal_length = 10;            % length of crystal (in mm)
    input.pw_mix_lp_deff = 1.5;                     % crystal nonlinear coefficient (in pm/V)
    input.pw_mix_lp_deltak = 0;                     % phase mismatch (in rad/mm)
    input.pw_mix_lp_nz = 50;                        % number of z integration steps (through the crystal)
    input.pw_mix_lp_nt = 200;                       % number of time points
    
    
    %% call pw mix lp, run model, load output
    output = cell(size(red1_crystal_loss_list)); % pre-allocate a cell array to store the contents of the output files in
    clear problem; % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    for K = 1:length(red1_crystal_loss_list)% loop over each element in the vector of crystal losses for red1 wave
        problem(K) = input;                 % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).pw_mix_lp_crystal_losses(1) = red1_crystal_loss_list(K); % modify this problem's red1 crystal loss to be those defined in the vector of crystal losses generated earlier
        fcn_handles = snlo_pw_mix_lp_func(problem(K)); % call the SNLO pw-mix-lp function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1};        % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        close_handle = fcn_handles{end};    % the function handle of the 'Close figure' callback is last returned
        run_handle();                       % call the 'Run' button callback - equivalent to clicking the run button
        output{K} = importdata('pw_mix_lp.dat'); % load the pw_mix_lp output file into memory, stick into the output cell array
%         close_handle();                     % call the 'Close figure' callback
    end
end

% process the model outputs
output_fluences = zeros(3,length(red1_crystal_loss_list));
for K = 1:length(red1_crystal_loss_list)
    times = output{K}(:,1); % the first column of the output file is the time in seconds
    % calculate the on-axis output fluences for the red1, red2, and blue waves
    output_red1_irrads = output{K}(:,2); % the second column is the red1 on-axis irradiance at that time
    output_red2_irrads = output{K}(:,3); % the third column is the red2 on-axis irradiance at that time
    output_blue_irrads = output{K}(:,4); % the fourth column is the blue on-axis irradiance at that time
    output_fluences(1,K) = trapz(times,output_red1_irrads); % use the trapezoidal method to integrate the red1 irradiance in time
    output_fluences(2,K) = trapz(times,output_red2_irrads); % use the trapezoidal method to integrate the red2 irradiance in time
    output_fluences(3,K) = trapz(times,output_blue_irrads); % use the trapezoidal method to integrate the blue irradiance in time
end

%% plot the processed model outputs
fh = figure;                    % make new figure
ax = axes('Parent',fh);         % make new axis in that figure
ph = plot(ax,red1_crystal_loss_list,output_fluences(1,:),'ro-',...     % plot the red1 fluence as fcn of blue input energy
    red1_crystal_loss_list,output_fluences(2,:),'bo-',...              % plot the red2 fluence as fcn of blue input energy
    red1_crystal_loss_list,output_fluences(3,:),'go-');                % plot the blue fluence as fcn of blue input energy
% label the plot appropriately
set(ax,'FontName','Times New Roman');
% xlh = xlabel(ax,'Input blue energy [J]');
xlh = xlabel(ax,'Red1 crystal loss [1/mm]');
ylh = ylabel(ax,'Output on-axis fluences [J/m^2]');
lh  = legend(ax,'Red1','Red2','Blue','Location','Best');