% Modified 2D-mix-LP example (with user-specified input fields)
%
% This is an example of two-stage third-harmonic generation (using two
% crystals in series) without propagation between them. The first crystal
% is Type I second-harmonic generation, and the second third-harmonic generation.
% This example is a simple batch of runs where the SHG red input energy is
% varied.
%
% Note: This example script doesn't generate any fields to use as input
% fields to the 2D-mix-LP model directly, but rather uses the outputs of
% the first 2D-mix-LP stage as inputs for the second. SHG runs, then the
% script loads the output fields from the SHG model run from disk and uses
% them as the input fields for the THG (after some manipulation.)
%
% If you wish to generate your own electric field arrays to use as inputs
% for the modified 2D-mix-LP function, the text file with the name
% 'notes_on_snlo_2d_mix_lp_with_electric_field_inputs.txt' contains further
% information.
%
% SHG (type I): 1064 nm -> 532 nm
%   * In SNLO we do type I SHG by splitting the red input into both
%   red waves with equal powers (with identical values for other properties).
%   * Blue wave (532 nm) input energy is 0.
% THG: 1064 nm + 532 nm -> 354.7 nm
%   * We will collect the residual of both red output waves of the SHG and 
%   combine them in one of the red THG waves, and use the frequency doubled 
%   output of the SHG stage as the other red wave.
%   * Blue wave (354.7 nm) input energy is 0.


c = 3e8;
epsilon_0 = 8.85e-12;

calculate_curves = true;
if calculate_curves
    % Construct the vector of red energies we'll consider (this energy goes
    % into each of the red waves, so the total input energy will be twice
    % this value)
    red1_red2_energy_vec = linspace(0.05,1.25,13)*1e-3; % in J

    % Construct the data structure of inputs for the modified 2D-mix-LP
    % function for SHG modeling
    % ------------------------------------------------------------
    clear shg_inputs
    shg_inputs.mix_2d_lp_wavelengths    = [1064,1064,532];              % wavelengths (in nm)
    shg_inputs.mix_2d_lp_ref_inds       = [1.655,1.655,1.655];          % refractive indices
    shg_inputs.mix_2d_lp_phase          = [0,0,0];                      % input phases (in radians)
    shg_inputs.mix_2d_lp_input_refl     = [0,0,0];                      % crystal input reflectivity
    shg_inputs.mix_2d_lp_output_refl    = [0,0,0];                      % crystal input reflectivity
    shg_inputs.mix_2d_lp_crystal_losses = [0,0,0];                      % linear absorption (in 1/mm)
    shg_inputs.mix_2d_lp_pulseenergy    = [10e-2,10e-2,0];              % input pulse energies (in joules)
    shg_inputs.mix_2d_lp_pulse_durations = [1,1;1,1;1,1];               % input pulse durations (in fwhm ns)
    shg_inputs.mix_2d_lp_beam_diameters = [1,1;1,1;1,1];                % input beam diameters (in fwhm mm) [red1 in walkoff dir, red2 in walkoff dir, blue in walkoff dir;red1 perp to walkoff, red2 perp to walkoff, blue perp to walkoff]
    shg_inputs.mix_2d_lp_supergaussian_coeff = [1,1,1];                 % supergaussian coefficient for spatial input beam profile
    shg_inputs.mix_2d_lp_pulse_delays   = [0,0,0];                      % Pulse delay (ns)
    shg_inputs.mix_2d_lp_n2_red1        = [0,0,0];                      % red1 nonlinear refractive index (in sq. cm/W)
    shg_inputs.mix_2d_lp_n2_red2        = [0,0,0];                      % red2 nonlinear refractive index (in sq. cm/W)
    shg_inputs.mix_2d_lp_n2_blue        = [0,0,0];                      % blue nonlinear refractive index (in sq. cm/W)
    shg_inputs.mix_2d_lp_beta_red1      = [0,0,0];                      % red1 2-photon absorption cofficient (in cm/W)
    shg_inputs.mix_2d_lp_beta_red2      = [0,0,0];                      % red2 2-photon absorption cofficient (in cm/W)
    shg_inputs.mix_2d_lp_beta_blue      = [0,0,0];                      % blue 2-photon absorption cofficient (in cm/W)
    shg_inputs.mix_2d_lp_wo_angles      = [0,0,0];                      % spatial walkoff angles (in milliradians)
    shg_inputs.mix_2d_lp_offset_wodir   = [0,0,0];                      % spatial offset in walkoff direction (in mm)
    shg_inputs.mix_2d_lp_rad_curv       = [1e9,1e9,1e9;1e9,1e9,1e9].';  % input beam radii of curvature for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction] (in mm in air)
    shg_inputs.mix_2d_lp_nz             = 64;                           % number of z integration steps
    shg_inputs.mix_2d_lp_nxny           = [64,64];                      % number of grid points in the transverse (x,y) directions
    shg_inputs.mix_2d_lp_crystal_length = [10,1];                       % crystal length (in mm); optionally, the second parameter is number of walkoff compensating sections
    shg_inputs.mix_2d_lp_lx_ly          = [3,3];                        % size of grid in walkoff direction, perpendicular to walkoff direction (in mm)
    shg_inputs.mix_2d_lp_deff           = 2;                            % crystal nonlinear coefficient (in pm/V)
    shg_inputs.mix_2d_lp_deltak         = 0;                            % phase mismatch (in rad/mm)
    shg_inputs.mix_2d_lp_nt             = 64;                           % number of time steps to include (Recommend starting with 32 and varying to check convergence.)    
    shg_inputs.mix_2d_lp_auto_analyze   = false;                        % analyze the output beam (calculate M-squared, tilt, second moment, and radius of curvature)
    shg_inputs.mix_2d_lp_save_movie     = true;                         % Write the file 'mix_2d_lp_movie.mat' which contains output electric fields (x,y,t); 
    
    % Construct the data structure of inputs for the modified 2D-mix-LP
    % function for THG modeling. In the loop where we consider each element
    % in the red1_red2_energy_vec, we'll load the output fields of the SHG
    % run and set them as the input fields for THG. Because the input
    % fields take precedence over the rest of the variables in the input
    % data structure, some of the values we set here will be ignored. These
    % ignored variables include:
    %   mix_2d_lp_pulseenergy, mix_2d_lp_pulse_durations, mix_2d_lp_pulse_delays,
    %   mix_2d_lp_beam_diameters, mix_2d_lp_supergaussian_coeff,
    %   mix_2d_lp_offset_wodir, mix_2d_lp_rad_curv, mix_2d_lp_nxny,
    %   mix_2d_lp_lx_ly, and mix_2d_lp_nt.
    % ------------------------------------------------------------
    clear thg_inputs
    thg_inputs.mix_2d_lp_wavelengths        = [1064, 532, 354.7];           % wavelengths (in nm)
    thg_inputs.mix_2d_lp_ref_inds           = [1.655,1.655,1.655];          % refractive indices
    thg_inputs.mix_2d_lp_phase              = [0,0,0];                      % input phases (in radians)
    thg_inputs.mix_2d_lp_input_refl         = [0,0,0];                      % crystal input reflectivity
    thg_inputs.mix_2d_lp_output_refl        = [0, 0, 0];                    % crystal output reflectivity
    thg_inputs.mix_2d_lp_crystal_losses     = [0,0,0];                      % linear absorption (in 1/mm)
    thg_inputs.mix_2d_lp_pulseenergy        = [0,0,0];                      % input pulse energies (in joules)
    thg_inputs.mix_2d_lp_pulse_durations    = [1,1,1];                      % input pulse durations (in fwhm ns)
    thg_inputs.mix_2d_lp_beam_diameters     = [1,1;1,1;1,1];                % input beam diameters (in fwhm mm) [red1 in walkoff dir, red2 in walkoff dir, blue in walkoff dir;red1 perp to walkoff, red2 perp to walkoff, blue perp to walkoff]
    thg_inputs.mix_2d_lp_supergaussian_coeff = [1,1,1];                     % supergaussian coefficient for spatial input beam profile
    thg_inputs.mix_2d_lp_pulse_delays       = [0,0,0];                      % Pulse delay (ns)
    thg_inputs.mix_2d_lp_n2_red1            = [0,0,0];                      % red1 nonlinear refractive index (in sq. cm/W)
    thg_inputs.mix_2d_lp_n2_red2            = [0,0,0];                      % red2 nonlinear refractive index (in sq. cm/W)
    thg_inputs.mix_2d_lp_n2_blue            = [0,0,0];                      % blue nonlinear refractive index (in sq. cm/W)
    thg_inputs.mix_2d_lp_beta_red1          = [0,0,0];                      % red1 2-photon absorption cofficient (in cm/W)
    thg_inputs.mix_2d_lp_beta_red2          = [0,0,0];                      % red2 2-photon absorption cofficient (in cm/W)
    thg_inputs.mix_2d_lp_beta_blue          = [0,0,0];                      % blue 2-photon absorption cofficient (in cm/W)
    thg_inputs.mix_2d_lp_wo_angles          = [0,0,0];                      % spatial walkoff angles (in milliradians)
    thg_inputs.mix_2d_lp_offset_wodir       = [0,0,0];                      % spatial offset in walkoff direction (in mm)
    thg_inputs.mix_2d_lp_rad_curv           = [1e9,1e9,1e9;1e9,1e9,1e9].';  % input beam radii of curvature for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction] (in mm in air)
    thg_inputs.mix_2d_lp_nz                 = 64;                           % number of z integration steps
    thg_inputs.mix_2d_lp_nxny               = [64,64];                      % number of grid points in the transverse (x,y) directions
    thg_inputs.mix_2d_lp_crystal_length     = [10,1];                       % crystal length (in mm); optionally, the second parameter is number of walkoff compensating sections
    thg_inputs.mix_2d_lp_lx_ly              = [3,3];                        % size of grid in walkoff direction, perpendicular to walkoff direction (in mm)
    thg_inputs.mix_2d_lp_deff               = 2;                            % crystal nonlinear coefficient (in pm/V)
    thg_inputs.mix_2d_lp_deltak             = 0;                            % phase mismatch (in rad/mm)
    thg_inputs.mix_2d_lp_nt                 = 64;                           % number of time steps to include (Recommend starting with 32 and varying to check convergence.)
    thg_inputs.mix_2d_lp_auto_analyze       = false;
    thg_inputs.mix_2d_lp_save_movie         = true;
    thg_inputs.input_fields = []; % Included to prevent an error because 
        % in the first iteration we add the input_fields field to the data structure, but
        % the THG problem set for the next iteration will have a different set of data
        % structure fields unless this is included.
    
    % Pre-allocate some cell arrays we'll stick the output SHG waves and THG waves in
    shg_output_red1     = cell(size(red1_red2_energy_vec)); % pre-allocate a cell array to store the contents of the red1 output files in
    shg_output_red2     = cell(size(red1_red2_energy_vec)); % pre-allocate a cell array to store the contents of the red2 output files in
    shg_output_blue     = cell(size(red1_red2_energy_vec)); % pre-allocate a cell array to store the contents of the blue output files in
    thg_output_red1     = cell(size(red1_red2_energy_vec)); % pre-allocate a cell array to store the contents of the red1 output files in
    thg_output_red2     = cell(size(red1_red2_energy_vec)); % pre-allocate a cell array to store the contents of the red2 output files in
    thg_output_blue     = cell(size(red1_red2_energy_vec)); % pre-allocate a cell array to store the contents of the blue output files in
    shg_red1_energy_out = zeros(size(red1_red2_energy_vec));% pre-allocate an array of zeros to store  1064 nm output energies
    shg_red2_energy_out = zeros(size(red1_red2_energy_vec));% pre-allocate an array of zeros to store  1064 nm output energies
    shg_blue_energy_out = zeros(size(red1_red2_energy_vec));% pre-allocate an array of zeros to store   532 nm output energies
    thg_red1_energy_out = zeros(size(red1_red2_energy_vec));% pre-allocate an array of zeros to store  1064 nm output energies
    thg_red2_energy_out = zeros(size(red1_red2_energy_vec));% pre-allocate an array of zeros to store   532 nm output energies
    thg_blue_energy_out = zeros(size(red1_red2_energy_vec));% pre-allocate an array of zeros to store 354.7 nm output energies
    clear shg_problem thg_problem % shg_problem and thg_problem will be vectors of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    for K = 1:length(red1_red2_energy_vec) % loop over each element in the vector of input pulse energies
        % I had some trouble keeping track of which part of the model was
        % running, so I print some strings to the SNLO console reflecting
        % the current step and its results
        fprintf(1,'Running SHG for %.6G mJ (red)\n', ...
            2*red1_red2_energy_vec(K)*1e3); % (Remember, red1_red2_energy_vec energy gets stuck in each of the red waves
        
        % Set up the SHG part of the problem for this iteration by copying
        % the inputs specified earlier
        shg_problem(K) = shg_inputs;
        shg_problem(K).mix_2d_lp_pulseenergy(1) = red1_red2_energy_vec(K); % But be sure to set the red pulse energies to the appropriate value for this iteration
        shg_problem(K).mix_2d_lp_pulseenergy(2) = red1_red2_energy_vec(K);
        
        % call the modified SNLO 2D-MIX-LP function with the problem set,
        % and assign the returned cell array of function handles which are
        % local to that file to 'fcn_handles'
        shg_fcn_handles   = snlo_2d_mix_lp_func(shg_problem(K)); 
        shg_run_handle    = shg_fcn_handles{1};     % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        shg_accept_handle = shg_fcn_handles{2};     % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        shg_close_handle  = shg_fcn_handles{end};   % the function handle of the 'Close figure' callback is last returned
        shg_accept_handle();                        % call the 'accept' button's callback - equivalent to clicking the accept button
        shg_run_handle();                           % call the 'Run' button callback - equivalent to clicking the run button
        % can pause to look things over here if you uncomment this:
        keyboard;
        
        % After running, 2D-mix-LP saves some output files to disk, with
        % one file per wave in ascii format. The first columns are
        % times, the second are the powers at these times, the third is the
        % phase record (summed up from the transverse field profile at this
        % time), fourth through ninth columns have to do with M-squared
        % calculations (M-squared, tilt, curvature for each direction)
        shg_output_red1{K} = importdata('mix_2d_lp_red1_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        shg_output_red2{K} = importdata('mix_2d_lp_red2_data.dat'); % load the red2 2d mix lp output file into memory, stick into the output cell array
        shg_output_blue{K} = importdata('mix_2d_lp_blue_data.dat'); % load the blue 2d mix lp output file into memory, stick into the output cell array

        % Integrate the power vector in time (using the trapezoidal method,
        % with the times given as the first argument and power as second)
        shg_red1_energy_out(K) = trapz(shg_output_red1{K}.data(:,1),shg_output_red1{K}.data(:,2))*1e3; % Multiply by 1000 to get in to mJ
        shg_red2_energy_out(K) = trapz(shg_output_red2{K}.data(:,1),shg_output_red2{K}.data(:,2))*1e3;
        shg_blue_energy_out(K) = trapz(shg_output_blue{K}.data(:,1),shg_output_blue{K}.data(:,2))*1e3;
        
        % Write some output to the command window showing our results
        fprintf(1,'SHG Output Energies: %.6G mJ @ 1064nm, %.6G mJ @ 532nm (total %.6G mJ)\n',...
            shg_red1_energy_out(K) + shg_red2_energy_out(K), shg_blue_energy_out(K),...
            sum([shg_red1_energy_out(K), shg_red2_energy_out(K), shg_blue_energy_out(K)]));
        
        % Now, load the output fields from the file saved by the SHG model
        shg_output_filename = 'mix_2d_lp_movie.mat';
        shg_output_fields = load(shg_output_filename);
        
        % Make a copy of the SHG output fields, generate the THG input fields
        thg_input_fields = shg_output_fields;
        
        % In the frequency doubling stage, we split the red input energy
        % equally between red1 and red2 (i.e. we multiply the true input by
        % 1/sqrt(2) ). Now when we combine the red1 and red2 fields to form
        % the new red1 input we sum the red1 and red2 output fields and
        % multiply by 1/sqrt(2). That's whats happening in this next line.
        % If the first stage is not type 1 SHG, all of this is irrelevant.
        thg_input_fields.field_red1_xyt = sqrt(0.5) * (shg_output_fields.field_red1_xyt + shg_output_fields.field_red2_xyt); % 1064 nm
        thg_input_fields.field_red2_xyt = shg_output_fields.field_blue_xyt;     % 532 nm
        % thg_input_fields.field_blue_xyt = 0*shg_output_fields.field_blue_xyt;   % 354.7nm % Note, the current version of mlSNLO's 2D-mix-LP model errors out if all of one of the input fields is 0, just use a very small portion of other wave
        thg_input_fields.field_blue_xyt = 1e-30*shg_output_fields.field_blue_xyt;

        % As a check, we calculate the THG input pulse energies by
        % integrating the THG input fields in x,y, and t
        dx = thg_input_fields.xgrid(2) - thg_input_fields.xgrid(1);     % x-step
        dy = thg_input_fields.ygrid(2) - thg_input_fields.ygrid(1);     % y-step
        dt = thg_input_fields.tvec(2)  - thg_input_fields.tvec(1);      % time-step
        thg_r1_energy_in = dx*dy*dt*sum(sum(sum( 0.5 * c * epsilon_0 .* abs(thg_input_fields.field_red1_xyt).^2)));
        thg_r2_energy_in = dx*dy*dt*sum(sum(sum( 0.5 * c * epsilon_0 .* abs(thg_input_fields.field_red2_xyt).^2)));
        thg_bl_energy_in = dx*dy*dt*sum(sum(sum( 0.5 * c * epsilon_0 .* abs(thg_input_fields.field_blue_xyt).^2)));
        fprintf(1,'Running THG (for SHG %.6G mJ red input)\n',2*red1_red2_energy_vec(K)*1e3);
        fprintf(1,'THG Input Energies: %.6G mJ @ 1064nm, %.6G mJ @ 532nm, %.6G mJ @ 354.7nm (total %.6G mJ)\n',thg_r1_energy_in*1e3, ...
            thg_r2_energy_in*1e3, thg_bl_energy_in*1e3, (thg_r1_energy_in + thg_r2_energy_in + thg_bl_energy_in)*1e3);
        % thg_r1_energy_in should be equal to
        % shg_red1_energy_out+shg_red2_energy_out, 
        % thg_r2_energy_in should be equal to shg_blue_energy_out,
        % and thg_bl_energy_in should be zero

        thg_problem(K) = thg_inputs;
        thg_problem(K).input_fields = thg_input_fields;
        thg_problem(K).mix_2d_lp_pulseenergy(1) = (shg_red1_energy_out(K) + shg_red2_energy_out(K)); % Since we specify field profiles, these values are ignored we use the values calculated from SHG output for clarity
        thg_problem(K).mix_2d_lp_pulseenergy(2) = (shg_blue_energy_out(K));
        thg_fcn_handles     = snlo_2d_mix_lp_func(thg_problem(K)); % call the SNLO 2D-MIX-LP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        thg_run_handle      = thg_fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        thg_accept_handle   = thg_fcn_handles{2}; % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        thg_close_handle    = thg_fcn_handles{end}; % the function handle of the 'Close figure' callback is last returned
        thg_accept_handle();    % call the 'accept' button's callback - equivalent to clicking the accept button
        thg_run_handle();       % call the 'Run' button callback - equivalent to clicking the run button
        
        % can pause to look things over here if you uncomment this:
        % keyboard

        thg_output_red1{K} = importdata('mix_2d_lp_red1_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        thg_output_red2{K} = importdata('mix_2d_lp_red2_data.dat'); % load the red2 2d mix lp output file into memory, stick into the output cell array
        thg_output_blue{K} = importdata('mix_2d_lp_blue_data.dat'); % load the blue 2d mix lp output file into memory, stick into the output cell array
        
        % integrate power vectors (using the trapz function) in time to get output energies
        thg_red1_energy_out(K) = trapz(thg_output_red1{K}.data(:,1),thg_output_red1{K}.data(:,2))*1e3; % in mJ
        thg_red2_energy_out(K) = trapz(thg_output_red2{K}.data(:,1),thg_output_red2{K}.data(:,2))*1e3; % in mJ
        thg_blue_energy_out(K) = trapz(thg_output_blue{K}.data(:,1),thg_output_blue{K}.data(:,2))*1e3; % in mJ
        
        fprintf(1,'THG Output Energies: %.6G mJ @ 1064nm, %.6G mJ @ 532nm, %.6G mJ @ 354.7nm (total %.6G mJ)\n',...
            thg_red1_energy_out(K), thg_red2_energy_out(K), thg_blue_energy_out(K), ...
            sum([thg_red1_energy_out(K),thg_red2_energy_out(K),thg_blue_energy_out(K)]));
%         close_handle(); % call the 'Close figure' callback - not usually
%         necessary, and it's often useful to keep it open to look at
%         outputs we haven't collected and saved
    end
    
else
    % If you don't want to run the batch of models again, but you're
    % running this script, try to load the results of previous runs and
    % plot them
    data = load('thg_results.mat');
    red1_red2_energy_vec = data.red1_red2_energy_vec;
    shg_red1_energy_out = data.shg_red1_energy_out;
    shg_red2_energy_out = data.shg_red2_energy_out;
    shg_blue_energy_out = data.shg_blue_energy_out;
    thg_red1_energy_out = data.thg_red1_energy_out;
    thg_red2_energy_out = data.thg_red2_energy_out;
    thg_blue_energy_out = data.thg_blue_energy_out;
end

% plot the output energies we calculated in the batch loop
fh = figure;            % make new outut figure

% first, plot the first stage SHG output energies
ax = axes('Parent',fh); % make new axes in that output figure
subplot(2,1,1,ax);      % set it as the first subplot
% plot the SHG red1, red2, and blue output energies as a function of total
% input pulse energy (red1+red2) in mJ
ph = plot(ax, red1_red2_energy_vec*2e3, shg_red1_energy_out, 'rs-',...
    red1_red2_energy_vec*2e3, shg_red2_energy_out, 'bo-',...
    red1_red2_energy_vec*2e3, shg_blue_energy_out, 'g*-');
% xlabel(ax, 'SHG input energy [mJ]');     % label x-axis 
xlabel(ax, sprintf('%.1f nm input energy [mJ]',shg_inputs.mix_2d_lp_wavelengths(1)));
ylabel(ax, 'SHG output energy [mJ]');    % label y-axis 
lh = legend(ax, 'Red1 (1064 nm)', 'Red2 (1064 nm)',...
    'Blue (532 nm)', 'Location','Best'); % make legend to describe each curve on the plot

% second, plot the second stage THG output energies
ax2 = axes('Parent',fh); % make new axes in that output figure
subplot(2,1,2,ax2); % set it as the 2nd subplot
% plot the THG red1, red2, and blue output energies as a function of total
% input pulse energy (red1+red2) in mJ
ph2 = plot(ax2,red1_red2_energy_vec*2e3,thg_red1_energy_out,'rs-',...
    red1_red2_energy_vec*2e3,thg_red2_energy_out,'bo-',...
    red1_red2_energy_vec*2e3,thg_blue_energy_out,'g*-');
% xlabel(ax2,'SHG input energy [mJ]');  % label x-axis 
xlabel(ax2, sprintf('%.1f nm input energy [mJ]',shg_inputs.mix_2d_lp_wavelengths(1)));
ylabel(ax2,'THG output energy [mJ]');    % label y-axis 
lh2 = legend(ax2,'1064 nm','532 nm',...
    '354.7 nm','Location','Best'); % make legend to describe each curve on the plot

if calculate_curves 
    % If you ran the batch set of models, save the results of the batch
    % script
    rmfield(thg_problem,'input_fields');
    saveas(fh,'thg_results.fig');
    save('thg_results.mat', 'c', 'epsilon_0', 'red1_red2_energy_vec', ...
        'shg_inputs', 'thg_inputs', 'shg_output_red1', 'shg_output_red2', ...
        'shg_output_blue', 'thg_output_red1', 'thg_output_red2', ...
        'thg_output_blue', 'shg_red1_energy_out', 'shg_red2_energy_out', ...
        'shg_blue_energy_out', 'thg_red1_energy_out', 'thg_red2_energy_out', ...
        'thg_blue_energy_out', 'shg_output_filename', ...
        'dx', 'dy', 'dt', 'thg_r1_energy_in', ...
        'thg_r2_energy_in', 'thg_bl_energy_in', 'shg_problem', 'thg_problem');
end