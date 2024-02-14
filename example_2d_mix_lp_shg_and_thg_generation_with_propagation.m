%% 2D-mix-LP: cascaded SHG to THG devices with propagation between crystals
%
% Jesse Smith [jesse.smith@as-photonics.com]
% last modified Feb 2, 2024
% 
% Example of how to use 2D-mix-LP (diffractive long pulse mixing) to model a two stage
% device, with the second using the outputs of the first.
%
% This script does not work with mlSNLO versions prior to 75.2.0
% (from April, 2021.)
% It will execute without throwing any errors but 2D-mix-LP will not 
% correctly use the input field profiles specified here and instead just 
% uses the input box values to create Gaussian or super Gaussian input fields.
% This is an example of two-stage third-harmonic generation (using two crystals in series)
% with simple, free-space propagation between them. The first crystal is second-harmonic
% generation, and the second third-harmonic generation. This example is a simple batch of
% runs where the SHG red input energy is varied.
%
% Note: This example script doesn't generate any fields to use as input fields to the
% modified 2D-mix-LP model directly, but rather uses the modified 2D-mix-LP model to
% generate them for the from first stage SHG, then loads the output fields from the SHG
% model run from disk and uses them as the input fields for the THG (after some
% manipulation.)
%
% If you wish to generate your own electric field arrays to use as inputs
% for the modified 2D-mix-LP function, the text file with the name
% 'notes_on_snlo_2d_mix_lp_with_electric_field_inputs.txt' contains further
% information.
%
% SHG: 1064 nm -> 532 nm
%   * In SNLO we do SHG by splitting the red input into both
%   red waves with equal powers (with identical values for other properties)
% THG: 1064 nm + 532 nm -> 354.7 nm
%   * We will collect both red output waves of the SHG and combine them for
%   one of the red THG waves, and use the 532 nm output of SHG as
%   the other red wave.
%
%   
%  Results and other thoughts
%  -Using 1064 nm input energies from 0.1 to 20 mJ, we find that the conversion efficiency
%   peaks around input energy of 1.5 mJ; peak efficiency is about 58%, with 870 mJ output at
%   354.7 nm.
% - Increasing output energy with more input energy still boosts the 354.7 nm
%   energy but with slightly decreased efficiency. 
% - If you look at the 1064 nm output from the THG crystal you can see there's been some 
%   substantial back-conversion at the center of the beam. This is a problem that will get
%   worse as you increase input energy. The beam quality of the THG light is still good,
%   but starting to deteriorate; at 2 mJ input energy M^2 = 1.10 
% - These are big beams. The actual example of propagating in free space between the crystals is
%   of dubious utility, but this example should help you start designing your own optimization.
%   If you care to, you could shrink the beams and use something other than collimated input
%   light. You could even add a thin lens between the two crystals if you propagate part way,
%   add phase curvature, and propagate the rest of the way.
% 


c = 3e8;
epsilon_0 = 8.85e-12;

calculate_curves = true; % generate new data? (after running through this once and saving the results to file, you can set this value to false to load the file and process & plot that data)
if calculate_curves
    %% Inputs
    % Construct the vector of red energies we'll consider (this energy goes
    % into each of the red waves, so the total input energy will be twice
    % this value)
    red1_red2_energy_vec = linspace(0.05,1.00,11)*1e-3; % pulse energies of red1 and red2 waves in joules (total red input energy is twice this value)
    
    % distance between first-stage shg crystal output face and second-stage thg crystal input face (we'll propagate
    % the waves this distance in free space)
    shg_thg_distance = 1.250; % distance [m] - propagate 1 meter between SHG crystal output and THG crystal input 
    
    %% Gui inputs
    % Construct the data structure of inputs for the modified 2D-mix-LP
    % function for SHG modeling
    % ------------------------------------------------------------
    clear shg_inputs
    shg_inputs.mix_2d_lp_wavelengths    = [1064,1064,532];              % wavelengths (nm)
    shg_inputs.mix_2d_lp_ref_inds       = [1.655,1.655,1.655];          % refractive indices
    shg_inputs.mix_2d_lp_phase          = [0,0,0];                      % input phases (radians)
    shg_inputs.mix_2d_lp_input_refl     = [0,0,0];                      % crystal input reflectivity
    shg_inputs.mix_2d_lp_output_refl    = [0,0,0];                      % crystal input reflectivity
    shg_inputs.mix_2d_lp_crystal_losses = [0,0,0];                      % linear absorption (1/mm)
    shg_inputs.mix_2d_lp_pulseenergy    = [10e-2,10e-2,0];              % input pulse energies (joules)
    shg_inputs.mix_2d_lp_pulse_durations = [1,1;1,1;1,1];               % input pulse durations (fwhm ns) -- note that this is a 3x2 array; the first column is fwhm pulse duration in ns, the second is the temporal supergaussian coefficient (use value of 1 for lowest order Gaussian, larger values are more flat-topped).
    shg_inputs.mix_2d_lp_beam_diameters = [1,1;1,1;1,1];                % input beam diameters in each transverse dimension (fwhm mm) [red1 in walkoff dir, red2 in walkoff dir, blue in walkoff dir;red1 perp to walkoff, red2 perp to walkoff, blue perp to walkoff]
    shg_inputs.mix_2d_lp_supergaussian_coeff = [1,1,1];                 % supergaussian coefficient for spatial input beam profile
    shg_inputs.mix_2d_lp_pulse_delays   = [0,0,0];                      % Pulse delay (ns)
    shg_inputs.mix_2d_lp_n2_red1        = [0,0,0];                      % red1 nonlinear refractive index (sq. cm/W)
    shg_inputs.mix_2d_lp_n2_red2        = [0,0,0];                      % red2 nonlinear refractive index (sq. cm/W)
    shg_inputs.mix_2d_lp_n2_blue        = [0,0,0];                      % blue nonlinear refractive index (sq. cm/W)
    shg_inputs.mix_2d_lp_beta_red1      = [0,0,0];                      % red1 2-photon absorption cofficient (cm/W)
    shg_inputs.mix_2d_lp_beta_red2      = [0,0,0];                      % red2 2-photon absorption cofficient (cm/W)
    shg_inputs.mix_2d_lp_beta_blue      = [0,0,0];                      % blue 2-photon absorption cofficient (cm/W)
    shg_inputs.mix_2d_lp_wo_angles      = [0,0,0];                      % spatial walkoff angles (milliradians)
    shg_inputs.mix_2d_lp_offset_wodir   = [0,0,0];                      % spatial offset in walkoff direction (mm)
    shg_inputs.mix_2d_lp_rad_curv       = [1e9,1e9,1e9;1e9,1e9,1e9].';  % input beam radii of curvature for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction] (mm in air)
    shg_inputs.mix_2d_lp_nz             = 64;                           % number of z integration steps
    shg_inputs.mix_2d_lp_nxny           = [64,64];                      % number of grid points in the transverse (x,y) directions
    shg_inputs.mix_2d_lp_crystal_length = [10,1];                       % crystal length (mm); optionally, the second parameter is number of walkoff compensating sections
    shg_inputs.mix_2d_lp_lx_ly          = [4,4];                        % size of grid in walkoff direction, perpendicular to walkoff direction (mm)
    shg_inputs.mix_2d_lp_deff           = 2;                            % crystal nonlinear coefficient (pm/V)
    shg_inputs.mix_2d_lp_deltak         = 0;                            % phase mismatch (rad/mm)
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
    thg_inputs.mix_2d_lp_wavelengths        = [1064, 532, 354.7];           % wavelengths (nm)
    thg_inputs.mix_2d_lp_ref_inds           = [1.655,1.655,1.655];          % refractive indices
    thg_inputs.mix_2d_lp_phase              = [0,0,0];                      % input phases (radians)
    thg_inputs.mix_2d_lp_input_refl         = [0,0,0];                      % crystal input reflectivity
    thg_inputs.mix_2d_lp_output_refl        = [0, 0, 0];                    % crystal output reflectivity
    thg_inputs.mix_2d_lp_crystal_losses     = [0,0,0];                      % linear absorption (1/mm)
    thg_inputs.mix_2d_lp_pulseenergy        = [0,0,0];                      % input pulse energies (joules)
    thg_inputs.mix_2d_lp_pulse_durations    = [1,1,1];                      % input pulse durations (fwhm ns)
    thg_inputs.mix_2d_lp_beam_diameters     = [1,1;1,1;1,1];                % input beam diameters (fwhm mm) [red1 in walkoff dir, red2 in walkoff dir, blue in walkoff dir;red1 perp to walkoff, red2 perp to walkoff, blue perp to walkoff]
    thg_inputs.mix_2d_lp_supergaussian_coeff = [1,1,1];                     % supergaussian coefficient for spatial input beam profile
    thg_inputs.mix_2d_lp_pulse_delays       = [0,0,0];                      % Pulse delay (ns)
    thg_inputs.mix_2d_lp_n2_red1            = [0,0,0];                      % red1 nonlinear refractive index (sq. cm/W)
    thg_inputs.mix_2d_lp_n2_red2            = [0,0,0];                      % red2 nonlinear refractive index (sq. cm/W)
    thg_inputs.mix_2d_lp_n2_blue            = [0,0,0];                      % blue nonlinear refractive index (sq. cm/W)
    thg_inputs.mix_2d_lp_beta_red1          = [0,0,0];                      % red1 2-photon absorption cofficient (cm/W)
    thg_inputs.mix_2d_lp_beta_red2          = [0,0,0];                      % red2 2-photon absorption cofficient (cm/W)
    thg_inputs.mix_2d_lp_beta_blue          = [0,0,0];                      % blue 2-photon absorption cofficient (cm/W)
    thg_inputs.mix_2d_lp_wo_angles          = [0,0,0];                      % spatial walkoff angles (milliradians)
    thg_inputs.mix_2d_lp_offset_wodir       = [0,0,0];                      % spatial offset in walkoff direction (mm)
    thg_inputs.mix_2d_lp_rad_curv           = [1e9,1e9,1e9;1e9,1e9,1e9].';  % input beam radii of curvature for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction] (mm in air)
    thg_inputs.mix_2d_lp_nz                 = 64;                           % number of z integration steps
    thg_inputs.mix_2d_lp_nxny               = [64,64];                      % number of grid points in the transverse (x,y) directions
    thg_inputs.mix_2d_lp_crystal_length     = [10,1];                       % crystal length (mm); optionally, the second parameter is number of walkoff compensating sections
    thg_inputs.mix_2d_lp_lx_ly              = [4,4];                        % size of grid in walkoff direction, perpendicular to walkoff direction (mm)
    thg_inputs.mix_2d_lp_deff               = 2;                            % crystal nonlinear coefficient (pm/V)
    thg_inputs.mix_2d_lp_deltak             = 0;                            % phase mismatch (rad/mm)
    thg_inputs.mix_2d_lp_nt                 = 64;                           % number of time steps to include (Recommend starting with 32 and varying to check convergence.)
    thg_inputs.mix_2d_lp_auto_analyze       = true;
    thg_inputs.mix_2d_lp_save_movie         = true;
    thg_inputs.input_fields = []; % Included to prevent an error because 
        % in the first iteration we add the input_fields field to the data structure, but
        % the THG problem set for the next iteration will have a different set of data
        % structure fields unless this is included.
    
    %% Pre-allocate some cell arrays we'll stick the output SHG waves and THG waves in
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

    thg_red1_pulse_m2s  = zeros(length(red1_red2_energy_vec),2);% pre-allocate an array of zeros to store pulse beam qualities (fluence-based) in x & y
    thg_red2_pulse_m2s  = zeros(length(red1_red2_energy_vec),2);% pre-allocate an array of zeros to store pulse beam qualities (fluence-based) in x & y
    thg_blue_pulse_m2s  = zeros(length(red1_red2_energy_vec),2);% pre-allocate an array of zeros to store pulse beam qualities (fluence-based) in x & y
    thg_red1_m2_vs_t    = zeros(length(red1_red2_energy_vec), thg_inputs.mix_2d_lp_nt+1, 2); % pre-allocate an array of zeros to store beam qualities vs t  in x & y
    thg_red2_m2_vs_t    = zeros(length(red1_red2_energy_vec), thg_inputs.mix_2d_lp_nt+1, 2); % pre-allocate an array of zeros to store beam qualities vs t  in x & y
    thg_blue_m2_vs_t    = zeros(length(red1_red2_energy_vec), thg_inputs.mix_2d_lp_nt+1, 2); % pre-allocate an array of zeros to store beam qualities vs t  in x & y

    clear shg_problem thg_problem % shg_problem and thg_problem will be vectors of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting

    %% Run
    % (For each energy in the list red1_red2_energy_vec run the following)
    for K = 1:length(red1_red2_energy_vec) % loop over each element in the vector of input pulse energies
        % Print some strings to the MATLAB commandline (console) to display the progress so far
        fprintf(1,'Run %i of %i: SHG stage for %.6G mJ red1 input energy\n', ...
            K, length(red1_red2_energy_vec), 2*red1_red2_energy_vec(K)*1e3); % (Remember, red1_red2_energy_vec energy gets stuck in each of the red waves
        
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
        
        %% SHG outputs
        % After running, 2D-mix-LP saves some output files to disk, with
        % one file per wave in ascii format. The first columns are
        % times, the second are the powers at these times, the third is the
        % phase record (summed up from the transverse field profile at this
        % time), fourth through ninth columns have to do with M-squared
        % calculations (M-squared, tilt, curvature for each direction)
        shg_output_red1{K} = importdata('mix_2d_lp_red1_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        shg_output_red2{K} = importdata('mix_2d_lp_red2_data.dat'); % load the red2 2d mix lp output file into memory, stick into the output cell array
        shg_output_blue{K} = importdata('mix_2d_lp_blue_data.dat'); % load the blue 2d mix lp output file into memory, stick into the output cell array

        % Integrate the power vector in time (using the trapezoidal method, with the times
        % given as the first argument and power as second); you could instead give trapz 1
        % input argument, the powers, and multiply by dt
        shg_red1_energy_out(K) = trapz(shg_output_red1{K}.data(:,1),shg_output_red1{K}.data(:,2));
        shg_red2_energy_out(K) = trapz(shg_output_red2{K}.data(:,1),shg_output_red2{K}.data(:,2));
        shg_blue_energy_out(K) = trapz(shg_output_blue{K}.data(:,1),shg_output_blue{K}.data(:,2));

        shg_red1_energy_out(K) = shg_red1_energy_out(K)*1e3; % Multiply by 1000 to convert J to mJ
        shg_red2_energy_out(K) = shg_red2_energy_out(K)*1e3; % Multiply by 1000 to convert J to mJ
        shg_blue_energy_out(K) = shg_blue_energy_out(K)*1e3; % Multiply by 1000 to convert J to mJ
        
        % Write some output to the command window showing our results
        fprintf(1,'SHG Output Energies: %.6G mJ @ 1064nm, %.6G mJ @ 532nm (total %.6G mJ)\n',...
            shg_red1_energy_out(K) + shg_red2_energy_out(K), shg_blue_energy_out(K),...
            sum([shg_red1_energy_out(K), shg_red2_energy_out(K), shg_blue_energy_out(K)]));
        
        % Now, load the output fields from the file saved by the SHG model
        shg_output_filename = 'mix_2d_lp_movie.mat';
        shg_output_fields = load(shg_output_filename);

        %% THG inputs
        % From the SHG output fields, generate the THG input fields
        thg_input_fields = shg_output_fields;
        
        % Because we have type 1 SHG in the first stage (the frequency doubling stage),
        % we split the red input energy equally between red1 and red2 (i.e. we multiply
        % the fields input by 1/sqrt(2) ). Now when we combine the red1 and red2 fields to
        % form the new red1 input we sum the red1 and red2 output fields and multiply by
        % 1/sqrt(2). If the first stage is not type 1 SHG, all of this is irrelevant.
        thg_input_fields.field_red1_xyt = sqrt(0.5) * (shg_output_fields.field_red1_xyt + shg_output_fields.field_red2_xyt); % 1064 nm
        thg_input_fields.field_red2_xyt = shg_output_fields.field_blue_xyt;     % 532 nm
        % Start with zero in the blue wave
        thg_input_fields.field_blue_xyt = zeros(size(shg_output_fields.field_blue_xyt));   % 354.7nm
        thg_input_fields.field_blue_xyt = thg_input_fields.field_blue_xyt + 1e-20; % sometimes values of 0 mess things up.

        %% Propagate between crystals
        %
        % We have a MATLAB function file simple_freespace_propagate.m which will propagate
        % a field (or series of fields, as in a pulsed wave) through freespace by some
        % distance specified using the FFT beam propagation method. This method of
        % propagating beams has the potential to run into trouble if your fields expand
        % (or propagate, in the case of tilted beams) to not be fully contained by the
        % spatial grid; the fft handles this by sticking the light on the face opposite where
        % it exited. This is a sort of aliasing effect, I think, and can produce some
        % unexpected interference effects. This propagator is "simple" in the sense that
        % it doesn't do anything to prevent this (for instance, you might pad the spatial
        % grid with more points out the faces, or do some resampling thing).
        % 
        % You can read the code in propagate_freespace_class.m for the mechanics about
        % what it actually does and what input variables you must provide in your function
        % call.
        % 

        % Propagate beams in free space the shg_thg_distance propagate the 1064 nm and the
        % 532 nm beams (which are called the red1 and red2 waves now in our THG modeling):

        % Old method:
        % red1_propagated_field = simple_freespace_propagate(thg_inputs.mix_2d_lp_wavelengths(1)*1e-9, thg_input_fields.xgridmat, thg_input_fields.ygridmat, thg_input_fields.field_red1_xyt, shg_thg_distance);
        % red2_propagated_field = simple_freespace_propagate(thg_inputs.mix_2d_lp_wavelengths(2)*1e-9, thg_input_fields.xgridmat, thg_input_fields.ygridmat, thg_input_fields.field_red2_xyt, shg_thg_distance);
        propgtobj = propagate_freespace_class;
        propgtobj.wavelength    = thg_inputs.mix_2d_lp_wavelengths(1)*1e-9;
        propgtobj.xgridmat      = thg_input_fields.xgridmat;
        propgtobj.ygridmat      = thg_input_fields.ygridmat;
        propgtobj.input_field   = thg_input_fields.field_red1_xyt;
        red1_propagated_field   = propgtobj.propagate(shg_thg_distance);
        % keyboard;
        % Redo for red2 wave
        propgtobj.wavelength    = thg_inputs.mix_2d_lp_wavelengths(2)*1e-9;
        propgtobj.input_field   = thg_input_fields.field_red2_xyt;
        red2_propagated_field   = propgtobj.propagate(shg_thg_distance);

        %% Testing propagated beams
        % The propagator might produce bad results; a couple of checks you should do are
        % check that the total energy or power is unchanged from propagation, and try to
        % check that the propagated fields haven't expanded enough to overflow the spatial
        % grid.
        test_beams = true;

        if test_beams
            % Keep in mind that the thg_input_fields.field_red1_xyt and
            % thg_input_fields.field_red2_xyt fields have not yet been propagated, but the
            % red1_propagated_field and red2_propagated_field have. Also, that for the THG
            % stage the red1 is all the 1064 nm light, red2 is all the 532 nm light, and
            % blue is all the 354.7 nm light

            % % calculate second moment widths before and after propagation:
            Nx = size(thg_input_fields.field_red1_xyt, 1);
            Ny = size(thg_input_fields.field_red1_xyt, 2);
            Nt = size(thg_input_fields.field_red1_xyt, 3);
            dx = thg_input_fields.xgrid(2) - thg_input_fields.xgrid(1);     % x-step
            dy = thg_input_fields.ygrid(2) - thg_input_fields.ygrid(1);     % y-step
            dt = thg_input_fields.tvec(2)  - thg_input_fields.tvec(1);      % time-step

            x2d     = thg_input_fields.xgridmat;
            y2d     = thg_input_fields.ygridmat;
            tvec    = thg_input_fields.tvec;
            % Replicate the 2D x and y matrices to have nt points in the third dimension for simplicity:
            x3 = repmat(thg_input_fields.xgridmat,[1,1,size(red1_propagated_field,3)]);
            y3 = repmat(thg_input_fields.ygridmat,[1,1,size(red1_propagated_field,3)]);

            % Here are some anonymous functions to calculate the beams' centroids, and the
            % beams' second moments: You call the xBar_f & yBar_f functions first to get centroids,
            % with 1 input argument (the input electric field array); Calling the second moment 
            % calculator next will need the input electric field array and the centroid value:
            xBar_f  = @(E) squeeze(sum(sum(x3.*abs(E).^2)) ./ sum(sum(abs(E).^2,1),2));
            yBar_f  = @(E) squeeze(sum(sum(y3.*abs(E).^2)) ./ sum(sum(abs(E).^2,1),2));
            wx_f    = @(E, xBar) sqrt( squeeze( 4*sum(sum((x3-repmat(permute(xBar,[3,2,1]),[Nx, Ny, 1])).^2 .* abs(E).^2))  ./ sum(sum(abs(E).^2,1),2)));
            wy_f    = @(E, yBar) sqrt( squeeze( 4*sum(sum((y3-repmat(permute(yBar,[3,2,1]),[Nx, Ny, 1])).^2 .* abs(E).^2))  ./ sum(sum(abs(E).^2,1),2)));

            output_red1_xBar = xBar_f(thg_input_fields.field_red1_xyt);
            output_red1_yBar = yBar_f(thg_input_fields.field_red1_xyt);
            output_red1_wx   = wx_f(thg_input_fields.field_red1_xyt, output_red1_xBar);
            output_red1_wy   = wy_f(thg_input_fields.field_red1_xyt, output_red1_yBar);

            input_red1_xBar  = xBar_f(red1_propagated_field);
            input_red1_yBar  = yBar_f(red1_propagated_field);
            input_red1_wx    = wx_f(red1_propagated_field, input_red1_xBar);
            input_red1_wy    = wy_f(red1_propagated_field, input_red1_yBar);
            
            output_red2_xBar =  xBar_f(thg_input_fields.field_red2_xyt);
            output_red2_yBar =  yBar_f(thg_input_fields.field_red2_xyt);
            output_red2_wx   =  wx_f(thg_input_fields.field_red2_xyt, output_red2_xBar);
            output_red2_wy   =  wy_f(thg_input_fields.field_red2_xyt, output_red2_yBar);

            input_red2_xBar  =  xBar_f(red2_propagated_field);
            input_red2_yBar  =  yBar_f(red2_propagated_field);
            input_red2_wx    =  wx_f(red2_propagated_field, input_red2_xBar);
            input_red2_wy    =  wy_f(red2_propagated_field, input_red2_yBar);
            
            % THG Red1 and Red2 (1064nm and 532 nm) beam powers and energies prior to propagating
            P_vs_t_R1_SHG_out = squeeze(dx*dy*0.5*epsilon_0*c*trapz(trapz(abs(thg_input_fields.field_red1_xyt).^2,1),2));
            E_R1_SHG_out = dt * trapz(P_vs_t_R1_SHG_out);
            P_vs_t_R2_SHG_out = squeeze(dx*dy*0.5*epsilon_0*c*trapz(trapz(abs(thg_input_fields.field_red2_xyt).^2,1),2));
            E_R2_SHG_out = dt * trapz(P_vs_t_R2_SHG_out);

            % After propagating:
            P_vs_t_R1_THG_in = squeeze(dx*dy*0.5*epsilon_0*c*trapz(trapz(abs(red1_propagated_field).^2,1),2));
            E_R1_THG_in = dt * trapz(P_vs_t_R1_THG_in);
            P_vs_t_R2_THG_in = squeeze(dx*dy*0.5*epsilon_0*c*trapz(trapz(abs(red2_propagated_field).^2,1),2));
            E_R2_THG_in = dt * trapz(P_vs_t_R2_THG_in);

            thg_input_fields.field_red1_xyt = red1_propagated_field;
            thg_input_fields.field_red2_xyt = red2_propagated_field;

            % Fluences for each field (1064 nm and 532 nm) at SHG output face
            F_or1 = dt*0.5*epsilon_0*c*trapz(abs(thg_input_fields.field_red1_xyt).^2, 3);
            F_or2 = dt*0.5*epsilon_0*c*trapz(abs(thg_input_fields.field_red2_xyt).^2, 3);

            % Fluences for each field (1064 nm and 532 nm) after propagating to TGH input face
            F_r1 = dt*0.5*epsilon_0*c*trapz(abs(red1_propagated_field).^2, 3);
            F_r2 = dt*0.5*epsilon_0*c*trapz(abs(red2_propagated_field).^2, 3);

        end
        % place the newly propagated fields into the thg_input_fields structure
        thg_input_fields.field_red1_xyt = red1_propagated_field;
        thg_input_fields.field_red2_xyt = red2_propagated_field;

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

        %% Run the THG stage

        thg_problem(K) = thg_inputs;
        thg_problem(K).input_fields = thg_input_fields;
        thg_problem(K).mix_2d_lp_pulseenergy(1) = (shg_red1_energy_out(K) + shg_red2_energy_out(K)); % these are ignored but set to the values they are supposed to be to avoid confusion
        thg_problem(K).mix_2d_lp_pulseenergy(2) = (shg_blue_energy_out(K));
        thg_fcn_handles = snlo_2d_mix_lp_func(thg_problem(K)); % call the SNLO 2D-MIX-LP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        thg_run_handle = thg_fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        thg_accept_handle = thg_fcn_handles{2}; % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        thg_close_handle = thg_fcn_handles{end}; % the function handle of the 'Close figure' callback is last returned
        thg_accept_handle();    % call the 'accept' button's callback - equivalent to clicking the accept button
        thg_run_handle();       % call the 'Run' button callback - equivalent to clicking the run button
        
        thg_output_red1{K} = importdata('mix_2d_lp_red1_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        thg_output_red2{K} = importdata('mix_2d_lp_red2_data.dat'); % load the red2 2d mix lp output file into memory, stick into the output cell array
        thg_output_blue{K} = importdata('mix_2d_lp_blue_data.dat'); % load the blue 2d mix lp output file into memory, stick into the output cell array
        
        moviedata = matfile('mix_2d_lp_movie.mat');
        
        M2_1064nm_pulse_xy = moviedata.red1_pulse_m2_output;
        M2_1064nm_pulse_xy = M2_1064nm_pulse_xy.crystal_output_face.pulse_msquaredxy;
        
        M2_512nm_pulse_xy = moviedata.red2_pulse_m2_output;
        M2_512nm_pulse_xy = M2_512nm_pulse_xy.crystal_output_face.pulse_msquaredxy;
        
        M2_354nm_pulse_xy = moviedata.blue_pulse_m2_output;
        M2_354nm_pulse_xy = M2_354nm_pulse_xy.crystal_output_face.pulse_msquaredxy;
        M2_vs_t = moviedata.xyt_beam_qual;
        
        thg_red1_pulse_m2s(K,:) = M2_1064nm_pulse_xy;
        thg_red2_pulse_m2s(K,:) = M2_512nm_pulse_xy;
        thg_blue_pulse_m2s(K,:) = M2_354nm_pulse_xy;
        thg_red1_m2_vs_t(K,:,:) = M2_vs_t.red1.crystal_output_face.msquaredxy.';
        thg_red2_m2_vs_t(K,:,:) = M2_vs_t.red2.crystal_output_face.msquaredxy.';
        thg_blue_m2_vs_t(K,:,:) = M2_vs_t.blue.crystal_output_face.msquaredxy.';

        % integrate power vectors (using the trapz function) in time to get output energies
        thg_red1_energy_out(K) = trapz(thg_output_red1{K}.data(:,1),thg_output_red1{K}.data(:,2))*1e3; % in mJ
        thg_red2_energy_out(K) = trapz(thg_output_red2{K}.data(:,1),thg_output_red2{K}.data(:,2))*1e3; % in mJ
        thg_blue_energy_out(K) = trapz(thg_output_blue{K}.data(:,1),thg_output_blue{K}.data(:,2))*1e3; % in mJ
        
        fprintf(1,'THG Output Energies: %.6G mJ @ 1064nm, %.6G mJ @ 532nm, %.6G mJ @ 354.7nm (total %.6f mJ)\n\n',...
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
subplot(2,2,1,ax);      % set it as the first subplot
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
subplot(2,2,2,ax2); % set it as the 2nd subplot
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

% third, plot the second stage THG output beam quality vs input energy
ax3 = axes('Parent',fh); % make new axes in that output figure
subplot(2,2,3,ax3); % set it as the 2nd subplot
% plot the THG red1, red2, and blue output energies as a function of total
% input pulse energy (red1+red2) in mJ
ph3 = plot(ax3, red1_red2_energy_vec*2e3, thg_red1_pulse_m2s,'rs-',...
    red1_red2_energy_vec*2e3, thg_red2_pulse_m2s,'bo-',...
    red1_red2_energy_vec*2e3, thg_blue_pulse_m2s,'g*-');
% xlabel(ax2,'SHG input energy [mJ]');  % label x-axis 
xlabel(ax3, sprintf('%.1f nm input energy [mJ]',shg_inputs.mix_2d_lp_wavelengths(1)));
ylabel(ax3,'THG output pulse beam quality');    % label y-axis 
lh3 = legend(ax3,'1064 nm','532 nm',...
    '354.7 nm','Location','Best'); % make legend to describe each curve on the plot

% plot efficiency: third-harmonic output energy / energy in
ax4 = axes('parent', fh);
subplot(2,2,4,ax4);
ph4 = plot(ax4, red1_red2_energy_vec*2e3, thg_blue_energy_out ./ (red1_red2_energy_vec*2e3)*100);
xlabel(ax4, sprintf('%.1f nm input energy [mJ]',shg_inputs.mix_2d_lp_wavelengths(1)));
ylabel(ax4, 'Conversion efficiency [%]');

if calculate_curves 
    % If you ran the batch set of models, save the results of the batch script 
    %  the figure saved to thg_results_xx.fig, 
    %  and .mat file of calculation results to the_results_xx.mat, 
    % where xx is the date & time in year, month, day, hour, minute order
    tstring = [datestr(now,'yymmdd'),'_',datestr(now,'hhMM')];
    figurefilename = ['thg_results_',tstring,'.fig'];
    rmfield(thg_problem,'input_fields');
    saveas(fh,figurefilename);
    savefilename = ['thg_results',tstring,'.mat'];
    save(savefilename, 'c', 'epsilon_0', 'red1_red2_energy_vec', ...
        'shg_inputs', 'thg_inputs', 'shg_output_red1', 'shg_output_red2', ...
        'shg_output_blue', 'thg_output_red1', 'thg_output_red2', ...
        'thg_output_blue', 'shg_red1_energy_out', 'shg_red2_energy_out', ...
        'shg_blue_energy_out', 'thg_red1_energy_out', 'thg_red2_energy_out', ...
        'thg_blue_energy_out', 'shg_output_filename', ...
        'dx', 'dy', 'dt', 'thg_r1_energy_in', ...
        'thg_r2_energy_in', 'thg_bl_energy_in', 'shg_problem', 'thg_problem');
end



