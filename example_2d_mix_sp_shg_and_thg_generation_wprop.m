% 2D-mix-SP example demonstrating how to use the output of one 2D-mix-SP run as the input
% of another, with freespace propagation between the two nonlinear crystals. 
%
% The 2D-mix-SP model can be fed user-defined custom input electric field distribution by
% including in the input data structure (passed as the single input argument when calling
% snlo_2d_mix_sp_func) an object of class mix_2d_sp_input_fields in the field named
% input_fields. (See the class definition file mix_2d_sp_input_fields_class.m, and its
% documentation summary shown for the command 'doc mix_2d_sp_input_fields' for some more
% details). Briefly, the mix_2d_sp_input_fields_class class has properties xgrid, ygrid,
% field_red1_xyt, field_red2_xyt, field_blue_xyt, and tgrid. xgrid and ygrid are
% meshgrids, two-dimensional position grids (in meters) where xgrid is increasing in the
% first dimension and constant in the second, and ygrid is constant in the first dimension
% and increasing in the second. The properties field_red1_xyt, field_red2_xyt, and
% field_blue_xyt are complex-valued 3-dimensional arrays of electric fields (sized nx x ny
% x nt) with units of V/m. The shape and phase of these fields are constant but amplitude
% is scaled to the input energy specified in the gui input boxes for input energy or
% through the field named mix_2d_sp_pulseenergy in the inputs data structure.
%
% This is an example of two-stage third-harmonic generation (using two crystals in series)
% without propagation between them. The first crystal is second-harmonic generation, and
% the second third-harmonic generation. This example is a simple batch of runs where the
% SHG red input energy is varied.
%
% Note: This example script doesn't generate any fields to use as input fields to the
% modified 2D-mix-SP model directly, but rather uses the 2D-mix-SP model to generate them
% for the from first stage SHG, then loads the output fields from the SHG model run from
% disk and uses them as the input fields for the THG after some manipulation.
%
%
% SHG: 1064 nm -> 532 nm
%   * In SNLO we do type I second harmonic generation by splitting the red input into 
%   both red waves with equal powers and identical values for other properties
% THG: 1064 nm + 532 nm -> 354.7 nm
%   * We will collect both red output waves of the SHG and combine them for
%   one of the red THG waves, and use the 532 nm output of SHG as
%   the other red wave.
%
%   

c = 3e8;
epsilon_0 = 8.85e-12;

% Make vector of red energies we'll consider (this energy goes
% into each of the red waves, so the total input energy will be twice
% this value)

red1_red2_energy_vec = linspace(0.1, 1, 4)*10e-6; % energy in Joules

%% Construct the data structure of inputs for the modified 2D-mix-SP function for SHG modeling. 
% Note the red energy is split between red1 and red2 waves. 
% 
% Most data structure fields are three valued vectors, with one value each for the red1,
% red2, and blue waves. There are also three fields that are 3x2 arrays:
% mix_2d_sp_pulse_durations, mix_2d_sp_beam_diameters, and mix_2d_sp_rad_curv. 
%  mix_2d_sp_pulse_durations are the pulse durations have full-width, half-max durations
%   in the first column, and the second column is the temporal supergaussisan coefficient.
%  mix_2d_sp_beam_diameters are full-width, half-max diameters for each wave; the first
%   column values are diameters in the direction of walk-off (sometimes called x elsewhere
%   in documentation), and the second column is perpendicular to that (also called y).
%  mix_2d_sp_rad_curv are radii of curvature (with positive values focusing); first column
%  in the x direction, second column y direction
% 
% Also the field mix_2d_sp_pulse_chirps is a 3x3 array, with the first column specifying
% linear chirp, second column quadratic chirp, and the third column cubic chirp. Units are THz/ps, THz/ps^2, and THz/ps^3.
%
% Some fields are vectors of length 2, and some scalars; individual comments by each field
% below describe their contents.
%  
% Note that if custom input electric fields are provided (which is done by placing an
%  object of class mix_2d_sp_input_fields_class as a field named .input_fields in the
%  inputs data structure) many of the values of in the inputs data structure are ignored,
%  such as: mix_2d_sp_pulse_durations, mix_2d_sp_pulse_delays, mix_2d_sp_pulse_chirps,
%  mix_2d_sp_beam_diameters, mix_2d_sp_supergaussian_coeff, mix_2d_sp_offset_wodir,
%  mix_2d_sp_rad_curv, mix_2d_sp_nt, mix_2d_sp_nxny, and mix_2d_sp_lx_ly. The size and
%  extent of the x&y transverse grids and temporal grid are set in the
%  mix_2d_sp_input_fields_class object, and the other structure fields are used in
%  constructing the transverse & temporal electric field profiles.

clear shg_inputs     % Clear in case shg_inputs is already defined
shg_inputs.mix_2d_sp_wavelengths    = [1064,1064,532];          % wavelengths (in nm)
shg_inputs.mix_2d_sp_ref_inds       = [2.142,2.142,2.184];      % refractive indices
shg_inputs.mix_2d_sp_gvi            = [2.1875,2.1875,2.2811];   % group velocity indices
shg_inputs.mix_2d_sp_gdd            = [88,88,384];              % group delay dispersion (in fs^2/mm)
shg_inputs.mix_2d_sp_phase          = [0,0,0];                  % phases at input (in rad)
shg_inputs.mix_2d_sp_input_refl     = [0,0,0];                  % input face power reflectivity (0-1)
shg_inputs.mix_2d_sp_output_refl    = [0,0,0];                  % output face power reflectivity (0-1)
shg_inputs.mix_2d_sp_crystal_losses = [0,0,0];                  % linear power absorption coefficient (in 1/mm)
shg_inputs.mix_2d_sp_n2_red1        = [0,0,0];                  % red1 nonlinear refractive index (in sq cm/W)
shg_inputs.mix_2d_sp_n2_red2        = [0,0,0];                  % red2 nonlinear refractive index (in sq cm/W)
shg_inputs.mix_2d_sp_n2_blue        = [0,0,0];                  % blue nonlinear refractive index (in sq cm/W)
shg_inputs.mix_2d_sp_beta_red1      = [0,0,0];                  % red1 two-photon absorption coefficient (in cm/W)
shg_inputs.mix_2d_sp_beta_red2      = [0,0,0];                  % red2 two-photon absorption coefficient (in cm/W)
shg_inputs.mix_2d_sp_beta_blue      = [0,0,0];                  % blue two-photon absorption coefficient (in cm/W)
shg_inputs.mix_2d_sp_pulseenergy    = [2e-6,2e-6,0];            % input pulse energies (in J)
shg_inputs.mix_2d_sp_pulse_durations = [5,5,5;1,1,1].';         % input pulse durations  (in fwhm ps), with optional (super) gaussian coefficient (one wave per row)
shg_inputs.mix_2d_sp_pulse_delays   = [-2.50,-2.50];            % pulse delay relative to the blue pulse (in ps)
shg_inputs.mix_2d_sp_pulse_chirps   = zeros(3,3);               % Linear chirp of input pulse (in THz/ps) (three orders of chirp, in THz/ps, THz/ps^2, and THz/ps^3; one wave for each row).
shg_inputs.mix_2d_sp_beam_diameters = ones(3,2)*0.850;          % input beam diameter (in fwhm mm), specified in walkoff direction, and perpendicular to walkoff direction (one wave per row)
shg_inputs.mix_2d_sp_supergaussian_coeff = [1,1,1];             % spatial super gaussian coefficients
shg_inputs.mix_2d_sp_wo_angles      = [0,0,0];                  % spatial walkoff angles (in mrad)
shg_inputs.mix_2d_sp_offset_wodir   = [0, 0];                   % spatial offset relative to blue beam center in walk off direction (in mm)
shg_inputs.mix_2d_sp_rad_curv       = ones(3,2)*1e9;            % input beam radii of curvature (in mm in air) for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction]
shg_inputs.mix_2d_sp_nt             = 64;                       % number of time points to consider
shg_inputs.mix_2d_sp_nxny           = [128,128];                  % number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
shg_inputs.mix_2d_sp_crystal_length = 7.5;                      % length of nonlinear crystal (in mm)
shg_inputs.mix_2d_sp_lx_ly          = [1,1]*10;                 % size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
shg_inputs.mix_2d_sp_deff           = 2;                        % crystal nonlinear coefficient (in pm/V)
shg_inputs.mix_2d_sp_deltak         = 0;                        % phase mismatch (in rad/mm)
shg_inputs.mix_2d_sp_nz             = 64;                       % number of (z) integration steps through the crystal

shg_thg_distance                    = 0.035; % propagate 3.5 cm in freespace between crystals

% Construct the data structure of inputs for the modified 2D-mix-SP function for THG
% modeling. In the loop where we consider each element in the red1_red2_energy_vec, we'll
% load the output fields of the SHG run and set them as the input fields for THG. We'll
% set the THG input energies to use the output of the SHG stage. As stated above, because
% the input fields take precedence over the rest of the variables in the input data
% structure, some of the values we set here will be ignored. These ignored variables
% include:
%   mix_2d_sp_pulse_durations, mix_2d_sp_pulse_delays,
%   mix_2d_sp_beam_diameters, mix_2d_sp_supergaussian_coeff,
%   mix_2d_sp_offset_wodir, mix_2d_sp_rad_curv, mix_2d_sp_nxny,
%   mix_2d_sp_lx_ly, and mix_2d_sp_nt.
% ------------------------------------------------------------
clear thg_inputs
thg_inputs.mix_2d_sp_wavelengths    = [1064,532,354.7];         % wavelengths (in nm)
thg_inputs.mix_2d_sp_ref_inds       = [2.142,2.184,2.190];      % refractive indices
thg_inputs.mix_2d_sp_gvi            = [2.1875,2.2811,2.282];    % group velocity indices
thg_inputs.mix_2d_sp_gdd            = [88,384,400];              % group delay dispersion (in fs^2/mm)
thg_inputs.mix_2d_sp_phase          = [0,0,0];                  % phases at input (in rad)
thg_inputs.mix_2d_sp_input_refl     = [0,0,0];                  % input face power reflectivity (0-1)
thg_inputs.mix_2d_sp_output_refl    = [0,0,0];                  % output face power reflectivity (0-1)
thg_inputs.mix_2d_sp_crystal_losses = [0,0,0];                  % linear power absorption coefficient (in 1/mm)
thg_inputs.mix_2d_sp_n2_red1        = [0,0,0];                  % red1 nonlinear refractive index (in sq cm/W)
thg_inputs.mix_2d_sp_n2_red2        = [0,0,0];                  % red2 nonlinear refractive index (in sq cm/W)
thg_inputs.mix_2d_sp_n2_blue        = [0,0,0];                  % blue nonlinear refractive index (in sq cm/W)
thg_inputs.mix_2d_sp_beta_red1      = [0,0,0];                  % red1 two-photon absorption coefficient (in cm/W)
thg_inputs.mix_2d_sp_beta_red2      = [0,0,0];                  % red2 two-photon absorption coefficient (in cm/W)
thg_inputs.mix_2d_sp_beta_blue      = [0,0,0];                  % blue two-photon absorption coefficient (in cm/W)
thg_inputs.mix_2d_sp_pulseenergy    = [2e-6,2e-6,0];            % input pulse energies (in J)
thg_inputs.mix_2d_sp_pulse_durations= [5,5,5;1,1,1].';          % input pulse durations  (in fwhm ps), with optional (super) gaussian coefficient (one wave per row)
thg_inputs.mix_2d_sp_pulse_delays   = [-2.50,-2.50];            % pulse delay relative to the blue pulse (in ps)
thg_inputs.mix_2d_sp_pulse_chirps   = [0,0,0;0,0,0;0,0,0];      % Linear chirp of input pulse (in THz/ps) (three orders of chirp, in THz/ps, THz/ps^2, and THz/ps^3; one wave for each row).
thg_inputs.mix_2d_sp_beam_diameters = shg_inputs.mix_2d_sp_beam_diameters;           % input beam diameter (in fwhm mm), specified in walkoff direction, and perpendicular to walkoff direction (one wave per row)
thg_inputs.mix_2d_sp_supergaussian_coeff = [1,1,1];                     % spatial super gaussian coefficients
thg_inputs.mix_2d_sp_wo_angles      = [0,0,0];                          % spatial walkoff angles (in mrad)
thg_inputs.mix_2d_sp_offset_wodir   = [0,0];                            % spatial offset relative to blue beam center in walk off direction (in mm)
thg_inputs.mix_2d_sp_rad_curv       = [1e9,1e9,1e9;1e9,1e9,1e9].';      % input beam radii of curvature (in mm in air) for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction]
thg_inputs.mix_2d_sp_nt             = shg_inputs.mix_2d_sp_nt;          % number of time points to consider
thg_inputs.mix_2d_sp_nxny           = shg_inputs.mix_2d_sp_nxny;        % number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
thg_inputs.mix_2d_sp_crystal_length = 7.5;                              % length of nonlinear crystal (in mm)
thg_inputs.mix_2d_sp_lx_ly          = shg_inputs.mix_2d_sp_lx_ly;       % size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
thg_inputs.mix_2d_sp_deff           = 2;                                % crystal nonlinear coefficient (in pm/V)
thg_inputs.mix_2d_sp_deltak         = 0;                                % phase mismatch (in rad/mm)
thg_inputs.mix_2d_sp_nz             = 64;                               % number of (z) integration steps through the crystal
thg_inputs.input_fields             = [];                               % make a new, empty field named input_fields; in each iteration looping through the input pulse energies we'll place an object of class mix_2d_sp_input_fields_class

% Pre-allocate some cell arrays we'll stick the output SHG waves and THG waves in
shg_output_field_red1   = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the 3d array of red1 electric field for shg output
shg_output_field_red2   = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the 3d array of red2 electric field for shg output
shg_output_field_blue   = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the 3d array of blue electric field for shg output

thg_output_field_red1   = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the 3d array of red1 electric field for thg output
thg_output_field_red2   = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the 3d array of red2 electric field for thg output
thg_output_field_blue   = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the 3d array of blue electric field for thg output

shg_red1_energy_out     = zeros(size(red1_red2_energy_vec));    % pre-allocate an array of zeros to store 1064 nm output energies for the shg stages
shg_red2_energy_out     = zeros(size(red1_red2_energy_vec));    % pre-allocate an array of zeros to store 1064 nm output energies for the shg stages
shg_blue_energy_out     = zeros(size(red1_red2_energy_vec));    % pre-allocate an array of zeros to store  532 nm output energies for the shg stages

thg_red1_energy_out     = zeros(size(red1_red2_energy_vec));    % pre-allocate an array of zeros to store  1064 nm output energies for the thg stages
thg_red2_energy_out     = zeros(size(red1_red2_energy_vec));    % pre-allocate an array of zeros to store   532 nm output energies for the thg stages
thg_blue_energy_out     = zeros(size(red1_red2_energy_vec));    % pre-allocate an array of zeros to store 354.7 nm output energies for the thg stages

shg_red1_power_out      = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the output power at each time point for red1 in the shg stage
shg_red2_power_out      = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the output power at each time point for red2 in the shg stage
shg_blue_power_out      = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the output power at each time point for blue in the shg stage

thg_red1_power_out      = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the output power at each time point for red1 in the thg stage
thg_red2_power_out      = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the output power at each time point for red1 in the thg stage
thg_blue_power_out      = cell(size(red1_red2_energy_vec));     % pre-allocate a cell array to store the output power at each time point for red1 in the thg stage

clear shg_problem thg_problem % shg_problem and thg_problem will be vectors of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting

%% Loop through each value in the input energy vector red1_red2_energy_vec
for K = 1:length(red1_red2_energy_vec) % loop over each element in the vector of input pulse energies
    % I had some trouble keeping track of which part of the model was
    % running, so I print some strings to the SNLO console reflecting
    % the current step and its results
    fprintf(1,'Running SHG for %.6G uJ (red) (run %i of %i)\n', ...
        2*red1_red2_energy_vec(K)*1e6, K, length(red1_red2_energy_vec)); % (Remember, red1_red2_energy_vec energy gets stuck in each of the red waves

    % Set up the SHG part of the problem for this iteration by copying a set of the inputs specified earlier
    shg_problem(K) = shg_inputs;
    shg_problem(K).mix_2d_sp_pulseenergy(1) = red1_red2_energy_vec(K); % But be sure to set the red pulse energies to the appropriate value for this iteration
    shg_problem(K).mix_2d_sp_pulseenergy(2) = red1_red2_energy_vec(K);

    % call the modified SNLO 2D-mix-SP function with the problem set, and assign
    % the returned cell array of function handles which are local to that file
    % to 'fcn_handles'
    shg_fcn_handles   = snlo_2d_mix_sp_func(shg_problem(K)); % Call the 2D-mix-SP function file; will be returned a cell array of the file's function handles, we'll use them to simualted pressing the 'Accept' and 'Run' buttons
    shg_run_handle    = shg_fcn_handles{1};     % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
    shg_accept_handle = shg_fcn_handles{2};     % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
    shg_close_handle  = shg_fcn_handles{end};   % the function handle of the 'Close figure' callback is last returned
    shg_accept_handle();                        % call the 'accept' button's callback - equivalent to clicking the accept button
    shg_run_handle();                           % call the 'Run' button callback - equivalent to clicking the run button

    % After running, 2D-mix-SP saves to disk in the file mix_2d_sp_output.mat,
    % containing the variables field_red1_xyt_crystaloutput, field_red2_xyt_crystaloutput, field_blue_xyt_crystaloutput,
    % tgrid, xgrid, ygrid, xmat, ymat, energy_vs_z, power, and fluences
    dat = matfile('mix_2d_sp_output.mat');
    shg_output_powers   = dat.power;                      % power is an nt x 4 array, with the 1st column the time, the 2nd column red1 power, 3rd red2 power, 4th blue power
    field_red1_xyt      = dat.field_red1_xyt_crystaloutput; % Keep a copy of the red1 1064 nm output electric fields (3d array)
    field_red2_xyt      = dat.field_red2_xyt_crystaloutput; % Keep a copy of the red2 1064 nm output electric fields (3d array)
    field_blue_xyt      = dat.field_blue_xyt_crystaloutput; % Keep a copy of the blue 532  nm output electric fields (3d array)

    % Propagate fields from SHG crystal output face to THG crystal input face
    % using the propagate_freespace_class. This class definition file is
    % propagate_freespace_class.m, and contains 4 class properties and 2 class
    % methods. It works similar to the simple_freespace_propagate function.
    fprintf(1, '  Propagating SHG to THG crystal (%g mm)\n', shg_thg_distance*1e3);
    p                   = propagate_freespace_class;
    p.input_field       = field_red1_xyt;
    p.xgridmat          = dat.xmat;
    p.ygridmat          = dat.ymat;
    p.wavelength        = thg_inputs.mix_2d_sp_wavelengths(1)*1e-9; % propagate_freespace_class requires wavelengths to be provided in meters
    
    prop_field_red1     = p.propagate(shg_thg_distance);
    
    p.input_field       = field_red2_xyt;
    p.wavelength        = thg_inputs.mix_2d_sp_wavelengths(2)*1e-9;
    prop_field_red2     = p.propagate(shg_thg_distance);
    
    p.input_field       = field_blue_xyt;
    p.wavelength        = thg_inputs.mix_2d_sp_wavelengths(3)*1e-9;
    prop_field_blue =    p.propagate(shg_thg_distance);

    xmat                = dat.xmat;                         % 
    ymat                = dat.ymat;                         % 
    tgrid               = dat.tgrid;                        % 

    shg_red1_power_out{K}    = shg_output_powers(:,2);  %
    shg_red2_power_out{K}    = shg_output_powers(:,3);  %
    shg_blue_power_out{K}    = shg_output_powers(:,4);  %
    shg_output_field_red1{K} = field_red1_xyt;          %
    shg_output_field_red2{K} = field_red2_xyt;          %
    shg_output_field_blue{K} = field_blue_xyt;          %
    
    % To get input energy values to use in THG stage, integrate the shg output power vectors in time (using the trapezoidal method,
    % with the times given as the first argument and power as second)
    shg_red1_energy_out(K)   = trapz(shg_output_powers(:,1), shg_output_powers(:,2)); % 
    shg_red2_energy_out(K)   = trapz(shg_output_powers(:,1), shg_output_powers(:,3)); % 
    shg_blue_energy_out(K)   = trapz(shg_output_powers(:,1), shg_output_powers(:,4)); % 

    % Write some output to the command window showing our results
    fprintf(1,'SHG Output Energies: %.6G uJ @ 1064nm, %.6G uJ @ 532nm (total %.6G uJ)\n',...
        (shg_red1_energy_out(K) + shg_red2_energy_out(K))*1e6, shg_blue_energy_out(K)*1e6, ...
        sum([shg_red1_energy_out(K), shg_red2_energy_out(K), shg_blue_energy_out(K)])*1e6 );

    % Make an object of class mix_2d_sp_input_fields_class:
    thg_input_fields = mix_2d_sp_input_fields_class;
    
    %  Set parameter values:
    % In the type I mixing shg stage, we split the red input energy equally
    % between red1 and red2. Now when we combine the red1 and red2 fields to
    % form the new red1 input we sum the red1 and red2 output fields and
    % multiply by 1/sqrt(2). If the first stage is not type 1 SHG, this is
    % irrelevant.
    thg_input_fields.field_red1_xyt = sqrt(0.5) * (prop_field_red1 + prop_field_red2);      % 1064 nm
    thg_input_fields.field_red2_xyt = prop_field_blue;                                      % 532 nm
    thg_input_fields.field_blue_xyt = prop_field_blue*1e-12;                                % 354.7nm starts with 0 energy, but setting its value to be just 0 can cause problems?

    thg_input_fields.xgrid = xmat;
    thg_input_fields.ygrid = ymat;
    thg_input_fields.tgrid = tgrid;

    % As a check, we calculate the THG input pulse energies by
    % integrating the THG input fields in x,y, and t
    dx = xmat(2,2) - xmat(1,1);     % x-step
    dy = ymat(2,2) - ymat(1,1);     % y-step
    dt = tgrid(2)  - tgrid(1);      % time-step

    % make a copy of the 'master' thg input set
    thg_problem(K) = thg_inputs;                    
    
    thg_problem(K).input_fields = thg_input_fields;

    % calculate the thg input energy in each wave; this must be placed in the field named mix_2d_sp_pulseenergy
    thg_r1_energy_in = dx*dy*dt*sum(sum(sum( 0.5 * c * epsilon_0 .* abs(thg_problem(K).input_fields.field_red1_xyt).^2)));
    thg_r2_energy_in = dx*dy*dt*sum(sum(sum( 0.5 * c * epsilon_0 .* abs(thg_problem(K).input_fields.field_red2_xyt).^2)));
    thg_bl_energy_in = dx*dy*dt*sum(sum(sum( 0.5 * c * epsilon_0 .* abs(thg_problem(K).input_fields.field_blue_xyt).^2)));

    thg_problem(K).mix_2d_sp_pulseenergy = [thg_r1_energy_in, thg_r2_energy_in, thg_bl_energy_in];

    % Write some output to the command window showing our progress:
    fprintf(1,'Running THG (for SHG %.6G uJ red input) - run %i of %i\n',2*red1_red2_energy_vec(K)*1e6,K,length(red1_red2_energy_vec));
    fprintf(1,'THG Input Energies: %.6G uJ @ 1064nm, %.6G uJ @ 532nm, %.6G uJ @ 354.7nm (total %.6G uJ)\n',thg_r1_energy_in*1e6, ...
        thg_r2_energy_in*1e6, thg_bl_energy_in*1e6, (thg_r1_energy_in + thg_r2_energy_in + thg_bl_energy_in)*1e6);

    thg_fcn_handles     = snlo_2d_mix_sp_func(thg_problem(K));  % call the SNLO 2D-mix-SP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
    thg_run_handle      = thg_fcn_handles{1};                   % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
    thg_accept_handle   = thg_fcn_handles{2};                   % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
    thg_close_handle    = thg_fcn_handles{end};                 % the function handle of the 'Close figure' callback is last returned
    thg_accept_handle();    % call the 'accept' button's callback - equivalent to clicking the accept button
    thg_run_handle();       % call the 'Run' button callback - equivalent to clicking the run button

    % After running this model load the results from file
    dat = matfile('mix_2d_sp_output.mat');
    thg_output_powers = dat.power;      % power is an nt x 4 array, with the 1 column the time, the 2nd column red1 power, 3rd red2 power, 4th blue power
    field_red1_xyt  = dat.field_red1_xyt_crystaloutput;
    field_red2_xyt  = dat.field_red2_xyt_crystaloutput;
    field_blue_xyt  = dat.field_blue_xyt_crystaloutput;
    xmat            = dat.xmat;
    ymat            = dat.ymat;
    tgrid           = dat.tgrid;

    thg_red1_power_out{K}    = thg_output_powers(:,2);  % 2nd column red1
    thg_red2_power_out{K}    = thg_output_powers(:,3);  % 3rd column red2
    thg_blue_power_out{K}    = thg_output_powers(:,4);  % 4th column blue
    thg_output_field_red1{K} = field_red1_xyt;          % keep a copy of thg output fields (red1 wave)
    thg_output_field_red2{K} = field_red2_xyt;          % keep a copy of thg output fields (red2 wave) 
    thg_output_field_blue{K} = field_blue_xyt;          % keep a copy of thg output fields (blue wave)
    thg_red1_energy_out(K)   = trapz(thg_output_powers(:,1), thg_output_powers(:,2));
    thg_red2_energy_out(K)   = trapz(thg_output_powers(:,1), thg_output_powers(:,3));
    thg_blue_energy_out(K)   = trapz(thg_output_powers(:,1), thg_output_powers(:,4));

    % Write some output to the command window indicating our results so far:
    fprintf(1,'THG Output Energies: %.6G uJ @ 1064nm, %.6G uJ @ 532nm, %.6G uJ @ 354.7nm (total %.6G uJ)\n',...
        thg_red1_energy_out(K)*1e6, thg_red2_energy_out(K)*1e6, thg_blue_energy_out(K)*1e6, ...
        sum([thg_red1_energy_out(K),thg_red2_energy_out(K),thg_blue_energy_out(K)])*1e6);
    
end

%% Plot results
fh = figure; % Make a new figure

% With 2 subplot axes:
ax1 = axes('parent', fh);
ax2 = axes('parent', fh);
subplot(2,1,1,ax1);
subplot(2,1,2,ax2);

% Add a plot of shg output energy to top plot;
ph = plot(ax1, 2*red1_red2_energy_vec*1e6, (shg_red1_energy_out+shg_red2_energy_out)*1e6, 'r-', ...
    2*red1_red2_energy_vec*1e6, shg_blue_energy_out*1e6,'g-');
xlabel(ax1, 'Input 1064 energy [{\mu}J]');
ylabel(ax1, 'SHG output energy [{\mu}J]');
legend(ax1, '1064 nm', '532 nm', 'location', 'nw');

% Add a plot of thg output energy to bottom plot:
ph2 = plot(ax2, 2*red1_red2_energy_vec*1e6, thg_red1_energy_out*1e6, 'r-', ...
    2*red1_red2_energy_vec*1e6, thg_red2_energy_out*1e6, 'g-', ...
    2*red1_red2_energy_vec*1e6, thg_blue_energy_out*1e6,'b-');
xlabel(ax2, 'Input 1064 energy [{\mu}J]');
ylabel(ax2, 'THG output energy [{\mu}J]');
legend(ax2, '1064 nm', '532 nm', '354.7 nm', 'location', 'nw');