% This is an example for the modified version of 2D mix LP which allows a user to specify
% arbitrary field profiles (in x,y and t). Calling the function snlo_2d_mix_lp_func with
% arbitrary field profiles will require the user to provide three electric fields (red1,
% red2, and blue) 2-dimensional arrays (x, y) of size nx x ny, as well as the spatial and
% grids which each of these arrays is defined on. When called with user-specified input
% fields, several of the input edit boxes in the main form of 2D mix LP will be ignored,
% and several of these should be filled with values appropriate to the user-specified
% fields (such as spatial and temporal grid extents and number of points). Included in
% this example script is a section generating the input electric fields for the three
% waves, a section constructing a data structure defining all of the function inputs, a
% section looping through set of input pulse energies which calls the mlSNLO function,
% loads & stores the results.
%
% We'll be modeling cw second-harmonic generation, doubling 1000 nm to 500 nm. We model
% SHG by including two identical red waves (total red energy split evenly between them).
% The 'Focus' function was used to calculate the radius of curvature and input beam diameter
% (FWHM) for the red wave that would create a 35 um waist 1/2 through the crystal

% Call the mlSNLO function for each input energy, run it, load the output from disk for
% processing.

c = 3e8; % speed of light in m/s
epsilon_0 = 8.85E-12; % vacuum permittivity in F/m

calculate_curves = true;
if calculate_curves
    %% Define some parameters to generate grids and user-specified input electric fields
    
    % specify parameters for generating the user-specified input fields and grids
    lambda = [1000, 1000, 500]*1e-9;    % red1, red2, and blue wavelengths (in meters)
    red1red2_power_list = linspace(10,250,15);       % Powers of each red1&red2 (in watts) - total power input is twice this value
    blue_power = 1e-40;                 % Blue wave input power (in watts); it's zero, really, but needs a small, finite value to avoid dividing by 0 in the mlSNLO function
    nx = 128;                           % number of spatial grid points in the x (walkoff) direction
    ny = 128;                           % number of spatial grid points in the y (non-walkoff) direction
    durations = [0,0,0];                % input pulse durations (fwhm ns)
    lx = 0.45e-3;                       % extent of spatial grid in the x (walkoff) direction (mm)
    ly = 0.45e-3;                       % extent of spatial grid in the y (non-walkoff) direction (mm)
    spatial_supergaussian_coeffs    = [1,1,1];
    input_field_diameters = 105.2e-6*[1, 1; 1, 1; 1, 1];             % input beam diameters (in fwhm mm) [red1 x-dir, red1 y-dir; red2 x-dir, red2 y-dir; blue x-dir, blue y-dir;] 
    position  = [0, 0; 0, 0; 0, 0]*1e-3;                            % offset of center of beam (in walkoff direction & perpendicular to it) for red1, red2, and blue fields (mm)
    beam_rocs = 8.85e-3 * [1,1;1,1;1,1];    % radius of curvature for (x&y planes) for red1, red2, and blue waves (mm)
    
    % Spatial grid:
    % first calculate x step size and y step size
    dx = lx / (nx - 1);
    dy = ly / (ny - 1);
    % now make row vector of spatial grid in x and y directions
    xgrid_vec = (-0.5*lx) + (0:(nx-1))*dx;
    ygrid_vec = (-0.5*ly) + (0:(ny-1))*dy;
    % make 2d-arrays of spatial grids; ygrid_mat is constant in the first dimension and
    % xgrid_mat is constant in the second dimensions
    [ygrid_mat,xgrid_mat] = meshgrid(ygrid_vec, xgrid_vec);
    
    % first generate the shapes of the input irradiance profiles, then we'll normalize so
    % they include the appropriate amount of power, take the square root and multiply by
    % arrays of phase (sizes nx x ny) to add input beam radius of curvature and/or tilt
    
    % Define some nx x 3 and ny x 3 arrays with Gaussian profiles
    exponent_x = zeros(nx,3);
    exponent_x(:,1) = ( (2./input_field_diameters(1,1)).*(xgrid_vec-position(1,1)) ).^2; % red 1 spatial profile in x (with proper diameter and center)
    exponent_x(:,2) = ( (2./input_field_diameters(2,1)).*(xgrid_vec-position(2,1)) ).^2; % red 2 spatial profile in x (with proper diameter and center)
    exponent_x(:,3) = ( (2./input_field_diameters(3,1)).*(xgrid_vec-position(3,1)) ).^2; % blue spatial profile in x (with proper diameter and center)
    exponent_y = zeros(ny,3);
    exponent_y(:,1) = ( (2./input_field_diameters(1,2)).*(ygrid_vec-position(1,2)) ).^2; % red 1 spatial profile in x (with proper diameter and center)
    exponent_y(:,2) = ( (2./input_field_diameters(2,2)).*(ygrid_vec-position(2,2)) ).^2; % red 2 spatial profile in x (with proper diameter and center)
    exponent_y(:,3) = ( (2./input_field_diameters(3,2)).*(ygrid_vec-position(3,2)) ).^2; % blue spatial profile in x (with proper diameter and center)

    % Add the nx x 3 and ny x 3 arrays into nx x ny x 3 array, stick em in exp
    irrad_red1_xy = exp( -log(2)*(repmat(exponent_x(:,1),[1,ny])+repmat(exponent_y(:,1).',[nx,1])).^(spatial_supergaussian_coeffs(1)) ); % red1 spatial 2d profile
    irrad_red2_xy = exp( -log(2)*(repmat(exponent_x(:,2),[1,ny])+repmat(exponent_y(:,2).',[nx,1])).^(spatial_supergaussian_coeffs(2)) ); % red2 spatial 2d profile
    irrad_blue_xy = exp( -log(2)*(repmat(exponent_x(:,3),[1,ny])+repmat(exponent_y(:,3).',[nx,1])).^(spatial_supergaussian_coeffs(3)) ); % blue spatial 2d profile

    % normalize so each irradiance profile has a total of 1 watt
    irrad_red1_xy_norm = irrad_red1_xy ./ (dx*dy*sum(sum(irrad_red1_xy)));
    irrad_red2_xy_norm = irrad_red2_xy ./ (dx*dy*sum(sum(irrad_red2_xy)));
    irrad_blue_xy_norm = irrad_blue_xy ./ (dx*dy*sum(sum(irrad_blue_xy)));
    
    % Convert irradiance into field
    field_red1_xy_norm = sqrt(irrad_red1_xy_norm / (0.5*epsilon_0*c));
    field_red2_xy_norm = sqrt(irrad_red2_xy_norm / (0.5*epsilon_0*c));
    field_blue_xy_norm = sqrt(irrad_blue_xy_norm / (0.5*epsilon_0*c));
    
    % Generate nx x ny arrays of phases for the E fields - they have radius of curvature and maybe tilt
    curvature_phase_red1 = exp( 1i*(-pi/lambda(1))*( repmat((ygrid_vec.^2),[nx,1])./beam_rocs(1,2) + repmat(((xgrid_vec-position(1)).^2).',[1,ny])./beam_rocs(1,1) ) );
    curvature_phase_red2 = exp( 1i*(-pi/lambda(2))*( repmat((ygrid_vec.^2),[nx,1])./beam_rocs(2,2) + repmat(((xgrid_vec-position(2)).^2).',[1,ny])./beam_rocs(2,1) ) );
    curvature_phase_blue = exp( 1i*(-pi/lambda(3))*( repmat((ygrid_vec.^2),[nx,1])./beam_rocs(3,2) + repmat(((xgrid_vec-position(3)).^2).',[1,ny])./beam_rocs(3,1) ) );

    % multiply fields by phase
    field_red1_xy_norm = field_red1_xy_norm.*curvature_phase_red1;
    field_red2_xy_norm = field_red2_xy_norm.*curvature_phase_red2;
    field_blue_xy_norm = field_blue_xy_norm.*curvature_phase_blue;

    % Now each nx x ny field should have 1 W of power and includes any phase from curvature
    
    %% Input parameter data structure
    % Define the full set of input parameters for 2D-MIX-LP; these parameters correspond to the input
    % boxes visible when you run the 2D-MIX-LP function. For most parameters, a vector with three
    % elements is defined and corresponds to the values for the red1, red2, and blue waves. Calling the
    % 2D-MIX-LP function with a set of input parameters requires specifying a single input argument to
    % the function snlo_2d_mix_lp_func, where the input parameter is a data structure with the fieldnames
    % specified below, and values in the form specified below. 2D-mix-LP will use default values for 
    % all fieldnames omitted; messages will be written in the command window indicating any fieldnames missing
    % and the default value used.
    clear inputs
    inputs.mix_2d_lp_wavelengths    = lambda*1e9;      % wavelengths (in nm)
    inputs.mix_2d_lp_ref_inds       = [1.655,1.655,1.655];  % refractive indices
    inputs.mix_2d_lp_phase          = [0,0,0];              % input phases (in radians)
    inputs.mix_2d_lp_input_refl     = [0,0,0];              % crystal input face reflectivity
    inputs.mix_2d_lp_output_refl    = [0,0,0];              % crystal output face reflectivity
    inputs.mix_2d_lp_crystal_losses = [0,0,0];              % linear absorption (in 1/mm)
    inputs.mix_2d_lp_pulseenergy    = [red1red2_power_list(1),red1red2_power_list(1),blue_power];        % input powers (in watts) or pulse energies (in joules)
    inputs.mix_2d_lp_pulse_durations = durations;           % input pulse durations (in fwhm ns)
    inputs.mix_2d_lp_beam_diameters = input_field_diameters*1e3;% input beam diameters (in fwhm mm) [red1 in walkoff dir, red2 in walkoff dir, blue in walkoff dir;red1 perp to walkoff, red2 perp to walkoff, blue perp to walkoff]
    inputs.mix_2d_lp_pulse_delays   = [0,0,0];              % delay (in ns)
    inputs.mix_2d_lp_supergaussian_coeff = [1,1,1];         % supergaussian coefficient for spatial input beam profile
    inputs.mix_2d_lp_n2_red1        = [0,0,0];              % red1 nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_n2_red2        = [0,0,0];              % red2 nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_n2_blue        = [0,0,0];              % blue nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_beta_red1      = [0,0,0];              % red1 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_beta_red2      = [0,0,0];              % red2 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_beta_blue      = [0,0,0];              % blue 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_wo_angles      = [0,0,0];              % spatial walkoff angles (in milliradians)
    inputs.mix_2d_lp_offset_wodir   = [0,0,0];              % offset in walk off direction (in mm)
    inputs.mix_2d_lp_rad_curv       = beam_rocs*1e3;        % input beam radii of curvature for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction] (in mm in air)
    inputs.mix_2d_lp_nxny           = [nx, ny];             % # pts in x, y (in walk off direction and perpendicular to it).
    inputs.mix_2d_lp_nz             = 100;                  % number of z integration steps
    inputs.mix_2d_lp_crystal_length = [26,1];               % crystal length (in mm); optionally, the second parameter is number of walkoff compensating sections
    inputs.mix_2d_lp_lx_ly          = [lx, ly]*1e3;         % Extent of spatial grid in (in mm) in walk off direction and perpendicular to it
    inputs.mix_2d_lp_deff           = 2;                    % crystal nonlinear coefficient (in pm/V)
    inputs.mix_2d_lp_deltak         = 0;                    % phase mismatch (in rad/mm)
    inputs.mix_2d_lp_dist_to_image  = 0;                    % Distance from the crystal exit face to image plane  (in mm) of the output beam profile displayed during the run and by the post-run fluence and movie buttons.
    inputs.mix_2d_lp_nt             = 32;                   % number of points in time (used for pulsed operation only)
    inputs.mix_2d_lp_save_movie     = false;                 % logical; if true, resulting electric fields at crystal output face are saved to disk; this can take tens of seconds in some cases with many points in x,y (and t, for pulsed operation)
    inputs.mix_2d_lp_auto_analyze   = false;                 % logical; if true, do beam analysis on output fields: M-squared, tilt, second moment, radius of curvature
    %% call 2d mix lp, run model, load output
    output_red1 = cell(size(red1red2_power_list)); % pre-allocate a cell array to store the contents of the red1 output files in
    output_red2 = cell(size(red1red2_power_list)); % pre-allocate a cell array to store the contents of the red2 output files in
    output_blue = cell(size(red1red2_power_list)); % pre-allocate a cell array to store the contents of the blue output files in
    red1_output_powers = zeros(size(red1red2_power_list)); % pre-allocate a vector for red1 output powers, same size as input power list
    red2_output_powers = zeros(size(red1red2_power_list)); % pre-allocate a vector for red2 output powers, same size as input power list
    blue_output_powers = zeros(size(red1red2_power_list)); % pre-allocate a vector for blue output powers, same size as input power list
    efficiencies = zeros(size(red1red2_power_list));
    clear problem % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    for K = 1:length(red1red2_power_list) % loop over each element in the vector of input pulse energies
        powers = [red1red2_power_list(K),red1red2_power_list(K),blue_power];
        problem{K} = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
        problem{K}.mix_2d_lp_pulseenergy = powers;  % even though this is cw, the snlo_2d_mix_lp_func expects a '_pulseenergy' data field
        input_fields.field_red1_xyt = sqrt(powers(1))*field_red1_xy_norm;
        input_fields.field_red2_xyt = sqrt(powers(2))*field_red2_xy_norm;
        input_fields.field_blue_xyt = sqrt(powers(3))*field_blue_xy_norm;
        input_fields.xgrid = xgrid_vec;
        input_fields.ygrid = ygrid_vec;
        problem{K}.input_fields = input_fields;
        
        fcn_handles = snlo_2d_mix_lp_func(problem{K}); % call the SNLO 2D-MIX-LP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        accept_handle = fcn_handles{2}; % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        close_handle = fcn_handles{end}; % the function handle of the 'Close figure' callback is last returned
        accept_handle();    % call the 'accept' button's callback - equivalent to clicking the accept button
        run_handle();       % call the 'Run' button callback - equivalent to clicking the run button
        output_red1{K} = importdata('mix_2d_lp_red1_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        output_red2{K} = importdata('mix_2d_lp_red2_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        output_blue{K} = importdata('mix_2d_lp_blue_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        red1_output_powers(K) = output_red1{K}.data(:,2);   % pick output power from output file
        red2_output_powers(K) = output_red2{K}.data(:,2);   % pick output power from output file
        blue_output_powers(K) = output_blue{K}.data(:,2);   % pick output power from output file
        efficiencies(K) = blue_output_powers(K)./(red2_output_powers(K)+red1_output_powers(K));
        
        % display results of this iteration in command window
        fprintf(1,'Run %i of %i: input red power = %.4g W; output red power = %.4g W; output blue power = %.4g W; efficiency = %.4f\n', ...
            K, length(red1red2_power_list), 2*red1red2_power_list(K), red1_output_powers(K)+red2_output_powers(K), ...
            blue_output_powers(K), efficiencies(K));
%         close_handle(); % call the 'Close figure' callback
    end
end


%% plot processed output
fh = figure;  % make new outut figure
ax = axes('Parent',fh); % make new axes in that output figure

% plot the efficiency
ph = plot(ax,2*red1red2_power_list,efficiencies);
xlabel(ax,'Input red power [W]');       % label x-axis 
ylabel(ax,'Conversion efficiency');     % label y-axis 