% Example of 2D-mix-SP custom input fields using Laguerre-Gaussian modes
%  We'll compare the efficiency of a frequency doubler with a purely Gaussian red field
%  input with one that has 20% of its power in LG_10 and the rest in LG_00.
% As with other 2D-mix-SP custom input field example scripts, we'll start by defining a
% set of model inputs, things like wavelength, refractive indices, and so forth. We'll
% then generate a shape to pass to the model as input field profiles (in x, y, and t).
% 

%  
  

%% Make inputs data structure specifying all of the input parameter values
clear inputs
inputs.mix_2d_sp_wavelengths = [1550,1550,775];     % wavelengths (in nm); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_ref_inds = [2.142,2.142,2.154];    % refractive indices; 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_gvi = [2.1875,2.1875,2.25];        % group velocity indices; 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_gdd = [1.04e2,1.04e2,3.99e2];      % group delay dispersion (in fs^2/mm); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_phase = [0,0,0];                   % phases at input (in rad); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_input_refl = [0,0,0];              % input face power reflectivity (0-1); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_output_refl = [0,0,0];             % output face power reflectivity (0-1); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_crystal_losses = [0,0,0];          % linear power absorption coefficient (in 1/mm); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_n2_red1 = [0,0,0];                 % red1 nonlinear refractive index (in sq cm/W); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_n2_red2 = [0,0,0];                 % red2 nonlinear refractive index (in sq cm/W); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_n2_blue = [0,0,0];                 % blue nonlinear refractive index (in sq cm/W); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_beta_red1 = [0,0,0];               % red1 two-photon absorption coefficient (in cm/W); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_beta_red2 = [0,0,0];               % red2 two-photon absorption coefficient (in cm/W); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_beta_blue = [0,0,0];               % blue two-photon absorption coefficient (in cm/W); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_pulseenergy = [5e-8,5e-8,0];       % input pulse energies (in J); 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_pulse_durations = [1,1,1;1,1,1].'; % ** input pulse durations  (in fwhm ps), with optional (super) gaussian coefficient; 3x2 array (with temporal supergaussian coefficients), or 3x1 array (without)
inputs.mix_2d_sp_pulse_delays = [0,0];              % ** pulse delay relative to the blue pulse (in ps); 2 element vector (for red1, red2)
inputs.mix_2d_sp_pulse_chirps = [0,0,0];            % ** Linear chirp of input pulse, for linear (THz/ps) and optionally quadratic (THz/ps^2) and cubic (THz/ps^3) chirps; requires a 3x1 array (linear chirps only), 3x2 array (linear and quadratic chirps), or 3x3 array (linear, quadratic, and cubic chirps); 1st dimension is (red1, red2, blue) order
inputs.mix_2d_sp_beam_diameters = [0.5,0.5,0.5;0.5,0.5,0.5].'; % ** input beam diameter (in fwhm mm), specified in walkoff direction, and perpendicular to walkoff direction
inputs.mix_2d_sp_supergaussian_coeff = [1,1,1];     % ** spatial super gaussian coefficients; 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_wo_angles = [0,0,0];               % spatial walkoff angles (in mrad);  3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_offset_wodir = [0,0];              % ** spatial offset relative to blue beam center in walk off direction (in mm);  2 element vector (for red1, red2)
inputs.mix_2d_sp_rad_curv = [1e7,1e7,1e7;1e7,1e7,1e7].';  % ** input beam radii of curvature (in mm in air) for [red1,red2,blue] and optionally 2x3 array for different curvatures in the walkoff direction and perpendicular to it: [in the walkoff direction;parallel to the walkoff direction]. 
inputs.mix_2d_sp_nt = 100;                           % ** number of time points to model
inputs.mix_2d_sp_nxny = [64,64];                    % ** number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
inputs.mix_2d_sp_crystal_length = 5;                % length of nonlinear crystal (in mm)
inputs.mix_2d_sp_lx_ly = [2,2];                     % ** size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
inputs.mix_2d_sp_deff = 30;                         % crystal nonlinear coefficient (in pm/V)
inputs.mix_2d_sp_deltak = 0;                        % phase mismatch (in rad/mm)
inputs.mix_2d_sp_nz = 30;                           % number of (z) integration steps through the crystal
% **: snlo_2d_mix_sp_func will ignore the values provided in these fields. They are
% defined here so we can use them in this script to generate electric field profiles.

%% Make a copy of some values specified above, and convert to mks units, for use in generating grids to feed to Laguerre-Gaussian calculator
wavelengths = inputs.mix_2d_sp_wavelengths*1e-9;        % (meters) wavelength of red1, red2, blue waves
nx = inputs.mix_2d_sp_nxny(1) ;                         % # pts in x
ny = inputs.mix_2d_sp_nxny(2) ;                         % # pts in y
nt = inputs.mix_2d_sp_nt;                               % # pts in time
Lx = inputs.mix_2d_sp_lx_ly(1)*1e-3;                    % (meters) xgrid span
Ly = inputs.mix_2d_sp_lx_ly(2)*1e-3;                    % (meters) ygrid span
diams_x = inputs.mix_2d_sp_beam_diameters(:,1)*1e-3;    % (meters) beam diameter in birefringent walkoff direction
diams_y = inputs.mix_2d_sp_beam_diameters(:,2)*1e-3;    % (meters) beam diameter in direction perpendicular to walkoff
offsets_x = [inputs.mix_2d_sp_offset_wodir*1e-3,0];     % (meters) beam offset in birefringent walkoff direction
durations = inputs.mix_2d_sp_pulse_durations(:,1)*1e-12;% (seconds) full-width half-max pulse duration
temporal_sgs = inputs.mix_2d_sp_pulse_durations(:,2);   % temporal super gaussian index (value of 1 is Gaussian, larger values are more flat-topped)
delays = [inputs.mix_2d_sp_pulse_delays*1e-12,0];       % (seconds) delay relative to blue pulse center
Rxs = inputs.mix_2d_sp_rad_curv(:,1)*1e-3;              % (meters) radius of curvature in direction of walk off
Rys = inputs.mix_2d_sp_rad_curv(:,2)*1e-3;              % (meters) radius of curvature in direction of walk off

%% Generate spatial (x,y) grids centered about 0, and temporal grid, for use in generating the input electric field profiles
xvec = linspace( -0.5, 0.5, nx)*Lx; % has nx points and span of lx, centered about 0
yvec = linspace( -0.5, 0.5, ny)*Ly; % has ny points and span of ly, centered about 0
tvec = linspace( -3, 3, nt)*max(inputs.mix_2d_sp_pulse_durations(:,1))*1e-12; % has nt points, centered about 0, spans 6 * largest pulse duration
[ygrid, xgrid] = meshgrid(yvec, xvec);

% pre-allocate array we'll fill with electric field profiles for the red1, red2, and blue
% waves (dimension order: x, y, t, and the last dimension is for which wave)
field_xyt_superposition  = zeros(nx,ny,nt,3);         % For field that is superposition of LG_00 and LG_10 modes;
field_xyt_gaussian       = zeros(nx,ny,nt,3);         % For field that is just LG_00 mode

LG_modes_00_and10 = zeros(nx, ny, 2, 3);              % Store the output of lg_calc, 3rd dimension is superposition/only LG_00 mode, 4th dimension is red1/red2/blue wave

red_ratio_lg00_over_lg10 = 4;           % 20% of input power in LG_10, the rest in LG_00
lg00_lg10_phase = 0;                    % relative phase of LG_10 from LG_00
rmat = sqrt(xgrid.^2 + ygrid.^2);
% Note that our Laguerre-Gaussian mode calculator function lg_calc requires r
% positions to be a vector, so we will reshape rmat to be a vector. We'll do the
% opposite with lg_calc output to give us LG mode amplitudes on the rmat array, which
% is nx*ny elements
rvec = reshape(rmat, 1, []);

% For each of the the red1, red2, and blue waves, calculate the LG_00 and LG_10 Laguerre-Gaussian modes of size diams_x,
% keep an array that is a superposition of them and an array that is just LG_00 mode
for wave = 1:3 % loop through 
    % Call lg_calc with input arguments pmax, radius, and w (width) - pmax is largest p
    % index of LG_p0 modes to return, radius is a vector in meters, w is a scalar width in
    % meters. For the LG_00 mode, irradiance full-width half-max is smaller than w by a factor
    % of sqrt(log(2)).
    %  We want LG_00 and LG_10, so pmax = 1. The lg_calc file optionally has 2 output
    %  arguments; the first is LG_(pmax,0) mode amplitude at positions rvec, while
    %  the second contains the mode amplitudes of each mode from LG_(00) to LG_(pmax,0).
    %  
    [~, LGpmat] = lg_calc(1, rvec, diams_x(wave)*sqrt(log(2))); % Using diams_x(wave) as Gaussian full-width half-max, calculate mode amplitudes

    % Undo the reshaping from 2d array to vector we did earlier
    LG_modes_00_and10(:,:,:,wave) = reshape(LGpmat.',[nx, ny, 2]); % array of size nx * ny * 2, 1st slice in 3rd dimension is LG_00, 2nd is LG_10
    
    % Superposition of LG_00 and LG_p0 modes, with ratio between modes set by red_ratio_lg00_over_lg10 and 
    % phase between them set by lg00_lg10_phase
    field_xy_superposition = LG_modes_00_and10(:,:,1,wave)*sqrt(red_ratio_lg00_over_lg10) + ...
        exp(1i*lg00_lg10_phase)*LG_modes_00_and10(:,:,2,wave);
    
    % Just the LG_00 mode
    field_xy_gaussian = LG_modes_00_and10(:,:,1,wave);

    % make pulse shape - gaussian or super-gaussian in time
    pulse = exp(-log(2)*(((2/durations(wave)).*(tvec-delays(wave))).^2).^(temporal_sgs(wave)));

    % Combine the xy and t components to have arrays of fields in x, y, and t for each red1/red2/blue wave
    % Replicate the xy matrices to have nt points in the 3rd dimension (for time),
    % multiply by a pulse profile vector replicated to be nx*ny size in 1st and 2nd
    % dimensions
    field_xyt_superposition(:,:,:,wave) = repmat(field_xy_superposition, [1, 1, nt]) .* repmat(reshape(sqrt(pulse),1,1,[]),[nx,ny,1]);
    field_xyt_gaussian(:,:,:,wave)      = repmat(field_xy_gaussian, [1, 1, nt])      .* repmat(reshape(sqrt(pulse),1,1,[]),[nx,ny,1]);
end


% 2D-mix-SP requires custom input field profiles to be provided as objects of class
% mix_2d_sp_input_fields_class. We'll set the input fields in the for loop below, and set
% the grid properties here before we get to the loop.
 
input_fields = mix_2d_sp_input_fields_class;    % Make such an object here:

% Set the object's properties for xgrid, ygrid, and tgrid. xgrid and ygrid are
% two-dimensional position arrays in meters and tgrid is a vector of times in seconds.
input_fields.xgrid = xgrid;
input_fields.ygrid = ygrid;
input_fields.tgrid = tvec;

% vector of energies to loop over (each red wave gets this value so total input energy is twice this value)
red1_red2_energy_input_vector = logspace(-10,-7.5,10);

% preallocate some arrays to store results in
output_energies = zeros(length(red1_red2_energy_input_vector), 3, 2); % for red1, red2, and blue energies, for both superposition of LG00 and LG10 modes and for only LG00 mode
efficiencies    = zeros(length(red1_red2_energy_input_vector), 1, 2); % for red1, red2, and blue doubling efficiencies, for both superposition of LG00 and LG10 modes and for only LG00 mode

clear problem

%% Run once using LG_00 and LG_10 input, and once using just LG_00 input
for J = 1:2

    % set the object's properties using the electric field distributions found above
    if J == 1
        input_fields.field_red1_xyt = field_xyt_superposition(:,:,:,1);
        input_fields.field_red2_xyt = field_xyt_superposition(:,:,:,2);
        input_fields.field_blue_xyt = field_xyt_superposition(:,:,:,3);
    else
        input_fields.field_red1_xyt = field_xyt_gaussian(:,:,:,1);
        input_fields.field_red2_xyt = field_xyt_gaussian(:,:,:,2);
        input_fields.field_blue_xyt = field_xyt_gaussian(:,:,:,3);
    end
    % add the input_fields object to inputs data structure
    inputs.input_fields = input_fields;

    %Loop through each input energy value stored in the vector red1_red2_energy_input_vector
    for K = 1:length(red1_red2_energy_input_vector)
        problem(K,J) = inputs; % make a copy of the 'master' inputs
        problem(K,J).mix_2d_sp_pulseenergy = [1,1,0]*red1_red2_energy_input_vector(K); % modify the problem structure to have energy from vector we're looping over

        % call 2D-mix-SP with input parameter values in problem(K), get back cell array of functions in the model we can call
        functions = snlo_2d_mix_sp_func(problem(K,J));
        
        % pick out the 'Run' button and 'Accept' button functions
        run_fcn     = functions{1};
        accept_fcn  = functions{2};

        % write to the console what we are doing
        if J == 1
            fprintf(1, 'Starting run %i of %i (LG00 and LG10 modes) input energy = %.4E J\n', K, 2*length(red1_red2_energy_input_vector), red1_red2_energy_input_vector(K)*2);
        else
            fprintf(1, 'Starting run %i of %i (LG00 mode          ) input energy = %.4E J\n', length(red1_red2_energy_input_vector)+K, 2*length(red1_red2_energy_input_vector), red1_red2_energy_input_vector(K)*2);
        end

        accept_fcn(); % simulate pressing the 'Accept' button
        run_fcn();    % simulate pressing the 'Run' button

        % After the 'Run' button function has concluded, load its output from file
        dat = load('mix_2d_sp_output.mat');

        % In dat, the variables we are interested in are called tgrid and power; tgrid is a
        % vector of times in seconds, and power is a 4xnt array which is comprised in the
        % second dimension of the times, output power for red1 wave, output power for red2
        % wave, and output power for blue wave.
        dt = dat.tgrid(2) - dat.tgrid(1);
        output_energies(K,1,J) = dt * sum( dat.power(:,2) ); % integrate output power vector in time for red1
        output_energies(K,2,J) = dt * sum( dat.power(:,3) ); % integrate output power vector in time for red1
        output_energies(K,3,J) = dt * sum( dat.power(:,4) ); % integrate output power vector in time for red1
        efficiencies(K,J)  = output_energies(K,3,J) ./ (2*red1_red2_energy_input_vector(K)); % doubling efficiency is second harmonic output energy / fundamental input energy

        % update the console on our progress and results
        if J == 1
            fprintf(1,' LG00 and LG10 : Output energies = %.4E, %.4E, %.4E J (doubling efficiency of %.2f%%)\n', output_energies(K,:,J), efficiencies(K,J)*100);
        else
            fprintf(1,' LG00 only     : Output energies = %.4E, %.4E, %.4E J (doubling efficiency of %.2f%%)\n', output_energies(K,:,J), efficiencies(K,J)*100);
        end
    end
end
%% Plot results

% Make or find figure to plot in
results_figtag = '2dmixsp_results';
fh = findobj('tag', results_figtag);
if isempty(fh)
    fh = figure('tag', results_figtag);
else
    clf(fh);
    figure(fh);
end

% Make new axes,
ax = axes('parent', fh);
% Add plots of efficiency vs red input energy for both the LG_00+LG_10 and LG_00 cases
ph = plot(ax, red1_red2_energy_input_vector*2*1e6, efficiencies(:,1)*100, 'r-', ...
    red1_red2_energy_input_vector*2*1e6, efficiencies(:,2)*100, 'b-');
xlabel(ax, 'Input energy [{\mu}J]');
ylabel(ax, 'Doubling efficiency [%]');
legend(ax, 'LG_{00} and LG_{10}', 'LG_{00} only', 'location', 'best');
% You'll see that the LG_00+LG_10 case has higher conversion efficiency than the LG_00
% case. The first case has narrower profile, and higher peak irradiance than second case.