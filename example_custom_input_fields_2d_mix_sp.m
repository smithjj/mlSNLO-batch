% 2D-mix-SP example with custom fields and looping through input energies
% Example script demonstrating how to construct and pass to the 2D-mix-SP (diffractive,
% short-pulse mixing) function user-specified electric field input profiles. The method
% used here to set the input values for 2D-mix-SP is based upon the same as shown in upon
% the example script test_calling_snlo_2d_mix_sp_func_deltak but here we add the electric
% fields. The electric field profiles have values for both the transverse (x,y) as well as
% temporal dimensions. The electric field profiles are specified by using an additional
% field in the inputs data structure with fieldname 'input_fields'; this field should be comprised
% of three fieldnames 'field_red1_xyt', 'field_red2_xyt', and 'field_blue_xyt'. Each of these
% must be an object of class mix_2d_sp_input_fields_class. (See the classdef file
% mix_2d_sp_input_fields_class for more information).
%  field_****_xyt are three dimensional (x, y, and t dimensions). 
%  The magnitudes of field_****_xyt are not important; they will be scaled to contain the 
% energy specified in the pulse energy input box (fieldname mix_2d_sp_pulseenergy)
%
%
% Since we will provide 2d-mix-sp with input fields, many of the values stored in the data
% structure are ignored which relate to the transverse grids, temporal grids, and electric
% field amplitude and phase distributions. The fields with following names are
% ignored by 2D-mix-SP if provided with input_fields:
%  mix_2d_sp_pulse_durations, mix_2d_sp_pulse_delays, mix_2d_sp_pulse_chirps, 
%  mix_2d_sp_beam_diameters, mix_2d_sp_supergaussian_coeff, mix_2d_sp_offset_wodir, 
%  mix_2d_sp_rad_curv, mix_2d_sp_nt, mix_2d_sp_nxny, and mix_2d_sp_lx_ly.

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
inputs.mix_2d_sp_pulse_delays = [1,1];              % ** pulse delay relative to the blue pulse (in ps); 2 element vector (for red1, red2)
inputs.mix_2d_sp_pulse_chirps = [0,0,0];            % ** Linear chirp of input pulse, for linear (THz/ps) and optionally quadratic (THz/ps^2) and cubic (THz/ps^3) chirps; requires a 3x1 array (linear chirps only), 3x2 array (linear and quadratic chirps), or 3x3 array (linear, quadratic, and cubic chirps); 1st dimension is (red1, red2, blue) order
inputs.mix_2d_sp_beam_diameters = [0.15,0.15,0.15;0.15,0.15,0.15].'; % ** input beam diameter (in fwhm mm), specified in walkoff direction, and perpendicular to walkoff direction
inputs.mix_2d_sp_supergaussian_coeff = [1,1,1];     % ** spatial super gaussian coefficients; 3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_wo_angles = [0,0,10];              % spatial walkoff angles (in mrad);  3 element vector (for red1, red2, blue)
inputs.mix_2d_sp_offset_wodir = [0,0];              % ** spatial offset relative to blue beam center in walk off direction (in mm);  2 element vector (for red1, red2)
inputs.mix_2d_sp_rad_curv = [1e7,1e7,1e7;1e7,1e7,1e7].';  % ** input beam radii of curvature (in mm in air) for [red1,red2,blue] and optionally 2x3 array for different curvatures in the walkoff direction and perpendicular to it: [in the walkoff direction;parallel to the walkoff direction]. 
% inputs.mix_2d_sp_rad_curv = [1e3,1e3,1e3;1e3,1e3,1e3].';  % ** input beam radii of curvature (in mm in air) for [red1,red2,blue] and optionally 2x3 array for different curvatures in the walkoff direction and perpendicular to it: [in the walkoff direction;parallel to the walkoff direction]. 
inputs.mix_2d_sp_nt = 100;                           % ** number of time points to model
inputs.mix_2d_sp_nxny = [64,50];                    % ** number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
inputs.mix_2d_sp_crystal_length = 5;                % length of nonlinear crystal (in mm)
inputs.mix_2d_sp_lx_ly = [1.2,0.8];                     % ** size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
inputs.mix_2d_sp_deff = 30;                         % crystal nonlinear coefficient (in pm/V)
inputs.mix_2d_sp_deltak = 0;                        % phase mismatch (in rad/mm)
inputs.mix_2d_sp_nz = 30;                           % number of (z) integration steps through the crystal
% **: snlo_2d_mix_sp_func will ignore the values provided in these fields. They are
% defined here so we can use them in this script to generate electric field profiles.

%% Make a copy of some values specified above, and convert to mks units
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
spatial_sgs = inputs.mix_2d_sp_supergaussian_coeff;     % spatial super gaussian coefficient (value of 1 is Gaussian, larger values are more flat-topped)
Rxs = inputs.mix_2d_sp_rad_curv(:,1)*1e-3;              % (meters) radius of curvature in direction of walk off
Rys = inputs.mix_2d_sp_rad_curv(:,2)*1e-3;              % (meters) radius of curvature in direction of walk off

%% Generate spatial (x,y) grids centered about 0, and temporal grid, for use in generating the input electric field profiles
xvec = linspace( -0.5, 0.5, nx)*Lx; % has nx points and span of lx, centered about 0
yvec = linspace( -0.5, 0.5, ny)*Ly; % has ny points and span of ly, centered about 0
tvec = linspace( -3, 3, nt)*max(inputs.mix_2d_sp_pulse_durations(:,1))*1e-12; % has nt points, centered about 0, spans 4 * largest pulse duration
[ygrid, xgrid] = meshgrid(yvec, xvec);

% pre-allocate array we'll fill with electric field profiles for the red1, red2, and blue waves
fields_xyt = zeros(nx,ny,nt,3); % last dimension is for which wave

for wave = 1:3 % loop through the red1, red2, and blue waves

    % make irradiance distribution - Gaussian or super-Gaussian in space
    irrad_x = (2./diams_x(wave) .* (xvec-offsets_x(wave))).^2;
    irrad_y = (2./diams_y(wave) .* yvec).^2;
    irrad_xy = exp(-log(2)*(repmat(irrad_x.', [1,ny]) + repmat(irrad_y, [nx, 1])).^(spatial_sgs(wave)));
    
    % define phase curvature using values from above; just using a spherical curvature
    % with radii of curvature set by Rxs and Rys, centered about the beam center
    curvature = exp(1i * (-pi/wavelengths(wave).*(repmat((xvec.' - offsets_x(wave)).^2,[1,ny])./Rxs(wave) + repmat(yvec.^2,[nx,1])./Rys(wave) )));
    
    % define field from irradiances and curvatures
    field_xy = sqrt(irrad_xy).*curvature;
    
    % make pulse shape - gaussian or super-gaussian in time
    pulse = exp(-log(2)*(((2/durations(wave)).*(tvec-delays(wave))).^2).^(temporal_sgs(wave)));
%     pulse = exp(-log(2)*((2/durations(wave)).*(tvec-delays(wave))).^(2*temporal_sgs(wave)));

    % add temporal dimension to 
    fields_xyt(:,:,:,wave) = repmat(field_xy,[1,1,nt]) .* repmat(reshape(sqrt(pulse),1,1,[]),[nx,ny,1]);
end

% Make an object of class mix_2d_sp_input_fields_class
input_fields = mix_2d_sp_input_fields_class;

% set the object's properties for xgrid, ygrid, and tgrid. xgrid and ygrid are
% two-dimensional arrays and tgrid is a vector.
input_fields.xgrid = xgrid;
input_fields.ygrid = ygrid;
input_fields.tgrid = tvec;
% set the object's properties using the electric field distributions found above
input_fields.field_red1_xyt = fields_xyt(:,:,:,1);
input_fields.field_red2_xyt = fields_xyt(:,:,:,2);
input_fields.field_blue_xyt = fields_xyt(:,:,:,3);

% add the input_fields object to inputs data structure
inputs.input_fields = input_fields;

% vector of energies to loop over (each red wave gets this value so total input energy is twice this value)
% red1_red2_energy_input_vector = linspace(1e-8,1e-6,4);
% red1_red2_energy_input_vector = logspace(-10,log10(1e-7),20);
% red1_red2_energy_input_vector = 1e-10;
red1_red2_energy_input_vector = logspace(-10,-7.5,5);

% preallocate some arrays to store results in
output_energies = zeros(length(red1_red2_energy_input_vector), 3); % for red1, red2, and blue energies
efficiencies    = zeros(length(red1_red2_energy_input_vector), 1); % for red1, red2, and blue doubling efficiencies
clear problem

%% Loop through each input energy value stored in the vector red1_red2_energy_input_vector
for K = 1:length(red1_red2_energy_input_vector)
    problem(K) = inputs; % make a copy of the 'master' inputs
    problem(K).mix_2d_sp_pulseenergy = [1,1,0]*red1_red2_energy_input_vector(K); % modify the problem structure to have energy from vector we're looping over

    % call 2D-mix-SP with input parameter values in problem(K), get back cell array of functions in the model we can call
    functions = snlo_2d_mix_sp_func(problem(K)); 
% return
    % pick out the 'Run' button and 'Accept' button functions
    run_fcn     = functions{1};
    accept_fcn  = functions{2};
    
    % write to the console what we are doing
    fprintf(1, 'Starting shg run %i of %i: input energy = %.4E J\n',K,length(red1_red2_energy_input_vector),red1_red2_energy_input_vector(K)*2);

    accept_fcn(); % simulate pressing the 'Accept' button
    run_fcn();    % simulate pressing the 'Run' button

    % After the 'Run' button function has concluded, load its output from file
    dat = load('mix_2d_sp_output.mat');
    % In dat, the variables we are interested in are called tgrid and power; tgrid is a
    % vector of times in seconds, and power is a 4xnt array which is comprised in the
    % second dimension of the times, output power for red1 wave, output power for red2
    % wave, and output power for blue wave.
    dt = dat.tgrid(2) - dat.tgrid(1);
    output_energies(K,1) = dt * sum( dat.power(:,2) ); % integrate output power vector in time for red1
    output_energies(K,2) = dt * sum( dat.power(:,3) ); % integrate output power vector in time for red1
    output_energies(K,3) = dt * sum( dat.power(:,4) ); % integrate output power vector in time for red1
    efficiencies(K)  = output_energies(K,3) ./ (2*red1_red2_energy_input_vector(K)); % doubling efficiency is second harmonic output energy / fundamental input energy
    
    % update the console on our progress and results
    fprintf(1,' Output energies = %.4E, %.4E, %.4E J (doubling efficiency of %.2f%%)\n', output_energies(K,:), efficiencies(K)*100);
end

%% 
results_figtag = '2dmixsp_results';
fh = findobj('tag', results_figtag);
if isempty(fh)
    fh = figure('tag', results_figtag);
else
    clf(fh);
    figure(fh);
end

ax = axes('parent', fh);
ph = plot(ax, red1_red2_energy_input_vector*2*1e6, efficiencies*100);
xlabel(ax, 'Input energy [{\mu}J]');
ylabel(ax, 'Doubling efficiency [%]');
