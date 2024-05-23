% Demonstration of how to use custom pulse profiles in 2D-mix-SP
%
% This example is the same as the script demonstrating batch calls of
% 2D-mix-SP with the capability of specifying custom temporal pulse
% profiles. It demonstrates how to generate a set of inputs, including a
% custom temporal pulse profile for each wave (in both power and phase). It
% shows how to set up the problem parameters for the 2d short-pulse mixing
% function 2D-MIX-SP (with custom pulse profiles), call the model with the
% parameters, run the model, load the model outputs, process and plot them.
%
% Custom pulse profiles require a time grid, a power at each time point
% for each wave, and a phase at each time point for each wave. The time
% grid must accomodate the full pulse profiles and also temporal walk off
% between the pulses plus temporal spreading due to group velocity
% dispersion. SNLO normally calculates a suitable time grid taking these
% into account. We include this calculation below ("Standard SNLO time grid
% calculation") but you are free to change it.
% 
% To use a custom pulse profile, 2D-mix-sp must be provided with the time vector (in
% seconds) upon which the values of power and phase are defined as a field named
% 'mix_2d_sp_timevector' in the inputs data structure. Also required are the powers and
% phases which are stored in the data structure mix_2d_sp_pulse_profiles with fieldnames
% red1.power, red1.phase, red2.power, red2.phase, blue.power, and blue.phase. See line
% 200.
%
% When you pass a custom pulse profile to 2D-mix-SP their magnitudes will
% be normalized to contain the pulse energy specified. You just specify a
% shape with no need to worry about the absolute magnitude. The phase you
% specify will not be modified.
%
% To avoid confusion when using custom temporal profiles the normal
% interface is changed so the edit boxes on the input form related to pulse
% profiles cannot be modified outside of this script. 
%
% If you already have a full set of output data and just want to process
% and plot it, you can elect to skip the calculation by setting the
% 'calculate_curves' variable to false.
calculate_curves = true;
if calculate_curves
    %% Generate all the inputs
    % In this example we will loop on a set of Delta k values.
    delta_k_list = linspace(-0.1,0.1,11); % make a vector of phase mismatches, spanning -.25 to .25 per mm with 15 elements.
   
    c = 3e8; % speed of light in m/s
        
    % Define the full set of input parameters for 2D-MIX-SP; these parameters correspond to the input
    % boxes visible when you run the 2D-MIX-SP function. For most parameters, a row vector with three
    % elements is defined and corresponds to the values for the red1, red2, and blue waves. calling the
    % 2D-MIX-SP function with a set of input parameters requires specifying a single input argument to
    % the function snlo_2d_mix_sp_func, where the input parameter is a data structure with the fieldnames
    % specified below, and values in the form specified below.
    
    % In this example we'll be frequency doubling a 1064 nm pulse using type I phase matching. In order
    % to model that, we'll split the red wave into two equal components
    % (which must be identical).
    
    clear inputs

    inputs.mix_2d_sp_wavelengths = [1064,1064,532];     % wavelengths (in nm) for each wave
    inputs.mix_2d_sp_ref_inds = [1.654,1.654,1.654];    % refractive indices for each wave
    inputs.mix_2d_sp_gvi = [1.674,1.674,1.774];         % group velocity indices for each wave
    inputs.mix_2d_sp_gdd = [44, 44, 192];               % group delay dispersion (in fs^2/mm) for each wave
    inputs.mix_2d_sp_phase = [0,0,0];                   % phases at input (in rad) for each wave
    inputs.mix_2d_sp_input_refl = [0,0,0];              % input face power reflectivity (0-1) for each wave
    inputs.mix_2d_sp_output_refl = [0,0,0];             % output face power reflectivity (0-1) for each wave
    inputs.mix_2d_sp_crystal_losses = [0,0,0];          % linear power absorption coefficient (in 1/mm) for each wave
    inputs.mix_2d_sp_n2_red1 = [0,0,0];                 % red1 nonlinear refractive index (in sq cm/W) for each wave
    inputs.mix_2d_sp_n2_red2 = [0,0,0];                 % red2 nonlinear refractive index (in sq cm/W) for each wave
    inputs.mix_2d_sp_n2_blue = [0,0,0];                 % blue nonlinear refractive index (in sq cm/W) for each wave
    inputs.mix_2d_sp_beta_red1 = [0,0,0];               % red1 two-photon absorption coefficient (in cm/W) for each wave
    inputs.mix_2d_sp_beta_red2 = [0,0,0];               % red2 two-photon absorption coefficient (in cm/W) for each wave
    inputs.mix_2d_sp_beta_blue = [0,0,0];               % blue two-photon absorption coefficient (in cm/W) for each wave
    inputs.mix_2d_sp_pulseenergy = [2e-6,2e-6,0];       % input pulse energies (in J) for each wave
    inputs.mix_2d_sp_pulse_durations = [5,5,5;1,1,1].'; % input pulse durations  (in fwhm ps), with optional (super) gaussian coefficient for each wave. One wave per column.
    inputs.mix_2d_sp_pulse_delays = [-2.5,-2.5,0];      % pulse delay for each red wave relative to the blue pulse (in ps). The third value by definition is 0 but we'll include three values for the sake of code readability.
    inputs.mix_2d_sp_pulse_chirps = [0,0,0; 0,0,0; 0,0,0].'; % First-, second-, and third-order chirp for each wave (in THz/ps, THz/ps^2, and THz/ps^3). One wave per row.
    inputs.mix_2d_sp_beam_diameters = [0.885,0.885,0.885;0.885,0.885,0.885].'; % input beam diameter (in fwhm mm), specified in walkoff direction, and perpendicular to walkoff direction for each wave. One wave per column.
    inputs.mix_2d_sp_supergaussian_coeff = [1,1,1];     % spatial super gaussian coefficients for each wave
    inputs.mix_2d_sp_wo_angles = [0,0,55.85];           % spatial walkoff angles (in mrad) for each wave
    inputs.mix_2d_sp_offset_wodir = [0,0];              % spatial offset relative to blue beam center in walk off direction (in mm)
    inputs.mix_2d_sp_rad_curv = [1e9,1e9,1e9;1e9,1e9,1e9].';  % input beam radii of curvature (in mm in air) for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction]
    inputs.mix_2d_sp_nt= 128;                           % number of time points to model
    inputs.mix_2d_sp_nxny = [96,64];                    % number of spatial grid points in the walkoff direction, and perpendicular to walkoff direction
    inputs.mix_2d_sp_crystal_length = 17.5;             % length of nonlinear crystal (in mm)
    inputs.mix_2d_sp_lx_ly = [4.25,3];                  % size of spatial grid (in mm) in walkoff direction, and perpendicular to walkoff direction
    inputs.mix_2d_sp_deff = 2.01;                       % crystal nonlinear coefficient (in pm/V)
    inputs.mix_2d_sp_deltak = 0;                        % phase mismatch (in rad/mm) (In this example we'll reset this value on each call to the 2D-mix-SP function)
    inputs.mix_2d_sp_nz = 50;                           % number of (z) integration steps through the crystal
    inputs.mix_2d_sp_dist_to_image = 0;                 % Distance from the crystal exit face to image plane (in mm) of the output beam profile displayed during the run and by the post-run fluence and movie buttons.
    
           
    %% Standard SNLO time grid calculation
    % We must also construct a vector of times. Unlike the free, old
    % versions of SNLO, nt isn't restricted to powers of 2 here. We'll use
    % some of the values we've assigned above in the inputs data structure
    % to create two new fields in the structure with fieldnames
    % 'mix_2d_sp_timevector' and 'mix_2d_sp_pulse_profiles'. The first of
    % these must be a vector of length nt with a range determined by the
    % 'inputs' assigned above. Pulse profiles must include the phase and
    % amplitude at each time point for each wave, so the variable
    % mix_2d_sp_pulse_profiles has size nt x 3. If the input phase of a
    % wave varies in time, then the pulse profile must be complex-valued
    % and include this phase.
    
    % You can replace this with your own time grid generating function.
    % Here we'll generate the vector of times. This vector will have length
    % nt. Our convention defines a value of t=0 at the center of the blue
    % input pulse. Importantly, our vector of times will not necessarily be
    % symmetric about t=0; the values should have a range sufficient to
    % contain the input pulses and any temporal walk off and any group
    % velocity dispersion.
    
    % First, assign some shorter variable names and convert values to MKS
    L_crystal = inputs.mix_2d_sp_crystal_length*1e-3;           % (mm->m) Length of crystal
    gvi = inputs.mix_2d_sp_gvi;                                 % Group velocity indices
    gdd = inputs.mix_2d_sp_gdd*1e-27;                           % (fs^2/mm -> s^2/m) Group delay dispersions
    
    % Convert the pulse chirps (specified above with units of THz/ps,
    % THz/ps^2 and THz/ps^3 for the three lowest degrees of chirp).
    pulse_chirps = inputs.mix_2d_sp_pulse_chirps .* ...         % (THz/ps->Hz/s, THz/ps^2->Hz/s^2, THz/ps^3->Hz/s^3) pulse chirps for 
        repmat([1e24, 1e36, 1e48], [3,1]);                      % 1st, 2nd, 3rd order for each wave
    nt = inputs.mix_2d_sp_nt;
    pulse_delays = inputs.mix_2d_sp_pulse_delays*1e-12;         % (ps->s) Delay of pulse of each red wave relative to blue
    sg_coeffs = inputs.mix_2d_sp_pulse_durations(:,2).';        % Temporal supergaussian coefficient for each wave
    durations = inputs.mix_2d_sp_pulse_durations(:,1).'*1e-12;  % (ps->s) Pulse duration for each wave (FWHM)
    
    t_0 = min(L_crystal*(gvi-gvi(3))/c + pulse_delays-2.5*durations); % Time farthest before middle of blue pulse at crystal input
    t_1 = max(L_crystal*(gvi-gvi(3))/c + pulse_delays+2.5*durations); % Time farthest after middle of blue pulse at crystal input
    t_gvd = max (abs (gdd*L_crystal*3./durations));                   
    t_0 = t_0 - t_gvd;
    t_1 = t_1 + t_gvd;
    dt = (t_1 - t_0)./(nt-1);
    t_vec = t_0 + ((1:nt)-1)*dt;

    %% Generate the pulse profile for each wave
    % In this example, we'll use pulse profiles with either hyperbolic
    % secant shape or (super) Gaussian shape for each wave at the crystal
    % input face. Hyperbolic secant shape corresponds to a sg_coeffs value
    % of 0. Otherwise, sg_coeffs will be the degree of the superGaussian
    % profiles.
    
    % Because we'll have to concatenate the profile of each wave, t_vec
    % must be a column vector. Make sure it is one
    t_vec = reshape(t_vec, [length(t_vec),1]);
    
    sech_fwhm_to_tau = 2*log(sqrt(2)+1); % A coefficient to convert from FWHM values to sech width.
    if sg_coeffs(1)~=0 % if not hyperbolic secant (then we have superGaussian profile)
        red1_pulse_power_vec = exp(-log(2)*( ( 2./durations(1).*(t_vec - pulse_delays(1)) ).^2).^sg_coeffs(1) );
    else
        red1_pulse_power_vec = sech(sech_fwhm_to_tau*(t_vec-pulse_delays(1))/durations(1)).^2;
    end
    if sg_coeffs(2)~=0 % if not hyperbolic secant (then we have superGaussian profile)
        red2_pulse_power_vec = exp(-log(2)*( ( 2./durations(2).*(t_vec - pulse_delays(2)) ).^2).^sg_coeffs(2) );
    else
        red2_pulse_power_vec = sech(sech_fwhm_to_tau*(t_vec-pulse_delays(2))/durations(2)).^2;
    end

    % Pulse delays are relative to the blue wave pulse center, so
    % pulse_delays(3) are by definition zero, but these are left in the
    % form used above for clarity's sake.
    if sg_coeffs(3)~=0 % if not hyperbolic secant (then we have superGaussian profile)
        blue_pulse_power_vec = exp(-log(2)*( ( 2./durations(3).*(t_vec - pulse_delays(3)) ).^2).^sg_coeffs(3) );
    else
        blue_pulse_power_vec = sech(sech_fwhm_to_tau*(t_vec-pulse_delays(3))/durations(3)).^2;
    end

    % Normalize each pulse's magnitude to include the energy specified in
    % inputs.mix_2d_sp_pulseenergy. Add a 1e-40 because some part of the
    % code goes nuts if we ever have exactly 0 in xxxx_pulse_power_vec
    red1_pulse_power_vec = red1_pulse_power_vec ./ ...
        (dt*trapz(red1_pulse_power_vec)) * (inputs.mix_2d_sp_pulseenergy(1) + 1e-40); 
    red2_pulse_power_vec = red2_pulse_power_vec ./ ...
        (dt*trapz(red2_pulse_power_vec)) * (inputs.mix_2d_sp_pulseenergy(2) + 1e-40);
    blue_pulse_power_vec = blue_pulse_power_vec ./ ...
        (dt*trapz(blue_pulse_power_vec)) * (inputs.mix_2d_sp_pulseenergy(3) + 1e-40);

    % Generate a vector of phase vs time for each wave using the 1st, 2nd
    % and 3rd order chirps specified in pulse_chirps.
    red1_pulse_phase_vec =  exp(1i*2*pi*(1/2 * pulse_chirps(1,1)*(t_vec-pulse_delays(1)).^2 + ...
        (1/3 * pulse_chirps(1,2)*(t_vec-pulse_delays(1)).^3) + ...
        (1/4 * pulse_chirps(1,3)*(t_vec-pulse_delays(1)).^4 )));
    red2_pulse_phase_vec =  exp(1i*2*pi*(1/2 * pulse_chirps(2,1)*(t_vec-pulse_delays(2)).^2 + ...
        (1/3 * pulse_chirps(2,2)*(t_vec-pulse_delays(2)).^3) + ...
        (1/4 * pulse_chirps(2,3)*(t_vec-pulse_delays(2)).^4 )));
    blue_pulse_phase_vec =  exp(1i*2*pi*(1/2 * pulse_chirps(3,1)*(t_vec-pulse_delays(3)).^2 + ...
        (1/3 * pulse_chirps(3,2)*(t_vec-pulse_delays(3)).^3) + ...
        (1/4 * pulse_chirps(3,3)*(t_vec-pulse_delays(3)).^4 )));

    % Custom pulses to 2d-mix-sp requires a vector of points in time the powers &
    % phases are defined, provided in the inputs data structure with fieldname
    % 'mix_2d_sp_timevector'.
    inputs.mix_2d_sp_timevector = t_vec;

    % Custom pulses also require the pulse profiles to be specified in the inputs data
    % structure with fieldname 'mix_2d_sp_pulse_profiles' which itself is a data structure
    % that contains the values for the red1, red2, and blue waves (with fieldnames that
    % match those labels)

    inputs.mix_2d_sp_pulse_profiles.red1.power = red1_pulse_power_vec;
    inputs.mix_2d_sp_pulse_profiles.red1.phase = red1_pulse_phase_vec;
    inputs.mix_2d_sp_pulse_profiles.red2.power = red2_pulse_power_vec;
    inputs.mix_2d_sp_pulse_profiles.red2.phase = red2_pulse_phase_vec;
    inputs.mix_2d_sp_pulse_profiles.blue.power = blue_pulse_power_vec;
    inputs.mix_2d_sp_pulse_profiles.blue.phase = blue_pulse_phase_vec;

    red1_energies = zeros(size(delta_k_list)); % pre-allocate a vector for red1 output energies as fcn of phase mismatch
    red2_energies = zeros(size(delta_k_list)); % pre-allocate a vector for red2 output energies as fcn of phase mismatch
    blue_energies = zeros(size(delta_k_list)); % pre-allocate a vector for blue output energies as fcn of phase mismatch
    
    %% call 2d mix sp, run model, load output
%     clear problem   % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    problem = cell(size(delta_k_list)); % pre-allocate a cell array to store complete sets of input values for each iteration in the loop 
    output  = cell(size(delta_k_list)); % pre-allocate a cell array to store the contents of the output files in
    for K = 1:length(delta_k_list) % loop over each element in the vector of phase mismatches
        problem{K} = inputs; % copy the inputs specified earlier as the starting point for this problem's inputs
        problem{K}.mix_2d_sp_deltak = delta_k_list(K);  % modify this problem's phase mismatch to be an element of the vector of phase mismatches generated earlier
        fcn_handles = snlo_2d_mix_sp_func(problem{K});  % call the SNLO 2D-MIX-SP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'        
        run_handle = fcn_handles{1};        % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        accept_handle = fcn_handles{2};     % the function handle of the 'Run' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        close_handle = fcn_handles{end};    % the function handle of the 'Close figure' callback is last returned
        accept_handle();                    % call the 'accept' button's callback - equivalent to clicking the accept button
        run_handle();                       % call the 'Run' button callback - equivalent to clicking the run button
        output{K} = load('mix_2d_sp_output.mat'); % load the 2d mix sp output file into memory, stick into the output cell array
%         close_handle(); % call the 'Close figure' callback. This isn't strictly necessary, and skipping it lets you print the results from the last run.
        red1_energies(K) = trapz(output{K}.power(:,1),output{K}.power(:,2))*1e3; % calculate red1 energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
        red2_energies(K) = trapz(output{K}.power(:,1),output{K}.power(:,3))*1e3; % calculate red2 energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
        blue_energies(K) = trapz(output{K}.power(:,1),output{K}.power(:,4))*1e3; % calculate blue energy by integrating the pulse power in time using the trapezoidal method, multiply by 1000 to get into millijoules
        fprintf(1,'Finished %i of %i runs (Delta k = %.3g /mm). Output energies %.4g, %.4g, %.4g uJ.\n', K, length(delta_k_list), delta_k_list(K), red1_energies(K)*1e3, red2_energies(K)*1e3, blue_energies(K)*1e3);
    end
end


%% plot processed output
fh = figure; % make new outut figure
ax = axes('Parent',fh); % make new axes in that output figure
% plot the red1, red2, and blue energies as a function of phase mismatch
ph = plot(ax,delta_k_list,red1_energies,'ro-',...
    delta_k_list,red2_energies,'bo-',...
    delta_k_list,blue_energies,'go-');
xlabel(ax,'\Delta k');                                  % label x-axis 
ylabel(ax,'Output energy [mJ]');                        % label y-axis 
lh = legend(ax,'Red1','Red2','Blue','Location','Best'); % make legend to describe each curve on the plot
