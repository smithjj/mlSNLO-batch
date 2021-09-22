% This example demonstrates how to set up the problem parameters for the 2d long-pulse
% mixing function 2D-MIX-LP, call the SNLO 2D-MIX-LP function with the parameters, run the
% model, load the model outputs, process the model outputs and plot them. This example
% demonstrates second-harmonic generation using continuous-wave beams. The red waves are
% small and focus to a waist of 25 micrometers halfway through the 10 mm long crystal. The
% red1 and red2 beam diameters and radii of curvature were calculated by the 'Focus'
% function of mlSNLO.

% Several sets of similar input parameters are used, with the input
% powers for the red1 and red2 waves varying from run to run. The output powers of the
% red1, red2, and blue waves are then plotted as a function of the sum of red1 and red2
% input power.

% Version 1.0 written 9/22/2021 by Jesse Smith (jesse.smith@as-photonics.com)
%  Last modified 9/22/2021

calculate_curves = true; % Set this to false if you wish to just alter the processing and plotting but continue using previously calculated data
if calculate_curves
    %% generate inputs
    
    % Make a vector of input powers between 100 W and 500 W with 5 elements. These
    % values are used in each of the red1 and red2 waves, so the total input powers are
    % twice these values
    red1red2_power_list = linspace(100,500,5); 
    
    % Pre-allocate some cell arrays to store the contents of the output files in
    output_red1 = cell(size(red1red2_power_list)); 
    output_red2 = cell(size(red1red2_power_list)); 
    output_blue = cell(size(red1red2_power_list)); 
    
    % Preallocate some arrays for storing the output powers
    red1_output_powers = zeros(size(red1red2_power_list));
    red2_output_powers = zeros(size(red1red2_power_list));
    blue_output_powers = zeros(size(red1red2_power_list));
    
    % define the full set of input parameters for 2D-MIX-LP; these parameters correspond to the input
    % boxes visible when you run the 2D-MIX-LP function. for most parameters, a row vector with three
    % elements is defined and corresponds to the values for the red1, red2, and blue waves. calling the
    % 2D-MIX-LP function with a set of input parameters requires specifying a single input argument to
    % the function snlo_2d_mix_lp_func, where the input parameter is a data structure with the fieldnames
    % specified below, and values in the form specified below.
    clear inputs
    inputs.mix_2d_lp_wavelengths= [1000,1000,500];          % wavelengths (in nm)
    inputs.mix_2d_lp_ref_inds   = [1.655,1.655,1.655];      % refractive indices
    inputs.mix_2d_lp_phase      = [0,0,0];                  % input phases (in radians)
    inputs.mix_2d_lp_input_refl = [0,0,0];                  % power reflectivity at crystal input face
    inputs.mix_2d_lp_output_refl= [0,0,0];                  % power reflectivity at crystal output face
    inputs.mix_2d_lp_crystal_losses = [0,0,0];              % linear absorption (in 1/mm)
    inputs.mix_2d_lp_pulseenergy    = [0,0,0];              % input powers (in watts)
    inputs.mix_2d_lp_pulse_durations= [0,0,0];              % input pulse durations of 0 for cw beams. Optionally, include a second row with temporal super gaussian coefficients of 1, but these are ignored for cw.
    inputs.mix_2d_lp_pulse_delays   = [0,0,0];              % input pulse delays relative to blue wave
    inputs.mix_2d_lp_beam_diameters = 58.91e-3*ones(3,2);   % input beam diameters (in fwhm mm) [red1 in walkoff dir, red2 in walkoff dir, blue in walkoff dir;red1 perp to walkoff, red2 perp to walkoff, blue perp to walkoff]
    inputs.mix_2d_lp_supergaussian_coeff = [1,1,1];         % supergaussian coefficient for spatial input beam profile
    inputs.mix_2d_lp_n2_red1    = [0,0,0];                  % red1 nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_n2_red2    = [0,0,0];                  % red2 nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_n2_blue    = [0,0,0];                  % blue nonlinear refractive index (in sq. cm/W)
    inputs.mix_2d_lp_beta_red1  = [0,0,0];                  % red1 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_beta_red2  = [0,0,0];                  % red2 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_beta_blue  = [0,0,0];                  % blue 2-photon absorption cofficient (in cm/W)
    inputs.mix_2d_lp_wo_angles  = [0,0,0];                  % spatial walkoff angles (in milliradians)
    inputs.mix_2d_lp_offset_wodir = [0,0,0];                % spatial offset in walkoff direction (in mm)
    inputs.mix_2d_lp_rad_curv   = 3.6859*[1,1,1]';          % input beam radii of curvature for [red1,red2,blue] and [in the walkoff direction;parallel to the walkoff direction] (in mm in air)
    inputs.mix_2d_lp_nz         = 100;                      % number of z integration steps
    inputs.mix_2d_lp_nxny       = [64,64];                  % number of grid points in the transverse (x,y) directions
    inputs.mix_2d_lp_crystal_length = [10,1];               % crystal length (in mm); optionally, the second parameter is number of walkoff compensating sections
    inputs.mix_2d_lp_lx_ly      = [0.40,0.40];              % size of grid in walkoff direction, perpendicular to walkoff direction (in mm)
    inputs.mix_2d_lp_deff       = 2.01;                     % crystal nonlinear coefficient (in pm/V)
    inputs.mix_2d_lp_deltak     = 0;                        % phase mismatch (in rad/mm)
    inputs.mix_2d_lp_dist_to_image = 0;                     % Distance from the crystal exit face to image plane  (in mm) of the output beam profile displayed during the run and by the post-run fluence and movie buttons.
    inputs.mix_2d_lp_nt         = 30;                       % number of time steps to include (Recommend starting with 32 and varying to check convergence.)
    inputs.mix_2d_lp_save_movie = false;                    % in order to save time set this to false; if set to true, a movie of each wave's irradiance at each point in z is saved to disk
    inputs.mix_2d_lp_auto_analyze = false;                  % in order to save time set this to false; if set to true, the model calculates the beam quality factor M-squared, radius of curvature, and tilt for each wave at each point in z
    
    
    clear problem % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    
    %% call 2d mix lp, run model, load output
    for K = 1:length(red1red2_power_list)               % loop over each element in the vector of input powers
        problem(K) = inputs;                            % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).mix_2d_lp_pulseenergy = [red1red2_power_list(K),red1red2_power_list(K),0]; % modify this problem's red1 and red2 input powers to be those defined in the vector of powers generated earlier
        fcn_handles = snlo_2d_mix_lp_func(problem(K));  % call the SNLO 2D-MIX-LP function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_handle = fcn_handles{1};                    % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        accept_handle = fcn_handles{2};                 % the function handle of the 'Accept' button callback is the first returned; the 'accept' button callback function is always the second element of the cell array of returned function handles, if the function has an 'accept' button
        close_handle = fcn_handles{end};                % the function handle of the 'Close figure' callback is last returned
        accept_handle();                                % call the 'Accept' button's callback - equivalent to clicking the accept button
        run_handle();                                   % call the 'Run' button callback - equivalent to clicking the run button
        
        output_red1{K} = importdata('mix_2d_lp_red1_data.dat'); % load the red1 2d mix lp output file into memory, stick into the output cell array
        output_red2{K} = importdata('mix_2d_lp_red2_data.dat'); % load the red2 2d mix lp output file into memory, stick into the output cell array
        output_blue{K} = importdata('mix_2d_lp_blue_data.dat'); % load the blue 2d mix lp output file into memory, stick into the output cell array
        
        % For cw continuous wave inputs, in the output files mix_2d_lp_xxx_data.dat there
        % are three numbers saved. The first is the time which can be ignored; the second
        % is the output power, and the third is the phase at the crystal output face. Pick
        % out the output powers:
        red1_output_powers(K) = output_red1{K}.data(2);
        red2_output_powers(K) = output_red2{K}.data(2);
        blue_output_powers(K) = output_blue{K}.data(2);
        
        % Print some information to the MATLAB command window about our progress so far
        fprintf(1,'%i of %i: input red power %.3g W; output powers red %.3g, blue %.3g W\n',...
            K,length(red1red2_power_list), 2*red1red2_power_list(K), ...
            red1_output_powers(K) + red2_output_powers(K), ...
            blue_output_powers(K));
    end
end

%% plot output
fh = figure;                    % make new outut figure
ax = axes('Parent',fh);         % make new axes in that output figure
% plot the red1, red2, and blue output powers as a function of total input powers (red1+red2)

ph = plot(ax,2*red1red2_power_list,red1_output_powers,'ro-',...
    2*red1red2_power_list,red2_output_powers,'bo-',...
    2*red1red2_power_list,blue_output_powers,'go-');
xlabel(ax,'Input red1+red2 power [W]');                 % label x-axis 
ylabel(ax,'Output power [W]');                          % label y-axis 
lh = legend(ax,'Red1','Red2','Blue','Location','Best'); % make legend to describe each curve on the plot
