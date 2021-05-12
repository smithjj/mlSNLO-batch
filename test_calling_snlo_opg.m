clear inputs
inputs.opg_2d_wavelengths = [1200,2400,800];    % wavelengths [nm]
inputs.opg_2d_ref_inds =    [2.14,2.11,2.17];% ref ind
inputs.opg_2d_gvi =         [2.16,2.17,2.18];  % gvi
inputs.opg_2d_gdd =         [0,0,0];        % gdd
inputs.opg_2d_phase =       [0,0,0];            % phases
inputs.opg_2d_input_refl =  [0,0,0];            % input reflectiviy
inputs.opg_2d_output_refl = [0,0,0];            % output reflectivity
inputs.opg_2d_crystal_losses=[0,0,0];           % crystal loss
inputs.opg_2d_n2_red1 =     [0,0,0];        % n2 (red1)
inputs.opg_2d_n2_red2 =     [0,0,0];        % n2 (red2)
inputs.opg_2d_n2_blue =     [0,0,0];        % n2 (blue)
inputs.opg_2d_beta_red1 =   [0,0,0];        % beta (red1)
inputs.opg_2d_beta_red2 =   [0,0,0];        % beta (red2)
inputs.opg_2d_beta_blue =   [0,0,0];        % beta (blue)
inputs.opg_2d_pulseenergy = [0,1e-17,1e-6];     % pulse energies
inputs.opg_2d_pulse_durations = [1,1,1;3,3,1].'; % pulse durations
inputs.opg_2d_pulse_delays =[0,0];          % pulse delays
inputs.opg_2d_pulse_chirps =[0,0,0];        % pulse chirps
inputs.opg_2d_beam_diameters = [0.0963,0.137,0.0781;0.0963,0.137,0.0781].';% beam diam
inputs.opg_2d_supergaussian_coeff = [1,1,1];% supergaussian coeff
inputs.opg_2d_wo_angles =   [0,0,0];        % wo angles
inputs.opg_2d_offset_wodir =[0,0,0];        % offset in wo dir
inputs.opg_2d_rad_curv =    [34.8,35.8,34.6;34.8,35.8,34.6].'; % rad curv
inputs.opg_2d_nz =          300;             % nz
inputs.opg_2d_nxny =        [128,128];        % nx_ny
inputs.opg_2d_grid_duration  = 6;           % grid duration [ps]
inputs.opg_2d_crystal_length = 30;           % crystal length
inputs.opg_2d_lx_ly =       [0.8,0.8];    % Lx, Ly
inputs.opg_2d_deff =        15;             % deff
inputs.opg_2d_deltak =      0;              % deltak
inputs.opg_2d_dist_to_image = 0;            % dist to image
inputs.opg_2d_nt =          256;             % nt

% blue_energy_vec = (5:5:20)*1e-7;
blue_energy_vec = (4:2:20)*1e-9;
n_runs_per_energy = 5;
% max_tilt = 6.5e4;


clear problem
red1_spectra = cell(length(blue_energy_vec),n_runs_per_energy);
red2_spectra = cell(length(blue_energy_vec),n_runs_per_energy);
blue_spectra = cell(length(blue_energy_vec),n_runs_per_energy);
output_energies = zeros(length(blue_energy_vec),3,n_runs_per_energy);
red1_spectral_widths = zeros(length(blue_energy_vec),n_runs_per_energy);
red1_spectral_centroid = zeros(length(blue_energy_vec),n_runs_per_energy);
red2_spectral_widths = zeros(length(blue_energy_vec),n_runs_per_energy);
red2_spectral_centroid = zeros(length(blue_energy_vec),n_runs_per_energy);
blue_spectral_widths = zeros(length(blue_energy_vec),n_runs_per_energy);
blue_spectral_centroid = zeros(length(blue_energy_vec),n_runs_per_energy);
avg_m2 = zeros(3,2,length(blue_energy_vec),n_runs_per_energy);
efficiencies = zeros(length(blue_energy_vec),n_runs_per_energy);
% low_tilt_spectra = cell(length(blue_energy_vec),n_runs_per_energy);
for K = 1:length(blue_energy_vec)
    %generate inputs
    problem(K) = inputs;
    problem(K).opg_2d_pulseenergy(3) = blue_energy_vec(K);
    for J = 1:n_runs_per_energy
        fcn_handles = snlo_opg_func(problem(K));
    %     keyboard;
        run_handle = fcn_handles{1};
        accept_handle = fcn_handles{2};
        spectra_handle = fcn_handles{3};
        analyze_handle = fcn_handles{4};
%         k_perp_filter_handle = fcn_handles{5};
%         k_perp_filter_accept_handle = fcn_handles{6};
        close_handle = fcn_handles{end};
%         keyboard;
        accept_handle();
        run_handle();
        spectra_handle();
        analyze_handle();
%         k_perp_filter_handle(max_tilt);
%         pause(1);
%         k_perp_filter_accept_handle();
%         k_perp_filter_accept_handle(max_tilt);
        output = load('opg_2d_output.mat');
        output_energies(K,1,J) = trapz(output.power(:,1),output.power(:,2));
        output_energies(K,2,J) = trapz(output.power(:,1),output.power(:,3));
        output_energies(K,3,J) = trapz(output.power(:,1),output.power(:,4));
        red1_spectrum = load('OPG_BEAM_3WS.DAT');
        red1_spectra{K,J} = red1_spectrum;
        red1_spectrum = red1_spectrum(:,1:2);
        freqs = red1_spectrum(:,1);
        amps = red1_spectrum(:,2);
        y = cumsum(amps)./sum(amps);
        ind1 = find(y>=0.1,1,'first');
        ind2 = find(y>=0.9,1,'first');
        red1_spectral_widths(K,J) = freqs(ind2) - freqs(ind1);
        red1_spectral_centroid(K,J) = trapz(freqs,freqs.*amps)./trapz(freqs,amps);

        red2_spectrum = load('OPG_BEAM_3WI.DAT');
        red2_spectra{K,J} = red2_spectrum;
        red2_spectrum = red2_spectrum(:,1:2);
        freqs = red2_spectrum(:,1);
        amps = red2_spectrum(:,2);
        y = cumsum(amps)./sum(amps);
        ind1 = find(y>=0.1,1,'first');
        ind2 = find(y>=0.9,1,'first');
        red2_spectral_widths(K,J) = freqs(ind2) - freqs(ind1);
        red2_spectral_centroid(K,J) = trapz(freqs,freqs.*amps)./trapz(freqs,amps);

        blue_spectrum = load('OPG_BEAM_3WP.DAT');
        blue_spectra{K,J} = blue_spectrum;
        blue_spectrum = blue_spectrum(:,1:2);
        freqs = blue_spectrum(:,1);
        amps = blue_spectrum(:,2);
        y = cumsum(amps)./sum(amps);
        ind1 = find(y>=0.1,1,'first');
        ind2 = find(y>=0.9,1,'first');
        blue_spectral_widths(K,J) = freqs(ind2) - freqs(ind1);
        blue_spectral_centroid(K,J) = trapz(freqs,freqs.*amps)./trapz(freqs,amps);

        m2s = load('opg_2d_beam_analysis.mat','avg_msquared');
        avg_m2(:,:,K,J) = m2s.avg_msquared;
        
%          (input_energy - data.output_energies(:,3).')./input_energy;
        efficiencies(K,J) = (blue_energy_vec(K) - output_energies(K,3,J))./blue_energy_vec(K);
%         low_tilt_spectra{K,J} = load('low_tilt_spectra.mat');
    end
    % run opg model
    
    % load output (specifically output energies)
    save('opg_batch.mat','inputs','problem','blue_energy_vec','output_energies',...
        'avg_m2','red1_spectral_widths','red1_spectral_centroid',...
        'red2_spectral_widths','red2_spectral_centroid',...
        'blue_spectral_widths','blue_spectral_centroid','red1_spectra','red2_spectra','efficiencies');
%         'blue_spectral_widths','blue_spectral_centroid','red1_spectra','red2_spectra','low_tilt_spectra');
end