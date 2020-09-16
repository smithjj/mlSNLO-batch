
% first define some values we'll use in specifying set of inputs for
% calling snlo_qmix_func
crystal_name = 'CLBO';           % name of crystal. string. must match an entry in the file crystal_list.txt.
temperature = 353;              % temperature in kelvin for use in temperature dependent Sellmeier calculations.
plane = 'YZ';                   % plane (ignored for uniaxial crystals)
wavelengths = [532,532,0];       % wavelengths (in nm) (here one must be zero)
type = 'Mix';                   % mixing type ('Mix' or 'OPO')

temperature_vec = 270:10:370;

% Qmix saves the contents of the large edit box at the bottom of the
% figure to disk (in the file "Qmix.dat"). This file contains a copy of the
% text strings populating the edit box. The following section loads the
% Qmix.dat file and parses is it
outputfilename = 'Qmix.dat';

% find which element in the list of crystals matches the crystal named
% above, which is what snlo_qmix_func needs to properly select the
% crystal in question
qmix_crystal_list = importdata('crystal_list.txt');
crystal_list_number = find_selection_number_from_string(qmix_crystal_list,crystal_name); % find the position in the list of crystals which corresponds to the name of the crystal you specified

% construct data structure to use as input when calling snlo_qmix_func
input_set.qmix_selected_crystal = crystal_list_number;
input_set.qmix_temperature = temperature;
input_set.qmix_principal_plane = plane;
input_set.qmix_wavelength_red1 = wavelengths(1); 
input_set.qmix_wavelength_red2 = wavelengths(2); 
input_set.qmix_wavelength_blue = wavelengths(3); 
input_set.qmix_type = type;

inputs = cell(size(temperature_vec));
output_wavelengths = cell(size(temperature_vec));
output_polarizations = cell(size(temperature_vec));
output_walkoffs = cell(size(temperature_vec));
output_phase_velocities = cell(size(temperature_vec));
output_group_velocities = cell(size(temperature_vec));
output_gdd = cell(size(temperature_vec));
output_theta = cell(size(temperature_vec));
output_d_eff = cell(size(temperature_vec));
output_SoL2 = cell(size(temperature_vec));
output_angle_tol = cell(size(temperature_vec));
output_temp_range = cell(size(temperature_vec));
output_mix_accpt_ang = cell(size(temperature_vec));
output_mix_accpt_bw = cell(size(temperature_vec));

for K = 1:length(temperature_vec)
    inputs{K} = input_set;
    inputs{K}.qmix_temperature = temperature_vec(K);

    % call snlo_qmix_func with set of inputs, get list of function handles
    % for run, print buttons and other functions we don't care about here
    function_handles = snlo_qmix_func(inputs{K});
    run_function = function_handles{1};

    % simulate pressing the 'Run' button
    run_function();

    % First, we have to find how many results for phasematching Qmix produces
    % (it's 0, 1, 2, or 3). To do this, we'll load all of the contents of
    % Qmix.dat as a single string and find how many times the word 'Walkoff'
    % appears (which should happen once for each phasematch). Then. we'll close
    % Qmix.dat and open it again
    fid = fopen(outputfilename, 'r');
    filestr = fscanf(fid,'%s');
    num_phasematches = length(strfind(filestr,'Walkoff'));
    frewind(fid);

    % If no phasematches were found, so there's nothing to read and parse
    if num_phasematches==0 
        return;
    end

    output_wavelengths{K} = zeros(3,num_phasematches);
    output_polarizations{K} = zeros(3,num_phasematches);
    output_walkoffs{K} = zeros(3,num_phasematches);
    output_phase_velocities{K} = zeros(3,num_phasematches);
    output_group_velocities{K} = zeros(3,num_phasematches);
    output_gdd{K} = zeros(3,num_phasematches);
    output_theta{K} = zeros(1,num_phasematches);
    output_d_eff{K} = zeros(1,num_phasematches);
    output_SoL2{K} = zeros(1,num_phasematches);
    output_angle_tol{K} = zeros(1,num_phasematches);
    output_temp_range{K} = zeros(1,num_phasematches);
    output_mix_accpt_ang{K} = zeros(2,num_phasematches);
    output_mix_accpt_bw{K} = zeros(2,num_phasematches);

    for J = 1:num_phasematches

        line_contents = fgetl(fid);
        parsed_contents = sscanf(line_contents, '%f(%c) +  %f(%c)  =  %f(%c)');
        output_wavelengths{K}(:,J) = [parsed_contents(1), parsed_contents(3), parsed_contents(5)];
        output_polarizations{K}(:,J) = char([(parsed_contents(2)),(parsed_contents(4)),...
            (parsed_contents(6))]);

        line_contents = fgetl(fid);
        output_walkoffs{K}(:,J) = sscanf(line_contents, 'Walkoff [mrad]      = %f %f %f');

        line_contents = fgetl(fid);
        output_phase_velocities{K}(:,J) = sscanf(line_contents, 'Phase velocities    = c/ %f %f %f');

        line_contents = fgetl(fid);
        output_group_velocities{K}(:,J) = sscanf(line_contents, 'Group velocities    = c/ %f %f %f');

        line_contents = fgetl(fid);
        output_gdd{K}(:,J) = sscanf(line_contents, 'GrpDelDisp(fs^2/mm) = %f %f %f');

        line_contents = fgetl(fid);
        output_theta{K}(:,J) = sscanf(line_contents, 'At theta            = %f deg.');

        line_contents = fgetl(fid);
        output_d_eff{K}(:,J) = sscanf(line_contents, 'd_eff               = %f pm/V');

        line_contents = fgetl(fid);
        output_SoL2{K}(:,J) = sscanf(line_contents, 'S_o * L^2           = %f watt');

        line_contents = fgetl(fid);
        output_angle_tol{K}(:,J) = sscanf(line_contents, 'Crystal ang. tol.   = %f mrad-cm');

        line_contents = fgetl(fid);
        output_temp_range{K}(:,J) = sscanf(line_contents, 'Temperature range   = %f K-cm');

        line_contents = fgetl(fid);
        output_mix_accpt_ang{K}(:,J) = sscanf(line_contents, 'Mix accpt ang   = %f %f mrad-cm');

        line_contents = fgetl(fid);
        output_mix_accpt_bw{K}(:,J) = sscanf(line_contents, 'Mix accpt bw    = %f %f cm^-1-cm');

        % the last line will be blank
        line_contents = fgetl(fid);
    end
    fclose(fid);

end

blue_walkoffs = cellfun(@(x)x(3), output_walkoffs);
fh = figure;
ax = axes('Parent',fh);
subplot(2,1,1,ax);
ph = plot(ax, temperature_vec, blue_walkoffs);
xlabel(ax, 'Temperature [K]');
ylabel(ax, 'Blue walkoff [mrad]');

thetas = cellfun(@(x)x, output_theta);
ax2 = axes('Parent',fh);
subplot(2,1,2,ax2);
ph = plot(ax2, temperature_vec, thetas);
xlabel(ax2, 'Temperature [K]');
ylabel(ax2, '{\Theta} [deg.]');

