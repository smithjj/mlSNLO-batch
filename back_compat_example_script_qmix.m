% This version of the example_script_qmix.m file is compatible with older versions of
% mlSNLO. They required the crystal to be specified as an integer matching the line number
% of crystal_list.txt matched the character array of the crystal (the shortened
% abbreviation, such as LBO or CLBO). This was an awkward way to operate, and would
% potentially mean any scripts you had set up would be using the wrong crystal if the
% crystal_list.txt file was changed (such as a crystal removed or added, or the list
% resorted). The preferred new way to specify crystal is simply by using the character
% array, but for backward compatibility snlo_qmix_func still handles the crystal being
% provided as an integer. 

% Example: using a batch script to call Qmix. In this example, as the temperature varies
% between 270 and 370 K, find the phasematch angle and blue walk off angle for frequency
% doubling 532 nm in CLBO in the YZ plane.

% This script will first create a data structure which will hold all the inputs for Qmix.
% Next, it will pass the inputs to the Qmix function. Doing so returns handles to the
% internal functions like qmix_run which is the callback function for the 'Run' button.
% Then the script executes the qmix_run function, loads the function's output from the
% ascii output file Qmix.dat, parses the output for important values, and stores it in
% variables.

% This script example is somewhat complicated due to the possibility of having multiple
% phase matches. In fact, the number of phase matches might change as the temperature
% varies. Because of this variable number of phase matches, the output of each Qmix run
% will be stuck into cell arrays. 

% Qmix needs a data structure with fields 'qmix_selected_crystal', 'qmix_temperature',
% 'qmix_principal_plane', 'qmix_wavelength_red1', 'qmix_wavelength_red2',
% 'qmix_wavelength_blue', and 'qmix_type'. The temperature and wavelengths are scalar
% floats with units in kelvin and nanometers respectively.

temperature_vec = 270:10:370; % vector of temperatures we'll loop through

% Define some values to use in specifying set of inputs for calling snlo_qmix_func.
% These values will be stuck in the data structure below.
desired_crystal_name = 'CLBO';      % name of crystal. string. must match an entry in the file crystal_list.txt.

temperature = 353;                  % temperature in kelvin for use in temperature dependent Sellmeier calculations.
plane = 'YZ';                       % plane (ignored for uniaxial crystals); should be a pair of characters 'XY', 'XZ', or 'YZ'
wavelengths = [532,532,0];          % wavelengths (in nm) (here one must be zero)
type = 'Mix';                       % mixing type ('Mix' or 'OPO')

% Note that new versions of mlSNLO handle the crystal being provided as a character array,
% but older versions cannot handle that. They require the qmix_selected_crystal field to
% be an integer which matches the line number of crystal_list.txt You can find that number
% by importing the file crystal_list.txt with the importdata function. The returns a cell
% vector with one element per line. Using the strcmpi function we can test whether each
% element of the cell array matches our character array, and using the find command will
% return which element of that array of logical values have true values in them.
crystal_list    = importdata('crystal_list.txt');
crystal_number  = find(strcmpi(crystal_list, desired_crystal_name));

% Qmix saves the contents of the large edit box at the bottom of the figure to disk (in
% the file "Qmix.dat"). This file contains a copy of the text strings populating the edit
% box. The following section loads the Qmix.dat file and parses is it
outputfilename = 'Qmix.dat';

% construct data structure to use as input when calling snlo_qmix_func
input_set.qmix_selected_crystal = crystal_number;
input_set.qmix_temperature = temperature;
input_set.qmix_principal_plane = plane;
input_set.qmix_wavelength_red1 = wavelengths(1); 
input_set.qmix_wavelength_red2 = wavelengths(2); 
input_set.qmix_wavelength_blue = wavelengths(3); 
input_set.qmix_type = type;

% run through a for loop for each value in temperature_vec.
%   pre-allocate some cell arrays for each iteration's set of inputs
%   and outputs.
%  Output cell arrays used because unknown number of phase matches for each iteration
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

% Loop through each value in temperature_vec
for K = 1:length(temperature_vec)
    % make a copy of the input_set data structure
    inputs{K} = input_set;
    % set the temperature value in qmix_temperature field to 
    inputs{K}.qmix_temperature = temperature_vec(K);

    % call snlo_qmix_func with set of inputs, get list of function handles
    % for run, print buttons and other functions unimportant here
    function_handles = g_snlo_qmix_func(inputs{K});
    run_function = function_handles{1};

    % simulate pressing the 'Run' button
    run_function();

    % First, find how many results for phasematching Qmix produces (it's 0, 1, 2, or 3).
    % To do this, load all of the contents of Qmix.dat as a single string and find how
    % many times the word 'Walkoff' appears (which should happen once for each
    % phasematch). Then rewind to the beginning of Qmix.dat
    fid = fopen(outputfilename, 'r');
    filestr = fscanf(fid,'%s'); % load entire contents into variable filestr 
    num_phasematches = length(strfind(filestr,'Walkoff')); % find length of matches for the word 'Walkoff' in contents of Qmix.dat
    frewind(fid); % go back to the beginning of 

    % If no phasematches were found, so there's nothing to read and parse
    if num_phasematches==0 
        continue;
    end

    % stick some appropriately sized arrays of 0 in the output cell arrays
    output_wavelengths{K} = zeros(3,num_phasematches);          % 3 entries per phase match: wavelength of red1, red2, blue
    output_polarizations{K} = zeros(3,num_phasematches);        % 3 entries: polarizations (o or e) for red1, red2, and blue
    output_walkoffs{K} = zeros(3,num_phasematches);             % 3 entries: walk off of red1, red2, blue
    output_phase_velocities{K} = zeros(3,num_phasematches);     % 3 entries: phase velocity indices of red1, red2, blue
    output_group_velocities{K} = zeros(3,num_phasematches);     % 3 entries: group velocity indices of red1, red2, blue
    output_gdd{K} = zeros(3,num_phasematches);                  % 3 entries: group delay dispersions of red1, red2, blue
    output_theta{K} = zeros(1,num_phasematches);                % 1 entry: phase match angle
    output_d_eff{K} = zeros(1,num_phasematches);                % 1 entry: d effective
    output_SoL2{K} = zeros(1,num_phasematches);                 % 1 entry: S_0 * L^2
    output_angle_tol{K} = zeros(1,num_phasematches);            % angle tolerance
    output_temp_range{K} = zeros(1,num_phasematches);           % temperature range
    output_mix_accpt_ang{K} = zeros(2,num_phasematches);        % acceptance angles
    output_mix_accpt_bw{K} = zeros(2,num_phasematches);         % acceptance bandwidths

    % loop through each phase match for this temperature
    for J = 1:num_phasematches

        line_contents = fgetl(fid);
        % first line of Qmix output: for each wave, wavelength and polarization listed in 
        parsed_contents = sscanf(line_contents, '%f(%c) +  %f(%c)  =  %f(%c)');
        output_wavelengths{K}(:,J) = [parsed_contents(1), parsed_contents(3), parsed_contents(5)];
        output_polarizations{K}(:,J) = char([(parsed_contents(2)),(parsed_contents(4)),...
            (parsed_contents(6))]);

        % second line: walk off angles
        line_contents = fgetl(fid);
        output_walkoffs{K}(:,J) = sscanf(line_contents, 'Walkoff [mrad]      = %f %f %f');

        % third line: phase velocity indices
        line_contents = fgetl(fid);
        output_phase_velocities{K}(:,J) = sscanf(line_contents, 'Phase velocities    = c/ %f %f %f');

        % third line: group velocity indices
        line_contents = fgetl(fid);
        output_group_velocities{K}(:,J) = sscanf(line_contents, 'Group velocities    = c/ %f %f %f');

        % fourth line: group delay dispersions
        line_contents = fgetl(fid);
        output_gdd{K}(:,J) = sscanf(line_contents, 'GrpDelDisp(fs^2/mm) = %f %f %f');

        % fifth line: phase match angle
        line_contents = fgetl(fid);
        output_theta{K}(:,J) = sscanf(line_contents, 'At theta            = %f deg.');

        % sixth line: d effective
        line_contents = fgetl(fid);
        output_d_eff{K}(:,J) = sscanf(line_contents, 'd_eff               = %f pm/V');

        % seventh line: product of characteristic irradiance S_o and square of crystal length
        line_contents = fgetl(fid);
        output_SoL2{K}(:,J) = sscanf(line_contents, 'S_o * L^2           = %f watt');

        % eigth line: crystal angle tolerance
        line_contents = fgetl(fid);
        output_angle_tol{K}(:,J) = sscanf(line_contents, 'Crystal ang. tol.   = %f mrad-cm');

        % ninth line: temperature range
        line_contents = fgetl(fid);
        output_temp_range{K}(:,J) = sscanf(line_contents, 'Temperature range   = %f K-cm');

        % mix accptance angles *** note: for OPO mix type rather than 'Mix' mix type, first 3 characters of sscanf are OPO rather than Mix
        line_contents = fgetl(fid);
        output_mix_accpt_ang{K}(:,J) = sscanf(line_contents, 'Mix accpt ang   = %f %f mrad-cm');

        % mix acceptance bandwidths *** note: for OPO mix type rather than 'Mix' mix type, first 3 characters of sscanf are OPO rather than Mix
        line_contents = fgetl(fid);
        output_mix_accpt_bw{K}(:,J) = sscanf(line_contents, 'Mix accpt bw    = %f %f cm^-1-cm');

        % the last line will be blank
        line_contents = fgetl(fid);
    end
    fclose(fid); % close file handle

end

% for each item in the cell array output_walkoffs, choose the 3rd element (walk off angle for blue wave);
% convert cell array output_walkoffs to vector blue_walkoffs
blue_walkoffs = cellfun(@(x)x(3), output_walkoffs);

% make new figure, populate it with 2 subplots, plot walk off angle vs temperature and phase match angle vs temperature
fh = figure;
ax = axes('Parent',fh);
subplot(2,1,1,ax);
ph = plot(ax, temperature_vec, blue_walkoffs);
xlabel(ax, 'Temperature [K]');
ylabel(ax, 'Blue walkoff [mrad]');

% convert cell array output_theta to vector thetas
thetas = cellfun(@(x)x, output_theta);
ax2 = axes('Parent',fh);
subplot(2,1,2,ax2);
ph2 = plot(ax2, temperature_vec, thetas);
xlabel(ax2, 'Temperature [K]');
ylabel(ax2, '{\theta} [deg.]');

