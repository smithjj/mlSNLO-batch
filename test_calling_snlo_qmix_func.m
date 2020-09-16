
% first define some values we'll use in specifying set of inputs for
% calling snlo_qmix_func
crystal_name = 'CLBO';           % name of crystal. string. must match an entry in the file crystal_list.txt.
temperature = 353;              % temperature in kelvin for use in temperature dependent Sellmeier calculations.
plane = 'YZ';                   % plane (ignored for uniaxial crystals)
wavelengths = [532,532,0];       % wavelengths (in nm) (here one must be zero)
type = 'Mix';                   % mixing type ('Mix' or 'OPO')

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

% call snlo_qmix_func with set of inputs, get list of function handles
% for run, print buttons and other functions we don't care about here
function_handles = snlo_qmix_func(input_set);
run_function = function_handles{1};

% simulate pressing the 'Run' button
run_function();

% Qmix saves the contents of the large edit box at the bottom of the
% figure to disk (in the file "Qmix.dat"). This file contains a copy of the
% text strings populating the edit box. The following section loads the
% Qmix.dat file and parses is it
outputfilename = 'Qmix.dat';

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

output_wavelengths = zeros(3,num_phasematches);
output_polarizations = zeros(3,num_phasematches);
output_walkoffs = zeros(3,num_phasematches);
output_phase_velocities = zeros(3,num_phasematches);
output_group_velocities = zeros(3,num_phasematches);
output_gdd = zeros(3,num_phasematches);
output_theta = zeros(1,num_phasematches);
output_d_eff = zeros(1,num_phasematches);
output_SoL2 = zeros(1,num_phasematches);
output_angle_tol = zeros(1,num_phasematches);
output_temp_range = zeros(1,num_phasematches);
output_mix_accpt_ang = zeros(2,num_phasematches);
output_mix_accpt_bw = zeros(2,num_phasematches);

for K = 1:num_phasematches

    line_contents = fgetl(fid);
    parsed_contents = sscanf(line_contents, '%f(%c) +  %f(%c)  =  %f(%c)');
    output_wavelengths(:,K) = [parsed_contents(1), parsed_contents(3), parsed_contents(5)];
    output_polarizations(:,K) = char([(parsed_contents(2)),(parsed_contents(4)),...
        (parsed_contents(6))]);

    line_contents = fgetl(fid);
    output_walkoffs(:,K) = sscanf(line_contents, 'Walkoff [mrad]      = %f %f %f');

    line_contents = fgetl(fid);
    output_phase_velocities(:,K) = sscanf(line_contents, 'Phase velocities    = c/ %f %f %f');

    line_contents = fgetl(fid);
    output_group_velocities(:,K) = sscanf(line_contents, 'Group velocities    = c/ %f %f %f');

    line_contents = fgetl(fid);
    output_gdd(:,K) = sscanf(line_contents, 'GrpDelDisp(fs^2/mm) = %f %f %f');

    line_contents = fgetl(fid);
    output_theta(:,K) = sscanf(line_contents, 'At theta            = %f deg.');

    line_contents = fgetl(fid);
    output_d_eff(:,K) = sscanf(line_contents, 'd_eff               = %f pm/V');

    line_contents = fgetl(fid);
    output_SoL2(:,K) = sscanf(line_contents, 'S_o * L^2           = %f watt');

    line_contents = fgetl(fid);
    output_angle_tol(:,K) = sscanf(line_contents, 'Crystal ang. tol.   = %f mrad-cm');

    line_contents = fgetl(fid);
    output_temp_range(:,K) = sscanf(line_contents, 'Temperature range   = %f K-cm');

    line_contents = fgetl(fid);
    output_mix_accpt_ang(:,K) = sscanf(line_contents, 'Mix accpt ang   = %f %f mrad-cm');

    line_contents = fgetl(fid);
    output_mix_accpt_bw(:,K) = sscanf(line_contents, 'Mix accpt bw    = %f %f cm^-1-cm');

    % the last line will be blank
    line_contents = fgetl(fid);
end
fclose(fid);