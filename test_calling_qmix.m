% from derrek wilson
%     qmix_crystal_list = importdata('crystal_list.txt');    
% 
%     inputs.qmix_wavelengths = [0,1700,1420];     % wavelengths (in nm)
%     inputs.qmix_temperature = 298;
%     inputs.currentCrystalPopupName = 'CSP';
%     inputs.qmix_currentCrystalPopupValue = find_selection_number_from_string(qmix_crystal_list,inputs.currentCrystalPopupName); % find the position in the list of crystals which corresponds to the name of the crystal you specified
%     inputs.qmix_type = 'mix';
%     
%     handles = snlo_qmix_func(inputs);

% my try:
qmix_crystal_list = importdata('crystal_list.txt');

% inputs.qmix_wavelengths = [0,1700,1420];     % wavelengths (in nm)
inputs.qmix_wavelength_red1 = 0;
inputs.qmix_wavelength_red2 = 1700;
inputs.qmix_wavelength_blue = 1420;
inputs.qmix_temperature = 298;
inputs.currentCrystalPopupName = 'CSP';
inputs.qmix_selected_crystal = find_selection_number_from_string(qmix_crystal_list,inputs.currentCrystalPopupName); % find the position in the list of crystals which corresponds to the name of the crystal you specified
inputs.qmix_type = 'Mix';

handles = snlo_qmix_func(inputs);