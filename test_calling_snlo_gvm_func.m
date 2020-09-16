% This example demonstrates how to set up the problem parameters for the group velocity mismatch mode,
% call the SNLO gvm function with the parameters, run the model, load and process the model outputs. We
% demonstrate several sets of similar input parameters, with a single input parameter varying from run to
% run. In this case, we'll vary the blue wave slant angle and find its effect on the apparent group
% velocities of the three pulses.

calculate_curves = true;
if calculate_curves
    slant_angle_list = linspace(0,5,6); % make vector of pummp slant angles between 0 and 5 degrees with 6 elements
    
    % based on gvm example 1 (#4 from snlo examples button)
    
    % then define the full set of input parameters for gvm; these correspond to the input boxes
    % visible when you run gvm. for most parameters, a row vector with three elements is defined
    % and corresponds to the values for the red1, red2, and blue waves. calling the gvm function
    % with a set of input parameters requires specifying a single input argument to the function
    % snlo_gvm_func, where the input parameter is a data structure with the fieldnames specified
    % below, and values in the form specified below. 
    
    % In this case, the crystal is passed to the snlo_gvm_func as its number in the list of crystals
    % (which changes if we add crystals). We use the function 'find_selection_number_from_string' to
    % determine which number corresponds to the string currentCrystalName.
    
    crystal_list = importdata('crystal_list.txt'); % load list of crystals
    inputs.currentCrystalName = 'BBO';
    inputs.gvm_currentCrystalValue = find_selection_number_from_string(crystal_list,inputs.currentCrystalName);
    inputs.gvm_temperature = 293;                   % temperature (in K)
    inputs.gvm_wavelengths = [0,825,532];           % wavelengths (in nm) (here one must be zero)
    inputs.gvm_slant_angle = 0;                     % pump slant angle (in degrees)
    inputs.gvm_polarization = 'ooe';                % red1,red2,blue wave polarization
    inputs.gvm_plane = 'XZ';                        % not actually used for the uniaxial crystal BBO, but would be used for a biaxial crystal
    
    %% call gvm, run model, load output
    clear problem; % problem will be a vector of data structures with appropriate fieldnames and values, but will change size in the loop so it should be undefined prior to starting
    outputs = cell(size(slant_angle_list)); % pre-allocate cell array to load gvm model outputs into
    for K = 1:length(slant_angle_list) % 
        problem(K) = inputs;        % copy the inputs specified earlier as the starting point for this problem's inputs
        problem(K).gvm_slant_angle = slant_angle_list(K); % modify this problem's slant angle to be an element of the vector of slant angles generated earlier
        fcn_handles = snlo_gvm_func(problem(K)); % call the SNLO gvm function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        run_fcn = fcn_handles{1};   % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        run_fcn();                  % call the 'Run' button callback - equivalent to clicking the run button
        outputs{K} = importdata('gvm.dat'); % load the gvm output from file and stick it in an output cell array element
    end
end

%% process and plot the output
fh = figure;                                % make new figu re
thepos = get(fh,'Position');                % get the figure position property
set(fh,'Position',[thepos(1),50,800,800]);  % set figure position property to be a large 800x800 pixel window
ax = axes('Parent',fh);                     % make new axes in the newly created figure
gvmat = cell(size(slant_angle_list));       % pre-allocate cell array to stick the red1,red2, and blue group velocity arrays in
anglemat = cell(size(slant_angle_list));    % pre-allocate cell array to stick the delta,gamma, and theta angle arrays in
pos_gv = cell(size(slant_angle_list));      % pre-allocate cell array to stick the group velocities associated with positive delta angles in
pos_ang = cell(size(slant_angle_list));     % pre-allocate cell array to stick the delta,gamma and theta angles associated with positive delta angles in
gvm_red1 = cell(size(slant_angle_list));    % pre-allocate cell array to stick the red1 group velocities converted to time relative to the blue wave in
gvm_red2 = cell(size(slant_angle_list));    % pre-allocate cell array to stick the red2 group velocities converted to time relative to the blue wave in

phs = zeros(length(slant_angle_list),2); % pre-allocate an array to stick plot handles in

for K = 1:length(slant_angle_list) % loop over each element in vector of slant angles 
    gvmat{K} = outputs{K}(:,[4,5,6]);   % assign red1,red2, and blue group velocities to cell array element 
    anglemat{K} = outputs{K}(:,[1,2,3]);% assign delta,gamma, and theta angles to cell array element 
    pos_gv{K} = gvmat{K}(anglemat{K}(:,1)>0,:); % take only subset of group velocities which correspond to delta > 0
    pos_ang{K} = anglemat{K}(anglemat{K}(:,1)>0,:); % take only subset of delta,gamma, and theta angles with delta > 0
    gvm_red1{K} = 3333*(pos_gv{K}(:,1)-pos_gv{K}(:,3)); % convert red1 group velocities to be time walk off relative to blue wave group velocities (in fs/mm)
    gvm_red2{K} = 3333*(pos_gv{K}(:,2)-pos_gv{K}(:,3)); % convert red2 group velocities to be time walk off relative to blue wave group velocities (in fs/mm)
    % consider only subset of angles and group velocities which correspond to theta between theta_limits
    theta_limits = [20,40];
    inds = find((pos_ang{K}(:,3)>=theta_limits(1))&(pos_ang{K}(:,3)<=theta_limits(2)));
    % plot this subset as a pair of curves in the axis
    phs(K,:) = plot(ax,pos_ang{K}(inds,3),gvm_red1{K}(inds),'r-',...
        pos_ang{K}(inds,3),gvm_red2{K}(inds),'b--');
    hold(ax,'on'); % tell axis the next plot commands add plot curves instead of replace existing ones
    
    % add a pair of text labels in the axis located on the right side next to the curves you just plotted
    % which indicate the pump slant angle used to calculate curves
    thr1(K) = text(theta_limits(2),gvm_red1{K}(inds(end)),sprintf('slant = %.2f deg.',slant_angle_list(K)),'FontSize',8);
    thr2(K) = text(theta_limits(2),gvm_red2{K}(inds(end)),sprintf('slant = %.2f deg.',slant_angle_list(K)),'FontSize',8);
end
legend(phs(1,:),'Red1','Red2','Location','Best'); % add legend to axis which labels red1 and red2 line styles
xlabel(ax,'Theta [deg.]');  % label x axis
ylabel(ax,'Delay [fs/mm]'); % label y axis
xlim(ax,[22,45]);           % set x axis limits