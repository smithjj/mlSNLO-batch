% example of how to call simple_freespace_propagate

% first construct x,y position matrices
% generate a Gaussian beam
% call simple_freespace_propagate with position wavelength, x matrix, y matrix, field, and distance to propagate

wavelength  = 1000e-9;  % 1000 m
prop_dist   = 15e-3;    % 15 mm
%% construct x,y position matrices
% Make 1mm x 1mm grid with 100x100 points centered around 0,0
Lx = 1e-3; 
Ly = 1e-3;
nx = 100;
ny = 100;
[ymat, xmat] = meshgrid(linspace(-Ly/2, Ly/2, ny), linspace(-Lx/2, Lx/2, nx));
dx = xmat(2,2) - xmat(1,1);
dy = ymat(2,2) - ymat(1,1);

%% generate electric field this 
% (symmetric) Gaussian beam width (fwhm) of 250 um with flat phase profile, centered about x,y=0
%  (or substitute your own electric field)
diameter = 250e-6;
Irrad = exp(-log(2) * ((2./diameter.*ymat).^2 + (2./diameter.*xmat).^2) );
field = sqrt(Irrad);

%% propagate field with call to simple_freespace_propagate function
propagated_field = simple_freespace_propagate(wavelength, xmat, ymat, field, prop_dist);


%% Visualize input and output irradiances

% make a new figure to plot our results (or re-use a previous made one)
fh = findobj('tag', 'ex_propagation_fig');
if isempty(fh)
    fh = figure('tag', 'ex_propagation_fig');
else
    clf(fh);
    figure(fh);
end

% add a couple of axes to the figure
ax1 = axes('parent', fh);
ax2 = axes('parent', fh);

% make side-by-side plots
subplot(1,2,1,ax1);
subplot(1,2,2,ax2);

% stick mesh plot of input irradiance in left plot
mh1 = mesh(ax1, xmat*1e3, ymat*1e3, Irrad);
xlabel(ax1, 'x [mm]');
ylabel(ax1, 'y [mm]');
zlabel(ax1, 'Input field [arb]');

% stick mesh plot of propagated irradiance in right plot
mh2 = mesh(ax2, xmat*1e3, ymat*1e3, abs(propagated_field).^2);
xlabel(ax2, 'x [mm]');
ylabel(ax2, 'y [mm]');
zlabel(ax2, 'Propagated field [arb]');

daspect(ax1, [1,1,1]);
daspect(ax2, [1,1,1]);
box(ax1,'on');
box(ax2,'on');
