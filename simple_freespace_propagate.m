function output_field = simple_freespace_propagate(wavelength, xgridmat, ygridmat, input_field, distance)
% Simple freespace propagation (in air) version 1.1
% Created by Jesse Smith (jesse.smith@as-photonics.com), last modified October 26, 2023
%
% Inputs:
%  wavelength [m]  (scalar)
%  xgridmat   [m]   (2D array)          constant value in 2nd dimension
%  ygridmat   [m]   (2D array)          constant value in 1st dimension
%  input_field [V/m] (2D or 3D array)    electric field in x,y, or x,y,time
%  distance   [m]   (scalar)            distance to propagate
%
% This function takes the input values, creates a 2D complex-valued array of phases to
% propagate the input field by distance specified, loops over 3rd dimension of input_field
% to take 2D fft of input_field, multiplies by the phase, and takes ifft.
%
% This method of beam propagation is simple and fast, but it is not a universally good
% method to use. For instance, if your beam overflows the spatial grid (if it expands a
% lot, or if is is tilted, for instance) the FFT BPM will place the light exiting a face
% at the opposite face. This is, I believe, sort of an aliasing effect. This could be
% handled by padding the spatial grids, through some sort of resampling method, or having
% some absorbing boundary at the transverse faces.

% Changelog:
% 1.0 -> 1.1 (Oct 26, 2023):
%  - Add 'arguments' section to validate shape and values of input variables,
%  - Change input array size checking to be 'asserts' and make error string display the actualy sizes of input variables that don't match
arguments
    wavelength  (1,1)   double {mustBePositive, mustBeFinite, mustBeReal}           = 1000e-9;                                      % Wavelength of light in field [m]
    xgridmat    (:,:)   double {mustBeFinite, mustBeReal}                           = repmat(linspace(-1e-3,1e-3,64).', [1, 64]);   % 2D array of x positions [meters]. Default values are equally spaced between -1 and 1 mm, and gridded such that the values are constant over the second dimension.
    ygridmat    (:,:)   double {mustBeFinite, mustBeReal}                           = repmat(linspace(-1e-3,1e-3,64), [64, 1]).';   % 2D array of y positions [meters]. Default values are equally spaced between -1 and 1 mm, and gridded such that the values are constant over the first dimension.
    input_field (:,:,:) double {mustBeFinite}                                       = zeros(size(xgridmat,1), size(xgridmat,2));    % 2D or 3D array of complex-valued electric field distribution [V/m]; default value is a 2D array of zeros with same dimensions as xgrid and ygrid
    distance    (1,1)   double {mustBeFinite, mustBeReal}                           = 1e-3;                                         % distance to propagate [meters].
end

% Check that same xgrid, ygridmat have same size, and match the size of the first two dimensions of input_field:
assert( all(size(xgridmat) == size(ygridmat)), 'The xgridmat and ygridmat input variables must have the same size, but xgridmat has size %ix%i while ygridmat has size %ix%i', size(xgridmat,1), size(xgridmat,2), size(ygridmat, 1), size(ygridmat, 2));
assert( (size(xgridmat, 1) == size(input_field, 1)) && (size(xgridmat, 2) == size(input_field, 2)), 'The xgridmat&ygridmat size must match the first two dimensions of input_field, but xgridmat has size %ix%i while input_field has size %ix%i', size(xgridmat,1), size(xgridmat,2), size(input_field, 1), size(input_field, 2));

nx = size(input_field,1);
ny = size(input_field,2);
nt = size(input_field,3);
dx = xgridmat(2,2) - xgridmat(1,1);
dy = ygridmat(2,2) - ygridmat(1,1);


phasex  = repmat(pi*(nx*dx).^(-2).*(wavelength),[1,nx]).*circshift(((1:nx)-0.5*nx).^2,[0,-(0.5*nx-1)]);
phasey  = repmat(pi*(ny*dy).^(-2).*(wavelength),[1,ny]).*circshift(((1:ny)-0.5*ny).^2,[0,-(0.5*ny-1)]);

% assemble the phasex and phasey vectors into 2D meshgrid
prop_phase_mat  = (repmat(phasex.',[1,ny])+repmat(phasey,[nx,1]));
% use the new 2D meshgrid variable to define anonymous function to multiply by -1i*distance, stick in exponent
prop_phase_func = @(distance) exp(1i*-distance*prop_phase_mat);

% anonymous function to propagate a field using phase array
propagate_func  = @(fieldmat, phasemat) ifftn(fftn(fieldmat).*phasemat);

% pre allocate an array to stick propagated fields into
output_field = zeros(size(input_field));

% loop through the input_field 3rd dimension ( probably time dimension, for electric field array of pulse wave )
for K = 1:nt
    output_field(:,:,K) = propagate_func(input_field(:,:,K), prop_phase_func(distance));
end
end
