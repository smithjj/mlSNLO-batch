function output_field = simple_freespace_propagate(wavelength, xgridmat, ygridmat, input_field, distance)
% Inputs:
%  wavelength [nm]  (scalar)
%  xgridmat   [m]   (2D array)          constant value in 2nd dimension
%  ygridmat   [m]   (2D array)          constant value in 1st dimension
%  field      [V/m] (2D or 3D array)    electric field in x,y, or x,y,time
%  distance   [m]   (scalar)            distance to propagate

% make sure 
assert( (ndims(input_field)==2) || (ndims(input_field)==3) );
assert( (ndims(xgridmat)==2) );
assert( (ndims(ygridmat)==2) );

nx = size(input_field,1);
ny = size(input_field,2);
nt = size(input_field,3);
dx = xgridmat(2,2) - xgridmat(1,1);
dy = ygridmat(2,2) - ygridmat(1,1);


phasex  = repmat(pi*(nx*dx).^(-2).*(wavelength),[1,nx]).*circshift(((1:nx)-0.5*nx).^2,[0,-(0.5*nx-1)]);
phasey  = repmat(pi*(ny*dy).^(-2).*(wavelength),[1,ny]).*circshift(((1:ny)-0.5*ny).^2,[0,-(0.5*ny-1)]);

% anonymous function to make 2d array of propagation phase from vectors we just make 
% (and multiply by propagation distance)
prop_phase_func = @(distance) exp(1i*-distance*(repmat(phasex.',[1,ny])+repmat(phasey,[nx,1])));

% anonymous function to propagate a field using phase array
propagate_func  = @(fieldmat, phasemat) ifftn(fftn(fieldmat).*phasemat);

% pre allocate an array to stick propagated fields into
output_field = zeros(size(input_field));

% loop through nt (the input_field 3rd dimension)
for K = 1:nt
    output_field(:,:,K) = propagate_func(input_field(:,:,K), prop_phase_func(distance));
end
end