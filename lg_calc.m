function varargout = lg_calc(p, r, w)
% Calculates Laguerre-Gaussian modes of order p,0 (just the modes that don't depend on
% angle). For m = 0 only, ignoring any angle dependence of the phase.
% Function inputs: 
%  - order p; scalar integer of value 0 or greater
%  - r radius vector to evaluate mode on, scalar or 1 dimensional vector
%  - w width of mode (same units as r); scalar
% 1, optionally 2 function output arguments:
%  - Output argument 1 is LG_p0 amplitude at vector of positions specified in input argument r
%  - (Optional) output arugment 2 is an array of LG amplitudes from LG_00 to LG_p0

% Block to declare input arguments, enforce shape, class (data type) and validation functions
arguments
    p (1,1) {mustBeInteger, mustBeNonnegative}  % integer p for LG_p0 calculation; scalar
    r (1,:) double      % radius; scalar or vector
    w (1,1) double      % beam width of gaussian part; scalar
end

x = 2*r.^2 ./ w.^2;

% Laguerre numbers of p indices > 2 build upon previous values
if p == 0
    mat = ones(size(r));
elseif p == 1
    mat = zeros(2, length(r));
    mat(1,:) = 1;
    mat(2,:) = -x + 1;
else
    mat = zeros(p+1, length(r));
    mat(1,:) = 1;                   % p=0
    mat(2,:) = -x + 1;              % p=1

    for n = 3:(p+1)
        mat(n,:) = (2*n-3)./(n-1).*mat(n-1,:) - x./(n-1).*mat(n-1,:) - (n-2)./(n-1).*mat(n-2,:);
    end
end

% add gaussian part
g = exp(-r.^2 ./ w.^2);

lg_p0_mat = zeros(p+1, length(r));
for n = 1:(p+1)
    lg_p0_mat(n,:) = mat(n,:).*g.*sqrt(2/pi); % combine Laguerre number polynomial and Gaussian profile (add sqrt(2/pi) normalization factor)
end

varargout{1} = lg_p0_mat(end,:); % LG_p0 amplitude (at vector of positions specified in input argument r)
if nargout == 2
    varargout{2} = lg_p0_mat;
end

end