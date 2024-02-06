% Calculate the M-squared beam quality of an input electric field (the electric field
% input argument should be a 2 or 3 dimensional array). This calculation is taken from the
% Crystal Nonlinear Optics: with SNLO examples (second edition) book, see page 299-303
% (chapter 7.7). See as-photonics.com/book for more details about the book.

% By Jesse Smith jesse.smith@as-photonics.com
% Last modified October 30, 2023


% changelog:
% 10/30/23
%  - clean more code
%  - add equation numbers, where they can be found in the Crystal Nonlinear Optics book by Arlee Smith 
% 8/30/23
%  - remove some old code and clean up comments, add more dsecriptions of calculation steps
% 8/26/22
%  - use loop thru nz/nt rather than big 3d arrays (if nt is large, you can
%    literally run out of memory) 
%  - add power_threshold property: if this is 
% 8/17/22
%  - add to results property a fieldname pulse_flattened_Exyt, and flattened_Exyz (which
%    are fields that have had their curvature removed three times)
classdef Msquared_class
    % Calculate beam profile quantities for 2d (x,y) or 3d (x,y,z) or (x,y,t) fields. Including 

    % Make handle? Copyable?
    properties
        %  x position matrix (2d) units m
        xgrid (:,:) double {mustBeFinite,mustBeReal} = zeros(32,32);
        % y position matrix (2d) units m
        ygrid (:,:) double {mustBeFinite,mustBeReal} = zeros(32,32);
        
        % electric field array (2d or 3d) units V/m (but can be arbitrary I think)
        field (:,:,:) double {mustBeFinite} = zeros(32,32);
        % wavelength units m
        wavelength (1,1) double {mustBePositive,mustBeReal} = 1000e-9;

        % analyze beam for only points in 3rd dimension that have at least this much power, relative to the peak power; set to 0 to analyze beam for all points in 3rd dimensions
        power_threshold (1,1) double {mustBeNonnegative,mustBeReal,mustBeLessThan(power_threshold, 1.0)} = 0;

        % empty property, we'll stick calculation results here
        results; 

        % for comparing slightly different algorithms to calculate Msquared (e.g. different 
        debugging (1,1) logical = false;
    end

    methods

        function varargout = calculate_pulse_flat(obj)
            %% Fluence-based beam analysis for pulsed waves
            %  This calculation uses the values stored in class properties xgrid, ygrid,
            %  field, and wavelength. The method is analogous to the instantaneous cw x&y
            %  calculation. See Chapter 7.7.3 of Crystal Nonlinear Optics: with SNLO
            %  examples (second edition). More details can be found at as-photonics.com/book

            assert(size(obj.field,3) ~= 1, 'field must have 3 dimensions longer than 1 (x,y,t)');

            Exyt = obj.field; % 3 dimensional array of electric field; first dimension is x, second is y, third is time

            % xgrid and ygrid must be either a vector or a matrix (increasing in one
            % dimension and constant in the other; for xgrid, second dimension is
            % constant, and for ygrid first is constant.)

            if isvector(obj.xgrid) || isvector(obj.ygrid)
                % if xgrid and ygrid are provided as vectors, make meshgrids out of them:
                [ymat,xmat]=meshgrid(obj.ygrid,obj.xgrid);
            end
            % obj.xgrid and obj.ygrid are not vectors; let's make sure they are oriented
            % correctly:
            if (max(abs(obj.xgrid-circshift(obj.xgrid,[0,1])),[], "all") == 0) && (max(abs(obj.ygrid-circshift(obj.ygrid,[1,0])),[], "all") == 0)
                % xgrid is constant in 2nd dimension, and ygrid in 1st
                xmat = obj.xgrid;
                ymat = obj.ygrid;
            elseif (max(abs(obj.ygrid-circshift(obj.ygrid,[0,1])),[], "all") == 0) && (max(abs(obj.xgrid-circshift(obj.xgrid,[1,0])),[], "all") == 0)
                % xgrid is constant in 1st dimension, and ygrid constant in 2nd; we want
                % the opposite, so swap them:
                xmat = obj.ygrid;
                ymat = obj.xgrid;
            else
                keyboard;
            end

            k0              = 2*pi/obj.wavelength;

            Nx              = size(Exyt, 1);
            Ny              = size(Exyt, 2);
            Nt              = size(Exyt, 3);
            dx              = xmat(2,2) - xmat(1,1);
            dy              = ymat(2,2) - ymat(1,1);

            delta_k_x       = 2*pi/dx/Nx;
            delta_k_y       = 2*pi/dy/Ny;
            if mod(Nx, 2)
            	kx_vec 		= ((0:(Nx-1))-(Nx-1)*0.5)*delta_k_x;
            else
            	kx_vec 		= ((0:(Nx-1))-Nx*0.5)*delta_k_x;
            end
			if mod(Ny, 2)
            	ky_vec 		= ((0:(Ny-1))-(Ny-1)*0.5)*delta_k_y;
            else
            	ky_vec 		= ((0:(Ny-1))-Ny*0.5)*delta_k_y;
            end
            [kymat, kxmat] 	= meshgrid(ky_vec, kx_vec);
            % [kymat, kxmat]  = meshgrid(((0:(Ny-1))-Ny*0.5)*delta_k_y, ((0:(Nx-1))-Nx*0.5)*delta_k_x);

            xmat3           = repmat(xmat,[1,1,Nt]);
            ymat3           = repmat(ymat,[1,1,Nt]);

            Exyt            = obj.field;

            % Repeatedly remove curvature, which improves accuracy of subsequent calculation
            N_removals      = 3;
            U               = zeros(N_removals,1);
            Uk              = zeros(N_removals,1);
            xBar_fl         = zeros(N_removals,1);
            yBar_fl         = zeros(N_removals,1);
            kxBar_fl        = zeros(N_removals,1);
            kyBar_fl        = zeros(N_removals,1);
            wxsquared_fl    = zeros(N_removals,1);
            wysquared_fl    = zeros(N_removals,1);
            wkxsquared_fl   = zeros(N_removals,1);
            wkysquared_fl   = zeros(N_removals,1);
            Ax              = zeros(N_removals,1);
            Ay              = zeros(N_removals,1);
            Rx              = zeros(N_removals,1);
            Ry              = zeros(N_removals,1);

            for P = 1:N_removals
                E_kxkyt         = fftshift(fftshift( fft(fft(Exyt, [], 1), [], 2), 1), 2);
                U(P)            = sum(abs(Exyt   ).^2, 'all'); % Denominator for wx/wy calculation (proportional to the energy, just left off dx*dy*dt*0.5*epsilon_0*c factor)
                Uk(P)           = sum(abs(E_kxkyt).^2, 'all');

                % calculate fluences (with arbitrary units)
                Fxy             = sum(abs(Exyt   ).^2, 3);  % proportional to the fluences in W/m^2, just missing factor of dt*0.5*epsilon_0*c)
                Fkxky           = sum(abs(E_kxkyt).^2, 3);

                % fluence-based xBar, yBar, kxBar, kyBar (centroids in xy-space and kxky-space) - see Equations 7.112 and 7.113
                xBar_fl(P)         = sum(xmat.*Fxy, 'all') ./ U(P);
                yBar_fl(P)         = sum(ymat.*Fxy, 'all') ./ U(P);
                kxBar_fl(P)        = sum(kxmat.*Fkxky, 'all') ./ Uk(P);
                kyBar_fl(P)        = sum(kymat.*Fkxky, 'all') ./ Uk(P);

                % Second moment widths in xy space and kxky space - see Equations 7.115 and 7.116
                wxsquared_fl(P)     = 4*sum( (xmat-xBar_fl(P)  ).^2 .* Fxy,   "all" )./U(P);     % scalar; should be unchanged when the phase corrections to remove curvature are added
                wysquared_fl(P)     = 4*sum( (ymat-yBar_fl(P)  ).^2 .* Fxy,   "all" )./U(P);     % scalar; should be unchanged when the phase corrections to remove curvature are added
                wkxsquared_fl(P)    = 4*sum( (kxmat-kxBar_fl(P)).^2 .* Fkxky, "all" )./Uk(P);    % scalar; 
                wkysquared_fl(P)   = 4*sum( (kymat-kyBar_fl(P)).^2 .* Fkxky, "all" )./Uk(P);    % scalar

                % R, A, M-squared, w0squared, and z0 in x direction

                % Find A - see Equation 7.118. Has a dE/dx term in it
                dE_dx               = (circshift(Exyt,[-1,0,0]) - circshift(Exyt,[1,0,0]))/(2*dx);  % 3d array; first derivative with respect to x (1st array dimension)
                % set end points to 0 because the circshift command rotates array and
                % wraps elements that are rotated off around to the other side; there
                % should be very little light at the edges anyway
                dE_dx([1,end],:,:)  = 0; 
                E_dE_dx_star        = Exyt.*conj(dE_dx);                                            % (3d array)

                Ax(P)               = 1i/k0 * sum(sum( (xmat-xBar_fl(P)).*sum(E_dE_dx_star - conj(E_dE_dx_star),3) )) ./ U(P); % (scalar) - Equation 7.118
                Rx(P)               = wxsquared_fl(P) ./ (2*Ax(P));                                 % generalized radius of curvature, scalar - Equation 74.117
                
                % Repeat everything for the y dimension
                % R, A, M-squared, w0squared, and z0 in y direction - repeat steps above, just swap to operating on the 2nd array dimension
                dE_dy               = (circshift(Exyt,[0,-1,0]) - circshift(Exyt,[0,1,0]))/(2*dy);  % (3d array)
                dE_dy(:,[1,end],:)  = 0;
                E_dE_dy_star    = Exyt.*conj(dE_dy);                                            % 3d
                Ay(P)           = 1i/k0 * sum(sum( (ymat-yBar_fl(P)).*sum(E_dE_dy_star - conj(E_dE_dy_star),3) )) ./ U(P); % scalar
                Ry(P)           = wysquared_fl(P) ./ (2*Ay(P));                                 % generalized radius of curvature, scalar

                if P ~= N_removals
                    % Generate phase arrays to remove curvature
                    xcurv   = k0*(xmat3-xBar_fl(P)).^2 ./ (2*Rx(P));
                    ycurv   = k0*(ymat3-yBar_fl(P)).^2 ./ (2*Ry(P));
                    xphase  = exp(-1i*xcurv);
                    yphase  = exp(-1i*ycurv);
                    Exyt    = Exyt.*xphase.*yphase; % remove the curvature
                end
            end

            % Harmonic sum of corrected terms
            Rx_corrected = 1./(sum(1./Rx));
            Ry_corrected = 1./(sum(1./Ry));

            wx0squared_fl_corrected = wxsquared_fl(1) - (1./wkxsquared_fl(1)).*(wxsquared_fl(1)*k0/Rx_corrected).^2; % Equation 7.119
            wy0squared_fl_corrected = wysquared_fl(1) - (1./wkysquared_fl(1)).*(wysquared_fl(1)*k0/Ry_corrected).^2; % Equation 7.119
            
            M2_x_corrected          = sqrt(wx0squared_fl_corrected)*sqrt(wkxsquared_fl(1))/2; % Equation 7.120
            M2_y_corrected          = sqrt(wy0squared_fl_corrected)*sqrt(wkysquared_fl(1))/2;  % Equation 7.120

            z0x                     = k0^2*wxsquared_fl./(Rx.*wkxsquared_fl); % distance of sample plane downstream from focus
            z0y                     = k0^2*wysquared_fl./(Ry.*wkysquared_fl); % distance of sample plane downstream from focus

            % Collect all the results as data struct fields in the 'results' class property
            obj.results.pulse_M2_x      = M2_x_corrected;
            obj.results.pulse_M2_y      = M2_y_corrected;
            obj.results.pulse_Rx        = Rx_corrected;
            obj.results.pulse_Ry        = Ry_corrected;
            obj.results.pulse_wx0       = sqrt(wx0squared_fl_corrected);
            obj.results.pulse_wy0       = sqrt(wy0squared_fl_corrected);
            obj.results.pulse_xBar      = xBar_fl(1);
            obj.results.pulse_yBar      = yBar_fl(1);
            obj.results.pulse_wx        = sqrt(wxsquared_fl(1));
            obj.results.pulse_wy        = sqrt(wysquared_fl(1));
            obj.results.pulse_kxBar     = kxBar_fl(1);
            obj.results.pulse_kyBar     = kyBar_fl(1);
            obj.results.pulse_wkx       = sqrt(wkxsquared_fl(1));
            obj.results.pulse_wky       = sqrt(wkysquared_fl(1));
            obj.results.pulse_z0x       = z0x;
            obj.results.pulse_z0y       = z0y;
            obj.results.pulse_flattened_Exyt  = Exyt;
            if nargout == 1
                varargout{1} = obj.results;
            else
                keyboard;
            end

        end

        function results = calculate(obj)
            % Calculate M-squared beam quality, as well as some other quantities like
            % centroid, second moment width, radius of curvature, angular width (second
            % moment in k-space), and the second moment at focus. M-squared is one half of
            % the sqrt of product of second moment at focus and angular width (Eq. 7.111
            % on page 302 of book). See Chapter 7.7.2 of Crystal Nonlinear Optics: with SNLO
            % examples (second edition). More details can be found at as-photonics.com/book
            %
            % The electric field array can be a 2d array or a 3d array and is provided to
            % this function as the as object property named 'field'. The 3rd dimension can
            % be time, if you are calculating M-squared vs t, or z, if you are calculating
            % M-squared through a crystal for instance. At each point in 3rd dimension,
            % calculate M-squared, which requires the beam size at focal waist, and the
            % second moment in k-space. We'll calculate the radius of curvature, beam size
            % at focal waist, distance to focal waist, spatial centroid (and spatial
            % angular centroid), spatial 2nd moment width (and spatial angular 2nd moment
            % width). Note that the M-squared calculator does not handle huge curvatures
            % that are non-spherical so we use a trick of repeatedly removing the
            % spherical curvatures

            k0 = 2*pi/obj.wavelength;
            if isvector(obj.xgrid) || isvector(obj.ygrid)
                % if xgrid and ygrid are provided as vectors, make meshgrids out of them:
                [ymat,xmat]=meshgrid(obj.ygrid,obj.xgrid);
            end

            % Make sure xgrid is constant in 2nd dimension and ygrid in 1st dimension
            if (max(abs(obj.xgrid-circshift(obj.xgrid,[0,1])),[], "all") == 0) && (max(abs(obj.ygrid-circshift(obj.ygrid,[1,0])),[], "all") == 0)
                xmat = obj.xgrid;
                ymat = obj.ygrid;
            elseif (max(abs(obj.ygrid-circshift(obj.ygrid,[0,1])),[], "all") == 0) && (max(abs(obj.xgrid-circshift(obj.xgrid,[1,0])),[], "all") == 0)
                xmat = obj.ygrid;
                ymat = obj.xgrid;
            else
                keyboard;
            end

            % Calculate spatial steps, create spatial frequency ("k-space") grids
            Nx                  = size(xmat,1);
            Ny                  = size(ymat,2);
            N3                  = size(obj.field,3);

            dx              = xmat(2,2) - xmat(1,1);
            dy              = ymat(2,2) - ymat(1,1);

            dkx             = 2*pi/dx/Nx;
            dky             = 2*pi/dy/Ny;
%             kxvec           = ((0:(nx-1))-0.5*nx)*dkx;
%             kyvec           = ((0:(ny-1))-0.5*ny)*dky;

            if mod(Nx, 2)
                kxvec 		= ((0:(Nx-1))-(Nx-1)*0.5)*dkx;
            else
                kxvec 		= ((0:(Nx-1))-Nx*0.5)*dkx;
            end
            if mod(Ny, 2)
                kyvec 		= ((0:(Ny-1))-(Ny-1)*0.5)*dky;
            else
                kyvec 		= ((0:(Ny-1))-Ny*0.5)*dky;
            end

                
            [kymat,kxmat]   = meshgrid(kyvec, kxvec);

                if obj.power_threshold ~= 0
                    norm_power_3rd = squeeze(sum(sum(abs(obj.field).^2, 1), 2)) ./ max(squeeze(sum(sum(abs(obj.field).^2, 1), 2)));
                    has_pwr = norm_power_3rd > obj.power_threshold;
                    t_start_ind = find(has_pwr, 1, 'first');
                    t_stop_ind  = find(has_pwr, 1, 'last');

                    if isempty(t_start_ind)
                        t_start_ind = 1;
                    end
                    if isempty(t_stop_ind)
                        t_stop_ind = size(obj.field,3);
                    end
                    if (t_start_ind~=1) || (t_stop_ind~=size(obj.field,3))
                        truncated = true;
                    else
                        truncated = false;
                    end
                    Exy3 = obj.field(:,:,t_start_ind:t_stop_ind);
                	N3              = size(Exy3,3);
                    original_nz = size(obj.field,3);
                    
                else
                    truncated = false;
                    Exy3 = obj.field;
                end

            % (Need to pre-allocate to store values at each index in 3rd dimension)
            xBar                = zeros(N3, 1);
            yBar                = zeros(N3, 1);
            wx                  = zeros(N3, 1);
            wy                  = zeros(N3, 1);
            kx_bar_original     = zeros(N3, 1);
            ky_bar_original     = zeros(N3, 1);
            wkx_original        = zeros(N3, 1);
            wky_original        = zeros(N3, 1);
            Rxmat               = zeros(N3, 3);
            Rymat               = zeros(N3, 3);
            Rx_corrected        = zeros(N3, 1);
            Ry_corrected        = zeros(N3, 1);
            wx0squared          = zeros(N3, 1);
            wy0squared          = zeros(N3, 1);
            M2_x                = zeros(N3, 1);
            M2_y                = zeros(N3, 1);
            Exyz_flat           = zeros(Nx, Ny, N3);

            for K = 1:N3 % Loop through third dimension (maybe time, or maybe z)
                    Exy = Exy3(:,:,K);
                    Ixy = abs(Exy).^2;
                    Ikxky   = abs( fftshift(fft2( Exy ))).^2;
                xBar(K)             = sum( xmat.*Ixy, "all" ) ./ sum(Ixy, "all");                                           % Centroid in x - see Equation 7.101
                yBar(K)             = sum( ymat.*Ixy, "all" ) ./ sum(Ixy, "all");                                           % Centroid in y - see Equation 7.101
                wx(K)               = sqrt( 4 * sum((xmat-xBar(K)).^2.*Ixy, "all") ./ sum(Ixy, "all") );                    % Second moment width in x - see Equation 7.100
                wy(K)               = sqrt( 4 * sum((ymat-yBar(K)).^2.*Ixy, "all") ./ sum(Ixy, "all") );                    % Second moment width in y - see Equation 7.100
                kx_bar_original(K)  = sum( kxmat.*Ikxky, "all") ./ sum(Ikxky, "all");                                       % Find centroid in kx (spatial frequency, the 2D fft of Exy) - see Equation 7.102
                ky_bar_original(K)  = sum( kymat.*Ikxky, "all") ./ sum(Ikxky, "all");                                       % Find centroid in ky (spatial frequency, the 2D fft of Exy) - see Equation 7.102
                wkx_original(K)     = sqrt( 4* sum((kxmat-kx_bar_original(K)).^2 .* Ikxky, "all") ./ sum(Ikxky, "all"));    % Find generalized angular width in x - see Equation 7.103
                wky_original(K)     = sqrt( 4* sum((kymat-ky_bar_original(K)).^2 .* Ikxky, "all") ./ sum(Ikxky, "all"));    % Find generalized angular width in y - see Equation 7.103

                % 'Remove curvature' three times, total curvature is harmonic sum,
                % with a big value with a couple of smaller correction terms.
                    for P = 1:3
                    % This is for Equation 7.108 - the first derivative of E with respect to x. Used in calculating A_x (used to define radius of curvature in x). 
                    % Derivative is using central differences.
                        dE_dx           = (circshift(Exy,[-1,0]) - circshift(Exy,[1,0])) / (2*dx);
                    % throw out the derivative at the end points because the circshift
                    % command will rotate to find the difference of the points on
                    % opposite ends, which is nonsense, and there should be
                    % essentially no light at the end points anyway.
                        dE_dx([1,end],:)= 0;
                        E_dE_dx_star    = Exy.*conj(dE_dx);
                    Ax              = 1i/k0 .* sum((xmat - xBar(K)) .* (E_dE_dx_star - conj(E_dE_dx_star)), "all") ./ sum(Ixy, "all");  % Equation 7.108
                    Rx              = wx(K).^2 ./ (2.*Ax);                                                                              % Equation 7.107                

                    % Repeat for y direction
                        dE_dy           = (circshift(Exy,[0,-1]) - circshift(Exy,[0,1])) / (2*dy);
                        dE_dy(:,[1,end])= 0;
                        E_dE_dy_star    = Exy.*conj(dE_dy);
                        Ay              = 1i/k0 .*  sum((ymat - yBar(K)) .* (E_dE_dy_star - conj(E_dE_dy_star)), "all") ./ sum(Ixy, "all");
                        Ry              = wy(K).^2 ./ (2.*Ay);

                        Rxmat(K,P)        = Rx;
                        Rymat(K,P)        = Ry;

                    % Generate the array of spherical phase front and apply opposite sign to electric field array
                        xcurv           = (xmat - xBar(K)).^2 ./ (2*Rx)*k0;
                        ycurv           = (ymat - yBar(K)).^2 ./ (2*Ry)*k0;
                        phase_xcurv     = exp(-1i.*xcurv);
                        phase_ycurv     = exp(-1i.*ycurv);
                    Exy             = Exy.*phase_xcurv.*phase_ycurv; % Remove the curvature
                end % end repeated removal of curvature

                % Harmonic sum of R which we'll use to calculate the beam size at focus
                    Rx_corrected(K) = 1./ (1/Rxmat(K,1) + 1/Rxmat(K,2) + 1/Rxmat(K,3));
                    Ry_corrected(K) = 1./ (1/Rymat(K,1) + 1/Rymat(K,2) + 1/Rymat(K,3));

                % Flattend electric field
                    Exyz_flat(:,:,K) = Exy;

                % Beam size at focal waist - Equation 7.110
                    wx0squared(K)  = wx(K).^2 - wx(K).^4.*k0^2 ./ (wkx_original(K).^2 .* Rx_corrected(K).^2);
                    wy0squared(K)  = wy(K).^2 - wy(K).^4.*k0^2 ./ (wky_original(K).^2 .* Ry_corrected(K).^2);

                % Beam quality - Equation 7.111
                    M2_x(K)        = sqrt(wx0squared(K)).*wkx_original(K)/2;
                    M2_y(K)        = sqrt(wy0squared(K)).*wky_original(K)/2;
            end % End loop through 3rd dimension of field

                if truncated
                % If this is a pulsed case, and there are many points at the start or end
                % of the array in 3rd dimension, replace the points where the is no (or
                % little) power by its nearest neighbor.
                    xBar            = cat(1, repmat(xBar(1),[t_start_ind-1,1]), xBar, repmat(xBar(end), [original_nz-t_stop_ind,1]));
                    yBar            = cat(1, repmat(yBar(1),[t_start_ind-1,1]), yBar, repmat(yBar(end), [original_nz-t_stop_ind,1]));
                    wx              = cat(1, repmat(wx(1),[t_start_ind-1,1]), wx, repmat(wx(end), [original_nz-t_stop_ind,1]));
                    wy              = cat(1, repmat(wy(1),[t_start_ind-1,1]), wy, repmat(wy(end), [original_nz-t_stop_ind,1]));
                    kx_bar_original = cat(1, repmat(kx_bar_original(1),[t_start_ind-1,1]), kx_bar_original, repmat(kx_bar_original(end), [original_nz-t_stop_ind,1]));
                    ky_bar_original = cat(1, repmat(ky_bar_original(1),[t_start_ind-1,1]), ky_bar_original, repmat(ky_bar_original(end), [original_nz-t_stop_ind,1]));
                    wkx_original    = cat(1, repmat(wkx_original(1),[t_start_ind-1,1]), wkx_original, repmat(wkx_original(end), [original_nz-t_stop_ind,1]));
                    wky_original    = cat(1, repmat(wky_original(1),[t_start_ind-1,1]), wky_original, repmat(wky_original(end), [original_nz-t_stop_ind,1]));
                    Rx_corrected    = cat(1, repmat(Rx_corrected(1),[t_start_ind-1,1]), Rx_corrected, repmat(Rx_corrected(end), [original_nz-t_stop_ind,1]));
                    Ry_corrected    = cat(1, repmat(Ry_corrected(1),[t_start_ind-1,1]), Ry_corrected, repmat(Ry_corrected(end), [original_nz-t_stop_ind,1]));
                    wx0squared      = cat(1, repmat(wx0squared(1),[t_start_ind-1,1]), wx0squared, repmat(wx0squared(end), [original_nz-t_stop_ind,1]));
                    wy0squared      = cat(1, repmat(wy0squared(1),[t_start_ind-1,1]), wy0squared, repmat(wy0squared(end), [original_nz-t_stop_ind,1]));
                    M2_x            = cat(1, repmat(M2_x(1),[t_start_ind-1,1]), M2_x, repmat(M2_x(end), [original_nz-t_stop_ind,1]));
                    M2_y            = cat(1, repmat(M2_y(1),[t_start_ind-1,1]), M2_y, repmat(M2_y(end), [original_nz-t_stop_ind,1]));

                end

%                z0x = k0^2 .* wx.^2 ./ (Rx_corrected.*wkx_original.^2);
                %z0y = k0^2 .* wy.^2 ./ (Ry_corrected.*wky_original.^2);

            % Collect into results class property
				z0x2 = -k0^2 .* wx.^2 ./ (Rx_corrected.*wkx_original.^2); % dist to focus
                z0y2 = -k0^2 .* wy.^2 ./ (Ry_corrected.*wky_original.^2); % dist to focus
                z0x  = pi.*sqrt(wx0squared)./(M2_x*obj.wavelength).*sqrt(wx.^2 - wx0squared); % dist to focus
                z0y  = pi.*sqrt(wy0squared)./(M2_y*obj.wavelength).*sqrt(wy.^2 - wy0squared); % dist to focus
                obj.results.z0x2        = z0x2; % 
                obj.results.z0y2        = z0y2; % 
				obj.results.z0x         = z0x; % dist to focus
                obj.results.z0y         = z0y; % dist to focus
                obj.results.M2_x        = M2_x;
                obj.results.M2_y        = M2_y;
                obj.results.Rx          = Rx_corrected;
                obj.results.Ry          = Ry_corrected;
                obj.results.wx0         = sqrt(wx0squared);
                obj.results.wy0         = sqrt(wy0squared);
                obj.results.xBar        = xBar;
                obj.results.yBar        = yBar;
                obj.results.wx          = wx;
                obj.results.wy          = wy;
                obj.results.kxBar       = kx_bar_original;
                obj.results.kyBar       = ky_bar_original;
                obj.results.wkx         = wkx_original;
                obj.results.wky         = wky_original;
                obj.results.flattened_Exyz = Exyz_flat;

            results = obj.results;

        end % end calculate method

        function w = calculate_second_moment(obj)
            % Find the second moment width of an electric field array; the calculation is
            % presented in Equations 7.100 and 7.101 of book Crystal Nonlinear Optics:
            % with SNLO examples (second edition) in chapter 7.7.2, page 300.
            arguments
                obj Msquared_class
            end
            if ~(isvector(obj.xgrid) || ismatrix(obj.xgrid))
%                 error('xgrid must be a vector or matrix; is has %i dimensions and size %s', ndims(obj.xgrid));
%                 sizestr = sprintf(join(repmat("%i ", ndims(obj.xgrid))), size(obj.xgrid));
                sizestr = sprintf(join(repmat("%i x", [1, ndims(obj.xgrid)-1]))+" %i", size(obj.xgrid));
                error('xgrid must be a vector or matrix; is has %i dimensions and size %s', ndims(obj.xgrid), sizestr);
            end
            if ~(isvector(obj.ygrid) || ismatrix(obj.ygrid))
                error('ygrid must be a vector or matrix; is has %i dimensions', ndims(obj.ygrid));
            end

            if ~isequal(ndims(obj.xgrid),ndims(obj.ygrid))
                error('xgrid and ygrid must have the same number of dimensions; xgrid has %i and ygrid has %i.', ndims(obj.xgrid), ndims(obj.ygrid));
            end

            if ismatrix(obj.xgrid)
                xmat = obj.xgrid;
                ymat = obj.ygrid;
            elseif isvector(obj.xgrid)
                [ymat, xmat] = meshgrid(obj.ygrid, obj.xgrid);
            end
            Nt = size(obj.field,3);
            xBar = zeros(Nt, 1);
            yBar = zeros(Nt, 1);
            wx   = zeros(Nt, 1);
            wy   = zeros(Nt, 1);
            % find centroid xBar, yBar
            for K = 1:Nt
                F = obj.field(:,:,K);
                I = abs(F).^2;
                denom = sum(sum(I));
                xBar(K) = sum(sum(xmat.*I)) ./ denom;
                yBar(K) = sum(sum(ymat.*I)) ./ denom;
                wx(K) = sqrt( 4 * sum(sum((xmat-xBar(K)).^2 .* I)) ./ denom );
                wy(K) = sqrt( 4 * sum(sum((ymat-yBar(K)).^2 .* I)) ./ denom );
            end
            w.wx = wx;
            w.wy = wy;
            % end "w = function calculate_second_moment";
        end

    end % end methods
end