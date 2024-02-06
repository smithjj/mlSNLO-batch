classdef propagate_freespace_class
    % Definition file for class which propagates a two- or three-dimensional electric
    % field in free space. See properties and methods.
    properties
        % input_field (in arbitrary units), a real- or complex-valued electric field
        % distribution. Two- or three-dimensions, the first two of which are transverse
        % distribution, and the third (provided optionally) is time. The length in the
        % first two dimensions must match those of xgridmat and ygridmat.
        input_field (:,:,:) double {mustBeFinite} = zeros(64,64);

        % Spatial grid (in meters) in x direction, real-valued. One- or Two-dimensional
        % positions. If 2D, is meshgrid with nx * ny points, and changes in the first
        % dimension and is constant in the second.
        xgridmat (:,:) double {mustBeReal, mustBeFinite} = repmat(linspace(-1,1,64).'*1e-3, [1, 64]);

        % Spatial grid (in meters) in y direction, real-valued. One- or Two-dimensional
        % positions. If 2D, is meshgrid with nx * ny points, and changes in the second
        % dimension and changes in the second.
        ygridmat (:,:) double {mustBeReal, mustBeFinite} = repmat(linspace(-1,1,64)  *1e-3, [64, 1]);

        % Wavelength (in nanometers), scalar.
        wavelength (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 1000e-9;

    end
    methods
        function out = propagate_demag(obj, dist_to_propagate, demag)
            % This is a quick and sloppy way to do a propagation on a sort of adaptive
            % spatial grid. It is useful in showing the shape, really, so the amplitude
            % should be determined by you. Multiply input field by phase from a lens with
            % focal length of (distance to propagate)/((demagnify value) - 1), do free
            % space propagation by 2D fft, multiply by phase from lens of opposite focal
            % length and it's now defined on an x, y grid that is (demagnify value) times
            % larger. 
            %  Two input arguments:
            %  1 : distance to propagate (in air), specified in meters. Can be a vector, in
            %  which case this method propagates the field from the start to each position
            %  listed in the vector.
            %
            %  2: demag: a coefficient of how to expand grid (useful for propagating long
            %  distances where beam would expand off grid). Must be greater than or equal to 1.
            arguments
                obj propagate_freespace_class
                dist_to_propagate (1,:) double {mustBeReal, mustBeFinite}
                demag (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(demag,1), mustBeFinite}
            end

            obj = validate_grids(obj);

            Nx = size(obj.input_field, 1);
            Ny = size(obj.input_field, 2);
            Nt = size(obj.input_field, 3);
            dx = obj.xgridmat(2,2) - obj.xgridmat(1,1);
            dy = obj.ygridmat(2,2) - obj.ygridmat(1,1);

            if mod(Nx,2)
                Nxvec = cat(2, 0:(Nx/2), (-floor(Nx/2)) : -1);
                Nxsq_vec = Nxvec.^2;
                phasex  = repmat(pi*(Nx*dx).^(-2).*(obj.wavelength),[1,Nx]).*Nxsq_vec;
            else
                phasex = repmat(pi*(Nx*dx).^(-2).*circshift(((1:Nx)-0.5*Nx).^2,[0,-(0.5*Nx-1)]),[Ny,1]).' .* obj.wavelength;
            end
            if mod(Ny,2)
                % Odd # pts
                Nyvec = cat(2, 0:(Ny/2), (-floor(Ny/2)) : -1);
                Nysq_vec = Nyvec.^2;
                phasey = repmat(pi*(Ny*dx).^(-2).*(obj.wavelength),[1,Ny]).*Nysq_vec;
            else
                phasey = repmat(pi*(Ny*dy).^(-2).*circshift(((1:Ny)-0.5*Ny).^2,[0,-(0.5*Ny-1)]),[Ny,1]) .* obj.wavelength;
            end
            N = demag;

            phasexy = (repmat(phasex.',[1,Ny])+repmat(phasey,[Nx,1]));
            
            focal_length = dist_to_propagate./(N-1);
            if N == 1 && any(dist_to_propagate == 0)
                focal_length(isnan(focal_length)) = inf;
            end
            
            % Radius of curvature
            lens_roc = 2*focal_length;

            % actually propagate by distance of D
            D = dist_to_propagate./N;

            % Define an anomyous function to calculate the complex phase over a propagation distance
            prop_phase = @(distance) exp(-1i*distance*phasexy);

            % Define an anonymous function for phase array for passing through lens of focal length focal_length = D/(N-1)
            real_lens_curvature_func            = @(roc) -2*pi.*(obj.xgridmat.^2 + obj.ygridmat.^2)./(obj.wavelength.*roc);
            complex_phase_lens_curvature_func   = @(roc) exp(1i*real_lens_curvature_func(roc));

            % Define anonymous function propagate with a phase - take 2D fft, multiply by 2D phase, take inverse 2D fft
            propagateFcn = @(fieldarray,phasearray) ifft2(fft2(fieldarray).*phasearray); % propagate a field with a phase

            % pre allocate an array to stick propagated fields into
            outfield = zeros(Nx, Ny, Nt, length(dist_to_propagate));
            for K = 1:length(dist_to_propagate) % Loop through each distance to propagate
                for J = 1:Nt % Loop through length of third dimension of input_field, probably time
                    % field_lensed      = exp(1i*real_lens_curvature_func(lens_roc(K))).*obj.input_field(:,:,J);
                    field_lensed      = complex_phase_lens_curvature_func(lens_roc(K)).*obj.input_field(:,:,J);
                    field_prop_lensed = propagateFcn(field_lensed, prop_phase(D(K)))/N;
                    outfield(:,:,J,K) = exp(-1i*real_lens_curvature_func(lens_roc(K))).*field_prop_lensed;
                end
            end
            
            out.field = outfield;
            out.xmat = obj.xgridmat*N;
            out.ymat = obj.ygridmat*N;
        end

        function outfield = propagate(obj, dist_to_propagate)
            % Input arguments:
            %  - distance to propagate (in air), specified in meters. Can be a vector, in
            %  which case this method propagates the field from the start to each position
            %  listed in the vector.
            arguments
                obj propagate_freespace_class
                dist_to_propagate (1,:) double {mustBeReal, mustBeNonnegative, mustBeFinite}
            end

            obj = validate_grids(obj);

            Nx = size(obj.input_field, 1);
            Ny = size(obj.input_field, 2);
            Nt = size(obj.input_field, 3);
            dx = obj.xgridmat(2,2) - obj.xgridmat(1,1);
            dy = obj.ygridmat(2,2) - obj.ygridmat(1,1);

            % The step size should also be Length / (Num_pts - 1), 
            if mod(Nx,2)
                % If #pts is odd value
                Nxvec = cat(2, 0:(Nx/2), (-floor(Nx/2)) : -1);
                Nxsq_vec = Nxvec.^2;
                phasex  = repmat(pi*(Nx*dx).^(-2).*(obj.wavelength),[1,Nx]).*Nxsq_vec;
            else % # pts is even value; easy-peasy:
            phasex  = repmat(pi*(Nx*dx).^(-2).*(obj.wavelength),[1,Nx]).*circshift(((1:Nx)-0.5*Nx).^2,[0,-(0.5*Nx-1)]);
            end

            if mod(Ny,2)
                % If #pts is odd value
                Nyvec = cat(2, 0:(Ny/2), (-floor(Ny/2)) : -1);
                Nysq_vec = Nyvec.^2;
                phasey  = repmat(pi*(Ny*dy).^(-2).*(obj.wavelength),[1,Ny]).*Nysq_vec;
            else
            phasey  = repmat(pi*(Ny*dy).^(-2).*(obj.wavelength),[1,Ny]).*circshift(((1:Ny)-0.5*Ny).^2,[0,-(0.5*Ny-1)]);
            end

            phasexy = (repmat(phasex.',[1,Ny])+repmat(phasey,[Nx,1]));
           
            % Define an anonymous function to calculate a 2d array of propagation phase
            % from the matrices we just made and the product of the propagation distance
            % (assumed to have refractive index of 1.0).
            prop_phase_func = @(distance) exp(1i*-distance*phasexy);

            % Define an anonymous function to propagate a field using phase array so we
            % have a shorter way of executing this repeatedly.
            propagate_func  = @(fieldmat, phasemat) ifft2(fft2(fieldmat).*phasemat);

            % pre allocate an array to stick propagated fields into
            outfield = zeros(Nx, Ny, Nt, length(dist_to_propagate));
            for K = 1:length(dist_to_propagate)
                for J = 1:Nt
                    outfield(:,:,J,K) = propagate_func(obj.input_field(:,:,J), prop_phase_func(dist_to_propagate(K)));
                end
            end
        end
    end
    
    methods (Access = protected)
        function obj = validate_grids(obj)
            % Given input of "obj", determine whether the properties "xgridmat", "ygridmat", and "input_field" are compatible.
            % The shape of xgridmat and ygridmat must match, and also matche the first two dimensions of input_field.
            % input_field should be a two- or three-dimensional array, sized nx*ny or nx*ny*nt.
            if isvector(obj.xgridmat) && isvector(obj.ygridmat)
                % Have provided only vectors for x and y positions
                [ygm, xgm] = meshgrid(obj.ygridmat, obj.xgridmat);

                if size(obj.input_field,3) ~= 1 % Nt > 1

                    if (size(xgm,1) ~= size(obj.input_field,1)) && (size(xgm,1) ~= size(obj.input_field,2))
                        obj.input_field = permute(obj.input_field,[2,1,3]);
                    elseif (size(xgm,1) ~= size(obj.input_field,1)) || (size(xgm,2) ~= size(obj.input_field,2))
                        shorten = @(chars, N) chars(1:(end-N));
                        xsizestr = shorten(sprintf('%i x ',size(xgm)),2);
                        ysizestr = shorten(sprintf('%i x ',size(ygm)),2);
                        fsizestr = shorten(sprintf('%i x ',size(obj.input_field),[1,2]),2);
                        error('Incompatible sizes of xgridmat, ygridmat, and input_field (sizes %s, %s, and %s, respectively)',...
                            xsizestr, ysizestr, fsizestr);
                    end
                else
                    if (size(xgm,1) ~= size(obj.input_field,1)) || (size(xgm,2) ~= size(obj.input_field,2))
                        shorten = @(chars, N) chars(1:(end-N));
                        xsizestr = shorten(sprintf('%i x ',size(xgm)),2);
                        ysizestr = shorten(sprintf('%i x ',size(ygm)),2);
                        fsizestr = shorten(sprintf('%i x ',size(obj.input_field),[1,2]),2);
                        error('Incompatible sizes of xgridmat, ygridmat, and input_field (sizes %s, %s, and %s, respectively)',...
                            xsizestr, ysizestr, fsizestr);
                    end
                end
                obj.xgridmat = xgm;
                obj.ygridmat = ygm;
            elseif ismatrix(obj.xgridmat) && ismatrix(obj.ygridmat) % Have provided 2D arrays for positions xgridmat and ygridmat
                % Make sure xgridmat, ygridmat same sizes
                if (size(obj.xgridmat,1) ~= size(obj.ygridmat,1)) || (size(obj.xgridmat,2) ~= size(obj.ygridmat,2))
                    shorten = @(chars, N) chars(1:(end-N));
                    xsizestr = shorten(sprintf('%i x ',size(obj.xgridmat)),2);
                    ysizestr = shorten(sprintf('%i x ',size(obj.ygridmat)),2);
                    error('Incompatible sizes of xgridmat and ygridmat (sizes %s and %s, respectively)',...
                        xsizestr, ysizestr);
                end
                % make sure xgridmat, input_field same sizes
                if (size(obj.xgridmat,1) ~= size(obj.input_field,1)) || (size(obj.xgridmat,2) ~= size(obj.input_field,2))
                    shorten = @(chars, N) chars(1:(end-N));
                    xsizestr = shorten(sprintf('%i x ',size(obj.xgridmat)),2);
                    ysizestr = shorten(sprintf('%i x ',size(obj.ygridmat)),2);
                    fsizestr = shorten(sprintf('%i x ',size(obj.input_field),[1,2]),2);
                    error('Incompatible sizes of xgridmat, ygridmat, and input_field (sizes %s, %s, and %s, respectively)',...
                        xsizestr, ysizestr, fsizestr);
                end
            else
                error('xgridmat and ygridmat must be vectors or matrices; xgridmat has %i dimensions, and ygridmat has %i dimensions', ndims(obj.xgridmat), ndims(obj.ygridmat));
            end
        end
    end
end
