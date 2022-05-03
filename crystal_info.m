% calculating crystal information (refractive index, group velocity index, group delay
% dispersion, third order dispersion, thermo-optic coefficient, and walkoff angles) using
% the functions ref_ind and N12 provided in mlSNLO as p-code files.

% 2022 May 3 by Jesse Smith (jesse.smith@as-photonics.com)

clear crystal
crystal.crystalname      = 'BBO';
crystal.wavelength       = 1000;
crystal.temperature      = 300;
crystal.prop_angle.theta = 80;
crystal.prop_angle.phi   = 20;

o = calc_crystal_info(crystal)

% o contains gvi, gdd, tod, dn_dT, n, walkoff angles, output_strings (a cell array of
% strings normally written to the output box in snlo_ref_ind_func)


function outputs = calc_crystal_info(inputs)
% input arguments: data structure 'inputs' with fieldnames crystalname, temperature, wavelength, and prop_angle.
% inputs.crystalname  must be a string matching one entry in the list of crystals (see crystals.txt)
% inputs.temperature  must be a scalar, units of kelvin
% inputs.wavelength   must be a scalar, units of nm
% inputs.prop_angle   must be a data structure
%    prop_angle must include fieldname prop_angle.theta for uniaxial crystals, 
%     and fieldnames prop_angle.theta and prop_angle.phi for biaxial crystals,
%     where theta is the 

% this function will do the following:
% read input variables:
% crystal name, temperature, wavelength, theta, phi, propagation axis
% convert wavelength to m units
%
% find refractive indices of crystal for wavelength and temperature
% compare it to result when using different temperature
%
% calculate group velocity index and group delay dispersion
% ( gvi is (d k)/(d omega) )
%  calculate omega (2*pi*c/wavelength)
%  set some omega increment d_omega
%  for omega-d_omega, omega, omega+d_omega
%    calculate ref inds for wavelengths corresponding to the three omega values above
% ( gdd is (d^2 k)/(d omega^2) )
% calculate group delay dispersion
% for omega-2*d_omega, omega+2*d_omega
%    calculate ref inds for wavelengths corresponding to the two omega values above

% remember ref_ind function behavior:
% ref_ind returns 0 if wavelength out of range for sellmeier
% ref_ind requires wavelengths to be in nanometers
% ref_ind returns 1 refractive index per wavelength for isotropic crystals, 2 for uniaxial crystals, and 3 for biaxial crystals

% constants
c = 3e8;                % speed of light (m/s)
epsilon_0 = 8.85e-12;   % vacuum permittivity (F/m)

% inputs data structure must include fieldnames
%  crystalname  (string corresponding to an entry in crystal_list.txt)
%  temperature  (scalar, in kelvin)
%  wavelength   (in nanometers)
%  prop_angle   (propagation angle, data structure with fieldnames theta and phi, scalars for angles in degrees)
%  axis_dir    (required for biaxial crystals; should be 'Z', 'X', or 'Y'
crystalname = inputs.crystalname;
temperature = inputs.temperature;
wavelength  = inputs.wavelength;
prop_angle  = inputs.prop_angle;

outputs.crystalname = crystalname;
outputs.temperature = temperature;
outputs.wavelength  = wavelength;
outputs.prop_angle  = prop_angle;

if isfield(prop_angle, 'phi')
    phi   = prop_angle.phi;
end

if isfield(prop_angle, 'theta')
    theta = prop_angle.theta;
end

% we try to keep all variables in mks units where where can
wavelength = wavelength*1e-9; % nm -> m

T_step = 25; % increment in temperature to use in determining whether crystal has a reliable thermo optic coefficient in ref_ind function

% calculate ref inds at temperature and at temperature + T_step:
% (note that ref_ind needs wavelength in nm)
ref_ind_T1 = ref_ind(crystalname, temperature, wavelength*1e9);
ref_ind_T2 = ref_ind(crystalname, temperature+T_step, wavelength*1e9);
temperature_dependent = any(ref_ind_T2 ~= ref_ind_T1);

% also note that for wavelengths out of range for sellmeier ref_ind returns 0
wavelength_out_of_range = any(ref_ind_T1 == 0);
if wavelength_out_of_range
    fprintf(1, 'wavelength out of range\n');
    return
end

% find group velocity index and group delay dispersion
% using their definitions (gvi = (d k) / (d omega)), (gdd = (d^2 k) / (d omega^2))
omega = 2*pi*c./wavelength; % rad / s?
% define increment d_omega for calculating discrete derivatives
d_omega     = 1e12; % rad / s
omega_p     = omega + d_omega;   % omega + 1*d_omega
omega_pp    = omega + 2*d_omega; % omega + 2*d_omega
omega_n     = omega - d_omega;   % omega - 1*d_omega
omega_nn    = omega - 2*d_omega; % omega - 2*d_omega

% convert to wavelengths
lambda_p    = 2*pi*c ./ omega_p;
lambda_pp   = 2*pi*c ./ omega_pp;
lambda_n    = 2*pi*c ./ omega_n;
lambda_nn   = 2*pi*c ./ omega_nn;

% find ref ind at wavelengths _p, _pp, _n, _nn, remembering to use nm units
n    = ref_ind(crystalname, temperature, wavelength*1e9);
rind_p  = ref_ind(crystalname, temperature, lambda_p*1e9);
rind_pp = ref_ind(crystalname, temperature, lambda_pp*1e9);
rind_n  = ref_ind(crystalname, temperature, lambda_n*1e9);
rind_nn = ref_ind(crystalname, temperature, lambda_nn*1e9);

% dn / dT
dT = 5;
rind_Tminus = ref_ind(crystalname, temperature-dT, wavelength*1e9);
rind_Tplus  = ref_ind(crystalname, temperature+dT, wavelength*1e9);


switch length(n)
    case 1
        % isoptropic crystal (ref_ind returns 1 refractive index for this crystal)
        k   = 2*pi*n     ./ wavelength;
        kp  = 2*pi*rind_p   ./ lambda_p;
        kpp = 2*pi*rind_pp  ./ lambda_pp;
        kn  = 2*pi*rind_n   ./ lambda_n;
        knn = 2*pi*rind_nn  ./ lambda_nn;

        % group velocity index: use finite difference for 1st derivative
        gvi = c * (kp - kn) ./ (d_omega*2);

        % group delay dispersion: use central 2nd order finite difference, for 2nd derivative
        gdd = c * (kp - kn - 2*k) / (d_omega.^2);  % units of s^2/m

        % third order dispersion: use central third order finite difference
        tod = (-0.5*knn + kn - kp + 0.5*kpp)./(d_omega.^3); % units of s^3 / m

        if temperature_dependent
            % first order central finite difference
            dn_dT = (rind_Tplus  - rind_Tminus) ./ (2*dT);
        else
            dn_dT = nan;
        end

        % strings
        angle_string = [];
        ref_ind_string  = sprintf('Ref. index       = %-8.5f',         n);
        gv_string       = sprintf('Grp vel ind      = c/%-8.5f',       gvi);
        gdd_string      = sprintf('GDD              = %-8.5G fs^2/mm', gdd*1e27); % 1e27 factor for s^2/m -> fs^2/mm
        tod_string      = sprintf('TOD              = %-8.5G fs^3/mm', tod*1e42);  % 1e42 factor for s^3/m -> fs^3/mm
        walkoff_string  = sprintf('Walkoff          = 0 mrad'); % its isotropic, no walkoff angle
        
        if temperature_dependent
            dn_dT = (rind_Tplus - rind_Tminus) / (2*dT);
            dn_dT_string = sprintf('dn/dT = %.3E /K', dn_dT);
            dn_dT_string = strrep(dn_dT_string, 'E-0', 'E-');
            dn_dT_string = strrep(dn_dT_string, 'E+', 'E');
            dn_dT_string = strrep(dn_dT_string, 'E0', 'E');
        end
        gdd_string = strrep(gdd_string, 'E-0', 'E-');
        gdd_string = strrep(gdd_string, 'E+', 'E');
        gdd_string = strrep(gdd_string, 'E0', 'E');
        tod_string = strrep(tod_string, 'E-0', 'E-');
        tod_string = strrep(tod_string, 'E+', 'E');
        tod_string = strrep(tod_string, 'E0', 'E');

        output_strings = {ref_ind_string; ...
            gv_string; ...
            gdd_string; ...
            tod_string; ...
            walkoff_string; ...
            dn_dT_string};

        outputs.gvi     = gvi;
        outputs.gdd     = gdd;
        outputs.tod     = tod;
        outputs.dn_dT   = dn_dT;
        outputs.n       = n;
        outputs.walkoff = 0;
        outputs.output_strings = output_strings;
    case 2
        % Uniaxial crystal (ref_ind returns 2 refractive indices for this crystal)
        angle_string = sprintf('Theta = %4.2f deg.',theta);

        % convert to radians
        theta = theta * pi/180;

        % e-polarization ref ind for propagation at angle theta
        rind = [n(1),...
            ((cos(theta)/n(1)).^2 + (sin(theta)./n(2)).^2).^(-0.5)];
        % and at slightly different angles (for use in calculating acceptance angles and walkoff)
        d_theta = 0.001;
        rind_a = [n(1),...
            ((cos(theta+d_theta)/n(1)).^2 + (sin(theta+d_theta)./n(2)).^2).^(-0.5)];
        rind_b = [n(1),...
            ((cos(theta-d_theta)/n(1)).^2 + (sin(theta-d_theta)./n(2)).^2).^(-0.5)];
        % and at slightly different wavelengths for use in calculating gv and gdd
        rind_p = [rind_p(1),...
            ((cos(theta)/rind_p(1)).^2 + (sin(theta)./rind_p(2)).^2).^(-0.5)];
        rind_n = [rind_n(1),...
            ((cos(theta)/rind_n(1)).^2 + (sin(theta)./rind_n(2)).^2).^(-0.5)];
        rind_pp = [rind_pp(1),...
            ((cos(theta)/rind_pp(1)).^2 + (sin(theta)./rind_pp(2)).^2).^(-0.5)];
        rind_nn = [rind_nn(1),...
            ((cos(theta)/rind_nn(1)).^2 + (sin(theta)./rind_nn(2)).^2).^(-0.5)];
        k   = 2*pi*rind    ./ wavelength;
        kp  = 2*pi*rind_p  ./ lambda_p;
        kpp = 2*pi*rind_pp ./ lambda_pp;
        kn  = 2*pi*rind_n  ./ lambda_n;
        knn = 2*pi*rind_nn ./ lambda_nn;
        % find group velocity index: (gvi = c * (d k)/(d omega))
        gvi = c * (kp - kn)./(d_omega*2);                       
        % find group delay dispersion: (gdd = (d^2 k)/(d omega^2))
        gdd  = ( kp + kn - 2*k ) / (d_omega^2);                 % s^2 / m units; use a factor of 1e27 to convert to fs^2/mm

        % find third order dispersion: (tod = (d^3 k)/(d omega^3));
        %  use -0.5*k_-2 + k_-1 - k_+1 + 0.5*k_+2
        tod = ( -0.5*knn + kn - kp + 0.5*kpp ) / (d_omega^3);   % s^3 / m units; use a factor of 1e42 to convert to fs^3/mm

        rho = [0,0.5E6*(rind_b(2)-rind_a(2))./rind(2)];         %walkoff angles in milliradians

        if temperature_dependent
            rind_Tplus = [rind_Tplus(1),...
                ((cos(theta)/rind_Tplus(1)).^2 + (sin(theta)./rind_Tplus(2)).^2).^(-0.5)];
            rind_Tminus = [rind_Tminus(1),...
                ((cos(theta)/rind_Tminus(1)).^2 + (sin(theta)./rind_Tminus(2)).^2).^(-0.5)];
            dn_dT = (rind_Tplus - rind_Tminus) / (2*dT);
        else
            dn_dT = nan;
        end

        ref_ind_string = sprintf('Ref. index (hi,lo)  = %-8.5f   %-8.5f',         rind);
        gv_string =      sprintf('Grp vel ind (hi,lo) = %-8.5f   %-8.5f',         gvi);
        gdd_string =     sprintf('GDD (hi,lo)         = %-8.5G   %-8.5G fs^2/mm', gdd*1e27);
        tod_string =     sprintf('TOD (hi,lo)         = %-8.5G   %-8.5G fs^3/mm', tod*1e42);
        walkoff_string = sprintf('Walkoff (hi,lo)     = %-8.3f   %-8.3f mrad',    abs(rho));
        dn_dT_string =   sprintf('dn/dT  (hi,lo)      = %-8.3E   %-8.3E /K', dn_dT);
        % remove extra +, or -0, or E0 from strings
        gdd_string = strrep(gdd_string, 'E-0', 'E-');
        gdd_string = strrep(gdd_string, 'E+',  'E');
        gdd_string = strrep(gdd_string, 'E0',  'E');
        tod_string = strrep(tod_string, 'E-0', 'E-');
        tod_string = strrep(tod_string, 'E+',  'E');
        tod_string = strrep(tod_string, 'E0',  'E');
        dn_dT_string = strrep(dn_dT_string, 'E-0', 'E-');
        dn_dT_string = strrep(dn_dT_string, 'E+', 'E');
        dn_dT_string = strrep(dn_dT_string, 'E0', 'E');

        output_strings = {angle_string;
            ref_ind_string;
            gv_string;
            gdd_string;
            tod_string;
            walkoff_string;
            dn_dT_string};

        outputs.gvi     = gvi;
        outputs.gdd     = gdd;
        outputs.tod     = tod;
        outputs.dn_dT   = dn_dT;
        outputs.n       = n;
        outputs.walkoff = rho;
        outputs.output_strings = output_strings;

    case 3
        % biaxial crystal (ref_ind returns 3 refractive indices for this crystal)
        theta = theta * pi/180;
        phi = phi * pi/180;
        if n(2)>n(3)
            % reorder n to be in increasing order
            n = n([1,3,2]);
            axis_dir = 'Y';
            
            % redefine angles to be relative to y axis?
            thetap = theta;
            phip = phi;
            theta = acos(sin(thetap).*sin(phip));
            phi = (-1.*sign(cos(thetap))).*acos( sin(thetap).*cos(phip)./sin(theta));
        else
            axis_dir = 'Z';
        end
        optic_axis = asin( sqrt(abs( (n(1).^(-2) - n(2).^(-2)) ./ (n(1).^(-2) - n(3).^(-2)) )) );
        optic_axis_string = sprintf('Optic axis          = %8.3f deg. from %s-axis', optic_axis*180/pi, axis_dir);

        % Compute delta (n rotate frame if optic axis is in xy-plane)
        delta = 0.5*atan( cos(theta).* sin(2*phi)./( sin(theta).^2./tan(optic_axis).^2  + sin(phi).^2 - cos(theta).^2 .* cos(phi).^2 ) );
        % delta is -45 to 45 deg. and is the angle to the nearest eigenpolarization, so must change it to be angle to hi
        % eigenpolarization
        signum = sign(sin(phi).*cos(phi).*cos(theta));
        if signum ~= sign(delta)
            delta = delta+signum*pi/2;
        end
        % delta doesn't always come out right for phi = 0 or
        % phi=180 deg, so fix that:
        if ((phi==pi)||(phi==0)) && (abs(theta - pi/2)<abs(optic_axis - pi/2))
            delta = 0;
        elseif ((phi==pi)||(phi==0)) && (abs(theta - pi/2)>abs(optic_axis - pi/2))
            delta = pi/2;
        end

        % Compute the walk off angles (in the rotated from if the
        % optic axis is in the xy-plane):
        % walk off doesn't come out right if phi=0 or 180 deg. so
        % add a small amount to phi in these cases:
        if (phi==0)||(phi==pi)
            phi = phi + 1E-6;
        end
        rho=[0;0];
        khat = [sin(theta).*cos(phi);...
            sin(theta).*sin(phi);...
            cos(theta)];
        rind = N12(n,[theta;phi].*180/pi,'Z');

        if (theta ~= pi/2)&&(theta ~= 0)
            Dprime =  sqrt(sum( (khat./(n.^(-2) - rind(1).^(-2)  ) ).^2 ));
            rho(1) = abs( 1E3*atan( rind(1).^2./Dprime ) );  % rho in milliradians
        end
        if (phi ~= pi/2) && (phi ~= 3*pi/2) && (theta ~= 0)
            Dprime =  sqrt(sum( (khat./(n.^(-2) - rind(2).^(-2)  ) ).^2 ));
            rho(2) = abs( 1E3*atan( rind(2).^2./Dprime ) );  % rho in milliradians
        end

        if axis_dir == 'Y' % return to original frame
            yp = -( cos(phi).*sin(delta) + cos(theta).*sin(phi).*cos(delta));
            delta = acos( -yp./sin(thetap) );
            if (delta>pi/2)
                delta = pi - delta;
            end
            signum = sin(phip).*cos(phip).*cos(thetap);
            delta = signum.*delta;
            if (phip == 0)
                delta = pi/2;
            end
            n=n([1,3,2]);
            theta = thetap;
            phi = phip;
        end
        angle_string = sprintf('theta,phi,delta     = %5.2f %5.2f %5.2f deg.',theta*180/pi,phi*180/pi,delta*180/pi);

        if (phi==0)||(phi==pi)
            phi = phi+1E-6;
        end

        rind    = N12(n,      [theta;phi] .*180/pi, 'Z');
        rind_p  = N12(rind_p, [theta;phi] .*180/pi, 'Z');
        rind_pp = N12(rind_pp,[theta;phi] .*180/pi, 'Z');
        rind_n  = N12(rind_n, [theta;phi] .*180/pi, 'Z');
        rind_nn = N12(rind_nn,[theta;phi] .*180/pi, 'Z');

        k   = 2*pi*rind    ./ wavelength;
        kp  = 2*pi*rind_p  ./ lambda_p;
        kpp = 2*pi*rind_pp ./ lambda_pp;
        kn  = 2*pi*rind_n  ./ lambda_n;
        knn = 2*pi*rind_nn ./ lambda_nn;
        gvi = c * (kp - kn)./(d_omega*2);

        % find group delay dispersion: (gdd = (d^2 k)/(d omega^2)) -- units s^2/m; a factor of 1e27 will change to fs^2/mm
        gdd  = ( kp + kn - 2*k ) ./ (d_omega^2); 

        % find third order dispersion: (tod = (d^3 k)/(d omega^3) -- units s^3/m; a factor of 1e42 will change to fs^3/mm
        tod = ( -0.5*knn + kn - kp + 0.5*kpp ) ./ (d_omega^3);

        if temperature_dependent
            % find dn / dT for both polarizations with propagation angles theta, phi
            rind_Tplus  = N12(rind_Tplus, [theta;phi].*180/pi,'Z');
            rind_Tminus = N12(rind_Tminus,[theta;phi].*180/pi,'Z');
            dn_dT = (rind_Tplus - rind_Tminus) ./ (2*dT);
            dn_dT_string = sprintf('dn/dT               = %-8.3E   %-8.3E /K', dn_dT);
            dn_dT_string = strrep(dn_dT_string, 'E-0', 'E-');
            dn_dT_string = strrep(dn_dT_string, 'E+', 'E');
            dn_dT_string = strrep(dn_dT_string, 'E0', 'E');
        else
            dn_dT = [0;0];
            dn_dT_string = sprintf('dn/dT               = %-8.3f   %-8.3f /K', 0, 0);
        end
        
        % These strings are written to the output box in snlo_ref_ind_func

        ref_ind_string = sprintf('Ref. index (hi,lo)  = %-8.5f   %-8.5f',       rind);
        gv_string =      sprintf('Grp vel ind (hi,lo) = %-8.5f   %-8.5f',       gvi);
        gdd_string =     sprintf('GDD (hi,lo)         = %-8.5G   %-8.5G fs^2/mm', gdd*1e27);
        tod_string =     sprintf('TOD (hi,lo)         = %-8.5G   %-8.5G fs^3/mm', tod*1e42);
        walkoff_string = sprintf('Walkoff (hi,lo)     = %-8.3f   %-8.3f mrad',  abs(rho));

        % remove extra +, or -0, or E0 from strings
        gdd_string = strrep(gdd_string, 'E-0', 'E-');
        gdd_string = strrep(gdd_string, 'E+',  'E');
        gdd_string = strrep(gdd_string, 'E0',  'E');
        tod_string = strrep(tod_string, 'E-0', 'E-');
        tod_string = strrep(tod_string, 'E+',  'E');
        tod_string = strrep(tod_string, 'E0',  'E');

        output_strings = {optic_axis_string; ...
            angle_string; ...
            ref_ind_string; ...
            gv_string; ...
            gdd_string; ...
            tod_string; ...
            walkoff_string;
            dn_dT_string};
        
        outputs.gvi     = gvi;
        outputs.gdd     = gdd;
        outputs.tod     = tod;
        outputs.dn_dT   = dn_dT;
        outputs.n       = rind;
        outputs.walkoff = rho;
        outputs.output_strings = output_strings;
%         fprintf(1,'%s',angle_string);
    otherwise
        keyboard;
end


end