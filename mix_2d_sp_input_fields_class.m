classdef mix_2d_sp_input_fields_class
    % Class to create object containing necessary data to pass snlo_2d_mix_sp_func for custom electric field input profiles. 
    
    % Required properties are field_red1_xyt, field_red2_xyt, field_blue_xyt, xgrid, and
    % ygrid; if cw, the field_ properties are 2-dimensional; if at least one field is
    % pulsed, the field_ properties are 3-dimensional and tgrid, upon which the fields are
    % defined, is also required
properties
    % 1- or 2-dimensional grid of spatial positions (in meters) on which the input electric field profiles are defined; centered about 0; if xgrid is a matrix, it must be increase in the first dimension and constant in the second
    xgrid (:,:) double {mustBeReal, mustBeFinite, mustBeNonempty} = linspace(-1,1,64);

    % 1- or 2-dimensional grid of spatial positions (in meters) on which the input electric field profiles are defined; centered about 0; if xgrid is a matrix, it must be increase in the second dimension and constant in the first
    ygrid (:,:) double {mustBeReal, mustBeFinite, mustBeNonempty} = linspace(-1,1,64); 

    % 2 or 3-dimensional grid of complex-valued electric field profiles (2d if cw, 3d if pulsed; 3rd dimension is time); Units are V/m, but everything is scaled to contain the energy stated in the field named mix_2d_sp_pulseenergy in the input data structure.
    field_red1_xyt (:,:,:) double {mustBeFinite, mustBeNonempty} = zeros(64,64,64);

    % 2 or 3-dimensional grid of complex-valued electric field profile (2d if cw, 3d if pulsed; 3rd dimension is time); Units are V/m, but everything is scaled to contain the energy stated in the field named mix_2d_sp_pulseenergy in the input data structure.
    field_red2_xyt (:,:,:) double {mustBeFinite, mustBeNonempty} = zeros(64,64,64);

    % 2 or 3-dimensional grid of complex-valued electric field profile (2d if cw, 3d if pulsed; 3rd dimension is time); Units are V/m, but everything is scaled to contain the energy stated in the field named mix_2d_sp_pulseenergy in the input data structure.
    field_blue_xyt (:,:,:) double {mustBeFinite, mustBeNonempty} = zeros(64,64,64);

    % 1-dimensional grid of time points (in seconds) on which 3rd dimension of field_red1_xyt, field_red2_xyt, and field_blue_xyt are defined (if pulsed)
    tgrid (1,:) double = linspace(-1.5, 1.5, 64)*1e-12;
end
% methods
% end
end