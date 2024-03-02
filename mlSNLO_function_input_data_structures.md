
# Ref ind
*snlo_ref_ind_func.m*
Ref.Ind. computes the refractive indexes, group velocities, group delay dispersions, and birefringent walkoff angles for uniaxial or biaxial crystals with any specified propagation direction. This function is useful if you wish to manually calculate phase matching or group velocities for situations not covered by SNLO functions Qmix, Bmix, QPM, Opoangles, Ncpm, or GVM.  

Optionally, this MATLAB function can be called with an input argument to pass the model input parameters. The input argument must be a data structure with fieldnames as follows:
| Fieldname | Shape | Units | Description | 
| - | - | - | - |
| ref_ind_selected_crystal  | scalar | integer | Of |
| ref_ind_temperature | scalar | K | Temperature |
| ref_ind_theta | scalar | degrees | Propagation angle |
| ref_ind_phi | scalar | degrees | Propagation angle |
| ref_ind_wavelength | scalar | nm | Wavelength |


# Qmix

*snlo_qmix_func.m*
QMIX helps you quickly select the best crystal for your application from a list of over 60 crystals.  It calculates all the possible phase-matching orientations for the chosen crystal at the wavelengths specified, returning phase-matching angles, polarizations, refractive indexes, group velocities, group delay dispersion, effective nonlinearity, crystal tilt tolerance, acceptance angles, acceptance bandwidths, and acceptance temperatures. The three beams are assumed to have collinear propagation vectors. If the nonlinearity, deff, is zero, the other properties are not calculated. (If you want to find noncritical wavelengths, use function Ncpm.) 
Optionally, this MATLAB function can be called with an input argument to pass the model input parameters. The input argument must be a data structure with fieldnames as follows:
| Fieldname | Shape | Units | Description | 
| - | - | - | - |
| qmix_selected_crystal | scalar | | Select your crystal by which line of crystal_list.txt it appears on
| qmix_temperature | scalar | kelvin | Temperature to calculate refractive indices at
| qmix_wavelength_red1 | scalar | nm | Wavelength of red2 wave |
| qmix_wavelength_red2 | scalar | nm | Wavelength of Red2 wave |
| qmix_wavelength_blue | scalar | nm | Wavelength of Blue  wave|
| qmix_principal_plane | 2 element character vector | Should be 'XY', 'XZ', or 'YZ'
| qmix_type | 3 element character vector | | Should be 'Mix' or 'OPO'

# Bmix
*snlo_bmix_func.m*
Bmix calculates the properties of biaxial crystals for propagation in any direction including out of the principal planes (in contrast to Qmix where propagation in a principal plane is assumed).  The calculations are based on the d-tensor and Sellmeier listed in Qmix. Phase matching curves are computed along with the associated effective nonlinearity, walkoff, and acceptance angles.  The purple circle on the phase matching plot indicates the location of the optic axes in the x-z plane. Parallel propagation directions (k vectors) of the three interacting waves red1, red2, and blue are assumed. 
For passing a set of inputs as an input argument to this function, the argument must be a data structure with fields named as follows:
|Fieldname | Shape | Units | Description |
| -- | -- | -- | -- |
| bmix_selected_crystal | String | | String containing abbreviation of the nonlinear crystal selection. See crystal_listb.txt for list of choices. |
| bmix_temperature | scalar | K | Temperature of crystal for calculations. | 
| bmix_wavelengths | 3 element vector | nm | Wavelengths for red1, red2, blue waves. | 
| bmix_rotation_axis | character | Member of set X, Y, Z | Rotation axis. |
| bmix_type | string | member of set "Type 1", "Type 2", or "Type 3" | Mixing type.|

# Cavity II
*snlo_brewster_displacement_cavity_func.m*
Updated cavity function. Can handle Brewster angle crystals, plots the beam width and radius of curvature through the cavity, also displays how a small displacement or tilt added at some specified location alters the cavity mode.
For passing a set of inputs as an input argument to this function, the argument must be a data structure with fields named as follows:
|Fieldname | Shape | Units | Description |
| -- | -- | -- | -- |
| brewster_dist_cav_wavelength | scalar | nm | Wavelength | 
| brewster_dist_cav_crystal_length  | scalar | mm | Length of the crystal |
| brewster_dist_cav_ref_ind         | scalar | | Crystal refractive index |
| brewster_dist_cav_mirror_roc      | scalar | mm | Radius of curvature of mirrors |
| brewster_dist_cav_mirror_angle    | scalar | degrees | Mirror angle
| brewster_dist_cav_leg1            | scalar | mm | Leg 1 length. The distance between the left mirror and right curved mirrors. | 
| brewster_dist_cav_leg2            | scalar | mm | Leg 2 length. Not used for linear cavities; for ring cavities, this is the length of the leg that bypasses the crystal. | 
| brewster_dist_cav_type            | string | 'Ring cavity' or 'Linear cavity' | Type of cavity - must be either 'Ring cavity' or 'Linear cavity'
| brewster_dist_cav_brewster_crystal | scalar | logical | Whether the crystal is tilted at the Brewster angle | 
| brewster_dist_cav_propagate_disturbance | scalar | logical | Whether to model and plot a disturbance, which can be either an tilt or a displacement added at some location specified. | 
| brewster_dist_cav_tilt_or_disp    | string | 'Tilt (urad)' or 'Displacement' | 
| brewster_dist_cav_dist_z          |  scalar | mm | Location of the displacement or tilt, measured from the left input mirror. | 
| brewster_dist_cav_dist_size       | scalar | um or urad | Magnitude of the displacement or tilt|

# Focus
*snlo_focus_func.m*
Calculate the width and curvature of a Gaussian beam at a specified location. These values are needed for the input parameters of 2D-mix-SP, 2D-mix-LP, and OPG. 

For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:

| Fieldname | shape | units | description|
| - | - | - |-|
| focus_wavelength | scalar | nm | Wavelength | 
| focus_ref_ind | scalar | | Crystal refractive index |
| focus_waist | scalar | mm | 1/e^2 radius of beam at waist |
| focus_face_to_focus | scalar | mm | Distance between face of crystal and focus waist.
| focus_dist | scalar | mm | Distance from the beam waist to calculate properties |

# GVM
Group velocity mismatch. 
*snlo_gvm_func.m*
GVM calculates group velocity (mis)match for noncollinear phase matching. The paper "Group-velocity-matched three-wave mixing in birefringent crystals" Optics Letters vol. 26, page 719 (2001) describes the method of group velocity matching in birefringent crystals implemented in this function. This paper is available online at http://www.as-photonics.com/Publications.html. The calculations are based on the d-tensor and Sellmeier equations displayed in Qmix.

For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
| Fieldname | shape | units | description |
|-|-|-|-|
|   gvm_temperature | scalar | K | Temperature |
|   gvm_currentCrystalValue | scalar |  Integer | Which index in cell array of crystal names matches the one we've selected. The cell array of crystal names is loaded from *crystals.txt* with one cell element per line. |
| gvm_wavelengths | 3 element vector | nm | Wavelengths. 1 element must be 0.
| gvm_slant_angle | scalar | degrees | Slant angle |
| gvm_polarization | character array, length 3 | | Red1, red2, blue polarizations (ordinary or extraordinary); choose from 'ooe', 'oee', 'eoe', 'eeo,' 'eoo', and 'oeo'. | 
|  gvm_plane | character array, length 2 | | The plane containing the optic axis (z axis) of the crystal; choose from XY, XZ, YZ |

# NCPM
*snlo_ncpm_func.m*
NCPM (noncritical phase matching) helps you quickly search the list of crystals for one that can noncritically phase match your wavelengths. Noncricitcal phase match refers to collinear phase matching with no birefringent walk off. It is desirable for efficient mixing of weak, tightly focused beams, and for maintaining the widest acceptance angle for all three beams. For uniaxial crystals the propagation direction is normal to the optic axis. For biaxial crystals the propagation direction is along one of the principal axes, X, Y, or Z.

For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
|ncpm_temperature|scalar|K|Temperature|
|ncpm_currentCrystalPopupValue | scalar | integer | index |
| ncpm_type | 6 elements |  character vector | Must be 'Type 1' or 'Type 2'. |
| ncpm_propagation_axis | 1 element | character | 'X' or 'Y' or 'Z' |
| ncpm_wavelengths | 2 element vector | nm | Red and blue wavelength target. |

# Opoangles
*snlo_opoangles_func.m*
Calculates red1 and red2 wavelengths in opo or opa as a function of propagation in birefringent crystals assuming the blue pump wavelength is fixed. For biaxial crystals, propagation is in one of the principal planes of the crystal so the Plane radio buttons are activated. Select the plane of interest. In biaxial crystals the angle is theta, measured from the Z axis for propagation in the XZ and YZ planes.  For propagation in the XY plane the angle is phi, measured from the X axis toward the Y axis. Calculations are based on the d-tensor and Sellmeier listed in Qmix.
For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
| opoangles_lambda | (3 element vector) | nm | Wavelengths for for red1, red2, and blue waves| 
| opoangles_temperature | scalar | kelvin | Temperature for refractive indices for phasematching | 
| opoangles_currentCrystalPopupValue | scalar | integer | Which line number of crystals.txt matches your crystal
| opoangles_type   | character vector | | Should be 'Type 1' or 'Type 2'
| opoangles_plane | 2 element character vector | | Should be 'XY' or 'XZ' or 'YZ'
| opoangles_delta | scalar | degrees | Pump tilt: the angle between signal wave vector and pump wave vector.
| opoangles_signal_pol | scalar character | | Should be 'o' or 'e' |

# PW-cav-LP
*snlo_pw_cav_lp_func.m*
PW-cav-LP models plane-wave mixing in cavities.  Input irradiances are those of the central rays of lowest order Gaussian beams with diameters specified by the Beam diam. input values.  You can model an OPO, cavity resonant SHG or any other cavity mixing process by specifying the input beams and the reflectivity of the mirrors at those wavelengths.  The crystal mixing is similar to PW-mix-LP. The cavity configuration can be standing wave or ring.
For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
| pw_cav_lp_wavelengths | (3 element vector) | nm | Wavelengths for red1, red2, and blue waves | 
| pw_cav_lp_ref_inds | (3 element vector) | |  Refractive indices for red1, red2, and blue waves (at wavelengths specified in pw_cav_lp_wavelengths) | 
| pw_cav_lp_left_refl_c | (3 element vector) | |  Left crystal reflectivity; this light is lost, not resonated. | 
| pw_cav_lp_right_refl_c | (3 element vector) | |  Right crystal reflectivity; this light is lost, not resonated. | 
| pw_cav_lp_crystal_loss | (3 element vector) | 1/mm | Crystal linear absorption for each wave | 
| pw_cav_lp_left_energy_pwr | (3 element vector) | J or W | For red1, red2, and blue waves this is the left input pulse energy (if wave is pulsed), or power (if wave is cw) 
| pw_cav_lp_right_energy_pwr | (3 element vector) | J or W | For red1, red2, and blue waves this is the right input pulse energy (if wave is pulsed), or power (if wave is cw) 
| pw_cav_lp_pulse_duration | (3 element vector)| FWHM ns | Red1, red2, and blue wave pulse durations (FWHM) | 
| pw_cav_lp_beam_diameters | (3 element vector) | mm | Red1, red2, and blue beam diameters (FWHM irrad.) at crystal entrance face. Lowest order Gaussian is assumed, and we model the ray at the beam center. To convert from the RADIUS of (1/e)^2 irradiance, multiply by 1.18.
| pw_cav_lp_left_refl_m | (3 element vector) | | Red1, red2, and blue wave power reflectivities at the left mirror.
| pw_cav_lp_right_refl_m | (3 element vector | | Red1, red2, and blue wave power reflectivities at the right mirror.
| pw_cav_lp_phase_left | (3 element vector) | radians |  Additional phase to add to the beam between the left mirror and crystal input face.
| pw_cav_lp_phase_right | (3 element vector) | radians |  Additional phase to add to the beam between the crystal output face and right mirror.
| pw_cav_lp_third_leg_phase | (3 element vector) | radians | For ring cavities, an additiona phase to add to the beam between the right mirror and left mirror | 
| pw_cav_lp_cavity_length | (scalar) | mm | Total cavity length | 
| pw_cav_lp_crystal_length | (scalar) | mm | Crystal length
| pw_cav_lp_deff | (scalar) | pm/V | Effective nonlinear coefficient Deff | 
| pw_cav_lp_cavity_type | (scalar) | integer | Sets the cavity type to linear cavity (value of 1) or ring cavity (value of 0). | 
| pw_cav_lp_deltak | (scalar) | 1/mm | Phase mismatch Delta k = k_blue - k_red1 - k_red2 | 
| pw_cav_lp_nz | (scalar) integer |  | # integration steps through the crystal | 
| pw_cav_lp_n_starts | (scalar) integer | | Number of starts; use a value of 1 for fasted computation; larger values provide better temporal resolution. Effective time step is (cavity length) / ((speed of light) * (# starts)). |

# PW-mix-BB
*snlo_pw_mix_bb_func.m*
PW-mix-BB computes three-wave mixing processes for plane waves with input irradiances equal to those of the central ray of lowest order Gaussian beams with energy and diameter specified on the input form. It ignores birefringent walkoff but includes group velocity effects. It is intended for modeling multi-longitudinal-mode pulses.  The number of modes populated in each model run is sufficient to cover the specified bandwidth for each beam. Quantum noise is automatically added to the red1 and red2 waves. The noise bandwidths of the red waves can be extended by specifying the bandwidth parameter. This inclusion of quantum noise is to approximate OPG modeling. However, to properly model OPG the spatial structure of the red beams should be noisy as well. In most other applications it will be invisible. The noise spectrum can be seen by setting d_eff to zero and setting the red energies or powers to zero, running and clicking
'Spectra'.
For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
| pw_mix_bb_wavelengths  		| (3 element vector) | nm |  Red1, red2, blue wave wavelengths |
| pw_mix_bb_ref_inds     		| (3 element vector) | |  Refractive indices | 
| pw_mix_bb_gvi          		| (3 element vector) | | Group velocity indices | 
| pw_mix_bb_gdd          		| (3 element vector) | fs^2/mm | gdd group delay dispersion | 
| pw_mix_bb_input_refl   		| (3 element vector) | | Crystal input power reflectivity| 
| pw_mix_bb_output_refl   		| (3 element vector) | | Crystal output power reflectivity |
| pw_mix_bb_crystal_loss       	| (3 element vector) | 1/mm | Linear absorption |
| pw_mix_bb_pulse_energy_power  | (3 element vector) | J or W | Input energy (for pulsed waves) or power (for cw waves) |
| pw_mix_bb_pulse_duration      | (3 element vector) | FWHM ns | Input pulse durations | 
| pw_mix_bb_beam_diameter       | (3 element vector) | FWHM mm | Input beam diameters at crystal input face. Lowest order Gaussian is assumed, and we model the ray at the beam center. To convert from the RADIUS of (1/e)^2 irradiance, multiply by 1.18.
| pw_mix_bb_bandwidth           | (3 element vector) | FWHM GHz | Bandwidth
| pw_mix_bb_mode_spacing	    | (3 element vector) | GHz | Mode spacing 
| pw_mix_bb_freq_mod            | (3 element vector, boolean) | | True: Frequency modulated; false: amplitude noise also.
| pw_mix_bb_crystal_length      | (scalar) | mm | Crystal length |
| pw_mix_bb_deff                | (scalar) | pm/V | Effective nonlinear coefficient deff
| pw_mix_bb_deltak              | (scalar) | 1/mm | Phase mismatch Delta k = (k_3-k_1-k_2)
| pw_mix_bb_nz | (scalar) | | Nnumber of points in z through the crystal for integration. | 



# PW-mix-BB
*snlo_pw_mix_bb_func.m*
PW-mix-LP computes three-wave mixing process for long-pulse plane- waves. The input irradiances are those of the central rays of spatial lowest-order Gaussian beams with diameters and energies specified in the input form.  This function ignores birefringent and group velocity walkoff.  It is intended for very quick evaluation of a proposed mixing process. For help in setting input values, right- click on the input edit box. Help text appears in the lower text box.
For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
pw_mix_lp_wavelengths | (3 element vector) | nm | Wavelengths
|   pw_mix_lp_ref_inds | (3 element vector) | | Refractive indices
|   pw_mix_lp_phases | (3 element vector) | radians | Input phases | 
|   pw_mix_lp_input_reflectivities | (3 element vector)  | | Power reflectivity at crystal input face
|   pw_mix_lp_output_reflectivities| (3 element vector) | | Power reflectivity at crystal output face
|   pw_mix_lp_crystal_losses 	| (3 element vector) | 1/mm | Crystal linear absorption | 
|   pw_mix_lp_energy_power 		| (3 element vector) | J or W | Input energy (if wave is pulsed) or power (if wave is cw)
|   pw_mix_lp_pulse_durations 	| (3 element vector) | fwhm ns | Pulse durations | 
|   pw_mix_lp_beam_diameters 	| (3 element vector) | fwhm mm | Beam diameters. Lowest order Gaussian is assumed, and we model the ray at the beam center. To convert from the RADIUS of (1/e)^2 irradiance, multiply by 1.18. |
|   pw_mix_lp_crystal_length 	| (scalar) | mm | crystal length |
|   pw_mix_lp_deff 				| (scalar) | pm/V | Effective nonlinear coefficient deff | 
|   pw_mix_lp_deltak | (scalar) | 1/mm | Phase mismatch Delta k = (k_3 - k_1 - k_2) | 
|   pw_mix_lp_nz 		| scalar integer | Number of integration z points through the crystal | 
|   pw_mix_lp_nt | (scalar integer | | Number of time steps (ignored if no waves are pulsed). |

The input values that have 3 elements are in order of red1, red2, and blue waves where blue has the shortest wavelength and the other two are interchangeable.| 

# PW-mix-SP
*snlo_pw_mix_sp_func.m*
PW-mix-SP computes three-wave mixing process for plane waves with irradiances equal to those of the central rays of spatial lowest order Gaussian beams with beam diameters and energies specified in the input form.  It ignores birefringent walkoff but includes group velocity effects, so it is useful for modeling very short pulses. cw light is not allowed in this function.  Pulses are Gaussian or Supergaussian in time at the input face. The three input pulses can have different durations and delays. For help in setting input values, right-click on the input edit box. Help text appears in the lower text box.
For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
|   pw_mix_sp_wavelengths | 3 element vector | nm | Wavelengths for red1, red2, blue waves.
|   pw_mix_sp_ref_inds | 		(3x1 array) | | ref ind |
|   pw_mix_sp_gvi |  		    (3x1 array) | | group velocity index |
|   pw_mix_sp_gdd | 			(3x1 array) | (fs^2/mm) |group delay dispersion gdd |
|   pw_mix_sp_phases | 			(3x1 array) | rad | Input phases |
|   pw_mix_sp_n2_red1| 			(3x1 array) | cm^2/W | nonlinear refractive index n2 (red1) |
|   pw_mix_sp_n2_red2| 			(3x1 array) nonlinear refractive index n2 (red2) (sq cm/W)
|   pw_mix_sp_n2_blue| 			(3x1 array) nonlinear refractive index n2 (blue) (sq cm/W)
|   pw_mix_sp_beta_red1 | 		(3x1 array) 2-photon absorption beta (red1) (cm/W)
|   pw_mix_sp_beta_red2 |		(3x1 array) 2-photon absorption beta (red2) (cm/W)
|   pw_mix_sp_beta_blue | 		(3x1 array) 2-photon absorption beta (blue) (cm/W)
|   pw_mix_sp_input_refl | 		(3x1 array) input power reflectivity
|   pw_mix_sp_output_refl |     (3x1 array) output power reflectivity
|   pw_mix_sp_crystal_losses|  	(3x1 array) crystal loss (1/mm)
|   pw_mix_sp_pulseenergy| 		(3x1 array) pulse energies (J)
|   pw_mix_sp_beam_diameter|	(3x1 array) beam diam (FWHM mm)
|   pw_mix_sp_pulse_durations| 	(3x2 array) pulse durations (FWHM ps), and temporal super gaussian index (0 = sech^2, 1 = Gaussian, >1 more tophat-like)
|   pw_mix_sp_pulse_chirp | 		(3x3 array) pulse chirp: red1, red2, blue in first dimension; linear, quadratic and cubic in second dimension (THz/ps, THz/ps^2, and THz/ps^3);
|   pw_mix_sp_pulse_delay | 		(2x1 array) pulse delay relative to center of blue pulse (ps)
|   pw_mix_sp_crystal_length | 	(scalar) crystal length (mm)
|   pw_mix_sp_deff |  			(scalar) effective nonlinear coefficient deff (pm/V)
|   pw_mix_sp_deltak | 			(scalar) phase mismatch (k_3 - k_1 - k_2) [1/mm]
|   pw_mix_sp_nz | 				(scalar) num z points 
|   pw_mix_sp_nt | 				(scalar) num time points

# PW-opo-BB
*snlo_pw_opo_bb_func.m*
PW-OPO-BB models broad-bandwidth, nanosecond, plane-wave (actually the central ray of a
 spatial Gaussian) OPO's.  It is similar to PW-cav-LP except it allows for differing
 group velocities for the three waves and broad bandwidths.  The mathematical methods are
 described in paper "Numerical models of broad bandwidth nanosecond optical parametric
 oscillators" JOSA B vol. 16 p. 609 (1999), available at
 http://www.as-photonics.com/Publications.html This function permits studies of unseeded
 OPO's and the transition to seeding. It also models double resonance effects and pumping
 by multi-mode lasers. I suggest you narrow the range of input parameters using the much
 faster function PW-cav-LP before you run this function because it is demanding of
 computer resources and user patience. For help in setting input values, right-click on
 the input edit box. Help text appears in the lower text box.
For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
| pw_opo_bb_wavelengths  |  3x1 array | nm | Wavelengths 
| pw_opo_bb_gvi |  3x1 array | | Group velocity indices | 
| pw_opo_bb_ref_inds |  3x1 array | |        Refractive indices | 
| pw_opo_bb_left_refl_c | 3x1 array |  | Left crystal reflectivities
| pw_opo_bb_right_refl_c |  3x1 array | | Right crystal reflectivity
| pw_opo_bb_crystal_loss | 3x1 array | 1/mm | Crystal linear absorption | 
| pw_opo_bb_left_energy_pwr | 3x1 array | J or W | Left mirror input power (if wave is cw) or energy (if wave is pulsed)
| pw_opo_bb_right_energy_pwr | 3x1 array | J or W | Right mirror input power (if wave is cw) or energy (if wave is pulsed)
| pw_opo_bb_pulse_duration | 3x1 array  | FWHM ns | Pulse durations (0 for cw).
| pw_opo_bb_beam_diameters | 3x1 array  | FWHM mm | Beam diameters. Lowest order Gaussian is assumed, and we model the ray at the beam center. To convert from the RADIUS of (1/e)^2 irradiance, multiply by 1.18. |
| pw_opo_bb_left_refl_m | 3x1 array | | Left cavity mirror reflectivity
| pw_opo_bb_right_refl_m | 3x1 array | | Right cavity mirror reflectivity
| pw_opo_bb_phase_left | 3x1 | radians | Additional phase to apply between the left cavity mirror and crystal input face | 
| pw_opo_bb_phase_right | 3x1 | radians | Additional phase to apply between the crystal output face and right cavity mirror | 
| pw_opo_bb_cavity_length |  scalar | mm | Cavity length | 
| pw_opo_bb_crystal_length |  scalar | mm  | Length of crystal | 
| pw_opo_bb_pump_bandwidth |  scalar | MHz | Pump bandwidth
| pw_opo_bb_pump_mode_spacing |  scalar | MHz | Pump mode spacing
| pw_opo_bb_pump_fm | scalar logical |  |  Whether the pump is just frequency modulated (value of 1) or amplitude modulated as well (value of 0) | 
| pw_opo_bb_deff |  scalar | pm/V | Effective nonlinear coefficient d|eff | 
| pw_opo_bb_cavity_type |  scalar logical | | Cavity type: 1 for linear cavity, 0 for ring cavity
| pw_opo_bb_deltak |  scalar | 1/mm | Phase mismatch Deltak = (k_blue - k_red1 - k_red2)
| pw_opo_bb_nz | scalar integer | | Number of integration steps through crystal
| pw_opo_bb_new_noise | scalar integer | |  Use new noise (value of 1) or try to re-use previous run's noise (value of 0)

For inputs with 3x1 arrays, first element is for red1 wave, second is for red2 wave, and third is for blue wave. The blue wave is the wave with the shortest wavelength, and red1 and red2 are interchangeable.

# PW-opo-SP
*snlo_pw_opo_sp_func.m*
PW-OPO-SP models a singly-resonant, plane-wave, synchronously-pumped OPO.  The
irradiance of the plane pump wave is set equal to the central ray of a lowest-order
spatial Gaussian whose beam diameter is specified in the input form under 'Beam
diameter'.  This function includes group velocity and group delay dispersion but not
higher order dispersion. The OPO cavity is a ring, with only the red1 wave resonated.
The duration of the blue (pump) pulses are specified. The shape of these pulses can be
changed by including an optional integer 0-10 following the duration, separated by a
space. For integer value of 0, the pulse is a hyperbolic secant squared; a value of 1 is
a typical Gaussian; values of 2+ are super Gaussians which increasingly resemble top-hat
profiles with increasing value. For help in setting input values, right-click on the
input edit box. Help text appears in the lower text box.

For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
| pw_opo_sp_wavelengths | 3x1 array | nm | Wavelengths | 
| pw_opo_sp_ref_inds | 3x1 array | | Refractive indices |
| pw_opo_sp_gvi | 3x1 array | | Group velocity indices  | 
| pw_opo_sp_gdd  | 3x1 array | fs^2/m | Group delay dispersion | 
| pw_opo_sp_crystal_loss | 3x1 array | 1/mm | Crystal linear absorption | 
| pw_opo_sp_n2_red1 | 3x1 array | cm^2/W | n2 from red1 light |
| pw_opo_sp_n2_red2 | 3x1 array | cm^2/W | n2 from red2 light |
| pw_opo_sp_n2_blue | 3x1 array | cm^2/W | n2 from blue light | 
| pw_opo_sp_beta_red1 | 3x1 array | cm/W | Two-photon absorption beta from red1 light  | 
| pw_opo_sp_beta_red2 | 3x1 array | cm/W | Two-photon absorption beta from red2 light  | 
| pw_opo_sp_beta_blue | 3x1 array | cm/W | Two-photon absorption beta from blue light  | 
| pw_opo_sp_pulse_energy | scalar | uJ | Input pulse energy | 
| pw_opo_sp_pulse_duration | scalar or 2x1 array | FWHM ps | Pulse duration (FWHM ps), optionally provide a second element for temporal Supergaussian index (value of 1 is Gaussian, larger values are more flat-topped) | 
| pw_opo_sp_chirp | scalar | THz/ps | Chirp of blue wave|
| pw_opo_sp_beam_diameter | 3x1 array | FWHM mm | Beam diameters. Lowest order Gaussian is assumed, and we model the ray at the beam center. To convert from the RADIUS of (1/e)^2 irradiance, multiply by 1.18. |
| pw_opo_sp_output_refl | scalar |  | Output coupler reflectivity for blue wave | 
| pw_opo_sp_loss | scalar | 1/mm | Other cavity loss | 
| pw_opo_sp_cavity_delay | scalar | ps | Cavity delay
| pw_opo_sp_cavity_gdd  | scalar | fs^2 | Cavity group delay dispersion | 
| pw_opo_sp_lcrystal | scalar | mm | Crystal length | 
| pw_opo_sp_deff  | scalar | pm/V | Effective nonlinear coefficient deff | 
| pw_opo_sp_deltak | scalar | 1/mm | Phase mismatch Delta k = k_blue - k_red1 - k_red2 | 
| pw_opo_sp_nz | scalar integer | | Number of z points through crystal for integration | 
| pw_opo_sp_nt | scalar, integer | | Number of time points | 
For inputs with 3x1 arrays, first element is for red1 wave, second is for red2 wave, and third is for blue wave. The blue wave is the wave with the shortest wavelength, and red1 and red2 are interchangeable.


# QPM
*snlo_qpm_func.m*
QPM calculates properties of quasiphasematched materials including the poling, or domain
reversal, period (always specified at room temperature), the temperature bandwidth, the
frequency bandwidth, the group velocities, the group delay dispersions.  In addition the
tuning of phase matched red1 and red2 wavelengths due to blue (pump) tuning or crystal
temperature tuning can be computed using the Pump Tune and Temp Tune buttons. The
Sellmeier equations are those referenced in Qmix when a crystal is selected. The
temperature range and acceptance bandwidths are calculated as described in Qmix help.

For passing a set of inputs as an input argument to this function, the following fields must be present in that data structure:
|Fieldname|Shape|Units|Description|
|-|-|-|-|
| qpm_temperature                | scalar | kelvin | Crystal temperature | 
| qpm_wavelengths                | 3x1 array | nm | Wavelengths of red1, red2, blue waves. | 
| qpm_currentCrystalPopupValue   | scalar integer |  | Set this to which line of crystals.txt your crystal appears on | 
| qpm_currentPolValue            | 
| qpm_temptune_temp_range        |  
| qpm_temptune_period            | 
| qpm_pumptune_wavelength_range  | 
| qpm_pumptune_period            | 

# Thermal
*snlo_thermal_func.m*
Simple 3D steady-state thermal solver. Cuboid volume temperature, volume heating. Input and output faces are insulated. The side faces can each be perfectly insulated or perfectly cooled with temperature of 0.
Uses spectral methods, cosine and sine transforms.

# 2D-cav-LP
Diffractive, long-pulse intracavity mixing 
*snlo_2d_cav_lp_func.m*

2D-cav-LP is similar to PW-cav-LP except it handles the full spatial profiles so it can realistically model beams with Gaussian or Supergaussian spatial profiles. It includes diffraction, birefringent walkoff, displaced beams, etc. for pulsed or cw light. The light is assumed to be monochromatic and the pulses are assumed to be long enough that group velocity effects are unimportant. The temporal shape of the pulse can be altered by including the super Gaussian index as an optional second argument in the duration box, separating the two arguments with a space; this second argument should be between 1 and 20, where a value of 1 is the typical Gaussian shape, and larger values are more flat-topped. For help in setting input values, right-click on the input edit box. Help text appears in the lower text box.

For passing a set of inputs as an input argument to this function, the argument must be a data structure with fields named as follows:
|Fieldname | Shape | Units | Description |
| -- | -- | -- | -- |
| cav_2d_lp_wavelengths    			| 3 element vector | nm | Wavelengths for red1, red2, and blue waves |
| cav_2d_lp_ref_inds 				| 3 element vector |    | Refractive indices for red1, red2, and blue waves |
| cav_2d_lp_left_refl_c 			| 3 element vector |    | Crystal input face power reflectivity for red1, red2, and blue waves |
| cav_2d_lp_right_refl_c 			| 3 element vector |    | Crystal output face power reflectivity for red1, red2, and blue waves |
| cav_2d_lp_crystal_loss 			| 3 element vector | 1 / mm | Linear absorption for red1, red2, and blue waves |
| cav_2d_lp_left_energy_pwr 		| 3 element vector | W (cw) or J (pulsed) | Left mirror input power in W or pulse energy in J for red1, red2, and blue waves |
| cav_2d_lp_right_energy_pwr 		| 2 element vector | W (cw) or J (pulsed) | Right mirror input power in W or pulse energy in J for red1, and red2 waves |
| cav_2d_lp_pulse_duration 			| 3 element vector or 2x3 array | fhwm ns | Input pulse duration in fwhm ns, optionally spatial Supergaussian coefficient for red1, red2, and blue waves. Supergaussian index should be between 1 and 20, 1 is Gaussian and larger values are increasingly flat-topped. |
| cav_2d_lp_pulse_delay 			| 2 element vector | ns   | Input pulse delay in ns for red1 and red2 wave relative to center of input blue wave |
| cav_2d_lp_beam_diameters 			| 3 element vector or 2x3 array | fwhm mm | Input beam diameters in fwhm mm for red1, red2, and blue waves. Optionally, use array sized 2x3 for elliptical beams where the 1st dimension is in the direction of walk off and perpendicular to it. |
| cav_2d_lp_spatial_supergaussian 	| 3 element vector |      | Super Gaussian coefficient for spatial profile of input beams for red1, red2, and blue waves. Values between 1 and 20, 1 is Gaussian and larger values are increasingly flat-topped. |
| cav_2d_lp_walkoff_angles 			| 3 element vector | mrad | birefringent walk off angles in mrad; for red1, red2, and blue waves |
| cav_2d_lp_beam_offset 			| 3 element vector | mm   | Input beam displacements in the walk off direction for red1, red2, and blue waves.|
| cav_2d_lp_beam_roc 				| 3 element vector or 2x3 array | mm | Input beam radius of curvature in mm for red1, red2, and blue waves. Optionally, array sized 2x3; first dimension is radius of curvature [in the dimension of birefringent walkoff, perpendicular to it]. |
| cav_2d_lp_left_refl_m 			| 3 element vector |      | Left mirror reflectivity for red1, red2, and blue waves. |
| cav_2d_lp_right_refl_m    		| 3 element vector | 	  | Right mirror reflectivity;  for red1, red2, and blue waves. |
| cav_2d_lp_phase_lc        		| 3 element vector | radians | Additional phase introduced between left mirror and crystal for red1, red2, and blue waves. |
| cav_2d_lp_phase_cr        		| 3 element vector | radians | Additional phase introduced between crystal and right mirror for red1, red2, and blue waves. |
| cav_2d_lp_phase_rl        		| 3 element vector | radians | Additional phase introduced between right mirror and left mirror for red1, red2, and blue waves. |
| cav_2d_lp_mirror_rocs     		| 2 element vector | mm      | Radius of curvature for mirrors (left, right order). |
| cav_2d_lp_leg_lengths     		| 3 element vector | mm      | Leg lengs in between left mirror and crystal, crystal and right mirror, and for ring cavitys between right mirror and left mirror. |
| cav_2d_lp_nz              		| (scalar) 		   | Integer | Number of steps through the crystal in z. |
| cav_2d_lp_nxny            		| (2 element vector) | Integer | Number of points in transverse dimensions, in the walk off direction and perpendicular to it.
| cav_2d_lp_crystal_length  		| scalar | mm | Length of opo crystal. |
| cav_2d_lp_lxly            		| 2 element vector | mm | Extent of transverse dimensions, in the direction of birefringent walkoff, and perpendicular to it. |
| cav_2d_lp_cavity_type     		| scalar | 0 or 1 | Cavity type: choose whether cavity is ring (value of 0) or linear (value of 1). |
| cav_2d_lp_cavity_inversion 		| scalar | 0 or 1 | Cavity inversion. 0 for no inversion between cavity passes; 1 inverts beam in walk off direction. |
| cav_2d_lp_deff            		| scalar | pm/V | Effective nonlinear coefficient. Find with Qmix.
| cav_2d_lp_deltak          		| scalar | 1/mm | Phase mismatch Delta k = k_blue - k_red1 - k_red2 |
| cav_2d_lp_nstarts 				| scalar | integer 1 or larger| For pulsed operation, the number of starts. Values of 1 give a time resolution of one round trip time, larger numbers subdivide the pulse into finer temporal resolution. |

# 2D-mix-LP
Long-pulse, diffractive mixing
*snlo_2d_mix_lp_func.m*
2D-mix-LP is similar to PW-mix-LP except it handles the full spatial profiles so it can realistically model circular or elliptical beams with Gaussian or Supergaussian transverse profiles.  It includes diffraction, birefringent walkoff, displaced beams, etc. for gaussian or super-guassian pulses or cw light. The light is assumed to be monochromatic and the pulses are assumed to be long enough that group velocity effects are unimportant. Maxwell's equations are integrated using split-step methods to give true diffractive calculations. For help in setting input values, right-click on the input edit box. Help text appears in the lower text box.

For passing a set of inputs as an input argument to this function, the argument must be a data structure with fields named as follows:
|Fieldname | Shape | Units | Description |
| -- | -- | -- | -- |
| mix_2d_lp_wavelengths | 3 element vector | nm      | Wavelengths for red1, red2, blue waves. |
| mix_2d_lp_ref_inds    | 3 element vector | 	       | Refractive indices in the crystal for the red1, red2, blue waves. |
| mix_2d_lp_phase       | 3 element vector | radians | Input phases for the red1, red2, blue waves. |
| mix_2d_lp_input_refl  | 3 element vector |         | Crystal input face (power) reflectivity for red1, red2, blue waves. |
| mix_2d_lp_output_refl | 3 element vector |         | Crystal output face (power) reflectivity for red1, red2, blue waves. |
| mix_2d_lp_crystal_losses | 3 element vector | 1/mm    | Linear absorption in crystal for red1, red2, blue waves. |
| mix_2d_lp_pulseenergy | 3 element vector | W (cw) or J (pulsed) | Input pulse energies/powers in J/W for red1, red2, blue waves. |
| mix_2d_lp_pulse_durations | 3 element vector or 2x3 array | ns | Pulse duration in ns for red1, red2, blue waves. Optionally, values in second row set temporal supergaussian; values of 1 are Gaussian, larger values are more flat-toppped. |
| mix_2d_lp_pulse_delays  | 3 element vector | ns | Pulse delay for red1, red2, blue waves. |
| mix_2d_lp_beam_diameters | 3 element vector or 2x3 array | fwhm mm | Beam diameter for red1, red2, blue waves. Optionally two rows; the first is in direction of birefringent walkoff, second is perpendicular to it.|
| mix_2d_lp_supergaussian_coeff | 3 element vector| | Spatial Supergaussian coefficient for red1, red2, blue waves. Values between 1 and 20; 1 is Gaussian, larger is more flat-topped.|
| mix_2d_lp_n2_red1 	| 3 element vector | m^2/W | Nonlinear index for red1, red2, blue waves caused by red1 light.|
| mix_2d_lp_n2_red2 	| 3 element vector | m^2/W | Nonlinear index for red1, red2, blue waves caused by red2 light.|
| mix_2d_lp_n2_blue 	| 3 element vector | m^2/W | Nonlinear index for red1, red2, blue waves caused by blue light.|
| mix_2d_lp_beta_red1 | 3 element vector | | beta. Two-photon absorption coefficient for red1, red2, blue waves caused by red1 light. | 
| mix_2d_lp_beta_red2 | 3 element vector | | beta. Two-photon absorption coefficient for red1, red2, blue waves caused by red2 light. | 
| mix_2d_lp_beta_blue | 3 element vector | | beta. Two-photon absorption coefficient for red1, red2, blue waves caused by blue light. | 
| mix_2d_lp_wo_angles | 3 element vector | mrad | Birefringent walk off angle for red1, red2, blue waves. |
| mix_2d_lp_offset_wodir | 3 element vector | mm | Input beam offset in the direction of walk off for red1, red2, blue waves. |
| mix_2d_lp_rad_curv | 3 element vector or 2x3 array | mm | Input beam radius of curvature for red1, red2, blue waves. Optionally two rows; the first is in direction of birefringent walkoff, second is perpendicular to it.|
| mix_2d_lp_nz | scalar | Integer | Number of z steps through crystal. | 
| mix_2d_lp_nxny | 2 element vector | mm | Number of transverse grid points, in the direction of walkoff and perpendicular to it. | 
| mix_2d_lp_crystal_length | scalar | mm | Crystal length. |
| mix_2d_lp_lx_ly | 2 element vector | mm | Transverse grid length, in the direction of walkoff and perpendicular to it. |
| mix_2d_lp_deff | scalar | pm/V | Effective nonlinear coefficient. |
| mix_2d_lp_deltak | scalar | 1/mm | Phase mismatch Deltak = k_blue - k_red1 - k_red2. | 
| mix_2d_lp_nt |scalar | Integer | Number of time points (for pulsed cases). |
| mix_2d_lp_save_movie | scalar | (logical) | Save movie data to file. |
| mix_2d_lp_auto_analyze | scalar | (logical) | Automatically analyze beam for each wave after model finishes (beam quality, tilt, radius of curvature). |

# 2D-mix-SP
Short-pulse, diffractive mixing
*snlo_2d_mix_sp_func.m*
% Units for chirps? Order for chirps? Does each wave get one column, or one row?

|Fieldname | Shape | Units | Description |
| -- | -- | -- | -- |
| mix_2d_sp_wavelengths | 3 element vector | nm | Wavelengths for red1, red2, blue waves. |
| mix_2d_sp_ref_inds    | 3 element vector |  | Crystal refractive indices for red1, red2, blue waves. |
| mix_2d_sp_gvi  		| 3 element vector |  | Crystal group velocity indices for red1, red2, blue waves. |
| mix_2d_sp_gdd 		| 3 element vector | fs^2/mm | Group delay dispersion for red1, red2, blue waves. |
| mix_2d_sp_phase 		| 3 element vector | radians | Input phase for red1, red2, blue waves. |
| mix_2d_sp_input_refl  | 3 element vector |  | Crystal input face  reflectiviy for red1, red2, blue waves. | 
| mix_2d_sp_output_refl | 3 element vector |  | Crystal output face reflectivity for red1, red2, blue waves. | 
| mix_2d_sp_crystal_losses | 3 element vector | % crystal loss (1/mm)
| mix_2d_sp_n2_red1 | 3 element vector |  | n2 (red1) sq cm/W
| mix_2d_sp_n2_red2 | 3 element vector |  | n2 (red2) sq cm/W
| mix_2d_sp_n2_blue | 3 element vector |  | n2 (blue) sq cm/W
| mix_2d_sp_beta_red1 | 3 element vector | | beta (red1) cm/W
| mix_2d_sp_beta_red2 | 3 element vector | | beta (red2) cm/W
| mix_2d_sp_beta_blue | 3 element vector | | beta (blue) cm/W
| mix_2d_sp_pulseenergy | 3 element vector | W or J | Input pulse energies (J) or powers (W) for red1, red2, blue waves.
| mix_2d_sp_pulse_durations |3 element vector or 2x3 array | ps |  Pulse durations (fwhm). Optionally, second row are temporal supergaussian indices -- values of 1 are Gaussian, larger values are more flat-topped. |
| mix_2d_sp_pulse_delays | 2 element vector | ps  | Pulse delays for red1, red2, blue waves. | 
| mix_2d_sp_pulse_chirps | 0, or 3 element vector, or 2x3 array, or 3x3 array | THz/ps, or THz/ps^2, or THz/ps^3 |  Pulse chirps for red1, red2, blue waves. You can enter 0, or a 3 element vector for linear chirps; optionally, add one row to specify the quadratic chirp, and another for cubic chirp. | 
| mix_2d_sp_beam_diameters | 3 element vector or 2x3 array | fwhm mm | Beam diam (fwhm), optionally second row (first is in the direction of birefringent walkoff, second is perpendicular to it). |
| mix_2d_sp_supergaussian_coeff | 3 element vector |  | Spatial supergaussian indices for red1, red2, blue waves. Value of 1 is Gaussian, larger values are more flat-topped. |
| mix_2d_sp_wo_angles | 3 element vector | mrad |  Birefringent walk off angles for red1, red2, blue waves. |
| mix_2d_sp_offset_wodir | 3 element vector | mm | Beam offset in the direction of birefringent walk off for red1, red2, blue waves. |
| mix_2d_sp_rad_curv | 3 element vector or 2x3 array | mm | Input beam radius of curvature (in air) for red1, red2, blue waves. With an optional second row, the first row is radius of curvature in the direction of birefringent walk off and the second is perpendicular to it. |
| mix_2d_sp_nz | scalar | Integer | Number of mixing & propagation z steps through crystal
| mix_2d_sp_nxny | 2 element vector |  Integer | Number of transverse grid points, in the direction of birefringent walkoff and perpendicular to it. | 
| mix_2d_sp_crystal_length | scalar | mm | Crystal length. |
| mix_2d_sp_lx_ly | 2 element vector | Integer | Crystal width and height. | 
| mix_2d_sp_deff | scalar | pm/V | Effective nonlinear coefficient. | 
| mix_2d_sp_deltak | scalar | 1/mm | Phase mismatch Deltak = k_blue - k_red1 - k_red2. | 
| mix_2d_sp_nt | scalar | Integer | Number of time points. | 

# OPG
OPG models single-pass mixing with diffraction, birefringent walkoff, group velocity walk off, and group velocity dispersion with simulated quantum noise in both red waves. The quantum noise is automatically computed based on the inputs. Optional red seed light can be added to the noise by specifying red pulse energies and pulse durations.
| Fieldname | Description |
| - | - |
|  opg_2d_wavelengths:      | (3 element vector) wavelengths in nm; for red1, red2, and blue waves
|  opg_2d_ref_inds:         | (3 element vector) refractive indices; for red1, red2, and blue waves
|  opg_2d_gvi:              | (3 element vector) group velocity indices; for red1, red2, and blue waves
|  opg_2d_gdd:              | (3 element vector) group delay dispersions in fs^2/mm; for red1, red2, and blue waves
|  opg_2d_phase:            | (3 element vector) input phase of waves red1, red2, blue; for red1, red2, and blue waves
|  opg_2d_input_refl:       | (3 element vector) input crystal face power reflectivity; for red1, red2, and blue waves
|  opg_2d_output_refl:      | (3 element vector) output crystal face power reflecitivity; for red1, red2, and blue waves
|  opg_2d_crystal_losses:   | (3 element vector) linear absorption in 1/mm; for red1, red2, and blue waves
|  opg_2d_n2_red1:          | (3 element vector) for red1, nonlinear refractive index in sq cm/W; for red1, red2, and blue waves
|  opg_2d_n2_red2:          | (3 element vector) for red2, nonlinear refractive index in sq cm/W; for red1, red2, and blue waves
|  opg_2d_n2_blue:          | (3 element vector) for blue, nonlinear refractive index in sq cm/W; for red1, red2, and blue waves
|  opg_2d_beta_red1:        | (3 element vector) for red1, two-photon absorption in cm/W; for red1, red2, and blue waves
|  opg_2d_beta_red2:        | (3 element vector) for red2, two-photon absorption in cm/W; for red1, red2, and blue waves
|  opg_2d_beta_blue:        | (3 element vector) for blue, two-photon absorption in cm/W; for red1, red2, and blue waves
|  opg_2d_pulseenergy:      | (3 element vector) input pulse energy in J; for red1, red2, and blue waves
|  opg_2d_pulse_durations:  | (3 element vector) input pulse durations power fwhm ps; for red1, red2, and blue waves; optionally array of 3x2 where 2nd dimension is integer 0-10 for sech^2 or (super) gaussian shape
|  opg_2d_pulse_delays:     | (3 element vector) input pulse delays relative to blue pulse center in ps; vector of 2 elements: for red1, red2
|  opg_2d_pulse_chirps:     | (3 element vector) chirp of input pulse in THz/ps; for red1, red2, and blue waves; optionally array of 3x2 or 3x3 for 2nd & 3rd order chirps (in THz/ps^2 and THz/ps^3); 
|  opg_2d_beam_diameters:   | (3 element vector) diameters of input beam fwhm mm; for red1, red2, and blue waves; optionally, 3x2 array with 1st column walk off direction diam, 2nd column perpendicular to that
|  opg_2d_supergaussian_coeff: | (3 element vector) spatial super gaussian coefficient (integer valued); for red1, red2, and blue waves
|  opg_2d_wo_angles:        | (3 element vector) birefringent walk off angle in mrad; for red1, red2, and blue waves
|  opg_2d_offset_wodir:     | (3 element vector) input beam displacement relative to center of grid in mm; ; for red1, red2, and blue waves
|  opg_2d_rad_curv:         | (3 element vector) radius of curvature for input beams in mm; for red1, red2, and blue waves; optionally, 3x2 array with 1st column walk off direction diam, 2nd column perpendicular to that
|  opg_2d_nz:               | (scalar) number of steps in z
|  opg_2d_nxny:             | (2 element vector) number of points in transverse dimensions (walk off direction and perpendicular to that)
|  opg_2d_grid_duration:    | (scalar) duration for simulation in ps
|  opg_2d_crystal_length:   | (scalar) crystal length in mm
|  opg_2d_lx_ly:            | (2 element vector) transverse spatial grid extend in mm; walk off direction and perpendicular to that
|  opg_2d_deff:             | (scalar) effective nonlinear coefficient in pm/V
|  opg_2d_deltak:           | (scalar) phase mismatch delta k in 1/mm; (k_blue - k_red1 - k_red2
|  opg_2d_nt:               | (scalar) number of points in time for simulation


