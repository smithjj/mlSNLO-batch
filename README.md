# mlSNLO-batch
Example batch scripts for mlSNLO nonlinear optics numerical modeling software. These examples work for the MATLAB App version of mlSNLO which requires a MATLAB license. For more information see [as-photonics.com/products/mlsnlo](http://www.as-photonics.com/products/mlsnlo).

Below are short descriptions of each mlSNLO model function, a link to the text help displayed in SNLO, and some example scripts which demonstrate how to call the model function from a MATLAB script for use in batch runs. 

### Table of contents
* [Introduction, conventions and syntax](#introduction)
* [SNLO functions](#snlo-functions)
 * [Qmix](#qmix)
 * [Ncpm](#ncpm)
 * [Bmix](#bmix)
 * [Focus](#focus)
 * [Cavity2](#cavity2)
 * [GVM](#gvm)
 * [OPO angles](#opo-angles)
 * [QPM](#qpm)
 * [PW-mix-LP](#pw-mix-lp)
 * [PW-mix-SP](#pw-mix-sp)
 * [PW-cav-LP](#pw-cav-lp)
 * [PW-mix-BB](#pw-mix-bb)
 * [PW-OPO-BB](#pw-opo-bb)
 * [PW-OPO-SP](#pw-opo-sp)
 * [2D-mix-LP](#2d-mix-lp)
 * [2D-mix-SP](#2d-mix-sp)
 * [2D-cav-LP](#2d-cav-lp)
 * [OPG](#opg)
 
### Introduction
These scripts are MATLAB scripts. They require the MATLAB App version of mlSNLO to run. These files are examples of how to call most SNLO models from a batch script. This is often useful to optimize the design of your nonlinear optical device; for example, you might want to choose the nonlinear crystal length with produces the largest output energy of a pulsed second harmonic generation device. Using a batch script like these examples could automate this optimization: set up the script, press 'Run' and view the results.

Each example script is comprised of a for loop in which iteratively creates a set of inputs for an SNLO model function, calls the model and passes it the set of inputs, simulates pressing the 'Run' (and sometimes 'Accept') buttons on the model, loads the model output from file, processes and stores the result.

The inputs for each SNLO function match the main input form input boxes. The inputs are passed to each SNLO function as an input argument consisting of a MATLAB data structure with fieldnames the model requires, each field populated with values with units that match the input form labels.

Each example script includes many comments explaining what they do. 

For each SNLO modeling function below, there is a link to the 'Help' text file and link to any example batch scripts.

The help text files here are the ones displayed after clicking the 'Help' button in the main menu of SNLO.

Additional information about how to specify the model's input parameters can be found at [mlSNLO_function_input_data_structures.md](mlSNLO_function_input_data_structures.md)

### SNLO functions
(in the order they appear in the main menu)

#### Qmix
QMIX helps you quickly select the best crystal for your application from a list of over 60 crystals. It calculates all the possible phase-matching orientations for the chosen crystal at the wavelengths specified, returning phase-matching angles, polarizations, refractive indexes, group velocities, group delay dispersion, effective nonlinearity, crystal tilt tolerance, acceptance angles, acceptance bandwidths, and acceptance temperatures. The three beams are assumed to have collinear propagation vectors. If the nonlinearity, deff, is zero, the other properties are not calculated. (If you want to find noncritical wavelengths, use function Ncpm.) 
* [helpfiles/Qmix_h.txt](helpfiles/Qmix_h.txt)
* [example_script_qmix.m](example_script_qmix.m)

### Ncpm
NCPM (noncritical phase matching) helps you quickly search the list of crystals for one that can noncritically phase match your wavelengths. Noncritical phasematching refers to collinear phase matching with no birefringent walk off. It is desirable for efficient mixing of weak, tightly focused beams, and for maintaining the widest acceptance angle for all three beams. 
* [helpfiles/Ncpm_h.txt](helpfiles/Ncpm_h.txt)
* [example_script_ncpm.m](example_script_ncpm.m)

#### Bmix
Bmix calculates the properties of biaxial crystals for propagation in any direction including out of the principal planes (in contrast to Qmix where propagation in a principal plane is assumed). The calculations are based on the d-tensor and Sellmeier listed in Qmix. Phase matching curves are computed along with the associated effective nonlinearity, walkoff, and acceptance angles. The purple circle on the phase matching plot indicates the location of the optic axes in the x-z plane. Parallel propagation directions (k vectors) of the three interacting waves red1, red2, and blue are assumed. 
* [helpfiles/Bmix_h.txt](helpfiles/Bmix_h.txt)
* [example_script_bmix.m](example_script_bmix.m)

#### Focus
FOCUS helps calculate the curvature of a Gaussian beam at a specific location that will give a certain waist size at a certain location. This information is needed in the two-dimensional mixing models 2D-mix-LP and 2D-mix-SP if the beams are focusing (have curved wavefronts).
* [helpfiles/Focus_h.txt](helpfiles/Focus_h.txt)
* [example_script_focus.m](example_script_focus.m)

#### Cavity2
This function is an extension of the original Cavity function. It adds the ability to model astigmatic cavities, with the astigmatism arising from tilted and curved mirrors and Brewster-cut crystals. The cavity is assumed to be a ring with two curved mirrors bracketing the nonlinear crystal. The ring is assumed to be planar so the in-plane and out-of-plane mode profiles are independent of one another and can be computed separately. The mode profiles for the two planes are both computed and displayed whenever the function is run. The function also allows user specification of a disturbance to the cavity. Disturbances can be either a tilt or a displacement of the mode in the plane of the cavity, and they can be located anywhere in the cavity. This calculation is useful in understanding how tilting a cavity mirror shifts the mode or how refraction in a tuning prism shifts the mode.
* [helpfiles/Brewster_cavity_help.txt](helpfiles/Brewster_cavity_help.txt)

#### GVM
GVM calculates group velocity (mis)match for noncollinear phase matching. The paper "Group-velocity-matched three-wave mixing in birefringent crystals" Optics Letters vol. 26, page 719 (2001) describes the method of group velocity matching in birefringent crystals implemented in this function. This paper is available online at http://www.as-photonics.com/Publications.html. The calculations are based on the d-tensor and Sellmeier equations displayed in Qmix. The book "Crystal Nonlinear Optics: with SNLO examples" (as-photonics.com/book) also describes the method. 
* [helpfiles/GVM_h.txt](helpfiles/GVM_h.txt)
* [example_script_gvm.m](example_script_gvm.m)

#### Opo angles
Opoangles computes OPO/OPA red1 and red2 wavelengths as a function of propagation angle in birefringent crystals assuming the blue (pump) wavelength is fixed. For biaxial crystals, propagation is in one of the principal planes of the crystal so the Plane radio buttons are activated. Select the plane of interest. In biaxial crystals the angle is theta, measured from the Z axis for propagation in the XZ and YZ planes. For propagation in the XY plane the angle is phi, measured from the X axis toward the Y axis. Calculations are based on the d-tensor and Sellmeier listed in Qmix. 
* [helpfiles/Opoangles_h.txt](helpfiles/Opoangles_h.txt)
* [example_script_opoangles.m](example_script_opoangles.m)

#### QPM
QPM calculates properties of quasiphasematched materials including the poling, or domain reversal, period (always specified at room temperature), the temperature bandwidth, the frequency bandwidth, the group velocities, the group delay dispersions.  In addition the tuning of phase matched red1 and red2 wavelengths due to blue (pump) tuning or crystal temperature tuning can be computed using the Pump Tune and Temp Tune buttons. The Sellmeier equations are those referenced in Qmix when a crystal is selected. The temperature range and acceptance bandwidths are calculated as described in Qmix help.
* [helpfiles/QPM_h.txt](helpfiles/QPM_h.txt)
* [example_script_qpm.m](example_script_qpm.m)

#### PW-mix-LP
PW-mix-LP computes three-wave mixing process for long-pulse plane- waves. The input irradiances are those of the central rays of spatial lowest-order Gaussian beams with diameters and energies specified in the input form. This function ignores birefringent and group velocity walkoff. It is intended for very quick evaluation of a proposed mixing process.
* [helpfiles/PWmixLP_h.txt](helpfiles/PWmixLP_h.txt)
* [example_script_pw_mix_lp.m](example_script_pw_mix_lp.m)

#### PW-mix-SP
PW-mix-SP computes three-wave mixing process for plane waves with irradiances equal to those of the central rays of spatial lowest order Gaussian beams with beam diameters and energies specified in the input form. It ignores birefringent walkoff but includes group velocity effects, so it is useful for modeling very short pulses. cw light is not allowed in this function. Pulses are Gaussian or Supergaussian in time at the input face. The three input pulses can have different durations and delays.
* [helpfiles/PWmixSP_h.txt](helpfiles/PWmixSP_h.txt)
* [example_script_pw_mix_sp.m](example_script_pw_mix_sp.m)

#### PW-cav-LP
PW-cav-LP models plane-wave mixing in cavities. Input irradiances are those of the central rays of lowest order Gaussian beams with diameters specified by the Beam diam. input values. You can model an OPO, cavity resonant SHG or any other cavity mixing process by specifying the input beams and the reflectivity of the mirrors at those wavelengths. The crystal mixing is similar to PW-mix-LP. The cavity configuration can be standing wave or ring. 
* [helpfiles/PWcavLP_h.txt](helpfiles/PWcavLP_h.txt)
* [example_script_pw_cav_lp.m](example_script_pw_cav_lp.m)

#### PW-mix-BB
PW-mix-BB computes three-wave mixing processes for plane waves with input irradiances equal to those of the central ray of lowest order Gaussian beams with energy and diameter specified on the input form. It ignores birefringent walkoff but includes group velocity effects. It is intended for modeling multi-longitudinal-mode pulses. The number of modes populated in each model run is sufficient to cover the specified bandwidth for each beam. Quantum noise is automatically added to the red1 and red2 waves. The noise bandwidths of the red waves can be extended by specifying the bandwidth parameter. This inclusion of quantum noise is to approximate OPG modeling. However, to properly model OPG the spatial structure of the red beams should be noisy as well.
* [helpfiles/PWmixBB_h.txt](helpfiles/PWmixBB_h.txt)
* [example_script_pw_mix_bb.m](example_script_pw_mix_bb.m)

#### PW-OPO-BB
PW-OPO-BB models broad-bandwidth, nanosecond, plane-wave (actually the central ray of a spatial Gaussian) OPO's. It is similar to PW-cav-LP except it allows for differing group velocities for the three waves and broad bandwidths. The mathematical methods are described in paper "Numerical models of broad bandwidth nanosecond optical parametric oscillators" JOSA B vol. 16 p. 609 (1999), available at http://www.as-photonics.com/Publications.html This function permits studies of unseeded OPO's and the transition to seeding. It also models double resonance effects and pumping by multi-mode lasers.
* [helpfiles/PWopoBB_h.txt](helpfiles/PWopoBB_h.txt)
* [example_script_pw_opo_bb.m](example_script_pw_opo_bb.m)

#### PW-OPO-SP
PW-OPO-SP models a singly-resonant, plane-wave, synchronously-pumped OPO.  The irradiance of the plane pump wave is set equal to the central ray of a lowest-order spatial Gaussian whose beam diameter is specified in the input form under 'Beam diameter'.  This function includes group velocity and group delay dispersion but not higher order dispersion. The OPO cavity is a ring, with only the red1 wave resonated. The duration of the blue (pump) pulses are specified. 
* [helpfiles/PWopoSP_h.txt](helpfiles/PWopoSP_h.txt)
* [example_script_pw_opo_sp.m](example_script_pw_opo_sp.m)

#### 2D-mix-LP
2D-mix-LP is similar to PW-mix-LP except it handles the full spatial profiles so it can realistically model circular or elliptical beams with Gaussian or Supergaussian transverse profiles.  It includes diffraction, birefringent walkoff, displaced beams, etc. for gaussian or super-guassian pulses or cw light. The light is assumed to be monochromatic and the pulses are assumed to be long enough that group velocity effects are unimportant. Maxwell's equations are integrated using split-step methods to give true diffractive calculations.
* 2D-mix-LP now allows arbitrary, user-specified spatial and temporal input electric field profiles. These examples show how to do this:
	 * [example_custom_input_fields_2d_mix_lp_cw.m](example_custom_input_fields_2d_mix_lp_cw.m): all-continuous wave case
	 * [example_custom_input_fields_2d_mix_lp_pulsed.m](example_custom_input_fields_2d_mix_lp_pulsed.m): all-pulsed case
	 * [example_custom_input_fields_2d_mix_lp_cw_and_pulse_mixed.m](example_custom_input_fields_2d_mix_lp_cw_and_pulse_mixed.m): a combination of pulsed and cw input beams
	 * [example_2d_mix_lp_shg_and_thg_generation.m](example_2d_mix_lp_shg_and_thg_generation.m): a two-crystal system using the output of an second-harmonic generator as the input for a third-harmonic generator
	 * [example_prop_2d_mix_lp_shg_and_thg_generation.m](example_prop_2d_mix_lp_shg_and_thg_generation.m): a two-crystal system using the output of an second-harmonic generator as the input for a third-harmonic generator with a free-space FFT propagation bwetween the two crystals
* [helpfiles/2DmixLP_h.txt](helpfiles/2DmixLP_h.txt)
* [example_script_2d_mix_lp.m](example_script_2d_mix_lp.m)
* [example_script_2d_mix_lp_vary_deltak.m](example_script_2d_mix_lp_vary_deltak.m)

#### 2D-mix-SP
2D-mix-SP models single-pass mixing with diffraction, birefringent walkoff, group velocity walk off, and group velocity dispersion. This is a slow memory hog so try to select the best values for the time and spatial grids using the faster functions PW-mix-SP, PW-mix-LP, and 2D-mix-LP before running it.
* 2D-mix-SP now allows arbitrary, user-specified input temporal pulse profiles (still limited to Gaussian or super-Gaussian spatial profiles). 1 example is included:
	 * [example_script_2d_mix_sp_custom_pulse_profile.m](example_script_2d_mix_sp_custom_pulse_profile.m): 
* [helpfiles/2DmixSP_h.txt](helpfiles/2DmixSP_h.txt)
* [example_script_2d_mix_sp_vary_deltak.m](example_script_2d_mix_sp_vary_deltak.m)
* [example_script_2d_mix_sp_vary_deltak_2.m](example_script_2d_mix_sp_vary_deltak_2.m)

#### 2D-cav-LP
2D-cav-LP is similar to PW-cav-LP except it handles the full spatial profiles so it can realistically model beams with Gaussian or Supergaussian transverse profiles. It includes diffraction, birefringent walkoff, displaced beams, etc. for pulsed or cw light. The light is assumed to be monochromatic and the pulses are assumed to be long enough that group velocity effects are unimportant. 
* [helpfiles/2DcavLP_h.txt](helpfiles/2DcavLP_h.txt)
* [example_script_2d_cav_lp.m](example_script_2d_cav_lp.m)

#### OPG
OPG models single-pass mixing with diffraction, birefringent walkoff, group velocity walk off, and group velocity dispersion with simulated quantum noise in both red waves. The quantum noise is automatically computed based on the inputs. Optional red seed light can be added to the noise by specifying red pulse energies and pulse durations. 
* [helpfiles/OPG_h.txt](helpfiles/OPG_h.txt)
* [example_script_opg.m](example_script_opg.m)
* [example_script_opg_vary_deltak.m](example_script_opg_vary_deltak.m)
* [example_script_opg_vary_pump_energy.m](example_script_opg_vary_pump_energy.m)

#### Other files:
* [simple_freespace_propagate.m](simple_freespace_propagate.m) A simple function to propagate a beam in free space using FFT
* [example_simple_freespace_propagate.m](example_simple_freespace_propagate.m) A simple example script to use the script (simple_freespace_propagate.m)
* [crystal_info.m](crystal_info.m) a script that does most of what snlo_ref_ind_func does but without the graphical user interface components; it uses ref_ind and N12 to find refractive indices for isotropic, uniaxial, and biaxial crystals, as well as group velocity indices, group delay dispersions, third order dispersions, walkoff angles, and thermo optic coefficients for a given crystal, temperature, wavelength, and propagation direction.
