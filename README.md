# mlSNLO-batch
Example batch scripts for mlSNLO nonlinear optics numerical modeling software. These examples work for the MATLAB App version of mlSNLO which requires a MATLAB license. For more information see [as-photonics.com/products/mlsnlo](http://www.as-photonics.com/products/mlsnlo).

## Table of contents
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
  
## Introduction
These scripts are MATLAB scripts. They require the MATLAB App version of mlSNLO to run. These files are examples of how to call most SNLO models from a batch script. This is often useful to optimize the design of your nonlinear optical device; for example, you might want to choose the nonlinear crystal length with produces the largest output energy of a pulsed second harmonic generation device. Using a batch script like these examples could automate this optimization: set up the script, press 'Run' and view the results.

Each example script is comprised of a for loop in which iteratively creates a set of inputs for an SNLO model function, calls the model and passes it the set of inputs, simulates pressing the 'Run' (and sometimes 'Accept') buttons on the model, loads the model output from file, processes and stores the result.

The inputs for each SNLO function match the main input form input boxes. The inputs are passed to each SNLO function as an input argument consisting of a MATLAB data structure with fieldnames the model requires, each field populated with values with units that match the input form labels.

Each example script includes many comments explaining what they do. 

For each SNLO modeling function below, there is a link to the 'Help' text file and link to any example batch scripts.

The help text files here are the ones displayed after clicking the 'Help' button in the main menu of SNLO.

## SNLO functions
(in the order they appear in the main menu)

### Qmix
* [helpfiles/Qmix_h.txt](helpfiles/Qmix_h.txt)
* [test_calling_snlo_qmix_func.m](test_calling_snlo_qmix_func.m)

### Ncpm
* [helpfiles/Ncpm_h.txt](helpfiles/Ncpm_h.txt)
* [test_calling_snlo_ncpm_func.m](test_calling_snlo_ncpm_func.m)

### Bmix
* [helpfiles/Bmix_h.txt](Bmix_h.txt)
* [test_calling_snlo_bmix_func.m](test_calling_snlo_bmix_func.m)

### Focus
* [helpfiles/Focus_h.txt](Focus_h.txt)
* [test_calling_snlo_focus_func.m](test_calling_snlo_focus_func.m)

### Cavity2
* [helpfiles/Brewster_cavity_help.txt](helpfiles/Brewster_cavity_help.txt)

### GVM
* [helpfiles/GVM_h.txt](GVM_h.txt)
* [test_calling_snlo_gvm_func.m](test_calling_snlo_gvm_func.m)

### Opo angles
* [helpfiles/Opoangles_h.txt](Opoangles_h.txt)
* [test_calling_snlo_opoangles_func.m](test_calling_snlo_opoangles_func.m)

### QPM
* [helpfiles/QPM_h.txt](QPM_h.txt)
* [test_calling_snlo_qpm_func.m](test_calling_snlo_qpm_func.m)

### PW-mix-LP
* [helpfiles/PWmixLP_h.txt](PWmixLP_h.txt)
* [test_calling_snlo_pw_mix_lp_func.m](test_calling_snlo_pw_mix_lp_func.m)

### PW-mix-SP
* [helpfiles/PWmixSP_h.txt](PWmixSP_h.txt)
* [test_calling_snlo_pw_mix_sp_func.m](test_calling_snlo_pw_mix_sp_func.m)

### PW-cav-LP
* [helpfiles/PWcavLP_h.txt](PWcavLP_h.txt)
* [test_calling_snlo_pw_cav_lp_func.m](test_calling_snlo_pw_cav_lp_func.m)

### PW-mix-BB
* [helpfiles/PWmixBB_h.txt](PWmixBB_h.txt)
* [test_calling_snlo_pw_mix_bb_func.m](test_calling_snlo_pw_mix_bb_func.m)

### PW-OPO-BB
* [helpfiles/PWopoBB_h.txt](PWopoBB_h.txt)
* [test_calling_snlo_pw_opo_bb_func.m](test_calling_snlo_pw_opo_bb_func.m)

### PW-OPO-SP
* [helpfiles/PWopoSP_h.txt](PWopoSP_h.txt)
* [test_calling_snlo_pw_opo_sp_func.m](test_calling_snlo_pw_opo_sp_func.m)

### 2D-mix-LP
* 2D-mix-LP now allows arbitrary, user-specified spatial and temporal input electric field profiles. These examples show how to do this:
	 *  [example_custom_input_fields_2d_mix_lp_cw.m](example_custom_input_fields_2d_mix_lp_cw.m): all-continuous wave case
	 *  [example_custom_input_fields_2d_mix_lp_pulsed.m](example_custom_input_fields_2d_mix_lp_pulsed.m): all-pulsed case
	 *  [example_custom_input_fields_2d_mix_lp_cw_and_pulse_mixed.m](example_custom_input_fields_2d_mix_lp_cw_and_pulse_mixed.m): a combination of pulsed and cw input beams
	 *  [example_2d_mix_lp_shg_and_thg_generation.m](example_2d_mix_lp_shg_and_thg_generation.m): a two-crystal system using the output of an second-harmonic generator as the input for a third-harmonic generator
	 *  [example_prop_2d_mix_lp_shg_and_thg_generation.m](example_prop_2d_mix_lp_shg_and_thg_generation.m): a two-crystal system using the output of an second-harmonic generator as the input for a third-harmonic generator with a free-space FFT propagation bwetween the two crystals
* [helpfiles/2DmixLP_h.txt](2DmixLP_h.txt)
* [test_calling_snlo_2d_mix_lp_func.m](test_calling_snlo_2d_mix_lp_func.m)
* [test_calling_snlo_2d_mix_lp_func_deltak.m](test_calling_snlo_2d_mix_lp_func_deltak.m)

### 2D-mix-SP
* 2D-mix-SP now allows arbitrary, user-specified input temporal pulse profiles (still limited to Gaussian or super-Gaussian spatial profiles). 1 example is included:
	 *  [test_calling_custom_pulse_2d_mix_sp.m](test_calling_custom_pulse_2d_mix_sp.m): 
* [helpfiles/2DmixSP_h.txt](2DmixSP_h.txt)
* [test_calling_snlo_2d_mix_sp_func_deltak.m](test_calling_snlo_2d_mix_sp_func_deltak.m)
* [test_calling_snlo_2d_mix_sp_func_deltak2.m](test_calling_snlo_2d_mix_sp_func_deltak2.m)

### 2D-cav-LP
* [helpfiles/2DcavLP_h.txt](2DcavLP_h.txt)
* [test_calling_snlo_2d_cav_lp_func.m](test_calling_snlo_2d_cav_lp_func)

### OPG
* [helpfiles/OPG_h.txt](OPG_h.txt)
* [test_calling_snlo_opg_func_deltak.m](test_calling_snlo_opg_func_deltak.m)
* [test_calling_snlo_opg_func_pump_energy.m](test_calling_snlo_opg_func_pump_energy.m)

### Other files:
* [simple_freespace_propagate.m](simple_freespace_propagate.m) A simple function to propagate a beam in free space using FFT
* [example_simple_freespace_propagate.m](example_simple_freespace_propagate.m) A simple example script to use the script (simple_freespace_propagate.m)
* [crystal_info.m](crystal_info.m) a script that does most of what snlo_ref_ind_func does but without the graphical user interface components; it uses ref_ind and N12 to find refractive indices for isotropic, uniaxial, and biaxial crystals, as well as group velocity indices, group delay dispersions, third order dispersions, walkoff angles, and thermo optic coefficients for a given crystal, temperature, wavelength, and propagation direction.
