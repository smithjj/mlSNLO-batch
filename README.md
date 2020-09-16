# mlSNLO-batch
Example batch scripts for mlSNLO nonlinear optics numerical modeling software. These examples work for the MATLAB App version of mlSNLO which requires a MATLAB license. For more information see [as-photonics.com/products/matlab-snlo](http://www.as-photonics.com/products/matlab-snlo).

## Table of contents
* [Introduction, conventions and syntax](#introduction)
* [SNLO functions](#snlo-functions)
  * [Qmix](#qmix)
  * [Bmix](#bmix)
  * [Focus](#focus)
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

Each example script is comprised of a for loop in which iteration creates a set of inputs for an SNLO model function, calls the model and passes it the set of inputs, simulates pressing the 'Run' (and sometimes 'Accept') buttons on the model, loads the model output from file, processes and stores the result.

The inputs for each SNLO function match the main input form input boxes. The inputs are passed to each SNLO function as an input argument consisting of a MATLAB data structure with fieldnames the model requires, each field populated with values with units that match the input form labels.

Each example script includes many comments explaining what they do. 

For each SNLO modeling function below, there is a link to the 'Help' text file and link to any example batch scripts.

The help text files here are the ones displayed after clicking the 'Help' button in the main menu of SNLO.

## SNLO functions
### Qmix
* <a href="/helpfiles/Qmix_h.txt">Qmix_h.txt</a>
* <a href="/test_calling_snlo_qmix_func.m">test_calling_snlo_qmix_func.m</a>
* <a href="/test_calling_snlo_qmix_func2.m">test_calling_snlo_qmix_func2.m</a>

### Bmix
* <a href="/helpfiles/Bmix_h.txt">Bmix_h.txt</a>
* <a href="/test_calling_snlo_bmix_func.m">test_calling_snlo_bmix_func.m</a>

### Focus
* <a href="/helpfiles/Focus_h.txt">Focus_h.txt</a>
* <a href="/test_calling_snlo_focus_func.m">test_calling_snlo_focus_func.m</a>

### GVM
* <a href="/helpfiles/GVM_h.txt">GVM_h.txt</a>
* <a href="/test_calling_snlo_gvm_func.m">test_calling_snlo_gvm_func.m</a>

### Opo angles
* <a href="/helpfiles/Opoangles_h.txt">Opoangles_h.txt</a>
* <a href="/test_calling_snlo_opoangles_func.m">test_calling_snlo_opoangles_func.m</a>

### QPM
* <a href="/helpfiles/QPM_h.txt">QPM_h.txt</a>
* <a href="/test_calling_snlo_qpm_func.m">test_calling_snlo_qpm_func.m</a>

### PW-mix-LP
* <a href="/helpfiles/PWmixLP_h.txt">PWmixLP_h.txt</a>
* <a href="/test_calling_snlo_pw_mix_lp_func.m">test_calling_snlo_pw_mix_lp_func.m</a>

### PW-mix-SP
* <a href="/helpfiles/PWmixSP_h.txt">PWmixSP_h.txt</a>
* <a href="/test_calling_snlo_pw_mix_sp_func.m">test_calling_snlo_pw_mix_sp_func.m</a>

### PW-cav-LP
* <a href="/helpfiles/PWcavLP_h.txt">PWcavLP_h.txt</a>
* <a href="/test_calling_snlo_pw_cav_lp_func.m">test_calling_snlo_pw_cav_lp_func.m</a>

### PW-mix-BB
* <a href="/helpfiles/PWmixBB_h.txt">PWmixBB_h.txt</a>
* <a href="/test_calling_snlo_pw_mix_bb_func.m">test_calling_snlo_pw_mix_bb_func.m</a>

### PW-OPO-BB
* <a href="/helpfiles/PWopoBB_h.txt">PWopoBB_h.txt</a>
* <a href="/test_calling_snlo_pw_opo_bb_func.m">test_calling_snlo_pw_opo_bb_func.m</a>

### PW-OPO-SP
* <a href="/helpfiles/PWopoSP_h.txt">PWopoSP_h.txt</a>
* <a href="/test_calling_snlo_pw_opo_sp_func.m">test_calling_snlo_pw_opo_sp_func.m</a>

### 2D-mix-LP
* <a href="/helpfiles/2DmixLP_h.txt">2DmixLP_h.txt</a>
* <a href="/test_calling_snlo_2d_mix_lp_func.m">test_calling_snlo_2d_mix_lp_func.m</a>
* <a href="/test_calling_snlo_2d_mix_lp_func_deltak.m">test_calling_snlo_2d_mix_lp_func_deltak.m</a>

### 2D-mix-SP
* <a href="/helpfiles/2DmixSP_h.txt">2DmixSP_h.txt</a>
* <a href="/test_calling_snlo_2d_mix_sp_func_deltak.m">test_calling_snlo_2d_mix_sp_func_deltak.m</a>
* <a href="/test_calling_snlo_2d_mix_sp_func_deltak2.m">test_calling_snlo_2d_mix_sp_func_deltak2.m</a>

### 2D-cav-LP
* <a href="/helpfiles/2DcavLP_h.txt">2DcavLP_h.txt</a>
* <a href="/test_calling_snlo_2d_cav_lp_func.m">test_calling_snlo_2d_cav_lp_func</a>

### OPG
* <a href="/helpfiles/OPG_h.txt">OPG_h.txt</a>
* <a href="/test_calling_snlo_opg_func_deltak.m">test_calling_snlo_opg_func_deltak.m</a>
* <a href="/test_calling_snlo_opg_func_pump_energy.m">test_calling_snlo_opg_func_pump_energy.m</a>
