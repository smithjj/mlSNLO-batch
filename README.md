# mlSNLO-batch
Example batch scripts for mlSNLO nonlinear optics numerical modeling software. These examples work for the MATLAB App version of mlSNLO which requires a MATLAB license. For more information see [as-photonics.com/products/matlab-snlo](http://www.as-photonics.com/products/matlab-snlo).

## Table of contents
* [Introduction, conventions and syntax](#introduction,-conventions-and-syntax)
* [SNLO functions](#funcs)
  * [Qmix](#qmix)
  * [Bmix](#bmix)
  * [Focus](#focus)
  * [GVM](#gvm)
  * [OPO angles](#opoangles)
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
  
## Introduction, conventions and syntax
These scripts are MATLAB scripts. They require the MATLAB App version of mlSNLO to run. These files are examples of how to call most mlSNLO models from a batch script. This is often useful to optimize the design of your nonlinear optical device; for example, you might want to choose the nonlinear crystal length with produces the largest output energy of a pulsed second harmonic generation device. Using a batch script like these examples could automate this optimization: set up the script, press 'Run' and view the results.

Each example script is comprised of a for loop in which iteration creates a set of inputs for an SNLO model function, calls the model and passes it the set of inputs, simulates pressing the 'Run' (and sometimes 'Accept') buttons on the model, loads the model output from file, processes and stores the result.

The inputs for each SNLO function match the main input form input boxes. The inputs are passed to each SNLO function as an input argument consisting of a MATLAB data structure with fieldnames the model requires, each field populated with values with units that match the input form labels.

Each example script includes many comments explaining what they do. 
