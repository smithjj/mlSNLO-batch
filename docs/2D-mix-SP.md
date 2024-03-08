# 2D-mix-SP
The units are usually kept in MKS for simplicity and consistency.

<a href="#inputs">Inputs</a> and <a href="#outputs">outputs</a> below.

# Inputs
<p id="inputs">2D-mix-SP *(snlo_2d_mix_sp_func)* accepts model parameters as a input argument to the MATLAB function file which is a data structure. The structure looks for matching fieldnames and tries to use the corresponding values to populate the edit boxes in the main GUI input form. Below is a list of the fieldnames with a remark about how the values should be formatted (and their units):</p>

|**Fieldname** | **Shape** | **Units** | **Description** |
| -- | -- | -- | -- |
| **mix_2d_sp_wavelengths** | 3 element vector | nm | Wavelengths for red1, red2, blue waves. |
| **mix_2d_sp_ref_inds**    | 3 element vector |  | Crystal refractive indices for red1, red2, blue waves. |
| **mix_2d_sp_gvi**  		    | 3 element vector |  | Crystal group velocity indices for red1, red2, blue waves. |
| **mix_2d_sp_gdd** 		    | 3 element vector | fs^2/mm | Group delay dispersion for red1, red2, blue waves. |
| **mix_2d_sp_phase** 		  | 3 element vector | radians | Input phase for red1, red2, blue waves. |
| **mix_2d_sp_input_refl**  | 3 element vector |  | Crystal input face  reflectiviy for red1, red2, blue waves. | 
| **mix_2d_sp_output_refl** | 3 element vector |  | Crystal output face reflectivity for red1, red2, blue waves. | 
| **mix_2d_sp_crystal_losses** | 3 element vector | Crystal linear absorption (1/mm) |
| **mix_2d_sp_n2_red1**     | 3 element vector |  | n2 (red1) sq cm/W
| **mix_2d_sp_n2_red2**     | 3 element vector |  | n2 (red2) sq cm/W
| **mix_2d_sp_n2_blue**     | 3 element vector |  | n2 (blue) sq cm/W
| **mix_2d_sp_beta_red1**   | 3 element vector | | beta (red1) cm/W
| **mix_2d_sp_beta_red2**   | 3 element vector | | beta (red2) cm/W
| **mix_2d_sp_beta_blue**   | 3 element vector | | beta (blue) cm/W
| **mix_2d_sp_pulseenergy** | 3 element vector | W or J | Input pulse energies (J) or powers (W) for red1, red2, blue waves.
| **mix_2d_sp_pulse_durations** |3 element vector or 2x3 array | ps |  Pulse durations (fwhm). Optionally, second row are temporal supergaussian indices -- values of 1 are Gaussian, larger values are more flat-topped. |
| **mix_2d_sp_pulse_delays** | 2 element vector | ps  | Pulse delays for red1, red2, blue waves. | 
| **mix_2d_sp_pulse_chirps** | 0, or 3 element vector, or 2x3 array, or 3x3 array | THz/ps, or THz/ps^2, or THz/ps^3 |  Pulse chirps for red1, red2, blue waves. You can enter 0, or a 3 element vector for linear chirps; optionally, add one row to specify the quadratic chirp, and another for cubic chirp. | 
| **mix_2d_sp_beam_diameters** | 3 element vector or 2x3 array | fwhm mm | Beam diam (fwhm), optionally second row (first is in the direction of birefringent walkoff, second is perpendicular to it). |
| **mix_2d_sp_supergaussian_coeff** | 3 element vector |  | Spatial supergaussian indices for red1, red2, blue waves. Value of 1 is Gaussian, larger values are more flat-topped. |
| **mix_2d_sp_wo_angles**   | 3 element vector | mrad |  Birefringent walk off angles for red1, red2, blue waves. |
| **mix_2d_sp_offset_wodir** | 3 element vector | mm | Beam offset in the direction of birefringent walk off for red1, red2, blue waves. |
| **mix_2d_sp_rad_curv** | 3 element vector or 2x3 array | mm | Input beam radius of curvature (in air) for red1, red2, blue waves. With an optional second row, the first row is radius of curvature in the direction of birefringent walk off and the second is perpendicular to it. |
| **mix_2d_sp_nz** | scalar | Integer | Number of mixing & propagation z steps through crystal
| **mix_2d_sp_nxny** | 2 element vector |  Integer | Number of transverse grid points, in the direction of birefringent walkoff and perpendicular to it. | 
| **mix_2d_sp_crystal_length** | scalar | mm | Crystal length. |
| **mix_2d_sp_lx_ly** | 2 element vector | Integer | Crystal width and height. | 
| **mix_2d_sp_deff** | scalar | pm/V | Effective nonlinear coefficient. | 
| **mix_2d_sp_deltak** | scalar | 1/mm | Phase mismatch Deltak = k_blue - k_red1 - k_red2. | 
| **mix_2d_sp_nt** | scalar | Integer | Number of time points. | 


# Outputs

<p id="outputs">Model outputs in the MATLAB file *mix_2d_sp_output.mat*</p>

| **Variable Name** | **Type** | **Description** |
| -- | -- | -- | 
|   field_red1_xyt_crystaloutput | (3D array) | Red1 wave output electric fields (units of V/m) at crystal output face (in air). NX x NY x NT sized array, specified at crystal output face in air. | 
|   field_red2_xyt_crystaloutput | (3D array) | Red2 wave output electric fields (units of V/m) at crystal output face (in air). NX x NY x NT sized array, specified at crystal output face in air. | 
|   field_blue_xyt_crystaloutput | (3D array) | Blue wave output electric fields (units of V/m) at crystal output face (in air). NX x NY x NT sized array, specified at crystal output face in air. | 
|   field_red1          | (3D array) | red1 wave electric fields (units of V/m) at plane of detector. NX x NY x NT sized array for each wave. This will be different from field_red1_xyt_crystaloutput if the 'Propagate' button has been used. | 
|   field_red2          | (3D array) | red2 wave electric fields (units of V/m) at plane of detector. NX x NY x NT sized array for each wave. This will be different from field_red2_xyt_crystaloutput if the 'Propagate' button has been used. | 
|   field_blue          | (3D array) | blue wave electric fields (units of V/m) at plane of detector. NX x NY x NT sized array for each wave. This will be different from field_blue_xyt_crystaloutput if the 'Propagate' button has been used. | 
|   xgrid               | (vector) | positions in x (meters), the direction of birefringent walk off. NX length vector. | 
|   ygrid               | (vector) | positions in y, the direction perpendicular to birefringent walk off. NY length vector. | 
|   xmat                | (2D meshgrid array) | positions in x (meters). NX x NY sized array, changing in 1st dimension and constant in the 2nd dimension. | 
|   ymat                | (2D meshgrid array) | positions in y (meters). NX x NY sized array, changing in 2nd dimension and constant in the 1st dimension. | 
|   tgrid               | (vector) | time (seconds). NT length vector. | 
|   zvec                | (vector) | positions in z (meters) through the crystal. NZ length vector. | 
|   fluences            | (3D array) |fluence (J/m^2) at crystal exit face. 3 x NX x NY sized array, 1st dimension is (red1,red2,blue), 2nd is x, 3rd is y. | 
|   farfield_fluences   | (3D array) | farfield fluences (arbitrary units). 3 x NX x NY sized array, 1st dimension is (red1,red2,blue), 2nd is kx, 3rd is ky. | 
|   fields_at_input     | (3D array) | data structure with fieldnames red1, red2, blue; each field is an array of electric field (V/m) in air prior to entering crystal. NX x NY x NT sized array, 1st dimension is x, 2nd is y, 3rd is t. | 
|   xBar_vs_t           | (2D array) | beam centroid (m) vs time in x direction (in birefringent walkoff direction) at crystal exit face. NT x 3 sized array, 1st dimension is time, 2nd dimension is (red1,red2,blue wave). | 
|   yBar_vs_t           | (2D array) | beam centroid (m) vs time in y direction (perpendicular to walkoff) at crystal exit face. NT x 3 sized array, 1st dimension is time, 2nd dimension is (red1,red2,blue wave). | 
|   wx_vs_t             | (2D array) | second moment width (m) in x direction at crystal exit face. NT x 3 sized array, 1st dimension is time, 2nd dimension is (red1,red2,blue wave). | 
|   wy_vs_t             | (2D array) | second moment width (m) in y direction at crystal exit face. NT x 3 sized array, 1st dimension is time, 2nd dimension is (red1,red2,blue wave). | 
|   chirp               | (2D array) | chirp (Hz) at crystal output. NT x 4 sized array, where the 1st dimension is time (seconds), the 2nd red1, red2, and blue wave chirp values. | 
|   xBar_vs_z           | (2D array) | centroid of fluence (meters) in x direction through the crystal. NZ x 3 sized array, 1st dimension is z position, 2nd dimension is (red1, red2, blue wave). | 
|   yBar_vs_z           | (2D array) | centroid of fluence (meters) in y direction through the crystal. NZ x 3 sized array, 1st dimension is z position, 2nd dimension is (red1, red2, blue wave). | 
|   wx_vs_z             | (2D array) | fluence 2nd moment width in x through the crystal (meters). NZ x 3 sized array, 1st dimension is z position, 2nd dimension is (red1, red2, blue wave). | 
|   wy_vs_z             | (2D array) | fluence 2nd moment width in y through the crystal (meters). NZ x 3 sized array, 1st dimension is z position, 2nd dimension is (red1, red2, blue wave). | 
|   fluence_xyz_red1    | (3D array) | fluences (J/m^2) for red1 wave at each x, y, z point in the crystal. NX x NY x NZ sized array, 1st dimension is x position, 2nd is y, 3rd is z. | 
|   fluence_xyz_red2    | (3D array) | fluences (J/m^2) for red2 wave at each x, y, z point in the crystal. NX x NY x NZ sized array, 1st dimension is x position, 2nd is y, 3rd is z. | 
|   fluence_xyz_blue    | (3D array) | fluences (J/m^2) for blue wave at each x, y, z point in the crystal. NX x NY x NZ sized array, 1st dimension is x position, 2nd is y, 3rd is z. | 
|   field_red1_spectrum | (vector) | spectrum (arbitrary units) for red1 wave at crystal output face. N_nu length vector, the variable freqs defines the grid for this spectra. N_nu might be larger than NT to pad fft and give smoother curves. | 
|   field_red2_spectrum | (vector) | spectrum (arbitrary units) for red2 wave at crystal output face. N_nu length vector, the variable freqs defines the grid for this spectra. N_nu might be larger than NT to pad fft and give smoother curves. | 
|   field_blue_spectrum | (vector) | spectrum (arbitrary units) for blue wave at crystal output face. N_nu length vector, the variable freqs defines the grid for this spectra. N_nu might be larger than NT to pad fft and give smoother curves. | 
|   freqs               | (vector) | Frequencies for the spectra (Hz). N_nu length vector. N_nu might be larger than NT to pad fft and give smoother curves. | 
|   kxarray             | (2D meshgrid array) |  angular freqeuency in x direction, for farfield irradiances. NX x NY sized array. | 
|   kyarray             | (2D meshgrid array) |  angular freqeuency in y direction, for farfield irradiances. NX x NY sized array. | 
|   energy_vs_z         | (2D array) | energy for each wave at each z position in crystal. NZ x 3 sized array, 1st dimension is z position (through crystal), 2nd dimension is (red1, red2, blue wave order.) | 
|   power               | (2D array) |  power (W) for each wave at crystal output face; NT x 4 sized array, 1st dimension is value at each time, and the first column in the 2nd dimension is time (seconds) and the second, third, and fourth red1, red2, and blue values. | 
