# 2D-mix-SP
The units are usually kept in MKS for simplicity and consistency.

[Inputs](#inputs) and [outputs](#outputs) below.

# Input parameters {#inputs}
2D-mix-SP *(snlo_2d_mix_sp_func)* accepts model parameters as a input argument to the MATLAB function file which is a data structure. The structure looks for matching fieldnames and tries to use the corresponding values to populate the edit boxes in the main GUI input form. Below is a list of the fieldnames with a remark about how the values should be formatted (and their units):
|**Fieldname** | **Shape** | **Units** | **Description** |
| -- | -- | -- | -- |
| **mix_2d_sp_wavelengths** | 3 element vector | nm | Wavelengths for red1, red2, blue waves. |
| **mix_2d_sp_ref_inds**    | 3 element vector |  | Crystal refractive indices for red1, red2, blue waves. |
| **mix_2d_sp_gvi**  		    | 3 element vector |  | Crystal group velocity indices for red1, red2, blue waves. |
| **mix_2d_sp_gdd** 		    | 3 element vector | fs^2/mm | Group delay dispersion for red1, red2, blue waves. |
| **mix_2d_sp_phase** 		  | 3 element vector | radians | Input phase for red1, red2, blue waves. |
| **mix_2d_sp_input_refl**  | 3 element vector |  | Crystal input face  reflectiviy for red1, red2, blue waves. | 
| **mix_2d_sp_output_refl** | 3 element vector |  | Crystal output face reflectivity for red1, red2, blue waves. | 
| **mix_2d_sp_crystal_losses** | 3 element vector | % crystal loss (1/mm)
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



# Outputs {#outputs}
Model outputs in the MATLAB file *mix_2d_sp_output.mat*
*   **field_red1_xyt_crystaloutput**, **field_red2_xyt_crystaloutput**, **field_blue_xyt_crystaloutput** 	: (3D arrays) Output electric fields (units of V/m), in x,y, and t dimensions, for each wave at the crystal output face after exiting the crystal (in air).
*   **field_red1**, **field_red2**, **field_blue** : (3D arrays) electric fields (units of V/m), in three dimensions (x,y, and t), for each wave. This will be different from field_red1_xyt_crystaloutput etc if the 'Propagate' button has been used.
*   **xgrid** 			: (vector) positions in x (meters), the direction of birefringent walk off. NX length vector.
*   **ygrid** 			: (vector) positions in y, the direction perpendicular to birefringent walk off. NY length vector.
*   **xmat** 				: (2D meshgrid array) positions in x (meters), NX x NY sized 2D array, changing in 1st dimension and constant in the 2nd dimension.
*   **ymat**				: (2D meshgrid array) positions in y (meters), NX x NY sized 2D array, changing in 2nd dimension and constant in the 1st dimension
*   **tgrid** 			: (vector) time (seconds), length NT  
*   **zvec** 				: (vector) positions in z (meters) through the crystal, NZ length vector
*   **fluences** 			: (3D array) fluence (J/m^2) at crystal exit face; array is sized 3 x NX x NY (1st dimension is red1,red2,blue, 2nd is x position, 3rd is y position).
*   **farfield_fluences** : (3D array) farfield fluences (arbitrary units); array is sized 3 x NX x NY (1st dimension is red1,red2,blue, 2nd is kx, 3rd is ky)
*   **fields_at_input** 	: (3D array) data structure with fieldnames red1, red2, blue; each field is an array of electric field (V/m) in air prior to entering crystal. NX x NY x NT sized 3D array
*   **xBar_vs_t** 		: (2D array) beam centroid (m) vs time in x direction (in birefringent walkoff direction) at crystal exit face; array is sized NT x 3 (1st dimension is time, 2nd dimension is red1,red2,blue wave)
*   **yBar_vs_t** 		: (2D array) beam centroid (m) vs time in y direction (perpendicular to walkoff) at crystal exit face; array is sized NT x 3 (1st dimension is time, 2nd dimension is red1,red2,blue wave)
*   **wx_vs_t**   		: (2D array) second moment width (m) in x direction at crystal exit face; array is sized NT x 3 (1st dimension is time, 2nd dimension is red1,red2,blue wave)
*   **wy_vs_t** 			: (2D array) second moment width (m) in y direction at crystal exit face; array is sized NT x 3 (1st dimension is time, 2nd dimension is red1,red2,blue wave)
*   **chirp**				: (2D array) chirp (Hz) at crystal output; array is sized NT x 4 where the 1st elements in the 2nd dimension are time (seconds), red1, red2, and blue wave values.
*   **xBar_vs_z**			: (2D array) centroid of fluence (meters) in x direction through the crystal, array is sized NZ x 3 (1st dimension is z position, 2nd dimension is red1, red2, blue wave)
*   **yBar_vs_z** 		: (2D array) centroid of fluence (meters) in y direction through the crystal, array is sized NZ x 3 (1st dimension is z position, 2nd dimension is red1, red2, blue wave)
*   **wx_vs_z** 			: (2D array) fluence 2nd moment width in x through the crystal (meters), array is sized NZ x 3 (1st dimension is z position, 2nd dimension is red1, red2, blue wave)
*   **wy_vs_z** 			: (2D array) fluence 2nd moment width in y through the crystal (meters), array is sized NZ x 3 (1st dimension is z position, 2nd dimension is red1, red2, blue wave)
*   **fluence_xyz_red1** 	: (3D array) fluences (J/m^2) for red1 wave at each x, y, z point in the crystal. Array is sized NX x NY x NZ.
*   **fluence_xyz_red2** 	: (3D array) fluences (J/m^2) for red2 wave at each x, y, z point in the crystal. Array is sized NX x NY x NZ.
*   **fluence_xyz_blue** 	: (3D array) fluences (J/m^2) for blue wave at each x, y, z point in the crystal. Array is sized NX x NY x NZ.
*   **field_red1_spectrum** 	: (vector) spectrum (arbitrary units) for red1 wave at crystal output face. Vector is length NT, the variable freqs defines the grid for this spectra.
*   **field_red2_spectrum**	: (vector) spectrum (arbitrary units) for red2 wave at crystal output face. Vector is length NT, the variable freqs defines the grid for this spectra.
*   **field_blue_spectrum**	: (vector) spectrum (arbitrary units) for blue wave at crystal output face. Vector is length NT, the variable freqs defines the grid for this spectra.
*   **freqs** 			: (vector) Frequencies for the spectra (Hz)
*   **kxarray** 			: (2D array) angular freqeuency in x direction, for farfield irradiances; NXxNY array
*   **kyarray** 			: (2D array) angular freqeuency in y direction, for farfield irradiances; NXxNY array
*   **energy_vs_z** 		: (2D array) energy for each wave at each z position in crystal; array is sized NZ x 3; 2nd dimension is red1, red2, blue wave order.
*   **power** 			: (2D array) power (W) for each wave at crystal output face; array is sized NT x 4, where the first elements in the 2nd dimension are time (seconds), red1, red2, and blue values.
