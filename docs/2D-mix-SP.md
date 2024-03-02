## 2D-mix-SP
The units are usually kept in MKS for simplicity and consistency.

## For inputs, use a data structure with fieldnames specified in [mlSNLO_function_input_data_structures.md](mlSNLO_function_input_data_structures.md)

## Model outputs in the MATLAB file *mix_2d_sp_output.mat*
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
