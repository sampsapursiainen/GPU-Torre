# GPU-ToRRe
GPU-accelerated Tomographic Radar Reconstruction Toolbox

GPU-ToRRe is a numerical toolbox for investigating mathematical non-linear tomographic radar inverse imaging problems. GPU-ToRRe aims at enabling fast and robust inversion of sparse full-wave tomography data in 2D. The forward approach utilizes a higher-order Born approximation (BA) to account for complex wave scattering inside the computational domain. The forward modelling approach uses a graphics processing unit (GPU) accelerated finite element time-domain (FETD) method. In the inversion stage, a GPU-computed multigrid-FETD deconvolution routine is applied to enhance computational performance. 

GPU-ToRRe toolbox has been developed for advancing the mathemati tomography of small solar system bodies (SSSB), and especially asteroid interiors which can contain internal details, for example layers, voids, and cracks, observable by a radar. The example 2D model included in this toolbox relates to a SSSB with three internal voids enclosed by a surface layer. 

To run the full simulation, follow these steps

1. Create the background and exact systems by running the scripts create_system_background.m and create_system_exact.m, respectively.

2. Compute the forward simulations for the background and the exact systems by running the scripts compute_field_data_background.m, and compute_field_data_exact.m. These are required for linear inversion. 

3. To run the non-linear inversion with higher-order BA, run the script compute_field_data_correction.m to update the wave-field corrections due to higher order scattering. 

4. Run the script load_field_data.m to load wave-field data for the inversion stage.

5. To invert the data, run the script inversion_procedure.m.

All the parameters for forward and inversion computations are stored in the file parameters.m.

A shortcut to the inversion stage: We have uploaded the system files and data from steps 1 and 2 to the folders data, field_data_background, and field_data_exact so it is possible to start directly from the step 3.
