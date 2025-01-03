# Frantz-Nodvik equation solver in 2+1D

Julia code capable of simulating the amplification of a laser pulse in a multi-pass amplifier. It is currently configured to simulate amplification in titanium-doped sapphire. The emission and absorption cross-sections of Ti:Sapphire, and the nonlinear refrative index of sapphire are provided. The 2D spatial profiles for both the seed and pump beams can be defined by the user, and the temporal profile of the seed pulse can also be defined. 

## To run: 
1) Compile ParamScan.jl
2) Run:
```julia
ParamScan(fn_params, <initial_seed_waist>, <adjusted_seed_waist>, <initial_pump_waist>, <initial_pump_energy>, <initial_seed_energy>, <number_of_steps_for_each_parameter>, <convert to LG Mode (true/false)>, <visualize during scan (true/false)>);
```

## For your information: 

- Provide values in SI units

- For given number of steps, each scan will increment waist sizes by 200 μm, pass number of seed waist adjustment by 1, pump energy by 5 mJ, and seed energy by 600 μJ. Each scan parameter will be saved in a vector in .csv format.

- Multi-dimensional arrays will be saved after each scan in .jld format. These include:
  - Seed energy for each pass
  - Extraction efficiency for each pass
  - B-integral for each pass
  - Effective area for each pass
  - Maximum energy after all passes
  - Pass number where maximum energy is reached
  - Maximum B-integral after all passes
  - Maximum effective area after all passes
 
# Vector Diffraction

This code is also capable of performing vector diffraction by solving the Richards-Wolf integral via FFTs. It currently supports polarization types: P, S, 45$$\deg$$, left-handed circular, right-handed circular, radial, and azimuthal. 

## To run:

1) Compile Testing.jl, you can adjust diffraction parameters in this file, such as focal length, f-number, initial beam waist, and refractive index, and Zernike polynomial coefficients to account for aberrations.
2) Run:
```julia
E, x, y = RichardsWolf(fn_params, diff_params, <Polarization_Type()>, <z_coordinate>, <OAM_index>, Z, aberration=(true/false), hole=(true/false));
```
- if you would like to see the 2D field distribution at a specific z-coordinate, e.g. z=0.

```julia
E, x, y, z = SpatioTemporalVectorDiffraction(fn_params, diff_params, <Polarization_Type()>, <z_min>, <z_max>, <z_steps>, <freq_steps>, <time_coordinate>, <OAM_index>, Z, aberration=(true/false), hole=(true/false), spectdata=(true/false))
```
- if you would like to get the 3D field distribution for a specific input pulse at a specific time, e.g. t=0. If spectdata=true, it will look for spectral data in the input_data folder. If not, it will create a Gaussian spectral profile with center wavelength defined in fn_params, and a total spectral phase which is defined by the user using
```julia
phi = SpectralPhase(<phi_0>, <phi_1>, <phi_2>, <phi_3>, <phi_4>, <freq_sample_vector>, <center_freq>)
```
- where spectral phases are given in units of fs, fs^2, fs^3, and fs^4.

There are visualization functions defined in Helpers.jl to look at slices of the 3D distribution. Using GLMakie.jl, it is also possible to look at isosurfaces of the 3D distributions. Some examples of volumetric plotting are also given in Helpers.jl.


## Data sources 
1) Cross-sections: Sorokin, E., "Solid-State Materials for Few-Cycle Pulse Generation and Amplification", 2004

2) Nonlinear refractive index: Major, A, et al., "Dispersion of the nonlinear refractive index in sapphire", 2004

3) Spectral data: DSCAN of post-compressed Salle Noire 2 beamline
