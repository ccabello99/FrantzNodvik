# Frantz-Nodvik equation solver in 2+1D

Julia code capable of simulating the amplification of a laser pulse in a multi-pass amplifier. It is currentlyconfigured to simulate amplification in titanium-doped sapphire. The emission and absorption cross-sections of Ti:Sapphire, and the nonlinear refrative index of sapphire are provided. The 2D spatial profiles for both the seed and pump beams can be defined by the user, and the temporal profile of the seed pulse can also be defined. 

## To run: 
1) Compile ParamScan.jl
2) Run:
```julia
ParamScan(fn_params, <initial_seed_waist>, <adjusted_seed_waist>, <initial_pump_waist>, <initial_pump_energy>, <initial_seed_energy>, <number_of_steps_for_each_parameter>, <convert to LG Mode (true/false)>, <visualize during scan (true/false)>)
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


## Data sources 
1) Cross-sections: Sorokin, E., "Solid-State Materials for Few-Cycle Pulse Generation and Amplification", 2004

2) Nonlinear refractive index: Major, A, et al., "Dispersion of the nonlinear refractive index in sapphire", 2004 
