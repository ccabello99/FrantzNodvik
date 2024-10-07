# Frantz-Nodvik equation solver in 2+1D

Julia code capable of simulating the amplification of a laser pulse in a multi-pass amplifier. It is currentlyconfigured to simulate amplification in titanium-doped sapphire. The emission and absorption cross-sections of Ti:Sapphire, and the nonlinear refrative index of sapphire are provided. The 2D spatial profiles for both the seed and pump beams can be defined by the user, and the temporal profile of the seed pulse can also be defined. 

## To run: 
1) Compile ParamScan.jl
2) Run ParamScan(fn_params, <initial_seed_waist>, <adjusted_seed_waist>, <initial_pump_waist>, <initial_pump_energy>, <initial_seed_energy>, <number_of_steps_for_each_parameters>, <convert to LG Mode (true/false)>, <visualize during scan (true/false)>)

## For your information: 

i) Provide values in SI units

ii) For given number of steps, each scan will increment waist sizes by 200 μm, pass number of seed waist adjustment by 1, pump energy by 5 mJ, and seed energy by 600 μJ. Each scan parameter will be saved in a vector in .csv format.

iii) Multi-dimensional arrays will be saved after each scan in .jld format. These include:
1) Seed energy for each pass
2) Extraction efficiency for each pass
3) B-integral for each pass
4) Effective area for each pass
5) Maximum energy after all passes
6) Pass number where maximum energy is reached
7) Maximum B-integral after all passes
8) Maximum effective area after all passes


## Data sources 
Cross-sections: Sorokin, E., "Solid-State Materials for Few-Cycle Pulse Generation and Amplification", 2004

Nonlinear refractive index: Major, A, et al., "Dispersion of the nonlinear refractive index in sapphire", 2004 
