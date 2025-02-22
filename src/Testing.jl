using CSV, Tables, DataFrames, Printf, JLD, Plots, LaTeXStrings, CairoMakie, GLMakie
using LinearAlgebra, FFTW, ForwardDiff, NumericalIntegration, DSP
using ThreadsX, Parameters, Dierckx, ClassicalOrthogonalPolynomials, SpecialFunctions

include("FN-Params.jl")
include("Gabor.jl")
include("CrystalProperties.jl")
include("EM-Field.jl")
include("Passes.jl")
include("Diffraction.jl")
include("Polarization.jl")
include("Helpers.jl")

@unpack h, c = fn_params

Z = ZernikeCoefficients(0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0)

#fn_params = FN_Params{Float64}(xmax=100e-3, ymax=100e-3, x0=25.5e-3, y0=0)

fn_params = FN_Params{Float64}(xmax=500e-3, ymax=500e-3, x0=200e-3, y0=0)

diff_params = Diffract{Float64}(fn_params, 26.795e-3, 10.7, 12.5e-3, 1.0)
#diff_params = Diffract{Float64}(fn_params, 51.25e-3, 1.3, 12.5e-3, 1.0)
#diff_params = Diffract{Float64}(fn_params, 26.975e-3, 1.3, 25e-3, 1.0)

Pol = P()
z0 = 0e-6
l = 0

aperture = elliptical_aperture(fn_params, 97e-3/2, 25e-3/2)

#Ef, xf, yf = RichardsWolf(fn_params, diff_params, Pol, z0, l, Z, aberration=false, hole=false, verbose=true);
#Hf, xf, yf = RichardsWolf(fn_params, diff_params, Pol, z0, l, Z, aberration=false, hole=false, verbose=false, magnetic=true);

#Ef, xf, yf, zf = FullSpatialProfile(fn_params, diff_params, Pol, -10e-6, 10e-6, 129, l, Z)
#Hf, xf, yf, zf = FullSpatialProfile(fn_params, diff_params, Pol, -10e-6, 10e-6, 129, l, Z; magnetic=true)

#fig, ax, hm = GLMakie.heatmap(xf.*1e6, yf.*1e6, abs2.(Ef[1]), colormap=:viridis)
#ax.xlabel="x (μm)"
#ax.ylabel="y (μm)"
#fig

#fig, ax, hm = GLMakie.heatmap(xf.*1e6, yf.*1e6, real.(Ef[1]), colormap=:thermometer)
#ax.xlabel="x (μm)"
#ax.ylabel="y (μm)"
#fig

#fig = getPolarizationEllipse2D(xf.*1e6, yf.*1e6, Ef[1], Ef[2]; amplification=0.75, num_ellipses=(11, 11), line_width=0.75, draw_arrow=true)

#DiffractionMovie(P(), "t", fn_params, diff_params, -15e-6, 15e-6, 129, 0, Z, aberration=false, intensity=true, hole=true, phase=false)


#E, x, y, z = FullSpatioTemporalVectorDiffraction(fn_params, diff_params, P(), -5e-6, 5e-6, 65, 65, 0, 0, Z, aberration=false, hole=false, spectdata=true, harmonic=false)
#E, x, y, z = SpatioTemporalVectorDiffraction(fn_params, diff_params, Radial(), -5e-6, 5e-6, 65, 65, 0, 0, Z, aberration=false, hole=false, spectdata=true)
#E, x, y, z = SpatioTemporalVectorDiffraction(fn_params, diff_params, Azimuthal(), -5e-6, 5e-6, 65, 65, 0, 1, "constant", Z; verbose=false, aberration=false, hole=false, spectdata=true)
