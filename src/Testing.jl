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

@unpack c = fn_params

Z = ZernikeCoefficients(0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0)

diff_params = Diffract{Float64}(fn_params, 54.4e-3, 1.36, 19e-3, 1.0)

Ef, xf, yf = RichardsWolf(fn_params, diff_params, RHC(), 0e-6, 1, Z, aberration=false, hole=false, verbose=true);

#fig, ax, hm = GLMakie.heatmap(xf.*1e6, yf.*1e6, abs2.(Ef[2]), colormap=:viridis)
#ax.xlabel="x (μm)"
#ax.ylabel="y (μm)"
#fig

#If = abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3])
#If ./= maximum(If)

#Plots.heatmap(If)

#Plots.heatmap(abs2.(Ef[2]))

#fig = getPolarizationEllipse2D(xf.*1e6, yf.*1e6, Ef[1], Ef[2]; amplification=0.75, num_ellipses=(11, 11), line_width=0.75, draw_arrow=true)


#println(e22D(xf, yf, If) .* 1e6)

#DiffractionMovie(P(), "t", fn_params, diff_params, -15e-6, 15e-6, 129, 0, Z, aberration=false, intensity=true, hole=true, phase=false)

#E, x, y, z = SpatioTemporalVectorDiffraction(fn_params, diff_params, Radial(), -5e-6, 5e-6, 65, 65, 0, 0, Z, aberration=false, hole=false, spectdata=true)


#E, x, y, z = SpatioTemporalVectorDiffraction(fn_params, diff_params, Azimuthal(), -5e-6, 5e-6, 65, 65, 0, 1, "constant", Z; verbose=false, aberration=false, hole=false, spectdata=true)
