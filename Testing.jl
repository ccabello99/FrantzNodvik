using CSV, Tables, DataFrames, Printf, JLD, Plots
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

Z = ZernikeCoefficients(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)

diff_params = Diffract{Float64}(fn_params, 54.4e-3, 1.36, 19e-3, 1.0)

#Ef, xf, yf = RichardsWolf(fn_params, diff_params, S(), 0e-6, 0, Z, aberration=true, hole=false);

#If = abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3])

#println(e22D(xf, yf, If) .* 1e6)

#DiffractionMovie(P(), "t", fn_params, diff_params, -15e-6, 15e-6, 129, 0, Z, aberration=true, intensity=true, hole=false)

E, x, y, z = SpatioTemporalVectorDiffraction(fn_params, diff_params, P(), -10e-6, 10e-6, 65, 65, 0, 0, Z, aberration=false, hole=false, spectdata=true)