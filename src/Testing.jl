using CSV, Tables, DataFrames, Printf, JLD, Plots, LaTeXStrings, CairoMakie, GLMakie
using LinearAlgebra, FFTW, ForwardDiff, NumericalIntegration, DSP, CUDA, BenchmarkTools, FourierTools
using ThreadsX, Parameters, Dierckx, ClassicalOrthogonalPolynomials, SpecialFunctions

include("FN-Params.jl")
include("Gabor.jl")
include("CrystalProperties.jl")
include("EM-Field.jl")
include("Passes.jl")
include("Diffraction.jl")
include("Polarization.jl")
include("Helpers.jl")
include("Utils.jl")

@unpack h, c = fn_params

Z = ZernikeCoefficients(0, 0, 0, 0, 0, 0, 0, 0.0, 0, 0, 0)

#x0=25.5e-3
fn_params = FN_Params{Float64}(xmax=30e-4, ymax=30e-4, x0=6e-4, y0=0e-4, N=2^8)

#fn_params = FN_Params{Float64}(xmax=100e-3, ymax=100e-3, N = 2^8+1)
#diff_params = Diffract{Float64}(fn_params, 400e-3, 16, 0.5e-3, 1.0)

Pol = Radial()
f = 51.25e-3
#f = 54.4e-3
#f = 400e-3
fnum = 1.3
#w = 19e-3
w = 5e-4
n = 1.0
z0 = 0e-6
l = 0

diff_params = Diffract{Float64}(fn_params, f, fnum, w, n)


for i in 0:120
    z = -0.5e-3*i

    Ex, Ey, Ez = TransmissionFunction(fn_params, diff_params, Pol, l, Z, aberration=false, hole=false, verbose=false, magnetic=false, OAP=true)

    Etx, Ety, Etz = AngularSpectrumPropagator(fn_params, diff_params, [Ex, Ey, Ez], z, lens=false, verbose=false, cuda=true, OAP=true)

    I = abs2.(Etx) .+ abs2.(Ety) .+ abs2.(Etz)
    
    
    #display(Plots.heatmap(fn_params.x.*1e3, fn_params.y.*1e3, I, title="z = $(z*1e3) mm", xlabel="x (mm)", ylabel="y (mm)", color=:inferno))
    #display(Plots.heatmap(fn_params.x.*1e3, fn_params.y.*1e3, real.(Etx), title="z = $(z*1e3) mm", xlabel="x (mm)", ylabel="y (mm)", color=:inferno))
end

#E, x, y = Bluestein(fn_params, diff_params, 100e-3, fn_params.x, fn_params.y, 0, focusing=false)

#XUVSTVD(fn_params, diff_params, Pol, -15e-6, 15e-6, 65, 64, l, Z, save=true, aberration=false, hole=false)
#E, x, y, z, ν_samples = XUVSTVD(fn_params, diff_params, Pol, -5e-6, 5e-6, 65, 64, l, Z, save=false)

#Ef, xf, yf = RichardsWolf(fn_params, diff_params, Pol, z0, l, Z, aberration=false, hole=false, verbose=false, OAP=true);
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

#fig = getPolarizationEllipse2D(xf.*1e6, yf.*1e6, Hf[1], Hf[2]; amplification=0.75, num_ellipses=(21, 21), line_width=0.75, draw_arrow=true)
#fig
#DiffractionMovie(P(), "t", fn_params, diff_params, -15e-6, 15e-6, 129, 0, Z, aberration=false, intensity=true, hole=true, phase=false, OAP=true)


#E, x, y, z = FullSpatioTemporalVectorDiffraction(fn_params, diff_params, P(), -5e-6, 5e-6, 65, 65, 0, 0, Z, aberration=false, hole=false, spectdata=true, harmonic=false)
#E, x, y, z = SpatioTemporalVectorDiffraction(fn_params, diff_params, Radial(), -5e-6, 5e-6, 65, 65, 0, 0, Z, aberration=false, hole=false, spectdata=true)
#E, x, y, z = SpatioTemporalVectorDiffraction(fn_params, diff_params, Azimuthal(), -5e-6, 5e-6, 65, 65, 0, 1, "constant", Z; verbose=false, aberration=false, hole=false, spectdata=true)

#@unpack ϕ = diff_params
#Ex = -sin.(ϕ) .* LaguerreGauss(fn_params, 0, 1, 1, w)
#Ey = cos.(ϕ) .* LaguerreGauss(fn_params, 0, 1, 1, w)

#Epx, Epy, Epz = Polarization(fn_params, diff_params, 1, Pol, Z)

#Etx, Ety, Etz = TransmissionFunction(fn_params, diff_params, Pol, 1, Z, aberration=false, hole=false, verbose=false, magnetic=false, OAP=true)

#fig = getPolarizationEllipse2D(x, y, Epx, Epy; amplification=0.75, num_ellipses=(21, 21), line_width=0.75, draw_arrow=true)
#display(fig)
