
function meshgrid(x::Vector, y::Vector)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function LaguerreGauss(params::FN_Params, P::Int, L::Int, A::Real, W::Real)
    
    # Laguerre-Gauss equation: 
    # (ref: N. Hodgson, 'Laser Resonators and Beam Propagation'.(Pg 222)) 

    @unpack N, dx, dy, x0, y0 = params
    
    t = zeros(N, N)
    Phi = zeros(N, N)
    Term1 = zeros(N, N)
    W2 = W^2
    sqrt2 = sqrt(2)
    
    
    ThreadsX.foreach(CartesianIndices((N, N)); simd=true) do I
        x = I[1] * dx
        y = I[2] * dy
        x_diff = x - x0
        y_diff = y - y0
        x_diff2 = x_diff^2
        y_diff2 = y_diff^2
        
        @inbounds t[I] = -(2 * x_diff2 / W2) - 2 * y_diff2 / W2
        @inbounds Phi[I] = L * atan(y_diff, x_diff)
        @inbounds Term1[I] = (sqrt2 * sqrt(x_diff2 + y_diff2))^L
    end
        
    C = A * sqrt(2*factorial(P)/(π*factorial(P+abs(L))))
    Term2 = laguerrel.(P, L, 2 .* t)
    Term3 = exp.(t)
    Term4 = exp.(1im .* Phi)
    Z = C .* (1 / W) .* Term1 .* Term2 .* Term3 .* Term4

    return Z
end


# wx and wy are the beam diameter along x and y direction
function Gaussian(params::FN_Params, wx::Real, wy::Real)

    @unpack N, dx, dy, x0, y0 = params

    gauss = zeros(N, N)
    wx2 = wx^2
    wy2 = wy^2

    ThreadsX.foreach(CartesianIndices((N, N)); simd=true) do I
        x = I[1] * dx
        y = I[2] * dy
        x_diff = x - x0
        y_diff = y - y0

        @inbounds gauss[I] = exp(-2 * (x_diff^2 / wx2) - 2 * (y_diff^2 / wy2))
    end

    return gauss

end


function SuperGaussian(params::FN_Params, w::Real, nsg::Int)

    @unpack N, dx, dy, x0, y0 = params

    super_gaussian = zeros(N, N)
    w2 = w^2

    ThreadsX.foreach(CartesianIndices((N, N)); simd=true) do I
        x = I[1] * dx
        y = I[2] * dy
        x_diff = x - x0
        y_diff = y - y0

        @inbounds super_gaussian[I] = exp(-2 * ((x_diff^2 + y_diff^2) / w2)^nsg)
    end


    return super_gaussian

end

function calcAeff(x::Vector, y::Vector, J::Matrix)

    Aeff = (NumericalIntegration.integrate((x,y), J)).^2 / NumericalIntegration.integrate((x,y), J.^2)

    return Aeff
end

# Electric field temporal profile  || if you want to input τ @ FWHM, scale by [sqrt(1/(-2*log(0.5)))*sqrt(1/2)]
function ComplexEnvelope(A0::Real, t0::Real, ϕ::Real, τ::Real, GDD::Real)
    C = sqrt(1 + (GDD / τ^2)^2)
    ϕ_σ = (1/2) * atan(GDD / τ^2)
    A(T) = (A0 / sqrt(C)) .* exp.(1im .* (ϕ .+ ϕ_σ)) .* exp.(-(T.-t0).^2 / (2 .* (C .* τ)^2)) .* exp.(-1im .* (GDD / τ.^2) .* ((T.-t0).^2 ./ (2 .* (C .* τ)^2)))
    return A
end


println("EM-Field.jl compiled")

#=

# To visualize LG-modes (spatio-temporally if wanted)

N = 100
xmax = 20e-4
ymax = xmax
dx = xmax / N
dy = dx

x = range(0, xmax, N) .* 1e3
y = range(0, ymax, N) .* 1e3
z = range(0, xmax, N) .* 1e3
x0 = xmax / 2
y0 = ymax / 2

p = 0
l = 1
a = 1
w = 2000e-6

# Laguerre-Gaussian spatial properties
Z = LaguerreGauss(fn_params, p, l, a, w)

E = real.(Z)
ϕ = angle.(Z)
Ixy = abs.(Z)

scale = sqrt(1/(-2*log(0.5)))*sqrt(1/2)
t = collect(range(-30e-15,30e-15,N))
At = ComplexEnvelope(1, 0, 0, 30e-15*scale, 0)
A_t = At.(t)
ω = 2.998e8*2π/800e-9

full_field = zeros(ComplexF64, N, N, N)

for I in CartesianIndices((N, N, N))
    i, j, k = (I[1], I[2], I[3])
    full_field[i, j, k] = Z[j, k] * exp(1im * ω * t[i]) * A_t[i]
end

Exyt = real.(full_field)
Ixyt = abs.(full_field)    
ϕxyt = angle.(full_field)

using GLMakie

scene = GLMakie.volume(t.*1e15, x.*20, y.*20, Exyt, algorithm = :iso, isorange = 0.1, isovalue = maximum(Exyt[50,:,:])*0.7)
GLMakie.volume!(t.*1e15, x.*20, y.*20, Exyt, algorithm = :iso, isorange = 0.1, isovalue = minimum(Exyt[50,:,:])*0.7)
#scene = GLMakie.volume(t.*1e15, x.*20, y.*20, Ixyt, colormap = :inferno, algorithm = :iso, isorange = 0.1, isovalue = maximum(Ixyt[50,:,:])*0.7)
display(scene)
=#


#=

# To visualize beam and pump spatial profiles

@unpack E0_p, A, ηc, ηq, w_xp, w_yp, w_xs, w_ys, N, x, y = fn_params
w_xp = 800e-6 # at tube window
supergauss = SuperGaussian(fn_params, w_xp, 6)
scale = NumericalIntegration.integrate((x,y), Gaussian(fn_params, w_xp, w_yp)) / NumericalIntegration.integrate((x,y), supergauss)
Jsto = scale .* supergauss
Eabs = ((E0_p) * (A * ηc * ηq)) * (2 - (A * ηc * ηq))
Aeff_p = calcAeff(x, y, Jsto)
Jsto0 = Eabs / Aeff_p
Jsto = Jsto0 .* Jsto
scaleJsto = Eabs / NumericalIntegration.integrate((x,y), Jsto)
Jsto = scaleJsto * Jsto

seed = Gaussian(fn_params, 1050e-6, 1050e-6)
Ein0 = 1.2e-3
Jin0 = Ein0 / calcAeff(x, y, seed)
seed .*= Jin0

using CairoMakie

fig = Figure(fontsize = 48, size=(1920, 1080))

#=
ax1 = Axis(fig[1, 1], 
	xlabel = L"\textbf{x (mm)}", ylabel = L"\textbf{y (mm)}",
	title = L"\textbf{LG}_{1, 0} \textbf{Intensity Profile}",
	ylabelpadding = 20)
=#


ax1 = Axis(fig[1, 1], limits=(2, 6, 2, 6),
	xlabel = L"\textbf{x (mm)}", ylabel = L"\textbf{y (mm)}",
	title = L"\textbf{Pump Beam Fluence Profile}",
	ylabelpadding = 20)
ax2 = Axis(fig[1, 3], limits=(2, 6, 2, 6),
	xlabel = L"\textbf{x (mm)}", ylabel = L"\textbf{y (mm)}",
	title = L"\textbf{Seed Beam Fluence Profile}",
	ylabelpadding = 20)

#=
ax1 = Axis(fig[1, 1], limits=(2, 6, 2, 6),
	xlabel = L"\textbf{x (mm)}", ylabel = L"\textbf{y (mm)}",
	title = L"\textbf{Seed Beam Fluence Profile}",
	ylabelpadding = 20)
=#
#=
ax2 = Axis(fig[1, 1],
	xlabel = L"\textbf{x (mm)}", ylabel = L"\textbf{y (mm)}",
	title = L"\textbf{LG}_{1, 0} \textbf{ Phase Profile}",
	ylabelpadding = 20)
=#

hm1 = CairoMakie.heatmap!(ax1, x .* 1e3, y .* 1e3, Jsto ./ maximum(Jsto), colormap = :inferno) #vikO100
hm2 = CairoMakie.heatmap!(ax2, x .* 1e3, y .* 1e3, seed ./ maximum(seed), colormap = :inferno) #vikO100
#hm2 = CairoMakie.heatmap!(ax2, x .* 1e3, y .* 1e3, (Jsto ./ maximum(Jsto)) .- (seed  ./ maximum(seed)), colormap = :inferno) #vikO100
cb = CairoMakie.Colorbar(fig[1, 2], hm1, size=20,)
cb = CairoMakie.Colorbar(fig[1, 4], hm2, size=20)

#hm = CairoMakie.heatmap!(ax1, x .* 1e3, y .* 1e3, Ixy, colormap = :vikO100)
#cb = CairoMakie.Colorbar(fig[1, 2], hm, size=20, label=L"\textbf{Normalized Intensity (arb. u)}")
#hm = CairoMakie.heatmap!(ax2, x, y, ϕ, colormap = :jet)
#cb = CairoMakie.Colorbar(fig[1, 2], hm, size=20, ticks=([-3.13, 0, 3.14], [L"-\pi", L"0", L"\pi"]))

fig

=#
