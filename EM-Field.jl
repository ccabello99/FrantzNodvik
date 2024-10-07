
function meshgrid(x::Vector, y::Vector)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function LaguerreGauss(params::FN_Params, P::Int, L::Int, A::Real, W::Real)

    # Define polarization (to do)
    #k = π/2  
    #ω = π/25  
    #E0 = 1 
    #Ex(z, t) = E0 .* real.(exp(im*(k .* z .- ω .* t))) .* exp(-(z).^2 ./ 6)
    #Ey1(z, t) = E0 .* real.(exp(im*(k .* z .- ω .* t - π/2))) .* exp(-(z).^2 ./ 6)
    #Ey2(z, t) = E0 .* real.(exp(im*(k .* z .- ω .* t + π/2))) .* exp(-(z).^2 ./ 6)

    
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

# Electric field temporal profile
function ComplexEnvelope(A0::Real, t0::Real, ϕ::Real, τ::Real, GDD::Real)
    C = sqrt(1 + (GDD / τ^2)^2)
    ϕ_σ = (1/2) * atan(GDD / τ^2)
    A(T) = (A0 / sqrt(C)) .* exp.(1im .* (ϕ .+ ϕ_σ)) .* exp.(-(T.-t0).^2 / (2 .* (C .* τ)^2)) .* exp.(-1im .* (GDD / τ.^2) .* ((T.-t0).^2 ./ (2 .* (C .* τ)^2)))
    return A
end


println("EM-Field.jl compiled")


#=

N = 250
xmax = 20e-4
ymax = xmax
dx = xmax / N
dy = dx

x = range(0, xmax, N) .* 1e3
y = range(0, ymax, N) .* 1e3
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

# Plots

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
