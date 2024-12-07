
function meshgrid(x::Vector, y::Vector)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function LaguerreGauss(params::FN_Params, P::Int, L::Int, A::Real, W::Real)
    
    # Laguerre-Gauss equation: 
    # (ref: N. Hodgson, 'Laser Resonators and Beam Propagation'.(Pg 222)) 

    @unpack N, x, y, x0, y0 = params
    
    t = zeros(Float64, N, N)
    Phi = zeros(Float64, N, N)
    Term1 = zeros(Float64, N, N)
    Z = zeros(ComplexF64, N, N)
    W2 = W^2
    
    X, Y = meshgrid(x, y)

    x_diff = X .- x0
    y_diff = Y .- y0
    x_diff2 = x_diff.^2
    y_diff2 = y_diff.^2

            
    t .= -(x_diff2 ./ (W2)) .- (y_diff2 ./ (W2))
    Phi .= L .* atan.(y_diff, x_diff)
    Term1 .= (sqrt(2) .* sqrt.(x_diff2 .+ y_diff2)).^L
    C = A * sqrt(2*factorial(P)/(π*factorial(P+abs(L))))
    Term2 = laguerrel.(P, L, 2 .* t)
    Term3 = exp.(t)
    Term4 = exp.(1im .* Phi)
    Z .= C .* (1 / W) .* Term1 .* Term2 .* Term3 .* Term4

    return Z
end


# wx and wy are the beam diameter along x and y direction
function Gaussian(params::FN_Params, wx::Real, wy::Real)

    @unpack N, x, y, x0, y0 = params

    gauss = zeros(N, N)
    wx2 = wx^2
    wy2 = wy^2

    X, Y = meshgrid(x, y)
    x_diff = X .- x0
    y_diff = Y .- y0
    x_diff2 = x_diff.^2
    y_diff2 = y_diff.^2

    gauss .= exp.(-(x_diff2 ./ (wx2)) .- (y_diff2 ./ (wy2)))

    return gauss

end


function SuperGaussian(params::FN_Params, w::Real, nsg::Int)

    @unpack N, dx, dy, x0, y0 = params

    super_gaussian = zeros(N, N)
    w2 = w^2

    X, Y = meshgrid(x, y)
    x_diff = X .- x0
    y_diff = Y .- y0
    x_diff2 = x_diff.^2
    y_diff2 = y_diff.^2

    super_gaussian .= exp.(-(x_diff2 .+ y_diff2) ./ (w2).^(2*nsg))

    return super_gaussian

end

function calcAeff(x::Vector, y::Vector, E::Matrix)

    Aeff = ((NumericalIntegration.integrate((x,y), abs2.(E))).^2) / (NumericalIntegration.integrate((x,y), abs2.(E).^2))

    return Aeff
end

# Electric field temporal profile  || if you want to input τ @ FWHM, scale by [sqrt(1/(-2*log(0.5)))*sqrt(1/2)]
function ComplexEnvelope(A0::Real, t0::Real, ϕ::Real, τ::Real, GDD::Real)
    C = sqrt(1 + (GDD / τ^2)^2)
    ϕ_σ = (1/2) * atan(GDD / τ^2)
    A(T) = (A0 / sqrt(C)) .* exp.(1im .* (ϕ .+ ϕ_σ)) .* exp.(-(T.-t0).^2 / (2 .* (C .* τ)^2)) .* exp.(-1im .* (GDD / τ.^2) .* ((T.-t0).^2 ./ (2 .* (C .* τ)^2)))
    return A
end

function Gauss3D(fn_params::FN_Params, w0::Real)

    @unpack x, y, z, x0, y0,dx, dy, N, λs = fn_params

    k = 2π/λs
    zR = π * w0^2 / λs
    w(zz) = w0 .* sqrt.(1 .+ (zz ./ zR).^2)
    R(zz) = zz .* (1 .+ (zR ./ zz)^2)
    E0 = zeros(ComplexF64, N, N, N) 

    lk = Threads.ReentrantLock()

    ThreadsX.foreach(CartesianIndices((N, N, N)); simd=true) do I
        xx = I[1] * dx
        yy = I[2] * dy
        x_diff = xx - x0
        y_diff = yy - y0
        x_diff2 = x_diff^2
        y_diff2 = y_diff^2
        ϕ = atan(z[I[3]]/zR)

        lock(lk) do
            E0[I] = Gaussian(fn_params, w(z[I[3]]), w(z[I[3]]))[I[1], I[2]] .* exp.(1im .* ((x_diff2 ./ (2 .* R(z[I[3]])))) .+ (y_diff2 ./ (2 .* R(z[I[3]])))) .* exp.(1im .* (k .* z[I[3]] - ϕ))
        end
    end

    return E0
end

function LG3D(fn_params::FN_Params, w0::Real, p::Int, l::Int)

    @unpack x, y, z, x0, y0, dx, dy, N, λs = fn_params

    k = 2π/λs
    zR = π * w0^2 / λs
    w(zz) = w0 .* sqrt.(1 .+ (zz ./ zR).^2)
    R(zz) = zz .* (1 .+ (zR ./ zz)^2)
    E0 = zeros(ComplexF64, N, N, N) 

    lk = Threads.ReentrantLock()

    ThreadsX.foreach(CartesianIndices((N, N, N)); simd=true) do I
        xx = I[1] * dx
        yy = I[2] * dy
        x_diff = xx - x0
        y_diff = yy - y0
        x_diff2 = x_diff^2
        y_diff2 = y_diff^2
        ϕ = (2*p + abs(l) + 1)*atan(z[I[3]]/zR)

        lock(lk) do
            E0[I] = w0 .* (sqrt(x_diff2 + y_diff2)*sqrt(2)/w(z[I[3]]))^(abs(l)) .* LaguerreGauss(fn_params, 0, 1, 1, w(z[I[3]]))[I[1], I[2]] .* exp.(1im .* ((x_diff2 ./ (2 .* R(z[I[3]])))) .+ (y_diff2 ./ (2 .* R(z[I[3]])))) .* exp.(1im .* (k .* z[I[3]] - ϕ))
        end
    end

    return E0

end


println("EM-Field.jl compiled")
