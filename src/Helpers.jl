
function zeropad(A, pad_size)
    M = size(A)[1]

    if M % 2 != 0
        m = Int((M-1)/2)
        pad_range = pad_size-m:pad_size+m
    else
        m = Int(M/2)
        pad_range = pad_size-m:pad_size+m-1
    end

    padded_vector = zeros(eltype(A), 2 * pad_size)
    padded_vector[pad_range] .= A
    
    return padded_vector
end

function zeropad_horizontal(A, pad_size)
    M, N = size(A)

    if M % 2 != 0
        m = Int((M-1)/2)
        pad_range = pad_size-m:pad_size+m
    else
        m = Int(M/2)
        pad_range = pad_size-m:pad_size+m-1
    end

    padded_matrix = zeros(eltype(A), M, 2 * pad_size)
    padded_matrix[:, pad_range] .= A
    
    return padded_matrix
end

function zeropad_vertical(A, pad_size)
    M, N = size(A)

    if N % 2 != 0
        n = Int((N-1)/2)
        pad_range = pad_size-n:pad_size+n
    else
        n = Int(N/2)
        pad_range = pad_size-n:pad_size+n-1
    end

    padded_matrix = zeros(eltype(A), 2 * pad_size, N)
    padded_matrix[pad_range, :] .= A
    
    return padded_matrix
end

function find_first(A::Vector, tolerance::Real, type::String)
    # type must be = "e2" or "fwhm"

    max = maximum(A)

    if type == "e2"
        threshold = max * exp(-2)
    elseif type == "fwhm"
        threshold = max * 0.5
    end

    for (index, value) in enumerate(A)
        if abs(value - threshold) <= tolerance
            return index
        end
    end

    # Return nothing if no index is found (should not occur if threshold is reasonable)
    return nothing
end

function find_last(A::Vector, tolerance::Real, type::String)
    # type must be = "e2" or "fwhm"

    max = maximum(A)

    if type == "e2"
        threshold = max * exp(-2)
    elseif type == "fwhm"
        threshold = max * 0.5
    end

    # Iterate in reverse to find the last index
    for index in length(A):-1:1
        if abs(A[index] - threshold) <= tolerance
            return index
        end
    end

    # Return nothing if no index is found (should not occur if threshold is reasonable)
    return nothing
end

function CartesiantoPolar(fn_params::FN_Params, A::Matrix{Float64})
    @unpack N, x, y, xmax, ymax = fn_params

    Apol = similar(A)

    spl2D = Spline2D(x, y, A)
    Rmax = sqrt(xmax^2 + ymax^2)
    R = collect(range(0, Rmax, N))
    Φ = collect(range(0, 2π, N))

    for i in 1:N, j in 1:N
        Apol[i, j] = spl2D(R[i]*cos(Φ[j]), R[i]*sin(Φ[j]))
    end

    return Apol

end

function CartesiantoPolar(fn_params::FN_Params, A::Matrix{ComplexF64})
    @unpack N, x, y, xmax, ymax = fn_params

    Apol = similar(A)

    Re_spl2D = Spline2D(x, y, real.(A))
    Im_spl2D = Spline2D(x, y, imag.(A))
    Rmax = sqrt(xmax^2 + ymax^2)
    R = collect(range(0, Rmax, N))
    Φ = collect(range(0, 2π, N))

    for i in 1:N, j in 1:N
        Apol[i, j] = Re_spl2D(R[i]*cos(Φ[j]), R[i]*sin(Φ[j])) .+ 1im .* Im_spl2D(R[i]*cos(Φ[j]), R[i]*sin(Φ[j]))
    end

    return Apol

end

function circular_aperture(fn_params::FN_Params, R::Real)
    @unpack N, x, y, x0, y0 = fn_params

    aperture = zeros(Float64, N, N)

    X, Y = meshgrid(x, y)
    x_diff = X .- x0
    y_diff = Y .- y0

    aperture[x_diff.^2 .+ y_diff.^2 .<= R^2] .= 1
    
    return aperture
end

function elliptical_aperture(fn_params::FN_Params, a::Real, b::Real)
    @unpack N, x, y, x0, y0 = fn_params

    aperture = zeros(Float64, N, N)

    X, Y = meshgrid(x, y)
    x_diff = X .- x0
    y_diff = Y .- y0

    aperture[(x_diff.^2 ./ a^2) .+ (y_diff.^2 ./ b^2) .<= 1] .= 1

    return aperture
end

function resize_symmetric(v::AbstractVector, new_size::Int)
    old_size = length(v)
    if old_size == new_size
        return v
    elseif old_size > new_size
        # Crop vector symmetrically
        start_idx = div(old_size - new_size, 2) + 1
        end_idx = start_idx + new_size - 1
        return v[start_idx:end_idx]
    else
        # Pad vector symmetrically with zeros
        pad_size = div(new_size - old_size, 2)
        padded = zeros(eltype(v), new_size)
        padded[pad_size+1:end-pad_size] .= v
        return padded
    end
end

function resize_symmetric(A::AbstractMatrix, new_size::Int)
    old_size = size(A, 1) # Assuming square matrix
    if old_size == new_size
        return A
    elseif old_size > new_size
        # Crop matrix symmetrically
        start_idx = div(old_size - new_size, 2) + 1
        end_idx = start_idx + new_size - 1
        return A[start_idx:end_idx, start_idx:end_idx]
    else
        # Pad matrix symmetrically with zeros
        pad_size = div(new_size - old_size, 2)
        padded = zeros(eltype(A), new_size, new_size)
        padded[pad_size+1:end-pad_size, pad_size+1:end-pad_size] .= A
        return padded
    end
end

function scaleField!(x::Vector, y::Vector, Ex::Matrix, Ey::Matrix, Ein::Real)

    current = calcEnergy(x, y, abs2.(Ex) .+ abs2.(Ey))
    
    scale = sqrt(Ein / current)
    
    Ex .*= scale
    Ey .*= scale
    
end

function calcEnergy(x::Vector, y::Vector, J::Matrix)

    En = NumericalIntegration.integrate((x, y), J)

    return En
end

function FWHM(X, Y)
    half_max = maximum(Y) / 2
    d = sign.(half_max .- Y[1:end-1]) .- sign.(half_max .- Y[2:end])
    left_idx = findfirst(d .> 0)
    right_idx = findlast(d .< 0)
    return (X[right_idx] - X[left_idx])
end

function e2(X, Y)
    threshold = maximum(Y) / exp(2)
    d = sign.(threshold .- Y[1:end-1]) .- sign.(threshold .- Y[2:end])
    left_idx = findfirst(d .> 0)
    right_idx = findlast(d .< 0)
    return (X[right_idx] - X[left_idx])
end

function FWHM2D(x, y, A)
    N = size(A)[1]

    if N % 2 != 0
        w0_x = FWHM(x, A[:, Int((N-1)/2+1)]) / 2
        w0_y = FWHM(y, A[Int((N-1)/2+1), :]) / 2
    else
        w0_x = FWHM(x, A[:, Int(N/2)]) / 2
        w0_y = FWHM(y, A[Int(N/2), :]) / 2
    end

    return w0_x, w0_y
end

function e22D(x, y, A)
    N = size(A)[1]

    if N % 2 != 0
        w0_x = e2(x, A[:, Int((N-1)/2+1)]) / 2
        w0_y = e2(y, A[Int((N-1)/2+1), :]) / 2
    else
        w0_x = e2(x, A[:, Int(N/2)]) / 2
        w0_y = e2(y, A[Int(N/2), :]) / 2 
    end

    return w0_x, w0_y
end

function FWHM2D(x, y, max_index1, max_index2, A)
    N = size(A)[1]

    if N % 2 != 0
        w0_x = FWHM(x, A[:, max_index2]) / 2
        w0_y = FWHM(y, A[max_index1, :]) / 2
    else
        w0_x = FWHM(x, A[:, max_index2]) / 2
        w0_y = FWHM(y, A[max_index1, :]) / 2
    end

    return w0_x, w0_y
end

function e22D(x, y, max_index1, max_index2, A)
    N = size(A)[1]

    if N % 2 != 0
        w0_x = e2(x, A[:, max_index2]) / 2
        w0_y = e2(y, A[max_index1, :]) / 2
    else
        w0_x = e2(x, A[:, max_index2]) / 2
        w0_y = e2(y, A[max_index1, :]) / 2
    end

    return w0_x, w0_y
end

function HoleyMirror!(fn_params::FN_Params, x0::Real, y0::Real, R::Real, E::Matrix)

    @unpack x, y = fn_params                        

    X, Y = meshgrid(x, y)
    X .-= x0
    Y .-= y0

    E[(X.^2 .+ Y.^2) .< R^2] .= 0
end

function readRefInd(RefInd::String)
    # Assuming data comes from refractiveindex.info, so always convert λ back to SI units

    data = CSV.read(RefInd, DataFrame)
    wl = data[!, 1] .* 1e-6
    n = data[!, 2]
    κ = data[!, 3]

    Respl = Spline1D(wl, n, k=3)
    Imspl = Spline1D(wl, κ, k=3) 

    𝑁(λ) = Respl(λ) + 1im * Imspl(λ)

    return 𝑁

end

# TODO include transmission coefficients
function FresnelCoefficients(θi::Real, λ0::Real)

    # Provide θi in degrees
    sinθi = sind(θi)
    cosθi = cosd(θi)

    # For silver mirror @ 785 nm
    𝑁 = readRefInd("input_data/Ag-RefInd.csv")(λ0)
    rp = (sqrt(𝑁^2 - sinθi^2) - 𝑁^2*cosθi) / (sqrt(𝑁^2 - sinθi^2) + 𝑁^2*cosθi)
    rs = (cosθi - sqrt(𝑁^2 - sinθi^2)) / (cosθi + sqrt(𝑁^2 - sinθi^2))

    return rp, rs

end

function SpectralPhase(ϕ0::Real, ϕ1::Real, ϕ2::Real, ϕ3::Real, ϕ4::Real, ν::Vector, ν0::Real)
    # Provide in units of fs, fs^2, fs^3, fs^4

    ϕ = (ϕ0 .+ (ϕ1 .* (ν .- ν0)) .+ (0.5 .* ϕ2 .* (ν .- ν0).^2) 
            .+ ((1/6) .* ϕ3 .* (ν .- ν0).^3) .+ ((1/24) .* ϕ4 .* (ν .- ν0).^4))

    return ϕ

end

function readSpect(fn_params::FN_Params, SpectData::String)
    @unpack c = fn_params

    data = CSV.read(SpectData, DataFrame)
    wl = data[!, 1] .* 1e-9
    I_ret = data[!, 3]
    phase = data[!, 4]

    spl1 = Spline1D(wl, I_ret, k=3)
    spl2 = Spline1D(wl, phase, k=3)

    I(λ) = spl1(λ)
    ϕ(λ) = spl2(λ)

    return wl, I, ϕ

end

function readHarmonicSpect(fn_params::FN_Params, SpectData::String)
    @unpack h, c = fn_params

    data = CSV.read(SpectData, DataFrame)
    conv = h * c / 1.602177e-19
    wl = (conv ./ data[!, 1])[:]
    I_ret = data[!, 2]

    spl1 = Spline1D(wl, I_ret, k=3)

    I(λ) = spl1(λ)
    ϕ(λ) = 0

    return wl, I, ϕ

end

function ZernikeCoefficients(Z1::Real, Z2::Real, Z3::Real, Z4::Real, 
                                Z5::Real, Z6::Real, Z7::Real, Z8::Real, 
                                    Z9::Real, Z10::Real, Z11::Real)

    return [Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11]

end

function Zernike(fn_params::FN_Params, E::Matrix, Z::Vector, l::Real)

    @unpack N, dx = fn_params

    if N % 2 != 0
        m = Int((N-1)/2)
    else
        m = Int(N/2)
    end

    E = abs.(E) ./ maximum(abs.(E))
    n = find_last(E[m, :], 1e-2, "e2") - find_first(E[m, :], 1e-2, "e2")
    n = n | 1

    if l != 0
        ϵ = find_first(E[m, m:end], 1e-2, "e2") * dx
    else
        ϵ = 0
    end

    x = collect(range(-1, 1, n))
    y = collect(range(-1, 1, n))
    X, Y = meshgrid(x, y)

    r = sqrt.(X.^2 .+ Y.^2)
    ϕ = atan.(Y, X)
    
    Zernike_terms = Dict(
        1 => () -> ones(size(r)),
        2 => () -> (2 .* r .* cos.(ϕ) ./ sqrt(ϵ^2 + 1)),
        3 => () -> (2 .* r .* sin.(ϕ) ./ sqrt(ϵ^2 + 1)),
        4 => () -> (sqrt(3) .* (-2 .* r.^2 .+ ϵ^2 .+ 1) ./ (ϵ^2 - 1)),
        5 => () -> (sqrt(6) .* r.^2 .* sin.(2 .* ϕ) ./ sqrt(ϵ^4 + ϵ^2 + 1)),
        6 => () -> (sqrt(6) .* r.^2 .* cos.(2 .* ϕ) ./ sqrt(ϵ^4 + ϵ^2 + 1)),
        7 => () -> (2 * sqrt(2) .* r .* (3 .* r.^2 .* (ϵ^2 + 1) .- 2 * (ϵ^4 + ϵ^2 + 1)) .* sin.(ϕ) ./ sqrt((ϵ^2 - 1)^2*(ϵ^2 + 1)*(ϵ^4 + 4*ϵ^2 + 1))),
        8 => () -> (2 * sqrt(2) .* r .* (3 .* r.^2 .* (ϵ^2 + 1) .- 2 * (ϵ^4 + ϵ^2 + 1)) .* cos.(ϕ) ./ sqrt((ϵ^2 - 1)^2*(ϵ^2 + 1)*(ϵ^4 + 4*ϵ^2 + 1))),
        9 => () -> (2 * sqrt(2) .* r.^3 .* sin.(3 .* ϕ) ./ sqrt(ϵ^6 + ϵ^4 + ϵ^2 + 1)),
        10 => () -> (2 * sqrt(2) .* r.^3 .* cos.(3 .* ϕ) ./ sqrt(ϵ^6 + ϵ^4 + ϵ^2 + 1)),
        11 => () -> (sqrt(5) * (6 .* r.^4 .- 6 * (ϵ^2 + 1) .* r.^2 .+ ϵ^4 .+ 4 * ϵ^2 .+ 1) ./ (ϵ^2 - 1)^2)
    )

    Z_tot = sum(Z[i] .* Zernike_terms[i]() for i in eachindex(Z) if Z[i] != 0)
    Z_tot[r .> 1] .= 0

    Ab = resize_symmetric(Z_tot, N)
    X, Y = meshgrid(fn_params.x, fn_params.y)
    Ab[(X.^2 .+ Y.^2) .< ϵ^2] .= 0

    return Ab

end

function getTotalSAMDensity(fn_params::FN_Params, diff_params::Diffract, Ef::Vector, Hf::Vector)

    @unpack ϵ0, λs, c, μ0 = fn_params
    @unpack ϕ = diff_params

    con = (λs) / (8π * c)

    Ex, Ey, Ez = Ef[1], Ef[2], Ef[3]
    Hx, Hy, Hz = Hf[1], Hf[2], Hf[3]

    sx = con .* imag((ϵ0.*conj(Ey) .* Ez - ϵ0.*conj(Ez) .* Ey) .+ (μ0.*conj(Hy) .* Hz - μ0.*conj(Hz) .* Hy))
    sy = con .* imag((ϵ0.*conj(Ez) .* Ex - ϵ0.*conj(Ex) .* Ez) .+ (μ0.*conj(Hz) .* Hx - μ0.*conj(Hx) .* Hz))

    st = sx .+ 1im .* sy
    sz = con .* imag((ϵ0.*conj(Ex) .* Ey - ϵ0.*conj(Ey) .* Ex) .+ (μ0.*conj(Hx) .* Hy - μ0.*conj(Hy) .* Hx))

    return abs.(st)./maximum(abs.(st)), real.(sz)./maximum(abs.(sz))
end

function partial_diff(E, dx, dy, dz)

    ∂E_∂x = diff(E, dims=1) / dx
    ∂E_∂y = diff(E, dims=2) / dy
    ∂E_∂z = diff(E, dims=3) / dz

    ∂E_∂x = vcat(∂E_∂x, reshape(∂E_∂x[end, :, :], 1, size(∂E_∂x, 2), size(∂E_∂x, 3)))
    ∂E_∂y = hcat(∂E_∂y, reshape(∂E_∂y[:, end, :], size(∂E_∂y, 1), 1, size(∂E_∂y, 3)))
    ∂E_∂z = cat(∂E_∂z, reshape(∂E_∂z[:, :, end], size(∂E_∂z, 1), size(∂E_∂z, 2), 1); dims=3)

    return ∂E_∂x, ∂E_∂y, ∂E_∂z
end


function getTotalOAMDensity(fn_params::FN_Params, diff_params::Diffract, Ef::Vector, Hf::Vector, x, y, z)

    @unpack ϵ0, λs, c, μ0 = fn_params
    @unpack ϕ = diff_params

    Ex, Ey, Ez = Ef[1], Ef[2], Ef[3]
    Hx, Hy, Hz = Hf[1], Hf[2], Hf[3]

    dx = x[2] - x[1]
    dy = y[2] - y[1]
    dz = z[2] - z[1]

    ∂Ex_∂x, ∂Ex_∂y, ∂Ex_∂z = partial_diff(Ex, dx, dy, dz)
    ∂Ey_∂x, ∂Ey_∂y, ∂Ey_∂z = partial_diff(Ey, dx, dy, dz)
    ∂Ez_∂x, ∂Ez_∂y, ∂Ez_∂z = partial_diff(Ez, dx, dy, dz)

    ∂Hx_∂x, ∂Hx_∂y, ∂Hx_∂z = partial_diff(Hx, dx, dy, dz)
    ∂Hy_∂x, ∂Hy_∂y, ∂Hy_∂z = partial_diff(Hy, dx, dy, dz)
    ∂Hz_∂x, ∂Hz_∂y, ∂Hz_∂z = partial_diff(Hz, dx, dy, dz)

    con = (λs) / (8π * c)
    px = con * imag((ϵ0 .* conj(Ex) .* ∂Ex_∂x .+ ϵ0 .* conj(Ey) .* ∂Ey_∂x .+ ϵ0 .* conj(Ez) .* ∂Ez_∂x) .+ 
                    (μ0 .* conj(Hx) .* ∂Hx_∂x .+ μ0 .* conj(Hy) .* ∂Hy_∂x .+ μ0 .* conj(Hz) .* ∂Hz_∂x))
    py = con * imag((ϵ0 .* conj(Ex) .* ∂Ex_∂y .+ ϵ0 .* conj(Ey) .* ∂Ey_∂y .+ ϵ0 .* conj(Ez) .* ∂Ez_∂y) .+
                    (μ0 .* conj(Hx) .* ∂Hx_∂y .+ μ0 .* conj(Hy) .* ∂Hy_∂y .+ μ0 .* conj(Hz) .* ∂Hz_∂y))
    pz = con * imag((ϵ0 .* conj(Ex) .* ∂Ex_∂z .+ ϵ0 .* conj(Ey) .* ∂Ey_∂z .+ ϵ0 .* conj(Ez) .* ∂Ez_∂z) .+
                    (μ0 .* conj(Hx) .* ∂Hx_∂z .+ μ0 .* conj(Hy) .* ∂Hy_∂z .+ μ0 .* conj(Hz) .* ∂Hz_∂z))                    

    x_grid = reshape(x, size(x, 1), 1, 1)
    y_grid = reshape(y, 1, size(y, 1), 1)
    z_grid = reshape(z, 1, 1, size(z, 1))

    Lx = y_grid .* pz .- z_grid .* py
    Ly = z_grid .* px .- x_grid .* pz
    Lz = x_grid .* py .- y_grid .* px

    Lt = Lx .+ 1im .* Ly

    return abs.(Lt)./maximum(abs.(Lt)), real.(Lz)./maximum(abs.(Lz))
end

function VisSpinDensity(s::Matrix, x::Vector, y::Vector; save=false)

    fig, ax, hm = GLMakie.heatmap(x, y, s, colormap=:thermometer)
    cb = CairoMakie.Colorbar(fig[1, 2], hm, size=30; label = L"\textbf{SAM Density (arb. u.)}")
    ax.xlabel = "x (μm)"
    ax.ylabel = "y (μm)"

    if save
        Makie.save("SpinAM-Density.png", fig)
    end

    fig

end

function VisOAMDensity(p::Matrix, x::Vector, y::Vector; save=false)

    fig, ax, hm = GLMakie.heatmap(x, y, p, colormap=:inferno)
    cb = CairoMakie.Colorbar(fig[1, 2], hm, size=30; label = L"\textbf{OAM Density (arb. u.)}")
    ax.xlabel = "x (μm)"
    ax.ylabel = "y (μm)"

    if save
        Makie.save("SpinAM-Density.png", fig)
    end
    
    fig

end

function getPolarizationEllipse2D(x, y, Ex, Ey; num_ellipses=(21, 21), draw_arrow=true,
                       amplification=0.75, color_line="white", line_width=0.5, save=true)
    # Approach taken from https://github.com/aocg-ucm/diffractio/blob/main/diffractio/vector_fields_XY.py

    intensity_max = maximum(abs.(Ex).^2 .+ abs.(Ey).^2)
    Dx = x[end] - x[1]
    Dy = y[end] - y[1]

    size_x = Dx / num_ellipses[1]
    size_y = Dy / num_ellipses[2]

    x_centers = size_x/2 .+ size_x .* collect(0:num_ellipses[1]-1)
    y_centers = size_y/2 .+ size_y .* collect(0:num_ellipses[2]-1)

    num_x, num_y = length(x), length(y)
    ix_centers = round.(Int, num_x / num_ellipses[1] / 2 .+ num_x / num_ellipses[1] .* collect(0:num_ellipses[1]-1))
    iy_centers = round.(Int, num_y / num_ellipses[2] / 2 .+ num_y / num_ellipses[2] .* collect(0:num_ellipses[2]-1))

    fig, ax1, hm1 = CairoMakie.heatmap(x, y, abs.(Ex).^2 .+ abs.(Ey).^2)
    ax1.xlabel = L"\textbf{x (μm)}"
    ax1.ylabel = L"\textbf{y (μm)}"
    cb = CairoMakie.Colorbar(fig[1, 2], hm1, size=30;
        label = L"\textbf{Peak Intensity (arb. u.)}")

    for (i, xi) in enumerate(ix_centers)
        for (j, yj) in enumerate(iy_centers)
            E0x = Ex[yj, xi]
            E0y = Ey[yj, xi]

            angles = LinRange(0, 2π, 64)
            Ex_real = real(E0x .* exp.(1im .* angles))
            Ey_real = real(E0y .* exp.(1im .* angles))

            max_r = maximum(sqrt.(Ex_real.^2 .+ Ey_real.^2))
            size_dim = min(size_x, size_y)

            if max_r > 0 && max_r^2 > exp(-2) * intensity_max
                Ex_real = Ex_real / max_r * size_dim * amplification / 2 .+ x[Int(xi)]
                Ey_real = Ey_real / max_r * size_dim * amplification / 2 .+ y[Int(yj)]

                lines!(ax1, Ex_real, Ey_real, color=color_line, linewidth=line_width, label="")

                if draw_arrow
                    arrows!(ax1, [Ex_real[1]], [Ey_real[1]], [Ex_real[1] - Ex_real[2]], 
                    [Ey_real[1] - Ey_real[2]], arrowsize = 11, color=color_line, linewidth=0)
                end
            end
        end
    end

    if save
        Makie.save("PolState_XY.png", fig)
    end

    return fig
end

function getPolarizationEllipse3D(x, y, z, Ex, Ey, Ez;
                                        num_ellipses=(10, 10, 10), draw_arrow=true, 
                                        amplification=0.75, color_line=:black, line_width=0.5, save=true)
    intensity_max = maximum(abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez))
    Dx, Dy, Dz = x[end] - x[1], y[end] - y[1], z[end] - z[1]

    size_x = Dx / num_ellipses[1]
    size_y = Dy / num_ellipses[2]
    size_z = Dz / num_ellipses[3]

    num_x, num_y, num_z = length(x), length(y), length(z)
    ix_centers = round.(Int, num_x / num_ellipses[1] / 2 .+ num_x / num_ellipses[1] .* collect(0:num_ellipses[1]-1))
    iy_centers = round.(Int, num_y / num_ellipses[2] / 2 .+ num_y / num_ellipses[2] .* collect(0:num_ellipses[2]-1))
    iz_centers = round.(Int, num_z / num_ellipses[3] / 2 .+ num_z / num_ellipses[3] .* collect(0:num_ellipses[3]-1))

    fig = Figure(resolution=(1280, 960))
    ax = Axis3(fig[1, 1], aspect = :equal, xlabel = "x (μm)", ylabel = "y (μm)", zlabel = "z (μm)")


    volume!(ax, x, y, z, abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez), algorithm = :iso, alpha = 0.5, transparency=true,
                                           isorange = maximum(It)*0.05, isovalue = maximum(It)*exp(-2))

    arrow_positions = Point3f[]  # Arrow start points
    arrow_directions = Vec3f[]  # Arrow direction vectors

    for (i, xi) in enumerate(ix_centers)
        for (j, yj) in enumerate(iy_centers)
            for (k, zk) in enumerate(iz_centers)
                E0x, E0y, E0z = Ex[xi, yj, zk], Ey[xi, yj, zk], Ez[xi, yj, zk]

                angles = LinRange(0, 2π, 64)
                Ex_real = real(E0x .* exp.(1im .* angles))
                Ey_real = real(E0y .* exp.(1im .* angles))
                Ez_real = real(E0z .* exp.(1im .* angles))

                max_r = maximum(sqrt.(Ex_real.^2 .+ Ey_real.^2 .+ Ez_real.^2))
                size_dim = min(size_x, size_y, size_z)

                if max_r > 0 && max_r^2 > exp(-2) * intensity_max
                    Ex_real = Ex_real / max_r * size_dim * amplification / 2 .+ x[Int(xi)]
                    Ey_real = Ey_real / max_r * size_dim * amplification / 2 .+ y[Int(yj)]
                    Ez_real = Ez_real / max_r * size_dim * amplification / 2 .+ z[Int(zk)]

                    GLMakie.lines!(ax, Ex_real, Ey_real, Ez_real, color=color_line, linewidth=line_width)

                    if draw_arrow
                        push!(arrow_positions, Point3f(Ex_real[1], Ey_real[1], Ez_real[1]))
                        push!(arrow_directions, Vec3f(Ex_real[1] - Ex_real[2], 
                                                      Ey_real[1] - Ey_real[2], 
                                                      Ez_real[1] - Ez_real[2]))
                    end
                end
            end
        end
    end

    if draw_arrow
        GLMakie.arrows!(ax, arrow_positions, arrow_directions, 
                linecolor=:gray, arrowcolor=color_line, 
                linewidth=0.0, arrowsize=Vec3f(0.1, 0.1, 0.1), align=:center)
    end

    if save
        Makie.save("PolState_XYZ.png", fig)
    end

    return fig
end

function Visualize3D(Pol::String, Comp::String, fn_params::FN_Params, diff_params::Diffract, 
                        zmin::Real, zmax::Real, zsteps::Int; slicex=false, slicey=false, l = 0, 
                            coeffs = 0, save=false, intensity=true, phase=false)
    # Provide Pol as "P", "S", "LHC", or "RHC"
    # Provide Comp as "t", "x", "y", or "z" for total intensity or x-, y-, or z-component of the intensity respectively
    @unpack N = fn_params
    
    if typeof(l) == Vector{Float64}
        Ef, x, y, z = FullSpatialProfile(fn_params, diff_params, Pol, zmin, zmax, zsteps, l, Z, coeffs=coeffs)
    else
        Ef, x, y, z = FullSpatialProfile(fn_params, diff_params, Pol, zmin, zmax, zsteps, l, Z)
    end

    if intensity

        It = zeros(Float64, N, N, zsteps)
        Ix = zeros(Float64, N, N, zsteps)
        Iy = zeros(Float64, N, N, zsteps)
        Iz = zeros(Float64, N, N, zsteps)

        It .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix .= abs2.(Ef[1])
        Iy .= abs2.(Ef[2])
        Iz .= abs2.(Ef[3])

    elseif phase

        ϕx = zeros(Float64, N, N, zsteps)
        ϕy = zeros(Float64, N, N, zsteps)
        ϕz = zeros(Float64, N, N, zsteps)

        magx = abs.(Ef[1])
        magy = abs.(Ef[2])
        magz = abs.(Ef[2])

        ϕx .= angle.(Ef[1])
        ϕx[magx .< exp(-2).*maximum(magx)] .= 0

        ϕy .= angle.(Ef[2])
        ϕy[magy .< exp(-2).*maximum(magy)] .= 0

        ϕz .= angle.(Ef[3])
        ϕz[magz .< exp(-2).*maximum(magz)] .= 0

    else

        Ex = zeros(Float64, N, N, zsteps)
        Ey = zeros(Float64, N, N, zsteps)
        Ez = zeros(Float64, N, N, zsteps)

        Ex .= real.(Ef[1])
        Ey .= real.(Ef[2])
        Ez .= real.(Ef[3])

    end

    if zsteps % 2 != 0
        z0 = Int((zsteps-1)/2)
    else
        z0 = Int(zsteps/2)
    end
    M = Int((fn_params.N-1) / 2)

    #Initialize fig
    if slicex == false && slicey == false
        fig = Figure(size=(940,940))
        ax1 = Axis3(fig[1, 1], aspect = :equal, xlabel = "x (μm)", ylabel = "y (μm)", zlabel = "z (μm)")
    end

    if Comp == "t"

        It ./= maximum(It[:, :, z0])

        if slicex
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(It[M, :, :]), colorscale=log10, colormap=:inferno, colorrange=(1e-5, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "x (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-It-xz.png", fig)
            end
            display(fig)
        elseif slicey
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, y.*1e6, transpose(It[:, M, :]), colorscale=log10, colormap=:inferno, colorrange=(1e-5, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "y (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-It-yz.png", fig)
            end
            display(fig)            
        else
            volume!(ax1, x.*1e6, y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.8,
                                    isorange = 0.1, isovalue = maximum(It[:, :, z0])*exp(-1))
            volume!(ax1, x.*1e6, y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.5,
                                    isorange = 0.05, isovalue = maximum(It[:, :, z0])*exp(-2))
            volume!(ax1, x.*1e6, y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.3,
                                    isorange = 0.01, isovalue = maximum(It[:, :, z0])*exp(-3))
            volume!(ax1, x.*1e6, y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.1,
                                    isorange = 0.001, isovalue = maximum(It[:, :, z0])*exp(-4))
            if save 
                Makie.save("FocusField-It-3D.png", fig)
            end
            
            display(fig)
        end
    
    elseif Comp == "x"

        if slicex
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Ix[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "x (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Ix-xz.png", fig)
            end
            display(fig)
        elseif slicex
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, y.*1e6, transpose(Ix[:, M, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "y (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Ix-yz.png", fig)
            end
            display(fig)
        else
            if intensity

                Ix ./= maximum(Ix[:, :, z0])

                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.8,
                                        isorange = 0.1, isovalue = maximum(Ix[:, :, z0])*exp(-1))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.5,
                                        isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-2))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.3,
                                        isorange = 0.01, isovalue = maximum(Ix[:, :, z0])*exp(-3))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.1,
                                        isorange = 0.001, isovalue = maximum(Ix[:, :, z0])*exp(-4))
                if save
                    Makie.save("FocusField-Ix-3D.png", fig)
                end

            elseif phase
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, ϕx, algorithm = :iso, colormap=:jet,
                                        isorange = 0.1, isovalue = 3.1)
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, ϕx, algorithm = :iso, colormap=:jet,
                                        isorange = 0.1, isovalue = -3.1)
                if save
                    Makie.save("FocusField-phix-3D.png", fig)
                end

            else

                Ex ./= maximum(abs.(Ex))

                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.8, colormap=:berlin,
                                        isorange = 0.1, isovalue = maximum(Ex[:, :, z0])*exp(-0.25))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.8, colormap=:berlin,
                                        isorange = 0.1, isovalue = minimum(Ex[:, :, z0])*exp(-0.25))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.5, colormap=:berlin,
                                        isorange = 0.1, isovalue = maximum(Ex[:, :, z0])*exp(-1))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.5, colormap=:berlin,
                                        isorange = 0.1, isovalue = minimum(Ex[:, :, z0])*exp(-1))
                if save
                    Makie.save("FocusField-Ex-3D.png", fig)
                end    

            end
            
            display(fig)
        end

    elseif Comp == "y"
        
        if slicex
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iy[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "x (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Iy-xz.png", fig)
            end
            display(fig)
        elseif slicex
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iy[:, M, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "y (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Iy-yz.png", fig)
            end
            display(fig)
        else
            if intensity

                Iy ./= maximum(Iy[:, :, z0])

                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.8,
                                        isorange = 0.1, isovalue = maximum(Iy[:, :, z0])*exp(-1))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.5,
                                        isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-2))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.3,
                                        isorange = 0.01, isovalue = maximum(Iy[:, :, z0])*exp(-3))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.1,
                                        isorange = 0.001, isovalue = maximum(Iy[:, :, z0])*exp(-4))
                if save
                    Makie.save("FocusField-Iy-3D.png", fig)
                end 

            else

                Ey ./= maximum(abs.(Ey))

                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ey, algorithm = :iso, alpha = 0.8, colormap=:berlin,
                                        isorange = 0.1, isovalue = maximum(Ey[:, :, z0])*exp(-0.25))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ey, algorithm = :iso, alpha = 0.8, colormap=:berlin,
                                        isorange = 0.1, isovalue = minimum(Ey[:, :, z0])*exp(-0.25))
                if save
                    Makie.save("FocusField-Ey-3D.png", fig)
                end    

            end
            display(fig)
        end

    elseif Comp == "z"
        
        if slicex
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iz[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "x (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Iz-xz.png", fig)
            end
            display(fig)
        elseif slicex
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iz[:, M, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "y (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Iz-yz.png", fig)
            end
            display(fig)
        else
            if intensity

                Iz ./= maximum(Iz[:, :, z0])

                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.8,
                                        isorange = 0.1, isovalue = maximum(Iz[:, :, z0])*exp(-1))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.5,
                                        isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-2))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.3,
                                        isorange = 0.01, isovalue = maximum(Iz[:, :, z0])*exp(-3))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.1,
                                        isorange = 0.001, isovalue = maximum(Iz[:, :, z0])*exp(-4))
                if save
                    Makie.save("FocusField-Iz-3D.png", fig)
                end      

            else

                Ez ./= maximum(abs.(Ez))

                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ez, algorithm = :iso, alpha = 0.8, colormap=:berlin,
                                        isorange = 0.1, isovalue = maximum(Ez[:, :, z0])*exp(-0.25))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ez, algorithm = :iso, alpha = 0.8, colormap=:berlin,
                                        isorange = 0.1, isovalue = minimum(Ez[:, :, z0])*exp(-0.25))
                if save
                    Makie.save("FocusField-Ez-3D.png", fig)
                end    

            end
            display(fig)
        end

    end
end

function DiffractionMovie(Pol, Comp::String, fn_params::FN_Params, diff_params::Diffract, 
                            zmin::Real, zmax::Real, zsteps::Int, l::Real, Z::Vector; save=false, 
                                intensity=true, phase=false, aberration=false, hole=false)
    @unpack kt, m, w = diff_params

    z = collect(range(zmin, zmax, zsteps))
    zR = π*w^2 / fn_params.λs

    Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, 0, l, Z, aberration=aberration, hole=hole)

    It_max = maximum(abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
    Ix_max = maximum(abs2.(Ef[1]))
    Iy_max = maximum(abs2.(Ef[2]))
    Iz_max = maximum(abs2.(Ef[3]))

    Ex_max = maximum(abs.(Ef[1]))
    Ey_max = maximum(abs.(Ef[2]))
    Ez_max = maximum(abs.(Ef[3]))

    if save
        if isdir("DiffractionMovie-Images")
            nothing
        else
            mkdir("DiffractionMovie-Images")
        end
    end

    for i in eachindex(z)
        Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, z[i], l, Z, aberration=aberration, hole=hole);

        ψg = (abs(l) + 1)*atan(z[i] / zR)
        Ef[1] .*= exp(1im * ψg)
        Ef[2] .*= exp(1im * ψg)
        Ef[3] .*= exp(1im * ψg)

        if Comp == "t"
            if (intensity == false && phase == false) || (intensity == false && phase == true)
                break
            end
            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
            ax.xlabel = L"\textbf{x (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Total intensity @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
            if save
                cd("DiffractionMovie-Images")
                Makie.save("Itot_z"*string(i)*".png", fig)
                cd("..")
            else
                display(fig)
                sleep(0.01)
            end
        elseif Comp == "x"
            if intensity
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, abs2.(Ef[1]))
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Intensity (x-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Itot_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            elseif phase
                ϕ = angle.(Ef[1])
                mag = abs.(Ef[1])
                ϕ[mag .< exp(-1).*maximum(mag)] .= NaN
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, ϕ, colorrange=(-π, π), colormap=:jet)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Phase (x-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Phase (rad)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("phix_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            else
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, real.(Ef[1]), colormap=:thermometer)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Electric field (x-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field Strength (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Ex_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            end

        elseif Comp == "y"
            if intensity
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, abs2.(Ef[2]))
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Intensity (y-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Itot_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            elseif phase
                ϕ = angle.(Ef[2])
                mag = abs.(Ef[2])
                ϕ[mag .< exp(-1).*maximum(mag)] .= NaN
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, ϕ, colorrange=(-π, π), colormap=:jet)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Phase (y-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Phase (rad)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("phiy_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            else
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, real.(Ef[2]), colormap=:thermometer)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Electric field (y-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field Strength (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Ey_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            end

        elseif Comp == "z"
            if intensity
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, abs2.(Ef[3]))
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Intensity (z-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Itot_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            elseif phase
                ϕ = angle.(Ef[3])
                mag = abs.(Ef[3])
                ϕ[mag .< exp(-1).*maximum(mag)] .= NaN
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, ϕ, colorrange=(-π, π), colormap=:jet)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Phase (z-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Phase (rad)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("phiz_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            else
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, real.(Ef[3]), colormap=:thermometer)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Electric field (z-comp) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field Strength (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Ez_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            end
            
        end
        
    end
end

function DiffractionMovie(Comp::String, Ex::Array, Ey::Array, Ez::Array, 
                            x::Vector, y::Vector, z::Vector; intensity=true, save=false)

    max_tot = maximum(abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez))
    max_Ex = maximum(abs.(Ex))
    max_Ey = maximum(abs.(Ey))
    max_Ez = maximum(abs.(Ez))
    max_Ix = maximum(abs2.(Ex))
    max_Iy = maximum(abs2.(Ey))
    max_Iz = maximum(abs2.(Ez))

    if save
        if isdir("DiffractionMovie-Images")
            nothing
        else
            mkdir("DiffractionMovie-Images")
        end
    end

    if Comp == "t"
        I_tot = abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez)
        I_tot ./= max_tot

        foreach(eachindex(z)) do i
            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, I_tot[:, :, i], colorrange=(0, 1))
            ax.xlabel = L"\textbf{x (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Total intensity @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
            if save
                cd("DiffractionMovie-Images")
                Makie.save("Itot_z"*string(i)*".png", fig)
                cd("..")
            else
                display(fig)
                sleep(0.01)
            end
        end

    elseif Comp == "x"
        if intensity
            Ix = abs2.(Ex)
            Ix ./= max_Ix

            foreach(eachindex(z)) do i
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, Ix[:, :, i], colorrange=(0, 1))
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Intensity (x-comp.) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Ix_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            end
        else

            foreach(eachindex(z)) do i
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, real.(Ex)[:, :, i], colormap=:thermometer)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Electric field (x-comp.) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Ex_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            end
        end

    elseif Comp == "y"
        if intensity
            Iy = abs2.(Ey)
            Iy ./= max_Iy

            foreach(eachindex(z)) do i
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, Iy[:, :, i], colorrange=(0,1))
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Intensity (y-comp.) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Iy_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            end
        else

            foreach(eachindex(z)) do i
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, real.(Ey)[:, :, i], colormap=:thermometer)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Electric field (y-comp.) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Ey_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            end
        end

    elseif Comp == "z"
        if intensity
            Iz = abs2.(Ez)
            Iz ./= max_Iz

            foreach(eachindex(z)) do i
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, Iz[:, :, i], colorrange=(0, 1))
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Intensity (z-comp.) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Iz_z"*string(i)*".png", fig)
                    cd("..")
                else
                    display(fig)
                    sleep(0.01)
                end
            end
        else

            foreach(eachindex(z)) do i
                fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, real.(Ez)[:, :, i], colormap=:thermometer)
                ax.xlabel = L"\textbf{x (μm)}"
                ax.ylabel = L"\textbf{y (μm)}"
                ax.title = "Electric field (z-comp.) @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
                cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
                if save
                    cd("DiffractionMovie-Images")
                    Makie.save("Ez_z"*string(i)*".png", fig)
                    cd("..0")
                else
                    display(fig)
                    sleep(0.01)
                end
            end
        end

    end

end

function DiffractionMovie(Comp::String, st::Array, sz::Array, Lt::Array, Lz::Array, 
                            x::Vector, y::Vector, z::Vector; save=false)


    if save
        if isdir("DiffractionMovie-Images")
            nothing
        else
            mkdir("DiffractionMovie-Images")
        end
    end

    if Comp == "St"

        foreach(eachindex(z)) do i
            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, st[:, :, i], colorrange=(0, 1), colormap=:thermometer)
            ax.xlabel = L"\textbf{x (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Transverse SAM Density @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{SAM Density (arb. u)}")
            if save
                cd("DiffractionMovie-Images")
                Makie.save("St_z"*string(i)*".png", fig)
                cd("..")
            else
                display(fig)
                sleep(0.01)
            end
        end

    elseif Comp == "Sz"

        foreach(eachindex(z)) do i
            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, sz[:, :, i], colormap=:thermometer)
            ax.xlabel = L"\textbf{x (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Longitudinal SAM Density @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{SAM Density (arb. u)}")
            if save
                cd("DiffractionMovie-Images")
                Makie.save("Sz_z"*string(i)*".png", fig)
                cd("..")
            else
                display(fig)
                sleep(0.01)
            end
        end

    elseif Comp == "Lt"

        foreach(eachindex(z)) do i
            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, Lt[:, :, i], colormap=:inferno)
            ax.xlabel = L"\textbf{x (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Transverse OAM Density @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
            if save
                cd("DiffractionMovie-Images")
                Makie.save("Lt_z"*string(i)*".png", fig)
                cd("..")
            else
                display(fig)
                sleep(0.01)
            end
        end

    elseif Comp == "Lz"

        foreach(eachindex(z)) do i
            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, Lz[:, :, i], colorrange=(0, 1), colormap=:inferno)
            ax.xlabel = L"\textbf{x (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Longitudinal OAM Density @ z = "*string(round(z[i]*1e6, digits=2))*" μm"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{OAM Density (arb. u)}")
            if save
                cd("DiffractionMovie-Images")
                Makie.save("Lz_z"*string(i)*".png", fig)
                cd("..")
            else
                display(fig)
                sleep(0.01)
            end
        end
    end

end

function XZSlices(Comp::String, Ex::Array, Ey::Array, Ez::Array, 
                            x::Vector, y::Vector, z::Vector; l=0, intensity=true, 
                            save=false, verbose=true)


    I = findmax(abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez))[2]
    maxx_index, maxy_index, maxz_index = I[1], I[2], I[3]

    if verbose
        open("output.txt", "w") do io
            NA = 1 / (2 * diff_params.fnum)
            println(io, "Numerical aperture of system = ", round(NA, digits=2))
    
            I_focus = abs2.(Ex[:, :, maxz_index]) .+ abs2.(Ey[:, :, maxz_index]) .+ abs2.(Ez[:, :, maxz_index])
            E_focus = calcEnergy(x, y, I_focus)
            println(io, "Energy @ focus = ", round(E_focus * 1e3, digits=3), " mJ")
    
            w0_x, w0_y = e22D(x, y, I_focus)
            println(io, "Beam waist (e2) @ focus =", round(w0_x*1e6, digits=2), " μm x ", round(w0_y*1e6, digits=2), " μm")
            w0_x, w0_y = FWHM2D(x, y, I_focus)
            println(io, "Beam spot size (FWHM) @ focus = ", round(w0_x*1e6, digits=2), " µm x ", round(w0_y*1e6, digits=2), " µm")


            w0_z = e2(z, abs2.(Ex[maxx_index, maxy_index, :]) .+ abs2.(Ey[maxx_index, maxy_index, :]) .+ abs2.(Ez[maxx_index, maxy_index, :]))
            τ = (w0_z / c) * 1e15
            println(io, "Pulse extent (1/e2) in z-direction = ", round(w0_z*1e6, digits=2), " μm")
            println(io, "Pulse duration (1/e2) = ", round(τ, digits=2), " fs")
            w0_z = FWHM(z, abs2.(Ex[maxx_index, maxy_index, :]) .+ abs2.(Ey[maxx_index, maxy_index, :]) .+ abs2.(Ez[maxx_index, maxy_index, :]))
            τ = (w0_z / c) * 1e15
            println(io, "Pulse extent (FWHM) in z-direction = ", round(w0_z*1e6, digits=2), " μm")
            println(io, "Pulse duration (FWHM) = ", round(τ, digits=2), " fs")
    
            Aeff = calcAeff(x, y, I_focus)
            println(io, "Effective area = ", round(Aeff*1e12, digits=2), " µm^2")
    
            Ppeak = 0.94 * E_focus / τ
            if l == 0
                I_target = 2 * Ppeak / Aeff
            else
                I_target = Ppeak / Aeff
            end
            println(io, "Peak intensity @ focus = ", round(I_target * 1e-4, digits=3), " W/cm^2")
    
            E_focus = calcEnergy(x, y, abs2.(Ex[:, :, maxz_index]))
            println(io, "Energy  in x-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
            E_focus = calcEnergy(x, y, abs2.(Ey[:, :, maxz_index]))
            println(io, "Energy  in y-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
            E_focus = calcEnergy(x, y, abs2.(Ez[:, :, maxz_index]))
            println(io, "Energy  in z-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
        end
    end

    max_tot = maximum(abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez))
    max_Ex = maximum(abs.(Ex))
    max_Ey = maximum(abs.(Ey))
    max_Ez = maximum(abs.(Ez))

    if Comp == "t"
        I_tot = abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez)
        I_tot ./= max_tot
        
        fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(I_tot[:, maxy_index, :]))
        ax.xlabel = L"\textbf{z (μm)}"
        ax.ylabel = L"\textbf{x (μm)}"
        ax.title = "Total intensity"
        cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
        display(fig)
        if save
            Makie.save("Itot_xz.png", fig)
        end

    elseif Comp == "x"
        if intensity
            Ix = abs2.(Ex)
            Ix ./= max_tot

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(Ix[:, maxy_index, :]))
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{x (μm)}"
            ax.title = "Intensity (x-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
            display(fig)
            if save
                Makie.save("Ix_xz.png", fig)
            end
        else

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(real.(Ex)[:, maxy_index, :]) ./ max_Ex, colormap=:RdBu)
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{x (μm)}"
            ax.title = "Electric field (x-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
            display(fig)
            if save
                Makie.save("Ex_xz.png", fig)
            end
        end

    elseif Comp == "y"
        if intensity
            Iy = abs2.(Ey)
            Iy ./= max_tot

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(Iy[:, maxy_index, :]))
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{x (μm)}"
            ax.title = "Intensity (y-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
            display(fig)
            if save
                Makie.save("Iy_xz.png", fig)
            end
        else

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(real.(Ey)[:, maxy_index, :]) ./ max_Ey, colormap=:RdBu)
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{x (μm)}"
            ax.title = "Electric field (y-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
            display(fig)
            if save
                Makie.save("Ey_xz.png", fig)
            end
        end

    elseif Comp == "z"
        if intensity
            Iz = abs2.(Ez)
            Iz ./= max_tot

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(Iz[:, maxy_index, :]))
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{x (μm)}"
            ax.title = "Intensity (z-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
            display(fig)
            if save
                Makie.save("Iz_xz.png", fig)
            end
        else

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(real.(Ez)[:, maxy_index, :]) ./ max_Ez, colormap=:RdBu)
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{x (μm)}"
            ax.title = "Electric field (z-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
            display(fig)
            if save
                Makie.save("Ez_xz.png", fig)
            end
        end

    end

end

function YZSlices(Comp::String, Ex::Array, Ey::Array, Ez::Array, 
                            x::Vector, y::Vector, z::Vector; l=0, intensity=true, 
                            save=false, verbose=true)


    I = findmax(abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez))[2]
    maxx_index, maxy_index, maxz_index = I[1], I[2], I[3]

    if verbose
        open("output.txt", "w") do io
            NA = 1 / (2 * diff_params.fnum)
            println(io, "Numerical aperture of system = ", round(NA, digits=2))
    
            maxz_index = findmax(abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez))[2][3]
            maxx_index = findmax(abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez))[2][1]
            I_focus = abs2.(Ex[:, :, maxz_index]) .+ abs2.(Ey[:, :, maxz_index]) .+ abs2.(Ez[:, :, maxz_index])
            E_focus = calcEnergy(x, y, I_focus)
            println(io, "Energy @ focus = ", round(E_focus * 1e3, digits=3), " mJ")
    
            w0_x, w0_y = FWHM2D(x, y, I_focus)
            println(io, "Beam spot size (FWHM) @ focus = ", round(w0_x*1e6, digits=2), " µm x ", round(w0_y*1e6, digits=2), " µm")


            w0_x, w0_y = e22D(x, y, I_focus)
            println(io, "Beam waist (e2) @ focus =", round(w0_x*1e6, digits=2), " μm x ", round(w0_y*1e6, digits=2), " μm")
            w0_x, w0_y = FWHM2D(x, y, I_focus)
            println(io, "Beam spot size (FWHM) @ focus = ", round(w0_x*1e6, digits=2), " µm x ", round(w0_y*1e6, digits=2), " µm")


            w0_z = e2(z, abs2.(Ex[maxx_index, maxy_index, :]) .+ abs2.(Ey[maxx_index, maxy_index, :]) .+ abs2.(Ez[maxx_index, maxy_index, :]))
            τ = (w0_z / c) * 1e15
            println(io, "Pulse extent (1/e2) in z-direction = ", round(w0_z*1e6, digits=2), " μm")
            println(io, "Pulse duration (1/e2) = ", round(τ, digits=2), " fs")
            w0_z = FWHM(z, abs2.(Ex[maxx_index, maxy_index, :]) .+ abs2.(Ey[maxx_index, maxy_index, :]) .+ abs2.(Ez[maxx_index, maxy_index, :]))
            τ = (w0_z / c) * 1e15
            println(io, "Pulse extent (FWHM) in z-direction = ", round(w0_z*1e6, digits=2), " μm")
            println(io, "Pulse duration (FWHM) = ", round(τ, digits=2), " fs")
    
            Aeff = calcAeff(x, y, I_focus)
            println(io, "Effective area = ", round(Aeff*1e12, digits=2), " µm^2")
    
            Ppeak = 0.94 * E_focus / τ
            if l == 0
                I_target = 2 * Ppeak / Aeff
            else
                I_target = Ppeak / Aeff
            end
            println(io, "Peak intensity @ focus = ", round(I_target * 1e-4, digits=3), " W/cm^2")
    
            E_focus = calcEnergy(x, y, abs2.(Ex[:, :, maxz_index]))
            println(io, "Energy  in x-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
            E_focus = calcEnergy(x, y, abs2.(Ey[:, :, maxz_index]))
            println(io, "Energy  in y-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
            E_focus = calcEnergy(x, y, abs2.(Ez[:, :, maxz_index]))
            println(io, "Energy  in z-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
        end
    end

    max_tot = maximum(abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez))
    max_Ex = maximum(abs.(Ex))
    max_Ey = maximum(abs.(Ey))
    max_Ez = maximum(abs.(Ez))

    if Comp == "t"
        I_tot = abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez)
        I_tot ./= max_tot
        
        fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(I_tot[maxx_index, :, :]))
        ax.xlabel = L"\textbf{z (μm)}"
        ax.ylabel = L"\textbf{y (μm)}"
        ax.title = "Total intensity"
        cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
        display(fig)
        if save
            Makie.save("Itot_yz.png", fig)
        end

    elseif Comp == "x"
        if intensity
            Ix = abs2.(Ex)
            Ix ./= max_tot

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(Ix[maxx_index, :, :]))
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Intensity (x-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
            display(fig)
            if save
                Makie.save("Ix_yz.png", fig)
            end
        else

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(real.(Ex)[maxx_index, :, :]), colormap=:RdBu)
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Electric field (x-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
            display(fig)
            if save
                Makie.save("Ex_yz.png", fig)
            end
        end

    elseif Comp == "y"
        if intensity
            Iy = abs2.(Ey)
            Iy ./= max_tot

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(Iy[maxx_index, :, :]))
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Intensity (y-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
            display(fig)
            if save
                Makie.save("Iy_yz.png", fig)
            end
        else

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(real.(Ey)[maxx_index, :, :]), colormap=:RdBu)
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Electric field (y-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
            display(fig)
            if save
                Makie.save("Ey_yz.png", fig)
            end
        end

    elseif Comp == "z"
        if intensity
            Iz = abs2.(Ez)
            Iz ./= max_tot

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(Iz[maxx_index, :, :]))
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Intensity (z-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Peak intensity (arb. u)}")
            display(fig)
            if save
                Makie.save("Iz_yz.png", fig)
            end
        else

            fig, ax, hm = CairoMakie.heatmap(x.*1e6, y.*1e6, transpose(real.(Ez)[maxx_index, :, :]), colormap=:RdBu)
            ax.xlabel = L"\textbf{z (μm)}"
            ax.ylabel = L"\textbf{y (μm)}"
            ax.title = "Electric field (z-comp.)"
            cbar = Colorbar(fig[1, 2], hm, label=L"\textbf{Field strength (arb. u)}")
            display(fig)
            if save
                Makie.save("Ez_yz.png", fig)
            end
        end

    end

end


function loadFieldData()

    Ex_real = load("Ex_real.jld", "Ex_real");

    Ey_real = load("Ey_real.jld", "Ey_real");

    Ez_real = load("Ez_real.jld", "Ez_real");

    Ex_imag = load("Ex_imag.jld", "Ex_imag");

    Ey_imag = load("Ey_imag.jld", "Ey_imag");

    Ez_imag = load("Ez_imag.jld", "Ez_imag");

    Ex = Ex_real .+ 1im .* Ex_imag;

    Ex = permutedims(Ex, (2, 1, 3));

    Ey = Ey_real .+ 1im .* Ey_imag;

    Ey = permutedims(Ey, (2, 1, 3));

    Ez = Ez_real .+ 1im .* Ez_imag;

    Ez = permutedims(Ez, (2, 1, 3));

    It = abs2.(Ex) .+ abs2.(Ey) .+ abs2.(Ez);

    x = CSV.read("x.csv", DataFrame)[!,1];

    y = CSV.read("y.csv", DataFrame)[!,1];

    z = CSV.read("z.csv", DataFrame)[!,1];

    return Ex, Ey, Ez, It, x, y, z

end

function normalize!(I::Array, νsteps::Int)

    for i in 1:νsteps

        I[i,:,:,:] ./= maximum(abs.(I[i,:,:,:]))

    end

end

println("Helpers.jl compiled")
