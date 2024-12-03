
function zeropad_horizontal(A, pad_size)
    M, N = size(A)
    
    padded_matrix = zeros(eltype(A), M, N + 2 * pad_size)  # 2 * pad_size for left and right padding
    
    padded_matrix[:, pad_size+1:pad_size+N] .= A
    
    return padded_matrix
end

function zeropad_vertical(A, pad_size)
    M, N = size(A)
    
    padded_matrix = zeros(eltype(A), M + 2 * pad_size, N)
    
    padded_matrix[pad_size+1:pad_size+M, :] .= A
    
    return padded_matrix
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

    for i in 1:N, j in 1:N
        x_diff = x[i] - x0
        y_diff = y[j] - y0
        if x_diff^2 + y_diff^2 <= R^2
            aperture[i, j] = 1
        end
    end
    
    return aperture
end

function smooth_circular_aperture(fn_params::FN_Params, R::Real)
    @unpack N, x, y, dx, dy, x0, y0 = fn_params

    aperture = zeros(Float64, N, N)
    ΔR = sqrt(dx^2 + dy^2)

    for i in 1:N, j in 1:N
        x_diff = x[i] - x0
        y_diff = y[j] - y0
        r = sqrt(x_diff^2 + y_diff^2)
        if x_diff^2 + y_diff^2 <= R^2
            aperture[i, j] = 0.5 * (1 + tanh((1.5/ΔR) * (R - r)))
        end
    end

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
        padded[pad_size+1:end-pad_size] = v
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
        padded[pad_size+1:end-pad_size, pad_size+1:end-pad_size] = A
        return padded
    end
end

function scaleField!(x::Vector, y::Vector, Ex::Matrix, Ey::Matrix, Ein::Real)

    Aeff = calcAeff(x, y, (abs2.(Ex) .+ abs2.(Ey))) # Initial effective area
    Jin = 2* Ein / Aeff # Initial pulse fluence
    Ex .*= sqrt(Jin)
    Ey .*= sqrt(Jin)
    
end

function calcEnergy(x::Vector, y::Vector, I_tot::Matrix)

    En = round(NumericalIntegration.integrate((x, y), I_tot), digits=3)

    return En
end

function FWHM(X, Y)
    half_max = maximum(Y) / 2
    d = sign.(half_max .- Y[1:end-1]) .- sign.(half_max .- Y[2:end])
    left_idx = findfirst(d .> 0)
    right_idx = findlast(d .< 0)
    return (X[right_idx] - X[left_idx]) / 2
end

function e2(X, Y)
    threshold = maximum(Y) / exp(2)
    d = sign.(threshold .- Y[1:end-1]) .- sign.(threshold .- Y[2:end])
    left_idx = findfirst(d .> 0)
    right_idx = findlast(d .< 0)
    return (X[right_idx] - X[left_idx]) / 2
end

function FWHM2D(x, y, A)
    N = size(A)[1]

    if N % 2 != 0
        w0_x = FWHM(x, A[Int((N-1)/2+1), :])
        w0_y = FWHM(y, A[:, Int((N-1)/2+1)])
    else
        w0_x = FWHM(x, A[Int(N/2), :])
        w0_y = FWHM(y, A[:, Int(N/2)])
    end

    return w0_x, w0_y
end

function e22D(x, y, A)
    N = size(A)[1]

    if N % 2 != 0
        w0_x = e2(x, A[Int((N-1)/2+1), :]) / 2
        w0_y = e2(y, A[:, Int((N-1)/2+1)]) / 2
    else
        w0_x = e2(x, A[Int(N/2), :]) / 2
        w0_y = e2(y, A[:, Int(N/2)]) / 2 
    end

    return w0_x, w0_y
end

function FresnelCoefficients(θi::Real)

    # Provide θi in degrees

    sinθi = sind(θi)
    cosθi = cosd(θi)

    # For silver mirror @ 785 nm
    𝑁 = 0.034455 + (1im * 5.4581)
    rp = (sqrt(𝑁^2 - sinθi^2) - 𝑁^2*cosθi) / (sqrt(𝑁^2 - sinθi^2) + 𝑁^2*cosθi)
    rs = (cosθi - sqrt(𝑁^2 - sinθi^2)) / (cosθi + sqrt(𝑁^2 - sinθi^2))

    return rp, rs

end

function FullSpatialProfile(::P, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real)
    @unpack N = fn_params

    diff_params = Diffract(f = f, fnum = fnum, w = w, nt = nt)

    E = Gaussian(fn_params, w, w)

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y, n = RichardsWolf(fn_params, diff_params, E, zeros(N, N), 0)


    It = zeros(Float64, n, n, zsteps)
    Ix = zeros(Float64, n, n, zsteps)
    Iy = zeros(Float64, n, n, zsteps)
    Iz = zeros(Float64, n, n, zsteps)

    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, E, zeros(N, N), z[I])

        It[:, :, I] .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix[:, :, I] .= abs2.(Ef[1])
        Iy[:, :, I] .= abs2.(Ef[2])
        Iz[:, :, I] .= abs2.(Ef[3])
    end

    I = [It, Ix, Iy, Iz]

    return I, x, y, z, n
end

function FullSpatialProfile(::S, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real)
    @unpack N = fn_params

    diff_params = Diffract(f = f, fnum = fnum, w = w, nt = nt)

    E = Gaussian(fn_params, w, w)

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y, n = RichardsWolf(fn_params, diff_params, zeros(N, N), E, 0)

    It = zeros(Float64, n, n, zsteps)
    Ix = zeros(Float64, n, n, zsteps)
    Iy = zeros(Float64, n, n, zsteps)
    Iz = zeros(Float64, n, n, zsteps)

    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, zeros(N, N), E, z[I])

        It[:, :, I] .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix[:, :, I] .= abs2.(Ef[1])
        Iy[:, :, I] .= abs2.(Ef[2])
        Iz[:, :, I] .= abs2.(Ef[3])
    end

    I = [It, Ix, Iy, Iz]

    return I, x, y, z, n
end

function FullSpatialProfile(::LHC, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real)
    @unpack N = fn_params

    diff_params = Diffract(f = f, fnum = fnum, w = w, nt = nt)

    E = (1/sqrt(2)) .* Gaussian(fn_params, w, w)

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y, n = RichardsWolf(fn_params, diff_params, E, 1im.*E, 0)

    It = zeros(Float64, n, n, zsteps)
    Ix = zeros(Float64, n, n, zsteps)
    Iy = zeros(Float64, n, n, zsteps)
    Iz = zeros(Float64, n, n, zsteps)

    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, E, 1im .* E, z[I])

        It[:, :, I] .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix[:, :, I] .= abs2.(Ef[1])
        Iy[:, :, I] .= abs2.(Ef[2])
        Iz[:, :, I] .= abs2.(Ef[3])
    end

    I = [It, Ix, Iy, Iz]

    return I, x, y, z, n
end

function FullSpatialProfile(::RHC, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real)
    @unpack N = fn_params

    diff_params = Diffract(f = f, fnum = fnum, w = w, nt = nt)

    E = (1/sqrt(2)) .* Gaussian(fn_params, w, w)

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y, n = RichardsWolf(fn_params, diff_params, E, -1im.*E, 0)

    It = zeros(Float64, n, n, zsteps)
    Ix = zeros(Float64, n, n, zsteps)
    Iy = zeros(Float64, n, n, zsteps)
    Iz = zeros(Float64, n, n, zsteps)


    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, E, -1im .* E, z[I])

        It[:, :, I] .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix[:, :, I] .= abs2.(Ef[1])
        Iy[:, :, I] .= abs2.(Ef[2])
        Iz[:, :, I] .= abs2.(Ef[3])
    end

    I = [It, Ix, Iy, Iz]

    return I, x, y, z, n
end

function FullSpatialProfile(::P, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real, l::Int)
    @unpack N = fn_params

    diff_params = Diffract(f = f, fnum = fnum, w = w, nt = nt)

    E = LaguerreGauss(fn_params, 0, l, 1, w)

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y, n = RichardsWolf(fn_params, diff_params, E, zeros(N, N), 0)

    It = zeros(Float64, n, n, zsteps)
    Ix = zeros(Float64, n, n, zsteps)
    Iy = zeros(Float64, n, n, zsteps)
    Iz = zeros(Float64, n, n, zsteps)

    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, E, zeros(N, N), z[I])

        It[:, :, I] .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix[:, :, I] .= abs2.(Ef[1])
        Iy[:, :, I] .= abs2.(Ef[2])
        Iz[:, :, I] .= abs2.(Ef[3])
    end

    I = [It, Ix, Iy, Iz]

    return I, x, y, z, n
end

function FullSpatialProfile(::S, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real, l::Int)
    @unpack N = fn_params

    diff_params = Diffract(f = f, fnum = fnum, w = w, nt = nt)

    E = LaguerreGauss(fn_params, 0, l, 1, w)

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y, n = RichardsWolf(fn_params, diff_params, zeros(N, N), E, 0)

    It = zeros(Float64, n, n, zsteps)
    Ix = zeros(Float64, n, n, zsteps)
    Iy = zeros(Float64, n, n, zsteps)
    Iz = zeros(Float64, n, n, zsteps)

    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, zeros(N, N), E, z[I])

        It[:, :, I] .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix[:, :, I] .= abs2.(Ef[1])
        Iy[:, :, I] .= abs2.(Ef[2])
        Iz[:, :, I] .= abs2.(Ef[3])
    end

    I = [It, Ix, Iy, Iz]

    return I, x, y, z, n
end

function FullSpatialProfile(::LHC, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real, l::Int)
    @unpack N = fn_params

    diff_params = Diffract(f = f, fnum = fnum, w = w, nt = nt)

    E = (1/sqrt(2)) .* LaguerreGauss(fn_params, 0, l, 1, w)

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y, n = RichardsWolf(fn_params, diff_params, E, 1im.*E, 0)

    It = zeros(Float64, n, n, zsteps)
    Ix = zeros(Float64, n, n, zsteps)
    Iy = zeros(Float64, n, n, zsteps)
    Iz = zeros(Float64, n, n, zsteps)

    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, E, 1im .* E, z[I])

        It[:, :, I] .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix[:, :, I] .= abs2.(Ef[1])
        Iy[:, :, I] .= abs2.(Ef[2])
        Iz[:, :, I] .= abs2.(Ef[3])
    end

    I = [It, Ix, Iy, Iz]

    return I, x, y, z, n
end

function FullSpatialProfile(::RHC, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real, l::Int)
    @unpack N = fn_params

    diff_params = Diffract(f = f, fnum = fnum, w = w, nt = nt)

    E = (1/sqrt(2)) .* LaguerreGauss(fn_params, 0, l, 1, w)

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y, n = RichardsWolf(fn_params, diff_params, E, -1im.*E, 0)

    It = zeros(Float64, n, n, zsteps)
    Ix = zeros(Float64, n, n, zsteps)
    Iy = zeros(Float64, n, n, zsteps)
    Iz = zeros(Float64, n, n, zsteps)

    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, E, -1im .* E, z[I])

        It[:, :, I] .= (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
        Ix[:, :, I] .= abs2.(Ef[1])
        Iy[:, :, I] .= abs2.(Ef[2])
        Iz[:, :, I] .= abs2.(Ef[3])
    end

    I = [It, Ix, Iy, Iz]

    return I, x, y, z, n
end

println("Helpers.jl compiled")