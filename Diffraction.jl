struct Gauss end
struct LG end


@with_kw mutable struct Diffract

    f::Real
    fnum::Real
    w::Real
    nt::Real = 1.0

end

struct P end
struct S end
struct LHC end
struct RHC end


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

function TransmissionFunction(fn_params::FN_Params, diff_params::Diffract, Ex::Matrix, Ey::Matrix)
    @unpack f, fnum, nt = diff_params
    @unpack N, x, y, λs = fn_params

    # Aperture
    R_aperture = f / (fnum * 2)
    #println(2*R_aperture)
    aperture = smooth_circular_aperture(fn_params, R_aperture)
    M = Int((count(!iszero, aperture[Int((N-1)/2 + 1), :]) - 1) / 2)
    #println(M)

    # Parameters
    NA = 1 / (2 * fnum)
    k = 2π/λs
    kt = k * nt 
    ΔK = k * NA / M
    n = Int(2*M + 1)

    # Fields
    Exa = zeros(eltype(Ex), n, n)
    Eya = zeros(eltype(Ey), n, n)
    Etx = zeros(ComplexF64, n, n)
    Ety = zeros(ComplexF64, n, n)
    Etz = zeros(ComplexF64, n, n)

    # Transmission coefficients (generally functions of θ and ϕ)
    tp = 1
    ts = 1

    # Grid
    kx = ΔK * collect(range(-M, M, n))
    ky = ΔK * collect(range(-M, M, n))
    kx_mesh, ky_mesh = meshgrid(kx, ky)

    θmn = zeros(Float64, n, n)
    ϕmn = zeros(Float64, n, n)
    #θmn .= asin.(sqrt.(kx_mesh.^2 + ky_mesh.^2) / kt)
    θmn .= asin.(clamp.(sqrt.(kx_mesh.^2 + ky_mesh.^2) / kt, -1, 1))
    ϕmn .= atan.(ky_mesh, kx_mesh)

    # Transmitted Field
    Exa .= resize_symmetric(aperture .* (Ex), n)
    Eya .= resize_symmetric(aperture .* (Ey), n)

    #display(Plots.heatmap(Exa))
    #println(e2(y, abs2.(Exa[M+1, :])))

    cosθ = cos.(θmn)
    sinθ = sin.(θmn)
    cosϕ = cos.(ϕmn)
    sinϕ = sin.(ϕmn)
    A = sqrt.(cosθ) # Apodization function (conservation of energy)

    Etx = Exa .* tp .* (cosϕ.^2 .* cosθ .+ sinϕ.^2) .+ Eya .* ts .* sinϕ .* cosϕ .* (cosθ .- 1)
    Ety = Exa .* tp .* sinϕ .* cosϕ .* (cosθ .- 1) .+ Eya .* ts .* (cosϕ.^2 .+ sinϕ.^2 .* cosθ)
    Etz = Exa .* tp .* cosϕ .* sinθ .+ Eya .* ts .* sinϕ .* sinθ

    Et = [Etx .* A, Ety .* A, Etz .* A]

    return Et, θmn, kx_mesh, ky_mesh

end

function RichardsWolf(fn_params::FN_Params, diff_params::Diffract, Ex::Matrix, Ey::Matrix, z::Real)
    @unpack f, fnum, nt, w = diff_params
    @unpack N, λs, y = fn_params

    R_aperture = f / (fnum * 2)
    aperture = smooth_circular_aperture(fn_params, R_aperture)
    M = Int((count(!iszero, aperture[Int((N-1)/2 + 1), :]) - 1) / 2)
    n = Int(2*M + 1)
    #println(M)

    NA = 1 / (2 * fnum)
    kt = (2π/λs) * nt 
    factor = -1im * R_aperture^2 / (λs * f * M^2)

    Exf = zeros(ComplexF64, n, n)
    Eyf = similar(Exf)
    Ezf = similar(Exf)

    Et, θmn, kx_mesh, ky_mesh = TransmissionFunction(fn_params, diff_params, Ex, Ey)

    cosθ = cos.(θmn)
    kz = sqrt.(max.(kt^2 .- (kx_mesh.^2 .+ ky_mesh.^2), 0.0))
    expz = exp.(1im * z .* kz)
    #expz = exp.(1im * z .* sqrt.(kt^2 .- (kx_mesh.^2 .+ ky_mesh.^2)))

    Etx = Et[1]
    Ety = Et[2]
    Etz = Et[3]

    Ŋ = 2^10
    dx = M * λs / (Ŋ * NA)
    #println("dx = ", dx)
    #println("NA = ", NA)
    xmax = dx * M * 1.4
    x = collect(range(-xmax/2, xmax/2, n))
    y = copy(x)
    pad_size = Int(Ŋ / 2)
 
    Etx_z = Etx .* expz ./ cosθ
    Ety_z = Ety .* expz ./ cosθ
    Etz_z = Etz .* expz ./ cosθ

    Etx_vertpad = zeropad_vertical(Etx_z, pad_size)
    Ety_vertpad = zeropad_vertical(Ety_z, pad_size)
    Etz_vertpad = zeropad_vertical(Etz_z, pad_size)

    tempx = fftshift(fft(Etx_vertpad, 1), 1)
    tempy = fftshift(fft(Ety_vertpad, 1), 1)
    tempz = fftshift(fft(Etz_vertpad, 1), 1)

    tempx_horpad = zeropad_horizontal(tempx[pad_size+1:pad_size+n, :], pad_size)
    tempy_horpad = zeropad_horizontal(tempy[pad_size+1:pad_size+n, :], pad_size)
    tempz_horpad = zeropad_horizontal(tempz[pad_size+1:pad_size+n, :], pad_size)

    Exf .= factor .* ((fftshift(fft(tempx_horpad, 2), 2)[:, pad_size+1:pad_size+n]))
    Eyf .= factor .* ((fftshift(fft(tempy_horpad, 2), 2)[:, pad_size+1:pad_size+n]))
    Ezf .= factor .* ((fftshift(fft(tempz_horpad, 2), 2)[:, pad_size+1:pad_size+n]))

    Ef = [Exf, Eyf, Ezf]

    return Ef, x, y, n

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

function Visualize3D(Pol::String, Comp::String, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real)
    # Provide Pol as "P", "S", "LHC", or "RHC"
    # Provide Comp as "t", "x", "y", or "z" for total intensity or x-, y-, or z-component of the intensity respectively

    if Pol == "P"
        I, x, y, z, n = FullSpatialProfile(P(), fn_params, zmin, zmax, zsteps, f, fnum, nt, w)
    elseif Pol == "S"
        I, x, y, z, n = FullSpatialProfile(S(), fn_params, zmin, zmax, zsteps, f, fnum, nt, w)
    elseif Pol == "LHC"
        I, x, y, z, n = FullSpatialProfile(LHC(), fn_params, zmin, zmax, zsteps, f, fnum, nt, w)
    elseif Pol == "RHC"
        I, x, y, z, n = FullSpatialProfile(RHC(), fn_params, zmin, zmax, zsteps, f, fnum, nt, w)
    end

    z0 = Int(zsteps/2)
    M = Int((n-1) / 2)

    scale = zmax / maximum(x) / 2

    if Comp == "t"

        It = I[1]
        It ./= maximum(It[:, :, z0])
        fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(It[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
        ax1.xlabel = "z (μm)"
        ax1.ylabel = "x (μm)"
        Colorbar(fig[1, 2], hm1)
        Makie.save("FocusField-It-xz.png", fig)
        display(fig)
        #= scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.8,
                                isorange = 0.05, isovalue = maximum(It[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.5,
                                isorange = 0.01, isovalue = maximum(It[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(It[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(It[:, :, z0])*exp(-4))
        display(scene)
 =#
    elseif Comp == "x"

        Ix = I[2]
        Ix ./= maximum(Ix[:, :, z0])
        fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Ix[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
        ax1.xlabel = "z (μm)"
        ax1.ylabel = "x (μm)"
        Colorbar(fig[1, 2], hm1)
        Makie.save("FocusField-Ix-xz.png", fig)
        display(fig)
#=         scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.8,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.3,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.1,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-4))
        display(scene) =#

    elseif Comp == "y"

        Iy = I[3]
        Iy ./= maximum(Iy[:, :, z0])
        fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iy[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
        ax1.xlabel = "z (μm)"
        ax1.ylabel = "x (μm)"
        Colorbar(fig[1, 2], hm1)
        Makie.save("FocusField-Iy-xz.png", fig)
        display(fig)
#=         scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.8,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.3,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.1,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-4))
        display(scene) =#

    elseif Comp == "z"

        Iz = I[4]
        Iz ./= maximum(Iz[:, :, z0])
        fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iz[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
        ax1.xlabel = "z (μm)"
        ax1.ylabel = "x (μm)"
        Colorbar(fig[1, 2], hm1)
        Makie.save("FocusField-Iz-xz.png", fig)
        display(fig)
#=         scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.8,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.3,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.1,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-4))
        display(scene) =#

    end
end

function Visualize3D(Pol::String, Comp::String, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, w::Real, l::Int)
    # Provide Pol as "P", "S", "LHC", or "RHC"
    # Provide Comp as "t", "x", "y", or "z" for total intensity or x-, y-, or z-component of the intensity respectively
    # Provide l as integer for LG mode

    if Pol == "P"
        I, x, y, z, n = FullSpatialProfile(P(), fn_params, zmin, zmax, zsteps, f, fnum, nt, w, l)
    elseif Pol == "S"
        I, x, y, z, n = FullSpatialProfile(S(), fn_params, zmin, zmax, zsteps, f, fnum, nt, w, l)
    elseif Pol == "LHC"
        I, x, y, z, n = FullSpatialProfile(LHC(), fn_params, zmin, zmax, zsteps, f, fnum, nt, w, l)
    elseif Pol == "RHC"
        I, x, y, z, n = FullSpatialProfile(RHC(), fn_params, zmin, zmax, zsteps, f, fnum, nt, w, l)
    end

    z0 = Int(zsteps/2)
    M = Int((n-1) / 2)
    scale = zmax / maximum(x) / 2

    if Comp == "t"

        It = I[1]
        It ./= maximum(It[:, :, z0])
        #display(CairoMakie.heatmap(z.*1e6, scale.*x.*1e6, transpose(It[:, 65, :]), colorscale=log10))
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.8,
                                isorange = 0.05, isovalue = maximum(It[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(It[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.3,
                                isorange = 0.05, isovalue = maximum(It[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.1,
                                isorange = 0.05, isovalue = maximum(It[:, :, z0])*exp(-4))
        display(scene)

    elseif Comp == "x"

        Ix = I[2]
        Ix ./= maximum(Ix[:, :, z0])
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.8,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.3,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.1,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-4))
        display(scene)

    elseif Comp == "y"

        Iy = I[3]
        Iy ./= maximum(Iy[:, :, z0])
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.8,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.3,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.1,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-4))
        display(scene)

    elseif Comp == "z"

        Iz = I[4]
        Iz ./= maximum(Iz[:, :, z0])
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.8,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.3,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.1,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-4))
        display(scene)

    end
end

function DiffractionMovie(Comp::String, fn_params::FN_Params, diff_params::Diffract, Ex::Matrix, Ey::Matrix, zmin::Real, zmax::Real, zsteps::Int)

    z = collect(range(zmin, zmax, zsteps))

    Ef, x, y, = RichardsWolf(fn_params, diff_params, Ex, Ey, 0e-6);
    Itmax = maximum(abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
    Ixmax = maximum(abs2.(Ef[1]))
    Iymax = maximum(abs2.(Ef[2]))
    Izmax = maximum(abs2.(Ef[3]))

    for i in eachindex(z)
        Ef, x, y = RichardsWolf(fn_params, diff_params, Ex, Ey, z[i]);
        if Comp == "t"
            display(Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))))
        elseif Comp == "x"
            display(Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[1]))))
        elseif Comp == "y"
            display(Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[2]))))
        elseif Comp == "z"
            display(Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[3]))))
        end
        
    end
end


println("Diffraction.jl compiled")


# To visualize LG-modes (spatio-temporally if wanted)
#=
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
t = collect(range(-10e-15,10e-15,N))
At = ComplexEnvelope(1, 0, 0, 3.5e-15*scale, 0)
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

scene = GLMakie.volume(t.*1e15, x.*10, y.*10, Exyt, algorithm = :iso, isorange = 0.1, isovalue = maximum(Exyt[50,:,:])*0.7)
GLMakie.volume!(t.*1e15, x.*10, y.*10, Exyt, algorithm = :iso, isorange = 0.1, isovalue = minimum(Exyt[50,:,:])*0.7)
#scene = GLMakie.volume(t.*1e15, x.*20, y.*20, Ixyt, colormap = :inferno, algorithm = :iso, isorange = 0.1, isovalue = maximum(Ixyt[50,:,:])*0.7)
display(scene)

=#

#=

# To visualize beam and pump spatial profiles

@unpack E0_p, A, ηc, ηq, w_xp, w_yp, w_xs, w_ys, N, x, y = fn_params
w_xp = 800e-6 # at tube window
supergauss = SuperGaussian(fn_params, w_xp, 6)
scale = integrate((x,y), Gaussian(fn_params, w_xp, w_yp)) / integrate((x,y), supergauss)
Jsto = scale .* supergauss
Eabs = ((E0_p) * (A * ηc * ηq) * (2 - (A * ηc * ηq))
Aeff_p = calcAeff(x, y, Jsto)
Jsto0 = Eabs / Aeff_p
Jsto = Jsto0 .* Jsto
scaleJsto = Eabs / integrate((x,y), Jsto)
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
