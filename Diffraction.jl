struct Gauss end
struct LG end


mutable struct Diffract{T}

    f::T
    fnum::T
    w::T
    R::T
    sinθmax::T
    nt::T
    kt::T
    Ein::T
    aperture::Matrix{T}
    θ::Matrix{T}
    ϕ::Matrix{T}



    function Diffract{T}(fn_params::FN_Params, f::Real, fnum::Real, w::Real, nt::Real) where T
        @unpack λs, x, y = fn_params

        # Parameters
        R = f / (fnum * 2)
        sinθmax = R ./ sqrt(R^2 + f^2)

        k = 2π/λs
        kt = k * nt 
        
        Ein = 2.5e-3

        aperture = smooth_circular_aperture(fn_params, R)

        # Grid
        X, Y = meshgrid(x, y)
        r = sqrt.(X.^2 .+ Y.^2)
        θ = r ./ f
        ϕ = atan.(Y, X)

        new{T}(f, fnum, w, R, sinθmax, nt, kt, Ein, aperture, θ, ϕ)
    end

end

struct P end
struct S end
struct LHC end
struct RHC end

function Polarization(fn_params::FN_Params, diff_params::Diffract, Pol::String)
    @unpack N, x, y = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(Float64, N, N)
    Ey = similar(Ex)

    @show Pol

    if Pol == "P"
        Ex .= Gaussian(fn_params, w, w)
        Ey .= zeros(N, N)
    elseif Pol == "S"
        Ex .= zeros(N, N)
        Ey .= Gaussian(fn_params, w, w)
    elseif Pol == "D"
        Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
        Ey .= Gaussian(fn_params, w, w) ./ sqrt(2)
    elseif Pol == "LHC"
        Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
        Ey .= -1im .* Gaussian(fn_params, w, w) ./ sqrt(2)
    elseif Pol == "RHC"
        Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
        Ey .= 1im .* Gaussian(fn_params, w, w) ./ sqrt(2)
    elseif Pol == "Radial"
        Ex .= cosϕ.* Gaussian(fn_params, w, w) ./ sqrt(2)
        Ey .= sinϕ .* Gaussian(fn_params, w, w) ./ sqrt(2)
    elseif Pol == "Azimuthal"
        Ex .= -sinϕ .* Gaussian(fn_params, w, w) ./ sqrt(2)
        Ey .= cosϕ .* Gaussian(fn_params, w, w) ./ sqrt(2)
    end

    scaleField!(x, y, Ex, Ey, Ein)

    rp, rs = FresnelCoefficients(28)

    M00 = rp .* cosϕ.^2 .* cosθ .+ rs .* sinϕ.^2
    M01 = sinϕ .* cosϕ .* (rp .* cosθ .- rs)

    M10 = sinϕ .* cosϕ .* (rp .* cosθ .- rs)
    M11 = rs .* cosϕ.^2 .+ rp .* sinϕ.^2 .* cosθ

    M20 = rp .* sinθ .* cosϕ
    M21 = rp .* sinθ .* sinϕ

    Epx = M00 .* Ex .+ M01 .* Ey
    Epy = M10 .* Ex .+ M11 .* Ey
    Epz = M20 .* Ex .+ M21 .* Ey

    return Epx, Epy, Epz
end



function TransmissionFunction(fn_params::FN_Params, diff_params::Diffract, Pol::String)
    @unpack sinθmax, R, aperture, θ, ϕ = diff_params
    @unpack N, x, y = fn_params
    
    # Fields
    Etx = zeros(ComplexF64, N, N)
    Ety = similar(Etx)
    Etz = similar(Etx)

    println("Effective aperture radius = ", round(R.*1e3, digits=2), " mm")

    # Apply polarization
    Epx, Epy, Epz = Polarization(fn_params, diff_params, Pol)

    # Apodization
    cosθ = cos.(θ)
    Apod = 1 ./ cosθ

    # Transmitted fields
    Etx .= Apod .* aperture .* Epx
    Ety .= Apod .* aperture .* Epy
    Etz .= Apod .* aperture .* Epz

    # Transmitted Energy
    It_trans = abs2.(Etx) + abs2.(Ety) + abs2.(Etz)
    E_trans = calcEnergy(x, y, It_trans)
    println("Energy after parabola = ", round(E_trans * 1e3, digits=3), " mJ")
    w0_x, w0_y = e22D(x, y, It_trans)
    println("Beam spot size (1/e^2) on parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
    #display(heatmap(x, y, It_trans))

    return Etx, Ety, Etz

end

function RichardsWolf(fn_params::FN_Params, diff_params::Diffract, Pol::String, z::Real)
    @unpack sinθmax, f, kt, θ = diff_params
    @unpack N, λs, dx, dy = fn_params

    factor = -(1im * sinθmax.^2 / (f * λs))

    Etx, Ety, Etz = TransmissionFunction(fn_params, diff_params, Pol)

    # Fields
    Efx = zeros(ComplexF64, N, N)
    Efy = similar(Efx)
    Efz = similar(Efx)

    # Propagation factor
    expz = exp.(1im .* kt .* z .* cos.(θ))

    M = 2^13
    pad_size = Int(M / 2)

    Etx_vertpad = zeropad_vertical(Etx .* expz, pad_size)
    Ety_vertpad = zeropad_vertical(Ety .* expz, pad_size)
    Etz_vertpad = zeropad_vertical(Etz .* expz, pad_size)

    tempx = fftshift(fft(Etx_vertpad, 1), 1)
    tempy = fftshift(fft(Ety_vertpad, 1), 1)
    tempz = fftshift(fft(Etz_vertpad, 1), 1)

    tempx_horpad = zeropad_horizontal(tempx[pad_size+1:pad_size+N, :], pad_size)
    tempy_horpad = zeropad_horizontal(tempy[pad_size+1:pad_size+N, :], pad_size)
    tempz_horpad = zeropad_horizontal(tempz[pad_size+1:pad_size+N, :], pad_size)

    Efx .= factor .* ((fftshift(fft(tempx_horpad, 2), 2)[:, pad_size+1:pad_size+N]))
    Efy .= factor .* ((fftshift(fft(tempy_horpad, 2), 2)[:, pad_size+1:pad_size+N]))
    Efz .= factor .* ((fftshift(fft(tempz_horpad, 2), 2)[:, pad_size+1:pad_size+N]))

    Efx .= Matrix(transpose(Efx))
    Efy .= Matrix(transpose(Efy))
    Efz .= Matrix(transpose(Efz))

    Ef = [Efx, Efy, Efz]

    freq_nyquist_x = 1 / (2 * dx)
    kx = collect(range(-freq_nyquist_x, freq_nyquist_x, N)) .* f
    freq_nyquist_y = 1 / (2 * dy)
    ky = collect(range(-freq_nyquist_y, freq_nyquist_y, N)) .* f

    xf = kx .* (λs / sinθmax * (N / M))
    yf = ky .* (λs / sinθmax * (N / M))

    I_focus = abs2.(Efx) .+ abs2.(Efy) .+ abs2.(Efz)
    E_focus = calcEnergy(xf, yf, I_focus)
    println("Energy @ focus = ", round(E_focus * 1e3, digits=3), " mJ")
    w0_x, w0_y = FWHM2D(xf, yf, I_focus)
    println("Beam spot size (FWHM) @ focus =", round(w0_x*1e6, digits=2), " μm x ", round(w0_y*1e6, digits=2), " μm")
    Aeff = calcAeff(xf, yf, I_focus)
    Ppeak = 0.94 * E_focus / 3.8e-15
    I_target = 2 * Ppeak / Aeff
    println("Peak intensity @ focus = ", round(I_target * 1e-4, digits=3), " W/cm^2")

    return Ef, xf, yf

end

function Polarization(fn_params::FN_Params, diff_params::Diffract, Pol::String, l::Int)
    @unpack N, x, y = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(Float64, N, N)
    Ey = similar(Ex)

    @show Pol

    if Pol == "P"
        Ex .= LaguerreGauss(fn_params, 0, l, 1, w)
        Ey .= zeros(N, N)
    elseif Pol == "S"
        Ex .= zeros(N, N)
        Ey .= LaguerreGauss(fn_params, 0, l, 1, w)
    elseif Pol == "D"
        Ex .= LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
        Ey .= LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
    elseif Pol == "LHC"
        Ex .= LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
        Ey .= -1im .* LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
    elseif Pol == "RHC"
        Ex .= LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
        Ey .= 1im .* LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
    elseif Pol == "Radial"
        Ex .= cosϕ .* LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
        Ey .= sinϕ .* LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
    elseif Pol == "Azimuthal"
        Ex .= -sinϕ .* LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
        Ey .= cosϕ .* LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
    end

    scaleField!(x, y, Ex, Ey, Ein)

    rp, rs = FresnelCoefficients(28)

    M00 = rp .* cosϕ.^2 .* cosθ .+ rs .* sinϕ.^2
    M01 = sinϕ .* cosϕ .* (rp .* cosθ .- rs)

    M10 = sinϕ .* cosϕ .* (rp .* cosθ .- rs)
    M11 = rs .* cosϕ.^2 .+ rp .* sinϕ.^2 .* cosθ

    M20 = rp .* sinθ .* cosϕ
    M21 = rp .* sinθ .* sinϕ

    Epx = M00 .* Ex .+ M01 .* Ey
    Epy = M10 .* Ex .+ M11 .* Ey
    Epz = M20 .* Ex .+ M21 .* Ey

    return Epx, Epy, Epz
end



function TransmissionFunction(fn_params::FN_Params, diff_params::Diffract, Pol::String, l::Int)
    @unpack sinθmax, R, aperture, θ, ϕ = diff_params
    @unpack N, x, y = fn_params
    
    # Fields
    Etx = zeros(ComplexF64, N, N)
    Ety = similar(Etx)
    Etz = similar(Etx)

    println("Effective aperture radius = ", round(R.*1e3, digits=2), " mm")

    # Apply polarization
    Epx, Epy, Epz = Polarization(fn_params, diff_params, Pol)

    # Apodization
    cosθ = cos.(θ)
    Apod = 1 ./ cosθ

    # Transmitted fields
    Etx .= Apod .* aperture .* Epx
    Ety .= Apod .* aperture .* Epy
    Etz .= Apod .* aperture .* Epz

    # Transmitted Energy
    It_trans = abs2.(Etx) + abs2.(Ety) + abs2.(Etz)
    E_trans = calcEnergy(x, y, It_trans)
    println("Energy after parabola = ", round(E_trans * 1e3, digits=3), " mJ")
    w0_x, w0_y = e22D(x, y, It_trans)
    println("Beam spot size (1/e^2) on parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
    #display(heatmap(x, y, It_trans))

    return Etx, Ety, Etz

end

function RichardsWolf(fn_params::FN_Params, diff_params::Diffract, Pol::String, z::Real, l::Int)
    @unpack sinθmax, f, kt, θ = diff_params
    @unpack N, λs, dx, dy = fn_params

    factor = -(1im * sinθmax.^2 / (f * λs))

    Etx, Ety, Etz = TransmissionFunction(fn_params, diff_params, Pol)

    # Fields
    Efx = zeros(ComplexF64, N, N)
    Efy = similar(Efx)
    Efz = similar(Efx)

    # Propagation factor
    expz = exp.(1im .* kt .* z .* cos.(θ))

    M = 2^13
    pad_size = Int(M / 2)

    Etx_vertpad = zeropad_vertical(Etx .* expz, pad_size)
    Ety_vertpad = zeropad_vertical(Ety .* expz, pad_size)
    Etz_vertpad = zeropad_vertical(Etz .* expz, pad_size)

    tempx = fftshift(fft(Etx_vertpad, 1), 1)
    tempy = fftshift(fft(Ety_vertpad, 1), 1)
    tempz = fftshift(fft(Etz_vertpad, 1), 1)

    tempx_horpad = zeropad_horizontal(tempx[pad_size+1:pad_size+N, :], pad_size)
    tempy_horpad = zeropad_horizontal(tempy[pad_size+1:pad_size+N, :], pad_size)
    tempz_horpad = zeropad_horizontal(tempz[pad_size+1:pad_size+N, :], pad_size)

    Efx .= factor .* ((fftshift(fft(tempx_horpad, 2), 2)[:, pad_size+1:pad_size+N]))
    Efy .= factor .* ((fftshift(fft(tempy_horpad, 2), 2)[:, pad_size+1:pad_size+N]))
    Efz .= factor .* ((fftshift(fft(tempz_horpad, 2), 2)[:, pad_size+1:pad_size+N]))

    Efx .= Matrix(transpose(Efx))
    Efy .= Matrix(transpose(Efy))
    Efz .= Matrix(transpose(Efz))

    Ef = [Efx, Efy, Efz]

    freq_nyquist_x = 1 / (2 * dx)
    kx = collect(range(-freq_nyquist_x, freq_nyquist_x, N)) .* f
    freq_nyquist_y = 1 / (2 * dy)
    ky = collect(range(-freq_nyquist_y, freq_nyquist_y, N)) .* f

    xf = kx .* (λs / sinθmax * (N / M))
    yf = ky .* (λs / sinθmax * (N / M))

    I_focus = abs2.(Efx) .+ abs2.(Efy) .+ abs2.(Efz)
    E_focus = calcEnergy(xf, yf, I_focus)
    println("Energy @ focus = ", round(E_focus * 1e3, digits=3), " mJ")
    w0_x, w0_y = FWHM2D(xf, yf, I_focus)
    println("Beam spot size (FWHM) @ focus =", round(w0_x*1e6, digits=2), " μm x ", round(w0_y*1e6, digits=2), " μm")
    Aeff = calcAeff(xf, yf, I_focus)
    Ppeak = 0.94 * E_focus / 3.8e-15
    I_target = 2 * Ppeak / Aeff
    println("Peak intensity @ focus = ", round(I_target * 1e-4, digits=3), " W/cm^2")

    return Ef, xf, yf

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

    scale = zmax / maximum(x) / 1.25

    if Comp == "t"

        It = I[1]
        It ./= maximum(It[:, :, z0])
#=         fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(It[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
        ax1.xlabel = "z (μm)"
        ax1.ylabel = "x (μm)"
        Colorbar(fig[1, 2], hm1)
        Makie.save("FocusField-It-xz.png", fig)
        display(fig) =#
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.8,
                                isorange = 0.1, isovalue = maximum(It[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(It[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(It[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(It[:, :, z0])*exp(-4))
        display(scene)
 
    elseif Comp == "x"

        Ix = I[2]
        Ix ./= maximum(Ix[:, :, z0])
#=         fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Ix[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
        ax1.xlabel = "z (μm)"
        ax1.ylabel = "x (μm)"
        Colorbar(fig[1, 2], hm1)
        Makie.save("FocusField-Ix-xz.png", fig)
        display(fig) =#
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.8,
                                isorange = 0.1, isovalue = maximum(Ix[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(Ix[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(Ix[:, :, z0])*exp(-4))
        display(scene)

    elseif Comp == "y"

        Iy = I[3]
        Iy ./= maximum(Iy[:, :, z0])
#=         fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iy[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
        ax1.xlabel = "z (μm)"
        ax1.ylabel = "x (μm)"
        Colorbar(fig[1, 2], hm1)
        Makie.save("FocusField-Iy-xz.png", fig)
        display(fig) =#
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.8,
                                isorange = 0.1, isovalue = maximum(Iy[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(Iy[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(Iy[:, :, z0])*exp(-4))
        display(scene)

    elseif Comp == "z"

        Iz = I[4]
        Iz ./= maximum(Iz[:, :, z0])
        #= fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iz[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
        ax1.xlabel = "z (μm)"
        ax1.ylabel = "x (μm)"
        Colorbar(fig[1, 2], hm1)
        Makie.save("FocusField-Iz-xz.png", fig)
        display(fig) =#
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.8,
                                isorange = 0.1, isovalue = maximum(Iz[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(Iz[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(Iz[:, :, z0])*exp(-4))

        display(scene)

    end
end

function Visualize3D(Pol::String, Comp::String, fn_params::FN_Params,  zmin::Real, zmax::Real, zsteps::Int, f::Real, fnum::Real, nt::Real, w::Real, l::Int)
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
                                isorange = 0.1, isovalue = maximum(It[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(It[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(It[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, It, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(It[:, :, z0])*exp(-4))
        Makie.save("It_LG_t0.png", scene)                                
        display(scene)

    elseif Comp == "x"

        Ix = I[2]
        Ix ./= maximum(Ix[:, :, z0])
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.8,
                                isorange = 0.1, isovalue = maximum(Ix[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Ix[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(Ix[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Ix, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(Ix[:, :, z0])*exp(-4))
        Makie.save("Ix_LG_t0.png", scene)
        display(scene)

    elseif Comp == "y"

        Iy = I[3]
        Iy ./= maximum(Iy[:, :, z0])
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.8,
                                isorange = 0.1, isovalue = maximum(Iy[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Iy[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(Iy[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iy, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(Iy[:, :, z0])*exp(-4))
        Makie.save("Iy_LG_t0.png", scene)
        display(scene)

    elseif Comp == "z"

        Iz = I[4]
        Iz ./= maximum(Iz[:, :, z0])
        scene = GLMakie.volume(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.8,
                                isorange = 0.1, isovalue = maximum(Iz[:, :, z0])*exp(-1))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.5,
                                isorange = 0.05, isovalue = maximum(Iz[:, :, z0])*exp(-2))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.3,
                                isorange = 0.01, isovalue = maximum(Iz[:, :, z0])*exp(-3))
        GLMakie.volume!(scale.*x.*1e6, scale.*y.*1e6, z.*1e6, Iz, algorithm = :iso, alpha = 0.1,
                                isorange = 0.001, isovalue = maximum(Iz[:, :, z0])*exp(-4))
        #Makie.save("Iz_LG_t0.png", scene)
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
            p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[1])), ylabel="y (μm)", xlabel="x (μm)", title="Ix @ z = "*string(round(z[i].*1e6, digits=2)))
            #savefig(p, "Ix_"*string(i)*".png")
            display(p)
        elseif Comp == "y"
            p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[2])), ylabel="y (μm)", xlabel="x (μm)", title="Iy @ z = "*string(round(z[i].*1e6, digits=2)))
            #savefig(p, "Iy_"*string(i)*".png")
            display(p)
        elseif Comp == "z"
            p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[3])), ylabel="y (μm)", xlabel="x (μm)", title="Iz @ z = "*string(round(z[i].*1e6, digits=2)))
            #savefig(p, "Iz_"*string(i)*".png")
            display(p)
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
