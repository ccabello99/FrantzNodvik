
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
    m::Int
    θ::Matrix{T}
    ϕ::Matrix{T}



    function Diffract{T}(fn_params::FN_Params, f::Real, fnum::Real, 
                            w::Real, nt::Real) where T
        @unpack λs, x, y, N = fn_params

        # Parameters
        R = f / (fnum * 2)
        sinθmax = R ./ sqrt(R^2 + f^2)

        k = 2π/λs
        kt = k * nt 
        
        Ein = 2.5e-3

        aperture = smooth_circular_aperture(fn_params, R)
        if N % 2 != 0
            n = Int((N-1)/2 + 1)
            m = Int((count(!iszero, aperture[n, :]) - 1) / 2)
        else
            n = Int(N/2)
            m = Int((count(!iszero, aperture[n, :])) / 2)
        end

        # Grid
        X, Y = meshgrid(x, y)
        r = sqrt.(X.^2 .+ Y.^2)
        θ = r ./ f
        ϕ = atan.(Y, X)

        new{T}(f, fnum, w, R, sinθmax, nt, kt, Ein, aperture, m, θ, ϕ)
    end

end


function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        Pol::String; verbose=false)
    @unpack N, x, y = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = similar(Ex)

    if Pol == "P"
        Ex .= Gaussian(fn_params, w, w)
        Ey .= zeros(N, N)
    elseif Pol == "S"
        Ex .= zeros(N, N)
        Ey .= Gaussian(fn_params, w, w)
    elseif Pol == "D"
        Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
        Ey .= Gaussian(fn_params, w, w) ./ sqrt(2)
    elseif Pol == "RHC"
        Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
        Ey .= -1im .* Gaussian(fn_params, w, w) ./ sqrt(2)
    elseif Pol == "LHC"
        Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
        Ey .= 1im .* Gaussian(fn_params, w, w) ./ sqrt(2)
    elseif Pol == "Radial"
        Ex .= cosϕ.* Gaussian(fn_params, w, w)
        Ey .= sinϕ .* Gaussian(fn_params, w, w)
    elseif Pol == "Azimuthal"
        Ex .= -sinϕ .* Gaussian(fn_params, w, w)
        Ey .= cosϕ .* Gaussian(fn_params, w, w)
    end

    scaleField!(x, y, Ex, Ey, Ein)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
    end

    rp, rs = FresnelCoefficients(28)

    M00 = rp .* cosϕ.^2 .* cosθ .+ rs .* sinϕ.^2
    M01 = sinϕ .* cosϕ .* (rp .* cosθ .- rs)

    M10 = sinϕ .* cosϕ .* (rp .* cosθ .- rs)
    M11 = rs .* cosϕ.^2 .+ rp .* sinϕ.^2 .* cosθ

    M20 = -rp .* sinθ .* cosϕ
    M21 = -rp .* sinθ .* sinϕ

    Epx = M00 .* Ex .+ M01 .* Ey
    Epy = M10 .* Ex .+ M11 .* Ey
    Epz = M20 .* Ex .+ M21 .* Ey

    return Epx, Epy, Epz
end

function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        Pol::String, l::Int; verbose=false)
    @unpack N, x, y = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = similar(Ex)

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
        Ex .= cosϕ.* LaguerreGauss(fn_params, 0, l, 1, w)
        Ey .= sinϕ .* LaguerreGauss(fn_params, 0, l, 1, w)
    elseif Pol == "Azimuthal"
        Ex .= -sinϕ .* LaguerreGauss(fn_params, 0, l, 1, w)
        Ey .= cosϕ .* LaguerreGauss(fn_params, 0, l, 1, w)
    end

    scaleField!(x, y, Ex, Ey, Ein)

    # Print some useful info about initial field
    if verbose
        @show Pol
        println("OAM l = ", l)
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
    end

    rp, rs = FresnelCoefficients(28)

    M00 = rp .* cosϕ.^2 .* cosθ .+ rs .* sinϕ.^2
    M01 = sinϕ .* cosϕ .* (rp .* cosθ .- rs)

    M10 = sinϕ .* cosϕ .* (rp .* cosθ .- rs)
    M11 = rs .* cosϕ.^2 .+ rp .* sinϕ.^2 .* cosθ

    M20 = -rp .* sinθ .* cosϕ
    M21 = -rp .* sinθ .* sinϕ

    Epx = M00 .* Ex .+ M01 .* Ey
    Epy = M10 .* Ex .+ M11 .* Ey
    Epz = M20 .* Ex .+ M21 .* Ey

    return Epx, Epy, Epz
end




function TransmissionFunction(fn_params::FN_Params, diff_params::Diffract, 
                                Pol::String, l::Int; verbose=false)
    @unpack sinθmax, R, aperture, θ, ϕ = diff_params
    @unpack N, x, y = fn_params
    
    # Fields
    Etx = zeros(ComplexF64, N, N)
    Ety = similar(Etx)
    Etz = similar(Etx)

    # Initialize fields and apply polarization matrix
    if l != 0
        Epx, Epy, Epz = Polarization(fn_params, diff_params, Pol, l)
    else
        Epx, Epy, Epz = Polarization(fn_params, diff_params, Pol)
    end

    # Apodization
    #Apod = sqrt.((cos.(θ)))
    Apod = 2 ./ (1 .+ cos.(θ))

    # Transmitted fields
    Etx .= Apod .* aperture .* Epx
    Ety .= Apod .* aperture .* Epy
    Etz .= Apod .* aperture .* Epz

    # Print some useful info about transmitted field
    if verbose
        println("Effective aperture radius = ", round(R.*1e3, digits=2), " mm")
        It_trans = abs2.(Etx) + abs2.(Ety) + abs2.(Etz)
        E_trans = calcEnergy(x, y, It_trans)
        println("Energy after parabola = ", round(E_trans * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, It_trans)
        println("Beam spot size (1/e2) on parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        display(heatmap(x, y, It_trans))
    end

    return Etx, Ety, Etz

end

function RichardsWolf(fn_params::FN_Params, diff_params::Diffract, 
                        Pol::String, z::Real, l::Int; verbose=false)
    @unpack sinθmax, f, R, kt, m, θ = diff_params
    @unpack N, λs, dx, dy = fn_params

    factor = -(1im * R^2 / (f * λs * m^2))

    Etx, Ety, Etz = TransmissionFunction(fn_params, diff_params, Pol, l)

    # Fields
    Efx = zeros(ComplexF64, N, N)
    Efy = similar(Efx)
    Efz = similar(Efx)

    # Propagation factor
    cosθ = cos.(θ)
    expz = exp.(1im .* kt .* z .* cosθ)

    # For zero-padding
    M = 2^12
    pad_size = Int(M/2)
    if N % 2 != 0
        n = Int((N-1)/2)
        pad_range = pad_size-n:pad_size+n
    else
        n = Int(N/2)
        pad_range = pad_size-n:pad_size+n-1
    end

    # Compute 2D FFT one dimension at a time
    Etx_vertpad = zeropad_vertical(Etx .* expz ./ cosθ, pad_size)
    Ety_vertpad = zeropad_vertical(Ety .* expz ./ cosθ, pad_size)
    Etz_vertpad = zeropad_vertical(Etz .* expz ./ cosθ, pad_size)

    tempx = fftshift(fft(fftshift(Etx_vertpad, 1), 1), 1)
    tempy = fftshift(fft(fftshift(Ety_vertpad, 1), 1), 1)
    tempz = fftshift(fft(fftshift(Etz_vertpad, 1), 1), 1)

    tempx_horpad = zeropad_horizontal(tempx[pad_range, :], pad_size)
    tempy_horpad = zeropad_horizontal(tempy[pad_range, :], pad_size)
    tempz_horpad = zeropad_horizontal(tempz[pad_range, :], pad_size)

    Efx .= factor .* fftshift(fft(fftshift(tempx_horpad, 2), 2), 2)[:, pad_range]
    Efy .= factor .* fftshift(fft(fftshift(tempy_horpad, 2), 2), 2)[:, pad_range]
    Efz .= factor .* fftshift(fft(fftshift(tempz_horpad, 2), 2), 2)[:, pad_range]

    Efx .= Matrix(transpose(Efx))
    Efy .= Matrix(transpose(Efy))
    Efz .= Matrix(transpose(Efz))

    Ef = [Efx, Efy, Efz]

    freq_nyquist_x = 1 / (2 * dx)
    kx = collect(range(-freq_nyquist_x, freq_nyquist_x, N+2)[2:end-1]) .* f
    freq_nyquist_y = 1 / (2 * dy)
    ky = collect(range(-freq_nyquist_y, freq_nyquist_y, N+2)[2:end-1]) .* f

    xf = kx .* (λs / sinθmax * (m / (M)))
    yf = ky .* (λs / sinθmax * (m / (M)))

    # Print some useful info about focus field
    if verbose
        NA = 1 / (2 * diff_params.fnum)
        absz = (m * λs * sqrt(1-NA^2) / (2*NA^2)) * 1e6
        println("Numerical aperture of system = ", round(NA, digits=2))
        println("Propagation is accurate over a distance |z| = ", round(absz, digits=2), " μm")
        I_focus = abs2.(Efx) .+ abs2.(Efy) .+ abs2.(Efz)
        E_focus = calcEnergy(xf, yf, I_focus)
        println("Energy @ focus = ", round(E_focus * 1e3, digits=3), " mJ")
        w0_x, w0_y = FWHM2D(xf, yf, I_focus)
        println("Beam spot size (FWHM) @ focus =", round(w0_x*1e6, digits=2), " μm x ", round(w0_y*1e6, digits=2), " μm")
        Aeff = calcAeff(xf, yf, I_focus)
        println("Effective area = ", round(Aeff*1e12, digits=2), " μm^2")
        Ppeak = 0.94 * E_focus / 3.8e-15
        I_target = 2 * Ppeak / Aeff
        println("Peak intensity @ focus = ", round(I_target * 1e-4, digits=3), " W/cm^2")
        println("Ratio of peak Ix to I_tot = ", maximum(abs2.(Efx))./maximum(I_focus))
        println("Ratio of peak Iy to I_tot = ", maximum(abs2.(Efy))./maximum(I_focus))
        println("Ratio of peak Iz to I_tot = ", maximum(abs2.(Efz))./maximum(I_focus))
    end

    return Ef, xf, yf

end

function FullSpatialProfile(fn_params::FN_Params, diff_params::Diffract, Pol::String, 
                                zmin::Real, zmax::Real, zsteps::Int, l::Int)
    @unpack N, λs = fn_params
    @unpack w, nt, kt = diff_params

    z = collect(range(zmin, zmax, zsteps))

    # Run once to compile and save x and y vectors
    Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, 0, l)

    Ex = zeros(ComplexF64, N, N, zsteps)
    Ey = zeros(ComplexF64, N, N, zsteps)
    Ez = zeros(ComplexF64, N, N, zsteps)

    zR = π * w^2 * nt / λs
    
    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, z[I], l)

        # Include Gouy phase
        ψg = (abs(l) + 1)*atan(z[I] / zR)

        # TODO create function to define wavefront curvature

        Ex[:, :, I] .= Ef[1] .* exp(1im * ψg)
        Ey[:, :, I] .= Ef[2] .* exp(1im * ψg)
        Ez[:, :, I] .= Ef[3] .* exp(1im * ψg)
    end

    E = [Ex, Ey, Ez]

    return E, x, y, z
end

function SpatioTemporalVectorDiffraction(fn_params::FN_Params, diff_params::Diffract, Pol::String, 
                                            zmin::Real, zmax::Real, zsteps::Int, fsteps::Int, t_now::Real, 
                                                l::Int; verbose=false)
    @unpack N, t0, ϕ0, τs, τ, ωs, nt, c = fn_params
    @unpack nt, w = diff_params

    # Define temporal profile
    scale = sqrt(1/(-2*log(0.5))) / sqrt(2)
    temp = ComplexEnvelope(1, t0, ϕ0, τs * scale, 0)
    t = collect(range(-3*τ, 3*τ, N))
    dt = t[2] - t[1]
    fs = 1/(dt)
    freq = (fs / N) .* range(-N/2, N/2, N)
    E_t = temp.(t) .* exp.(1im .* ωs .* t)
    E_ω = fftshift(fft(fftshift(E_t)))

    # Retrieve spectral profile for sampling
    spl_re = Spline1D(freq, real.(E_ω), k=3)
    spl_im = Spline1D(freq, imag.(E_ω), k=3)
    Eω(ω) = spl_re(ω) + 1im .* spl_im(ω)
    E_ω = Eω.(freq) ./ maximum(abs.(Eω.(freq)))
    E_ω = abs.(E_ω)

    # Define grid for wavelength sampling
    fmin = freq[find_first_e2(E_ω, 1e-1)]
    fmax = freq[find_last_e2(E_ω, 1e-1)]
    f_samples = collect(range(fmin, fmax, fsteps))
    dω = (f_samples[2] - f_samples[1]) * 2π
    λ_samples = collect(c ./ reverse(f_samples))

    # Run once to compile and save x and y vectors
    Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, 0, l)

    Ex = zeros(ComplexF64, N, N, zsteps, )
    Ey = zeros(ComplexF64, N, N, zsteps)
    Ez = zeros(ComplexF64, N, N, zsteps)

    z = collect(range(zmin, zmax, zsteps))

    foreach(eachindex(λ_samples)) do i

        fn_params.λs = λ_samples[i]
        k = 2π/fn_params.λs
        diff_params.kt = k * nt

        zR = π * w^2 * nt / fn_params.λs
        
        foreach(eachindex(z)) do I
            Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, z[I], l)

            # Include Gouy phase
            ψg = (abs(l) + 1)*atan(z[I] / zR)

            # Spectral Contribution
            ϕ = unwrap(angle.(Eω(f_samples[i])))
            expω = Eω(f_samples[i]) .* exp(1im * ψg) .* exp.(-1im * 2π * f_samples[i] * t_now + ϕ) * dω

            # TODO create function to define wavefront curvature

            Ex[:, :, I] .+= Ef[1] .* expω
            Ey[:, :, I] .+= Ef[2] .* expω
            Ez[:, :, I] .+= Ef[3] .* expω
        end

    end

    E = [Ex, Ey, Ez]

    return E, x, y, z

end
println("Diffraction.jl compiled")

