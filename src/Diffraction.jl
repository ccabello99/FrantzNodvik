
mutable struct Diffract{T}

    f::T
    fnum::T
    w::T
    R::T
    sinθmax::T
    nt::T
    kt::T
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

        aperture = circular_aperture(fn_params, R)
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

        new{T}(f, fnum, w, R, sinθmax, nt, kt, aperture, m, θ, ϕ)
    end

end


function TransmissionFunction(fn_params::FN_Params, diff_params::Diffract, 
                                Pol, l::Real, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)
    @unpack sinθmax, R, aperture, θ, ϕ = diff_params
    @unpack N, x, y = fn_params
    
    # Fields
    Etx = zeros(ComplexF64, N, N)
    Ety = similar(Etx)
    Etz = similar(Etx)

    # Initialize fields and apply polarization matrix
    if l != 0
        Epx, Epy, Epz = Polarization(fn_params, diff_params, l, Pol, Z, aberration=aberration, hole=hole, verbose=verbose, magnetic=magnetic)
    else
        Epx, Epy, Epz = Polarization(fn_params, diff_params, Pol, Z, aberration=aberration, hole=hole, verbose=verbose, magnetic=magnetic)
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
        println("Beam waist (1/e2) on parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
    end

    return Etx, Ety, Etz

end

function RichardsWolf(fn_params::FN_Params, diff_params::Diffract, 
                        Pol, z::Real, l::Real, Z::Vector; fft_plan=false, 
                        verbose=false, aberration=false, hole=false,
                        fft_plan_vert=nothing, fft_plan_hor=nothing, magnetic=false)
    @unpack sinθmax, f, R, kt, m, θ = diff_params
    @unpack N, λs, dx, dy = fn_params

    factor = -(1im * R^2 / (f * λs * m^2))

    Etx, Ety, Etz = TransmissionFunction(fn_params, diff_params, Pol, l, Z, aberration=aberration, hole=hole, verbose=verbose, magnetic=magnetic)

    # Fields
    Efx = zeros(ComplexF64, N, N)
    Efy = similar(Efx)
    Efz = similar(Efx)

    # Propagation factor
    cosθ = cos.(θ)
    expz = exp.(1im .* kt .* z .* cosθ)

    # For zero-padding
    M = 2^13
    pad_size = Int(M/2)
    if N % 2 != 0
        n = Int((N-1)/2)
        pad_range = pad_size-n:pad_size+n
    else
        n = Int(N/2)
        pad_range = pad_size-n:pad_size+n-1
    end

    # Compute 2D FFT one dimension at a time
    if fft_plan

        Etx_vertpad = zeropad_vertical(Etx .* expz ./ cosθ, pad_size)
        Ety_vertpad = zeropad_vertical(Ety .* expz ./ cosθ, pad_size)
        Etz_vertpad = zeropad_vertical(Etz .* expz ./ cosθ, pad_size)

        tempx = fftshift(fft_plan_vert * fftshift(Etx_vertpad, 1), 1)
        tempy = fftshift(fft_plan_vert * fftshift(Ety_vertpad, 1), 1)
        tempz = fftshift(fft_plan_vert * fftshift(Etz_vertpad, 1), 1)

        tempx_horpad = zeropad_horizontal(tempx[pad_range, :], pad_size)
        tempy_horpad = zeropad_horizontal(tempy[pad_range, :], pad_size)
        tempz_horpad = zeropad_horizontal(tempz[pad_range, :], pad_size)

        Efx .= factor .* fftshift(fft_plan_hor * fftshift(tempx_horpad, 2), 2)[:, pad_range]
        Efy .= factor .* fftshift(fft_plan_hor * fftshift(tempy_horpad, 2), 2)[:, pad_range]
        Efz .= factor .* fftshift(fft_plan_hor * fftshift(tempz_horpad, 2), 2)[:, pad_range]

    else

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

    end

    #Efx .= Matrix(transpose(Efx))
    #Efy .= Matrix(transpose(Efy))
    #Efz .= Matrix(transpose(Efz))

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
        w0_x, w0_y = e22D(xf, yf, I_focus)
        println("Beam waist (e2) @ focus =", round(w0_x*1e6, digits=2), " μm x ", round(w0_y*1e6, digits=2), " μm")
        w0_x, w0_y = FWHM2D(xf, yf, I_focus)
        println("Beam waist (FWHM) @ focus =", round(w0_x*1e6, digits=2), " μm x ", round(w0_y*1e6, digits=2), " μm")
        Aeff = calcAeff(xf, yf, I_focus)
        println("Effective area = ", round(Aeff*1e12, digits=2), " μm^2")
        Ppeak = 0.94 * E_focus / 3.8e-15
        I_target = Ppeak / Aeff
        println("Peak intensity @ focus = ", round(I_target * 1e-4, digits=3), " W/cm^2")
        E_focus = calcEnergy(xf, yf, abs2.(Efx))
        println("Energy  in x-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
        E_focus = calcEnergy(xf, yf, abs2.(Efy))
        println("Energy  in y-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
        E_focus = calcEnergy(xf, yf, abs2.(Efz))
        println("Energy  in z-component @ focus = ", round(E_focus * 1e3, digits=5), " mJ")
    end

    return Ef, xf, yf

end

function FullSpatialProfile(fn_params::FN_Params, diff_params::Diffract, Pol, 
                                zmin::Real, zmax::Real, zsteps::Int, l::Real, Z::Vector; 
                                    coeffs = 0, aberration=false, hole=false, magnetic=false)
    @unpack N, λs = fn_params
    @unpack w, nt, kt = diff_params

    z = collect(range(zmin, zmax, zsteps))
    zR = π * w^2 * nt / λs

    # Fields
    Ex = zeros(ComplexF64, N, N, zsteps)
    Ey = zeros(ComplexF64, N, N, zsteps)
    Ez = zeros(ComplexF64, N, N, zsteps)

    # Create fft plan
    M = 2^13
    pad_size = Int(M/2)
    fft_plan_vert = plan_fft(zeropad_vertical(Ex[:, :, 33], pad_size), 1; flags=FFTW.MEASURE)
    fft_plan_hor = plan_fft(zeropad_horizontal(Ex[:, :, 33], pad_size), 2; flags=FFTW.MEASURE)

    # Run once to compile and save x and y vectors
    Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, 0, l, Z, 
                aberration=aberration, hole=hole, fft_plan=true,
                fft_plan_vert=fft_plan_vert, fft_plan_hor=fft_plan_hor, magnetic=magnetic)

    foreach(eachindex(z)) do I
        Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, z[I], l, Z, aberration=aberration, hole=hole, 
                                fft_plan=true, fft_plan_vert=fft_plan_vert, fft_plan_hor=fft_plan_hor, magnetic=magnetic)
    
        # Include Gouy phase
        ψg = (abs(l) + 1)*atan(z[I] / zR)

        Ex[:, :, I] .= Ef[1] .* exp(1im * ψg)
        Ey[:, :, I] .= Ef[2] .* exp(1im * ψg)
        Ez[:, :, I] .= Ef[3] .* exp(1im * ψg)
    end

    E = [Ex, Ey, Ez]

    return E, x, y, z
end

function SpatioTemporalVectorDiffraction(fn_params::FN_Params, diff_params::Diffract, Pol, 
                zmin::Real, zmax::Real, zsteps::Int, νsteps::Int, t_now::Real, l::Real, Z::Vector; 
                verbose=false, aberration=false, hole=false, spectdata=false)
    @unpack N, t0, ϕ0, τs, τ, ωs, nt, c = fn_params
    @unpack nt, w = diff_params

    # Define spectral profile
    if spectdata
        cd("input_data")
        λ, Iλ, ϕ = readSpect(fn_params, "sample-spect.txt")
        λ_samples = collect(range(500e-9, λ[end], 65))
        ν_samples = c ./ reverse(λ_samples)
        dν = ν_samples[2] - ν_samples[1]
        Iν = reverse(Iλ.(λ_samples))
        ϕν = reverse(ϕ.(λ_samples))
        
        norm = NumericalIntegration.integrate(ν_samples, sqrt.(Iν))
        Eν_samples = sqrt.(Iν) .* exp.(1im .* ϕν) ./ norm
        cd("..")
    else
            
        Δν = 2 * log(2) / (π * τs)
        νs = ωs/2π
        ν = collect(range((νs - 2*Δν), νs + 2*Δν, N))
        Eν(ν) = exp(-4 * log(2) * (ν - νs)^2 / (Δν^2))
        E_ν = Eν.(ν)

        # Define grid for wavelength sampling
        νmin = ν[find_first(E_ν ./ maximum(E_ν), 1e-1, "e2")]
        νmax = ν[find_last(E_ν ./ maximum(E_ν), 1e-1, "e2")]
        ν_samples = collect(range(νmin, νmax, νsteps))
        dν = ν_samples[2] - ν_samples[1]
        λ_samples = collect(c ./ reverse(ν_samples))

        # Sampled spectrum + define spectral phase
        norm = NumericalIntegration.integrate(ν_samples, Eν.(ν_samples))
        ϕ = SpectralPhase(0, 0, 15, 0, 0, ν_samples, νs)
        Eν_samples = Eν.(ν_samples) .* exp.(1im .* ϕ) ./ norm

    end

    Ex = zeros(ComplexF64, N, N, zsteps)
    Ey = zeros(ComplexF64, N, N, zsteps)
    Ez = zeros(ComplexF64, N, N, zsteps)

    # Create fft plan
    M = 2^13
    pad_size = Int(M/2)
    fft_plan_vert = plan_fft(zeropad_vertical(Ex[:, :, 33], pad_size), 1; flags=FFTW.MEASURE)
    fft_plan_hor = plan_fft(zeropad_horizontal(Ex[:, :, 33], pad_size), 2; flags=FFTW.MEASURE)

    # Run once to compile and save x and y vectors
    Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, 0, l, Z, 
                aberration=aberration, hole=hole, fft_plan=true,
                fft_plan_vert=fft_plan_vert, fft_plan_hor=fft_plan_hor)

    z = collect(range(zmin, zmax, zsteps))

    foreach(CartesianIndices((νsteps, zsteps))) do I

        # Spatial contribution
        fn_params.λs = λ_samples[I[1]]
        k = 2π/fn_params.λs
        diff_params.kt = k * nt
        zR = π * w^2 * nt / fn_params.λs
        Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, z[I[2]], l, Z, 
                    aberration=aberration, hole=hole, fft_plan=true,
                    fft_plan_vert=fft_plan_vert, fft_plan_hor=fft_plan_hor)

        # Include Gouy phase
        ψg = (abs(l) + 1)*atan(z[I[2]] / zR)

        # Spectral Contribution (normalized)
        E_ν = Eν_samples[I[1]] .* exp(1im * ψg) .* exp.(-1im * 2π * ν_samples[I[1]] * t_now) * dν

        Ex[:, :, I[2]] .+= Ef[1] .* E_ν
        Ey[:, :, I[2]] .+= Ef[2] .* E_ν
        Ez[:, :, I[2]] .+= Ef[3] .* E_ν


    end

    E = [Ex, Ey, Ez]

    return E, x, y, z

end


function SpatioTemporalVectorDiffraction(fn_params::FN_Params, diff_params::Diffract, Pol, 
                zmin::Real, zmax::Real, zsteps::Int, νsteps::Int, t_now::Real, l0::Real, l_var::String, Z::Vector; 
                verbose=false, aberration=false, hole=false, spectdata=false)

    @unpack N, t0, ϕ0, τs, τ, ωs, nt, λs, c = fn_params
    @unpack nt, w = diff_params

    # Define spectral profile
    νs = ωs/2π
    if spectdata
        cd("input_data")
        λ, Iλ, ϕ = readSpect(fn_params, "sample-spect.txt")
        λ_samples = collect(range(500e-9, λ[end], 65))
        ν_samples = c ./ reverse(λ_samples)
        dν = ν_samples[2] - ν_samples[1]
        Iν = reverse(Iλ.(λ_samples))
        ϕν = reverse(ϕ.(λ_samples))
        
        norm = NumericalIntegration.integrate(ν_samples, sqrt.(Iν))
        Eν_samples = sqrt.(Iν) .* exp.(1im .* ϕν) ./ norm
        νmin = ν_samples[find_first(abs.(Eν_samples) ./ maximum(abs.(Eν_samples)), 1e-1, "e2")]
        νmax = ν_samples[find_last(abs.(Eν_samples) ./ maximum(abs.(Eν_samples)), 1e-1, "e2")]
        cd("..")

    else
            
        Δν = 2 * log(2) / (π * τs)
        ν = collect(range((νs - 2*Δν), νs + 2*Δν, N))
        Eν(ν) = exp(-4 * log(2) * (ν - νs)^2 / (Δν^2))
        E_ν = Eν.(ν)

        # Define grid for wavelength sampling
        νmin = ν[find_first(E_ν ./ maximum(E_ν), 1e-1, "e2")]
        νmax = ν[find_last(E_ν ./ maximum(E_ν), 1e-1, "e2")]
        ν_samples = collect(range(νmin, νmax, νsteps))
        dν = ν_samples[2] - ν_samples[1]
        λ_samples = collect(c ./ reverse(ν_samples))

        # Sampled spectrum + define spectral phase
        norm = NumericalIntegration.integrate(ν_samples, Eν.(ν_samples))
        ϕ = SpectralPhase(0, 0, 15, 0, 0, ν_samples, νs)
        Eν_samples = Eν.(ν_samples) .* exp.(1im .* ϕ) ./ norm

    end
    
    if l_var == "constant"
        # Case where OAM = l0 across entire bandwidth
        l = ν -> l0
    elseif l_var == "harmonic"
        # Case where OAM = n*l0 for n ~ harmonic order
        l = ν -> l0 * (ν / νs)
    else
        error("Unknown OAM variation: $l_var")
    end

    l_samples = round.(l.(ν_samples), digits=3)

    if verbose
        println("Spectral width = ", round(FWHM(ν_samples, abs.(Eν_samples)) * λs^2 / c * 1e9, digits=2), " nm")
        println("OAM @ λmin = ", round(c/νmax*1e9, digits=2), " nm : ", l_samples[end])
        println("OAM @ λ0 = ", c/νs*1e9, " nm : ", l(νs))
        println("OAM @ λmax = ", round(c/νmin*1e9, digits=2), " nm : ", l_samples[1])
        println("Mean OAM = ", round(NumericalIntegration.integrate(ν_samples, abs.(Eν_samples) .* l_samples), digits=3))


        p = Plots.plot(ν_samples.*1e-15, abs.(Eν_samples), 
                    title="Sampled frequencies and OAM", color="blue", 
                    ylabel="Spectral Amp. (a.u.)", xlabel="Freq. (PHz)")
        ax2 = twinx()
        Plots.scatter!(ax2, ν_samples.*1e-15, l_samples, color="red")
        display(p)
    end


    Ex = zeros(ComplexF64, N, N, zsteps)
    Ey = zeros(ComplexF64, N, N, zsteps)
    Ez = zeros(ComplexF64, N, N, zsteps)

    # Create fft plan
    M = 2^13
    pad_size = Int(M/2)
    fft_plan_vert = plan_fft(zeropad_vertical(Ex[:, :, 33], pad_size), 1; flags=FFTW.MEASURE)
    fft_plan_hor = plan_fft(zeropad_horizontal(Ex[:, :, 33], pad_size), 2; flags=FFTW.MEASURE)

    # Run once to compile and save x and y vectors
    Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, 0, l0, Z, 
                aberration=aberration, hole=hole, fft_plan=true,
                fft_plan_vert=fft_plan_vert, fft_plan_hor=fft_plan_hor)

    z = collect(range(zmin, zmax, zsteps))

    foreach(CartesianIndices((νsteps, zsteps))) do I

        # Spatial contribution
        fn_params.λs = λ_samples[I[1]]
        k = 2π/fn_params.λs
        diff_params.kt = k * nt
        zR = π * w^2 * nt / fn_params.λs
        Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, z[I[2]], l_samples[I[1]], Z, 
                    aberration=aberration, hole=hole, fft_plan=true,
                    fft_plan_vert=fft_plan_vert, fft_plan_hor=fft_plan_hor)

        # Include Gouy phase
        ψg = (abs(l_samples[I[1]]) + 1)*atan(z[I[2]] / zR)

        # Spectral Contribution
        E_ν = Eν_samples[I[1]] .* exp(1im * ψg) .* exp.(-1im * 2π * ν_samples[I[1]] * t_now) * dν

        Ex[:, :, I[2]] .+= Ef[1] .* E_ν
        Ey[:, :, I[2]] .+= Ef[2] .* E_ν
        Ez[:, :, I[2]] .+= Ef[3] .* E_ν


    end

    E = [Ex, Ey, Ez]

    return E, x, y, z

end

println("Diffraction.jl compiled")

