
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

    Aeff = calcAeff(x, y, Ex .+ Ey) # Initial effective area
    Jin = 2 * Ein / Aeff # Initial pulse fluence
    Ex .*= sqrt(Jin)
    Ey .*= sqrt(Jin)
    
end

function calcEnergy(x::Vector, y::Vector, I_tot::Matrix)

    En = round(NumericalIntegration.integrate((x, y), I_tot), digits=5)

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
        w0_x = FWHM(x, A[Int((N-1)/2+1), :]) / 2
        w0_y = FWHM(y, A[:, Int((N-1)/2+1)]) / 2
    else
        w0_x = FWHM(x, A[Int(N/2), :]) / 2
        w0_y = FWHM(y, A[:, Int(N/2)]) / 2
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

function HoleyMirror(fn_params::FN_Params, x0::Real, y0::Real, R::Real, E::Matrix)

    @unpack x, y = fn_params                        

    X, Y = meshgrid(x, y)
    X .-= x0
    Y .-= y0

    E[(X.^2 .+ Y.^2) .< R^2] .= 0

    return Matrix(transpose(E))
end

function readRefInd(RefInd::String)

    data = CSV.read(RefInd, DataFrame)
    wl = data[!, 1]
    n = data[!, 2]
    κ = data[!, 3]

    Respl = Spline1D(wl, n, k=3)
    Imspl = Spline1D(wl, κ, k=3) 

    𝑁(λ) = Respl(λ) + 1im * Imspl(λ)

    return 𝑁

end

function FresnelCoefficients(θi::Real, λ0::Real)

    # Provide θi in degrees
    sinθi = sind(θi)
    cosθi = cosd(θi)

    # For silver mirror @ 785 nm
    𝑁 = readRefInd("Ag-RefInd.csv")(λ0)
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
        7 => () -> (2 * sqrt(2) .* r .* (3 .* r.^2 .* (ϵ^2 + 1) - 2 * (ϵ^4 + ϵ^2 + 1)) .* sin.(ϕ) ./ sqrt((ϵ^2 - 1)^2*(ϵ^2 + 1)*(ϵ^4 + 4*ϵ^2 + 1))),
        8 => () -> (2 * sqrt(2) .* r .* (3 .* r.^2 .* (ϵ^2 + 1) - 2 * (ϵ^4 + ϵ^2 + 1)) .* cos.(ϕ) ./ sqrt((ϵ^2 - 1)^2*(ϵ^2 + 1)*(ϵ^4 + 4*ϵ^2 + 1))),
        9 => () -> (2 * sqrt(2) .* r.^3 .* sin.(3 .* ϕ) ./ sqrt(ϵ^6 + ϵ^4 + ϵ^2 + 1)),
        10 => () -> (2 * sqrt(2) .* r.^3 .* cos.(3 .* ϕ) ./ sqrt(ϵ^6 + ϵ^4 + ϵ^2 + 1)),
        11 => () -> (sqrt(5) * (6 .* r.^4 .- 6 * (ϵ^2 + 1) .* r.^2 .+ ϵ^4 + 4 * ϵ^2 + 1) ./ (ϵ^2 - 1)^2)
    )

    Z_tot = sum(Z[i] .* Zernike_terms[i]() for i in eachindex(Z) if Z[i] != 0)
    Z_tot[r .> 1] .= 0

    Ab = resize_symmetric(Z_tot, N)
    X, Y = meshgrid(fn_params.x, fn_params.y)
    Ab[(X.^2 .+ Y.^2) .< ϵ^2] .= 0

    return Matrix(transpose(Ab))

end

function Visualize3D(Pol::String, Comp::String, fn_params::FN_Params, diff_params::Diffract, 
                        zmin::Real, zmax::Real, zsteps::Int; slicex=false, slicey=false, l = 0, 
                            coeffs = 0, save=false, intensity=true, phase=false)
    # Provide Pol as "P", "S", "LHC", or "RHC"
    # Provide Comp as "t", "x", "y", or "z" for total intensity or x-, y-, or z-component of the intensity respectively
    @unpack N = fn_params
    
    if typeof(l) == Vector{Float64}
        Ef, x, y, z = FullSpatialProfile(fn_params, diff_params, Pol, zmin, zmax, zsteps, l=l, coeffs=coeffs)
    else
        Ef, x, y, z = FullSpatialProfile(fn_params, diff_params, Pol, zmin, zmax, zsteps, l=l)
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

function DiffractionMovie(Pol::String, Comp::String, fn_params::FN_Params, diff_params::Diffract, 
                            zmin::Real, zmax::Real, zsteps::Int, l::Real; save=false, intensity=true, phase=false)
    @unpack kt, m, w = diff_params

    z = collect(range(zmin, zmax, zsteps))
    zR = π*w^2 / fn_params.λs

    Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, 0, l)

    It_max = maximum(abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
    Ix_max = maximum(abs2.(Ef[1]))
    Iy_max = maximum(abs2.(Ef[2]))
    Iz_max = maximum(abs2.(Ef[3]))

    Ex_max = maximum(abs.(Ef[1]))
    Ey_max = maximum(abs.(Ef[2]))
    Ez_max = maximum(abs.(Ef[3]))

    for i in eachindex(z)
        Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, z[i], l);

        ψg = (abs(l) + 1)*atan(z[i] / zR)
        Ef[1] .*= exp(1im * ψg)
        Ef[2] .*= exp(1im * ψg)
        Ef[3] .*= exp(1im * ψg)

        if Comp == "t"
            if (intensity == false && phase == false) || (intensity == false && phase == true)
                break
            end
            p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3])), ylabel="y (μm)", xlabel="x (μm)", 
            title="It @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", clims=(0, It_max))
            if save
                savefig(p, "It_"*string(i)*".png")
            end
            display(p)
        elseif Comp == "x"
            if intensity
                p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[1])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Ix @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", clims=(0, Ix_max))
                if save
                    savefig(p, "Ix_"*string(i)*".png")
                end
            elseif phase
                ϕ = angle.(Ef[1])
                mag = abs.(Ef[1])
                ϕ[mag .< exp(-2).*maximum(mag)] .= 0
                p = Plots.heatmap(x.*1e6, y.*1e6, ϕ, ylabel="y (μm)", xlabel="x (μm)", title="ϕx @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:jet, clims=(-π, π))
                if save
                    savefig(p, "phix_"*string(i)*".png")
                end
            else
                p = Plots.heatmap(x.*1e6, y.*1e6, (real.(Ef[1])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Ex @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:berlin, clims=(-Ex_max, Ex_max))
                if save
                    savefig(p, "Ex_"*string(i)*".png")
                end
            end
            display(p)

        elseif Comp == "y"
            if intensity
                p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[2])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Iy @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", clims=(0, Iy_max))
                if save
                    savefig(p, "Iy_"*string(i)*".png")
                end
            elseif phase
                ϕ = angle.(Ef[2])
                mag = abs.(Ef[2])
                ϕ[mag .< exp(-2).*maximum(mag)] .= 0
                p = Plots.heatmap(x.*1e6, y.*1e6, ϕ, ylabel="y (μm)", xlabel="x (μm)", title="ϕy @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:jet, clims=(-π, π))
                if save
                    savefig(p, "phiy_"*string(i)*".png")
                end                
            else
                p = Plots.heatmap(x.*1e6, y.*1e6, (real.(Ef[2])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Ey @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:berlin, clims=(-Ey_max, Ey_max))
                if save
                    savefig(p, "Ey_"*string(i)*".png")
                end
            end
            display(p)

        elseif Comp == "z"
            if intensity
                p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[3])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Iz @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", clims=(0, Iz_max))
                if save
                    savefig(p, "Iz_"*string(i)*".png")
                end
            elseif phase
                ϕ = angle.(Ef[3])
                mag = abs.(Ef[3])
                ϕ[mag .< exp(-2).*maximum(mag)] .= 0
                p = Plots.heatmap(x.*1e6, y.*1e6, ϕ, ylabel="y (μm)", xlabel="x (μm)", title="ϕz @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:jet, clims=(-π, π))
                if save
                    savefig(p, "phiz_"*string(i)*".png")
                end                
            else
                p = Plots.heatmap(x.*1e6, y.*1e6, (real.(Ef[3])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Ez @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:berlin, clims=(-Ez_max, Ez_max))
                if save
                    savefig(p, "Ez_"*string(i)*".png")
                end
            end
            display(p)
            
        end
        
    end
end


function DiffractionMovie(Pol::String, Comp::String, fn_params::FN_Params, diff_params::Diffract, 
                            zmin::Real, zmax::Real, zsteps::Int, l::Vector, coeffs::Vector; save=false, intensity=true, phase=false)
    @unpack kt, m, w = diff_params

    z = collect(range(zmin, zmax, zsteps))
    zR = π*w^2 / fn_params.λs

    Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, 0, l, coeffs)

    It_max = maximum(abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3]))
    Ix_max = maximum(abs2.(Ef[1]))
    Iy_max = maximum(abs2.(Ef[2]))
    Iz_max = maximum(abs2.(Ef[3]))

    Ex_max = maximum(abs.(Ef[1]))
    Ey_max = maximum(abs.(Ef[2]))
    Ez_max = maximum(abs.(Ef[3]))

    for i in eachindex(z)
        Ef, x, y = RichardsWolf(fn_params, diff_params, Pol, z[i], l, coeffs);

        #ψg = (abs(l) + 1)*atan(z[i] / zR)
        ψg = 0
        Ef[1] .*= exp(1im * ψg)
        Ef[2] .*= exp(1im * ψg)
        Ef[3] .*= exp(1im * ψg)

        if Comp == "t"
            if (intensity == false && phase == false) || (intensity == false && phase == true)
                break
            end
            p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[1]) .+ abs2.(Ef[2]) .+ abs2.(Ef[3])), ylabel="y (μm)", xlabel="x (μm)", 
            title="It @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", clims=(0, It_max))
            if save
                savefig(p, "It_"*string(i)*".png")
            end
            display(p)
        elseif Comp == "x"
            if intensity
                p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[1])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Ix @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", clims=(0, Ix_max))
                if save
                    savefig(p, "Ix_"*string(i)*".png")
                end
            elseif phase
                ϕ = angle.(Ef[1])
                mag = abs.(Ef[1])
                ϕ[mag .< exp(-2).*maximum(mag)] .= 0
                p = Plots.heatmap(x.*1e6, y.*1e6, ϕ, ylabel="y (μm)", xlabel="x (μm)", title="ϕx @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:jet, clims=(-π, π))
                if save
                    savefig(p, "phix_"*string(i)*".png")
                end
            else
                p = Plots.heatmap(x.*1e6, y.*1e6, (real.(Ef[1])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Ex @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:berlin, clims=(-Ex_max, Ex_max))
                if save
                    savefig(p, "Ex_"*string(i)*".png")
                end
            end
            display(p)

        elseif Comp == "y"
            if intensity
                p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[2])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Iy @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", clims=(0, Iy_max))
                if save
                    savefig(p, "Iy_"*string(i)*".png")
                end
            elseif phase
                ϕ = angle.(Ef[2])
                mag = abs.(Ef[2])
                ϕ[mag .< exp(-2).*maximum(mag)] .= 0
                p = Plots.heatmap(x.*1e6, y.*1e6, ϕ, ylabel="y (μm)", xlabel="x (μm)", title="ϕy @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:jet, clims=(-π, π))
                if save
                    savefig(p, "phiy_"*string(i)*".png")
                end                
            else
                p = Plots.heatmap(x.*1e6, y.*1e6, (real.(Ef[2])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Ey @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:berlin, clims=(-Ey_max, Ey_max))
                if save
                    savefig(p, "Ey_"*string(i)*".png")
                end
            end
            display(p)

        elseif Comp == "z"
            if intensity
                p = Plots.heatmap(x.*1e6, y.*1e6, (abs2.(Ef[3])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Iz @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", clims=(0, Iz_max))
                if save
                    savefig(p, "Iz_"*string(i)*".png")
                end
            elseif phase
                ϕ = angle.(Ef[3])
                mag = abs.(Ef[3])
                ϕ[mag .< exp(-2).*maximum(mag)] .= 0
                p = Plots.heatmap(x.*1e6, y.*1e6, ϕ, ylabel="y (μm)", xlabel="x (μm)", title="ϕz @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:jet, clims=(-π, π))
                if save
                    savefig(p, "phiz_"*string(i)*".png")
                end                
            else
                p = Plots.heatmap(x.*1e6, y.*1e6, (real.(Ef[3])), ylabel="y (μm)", xlabel="x (μm)", 
                title="Ez @ z = "*string(round(z[i].*1e6, digits=2))*" μm from focal plane", cmap=:berlin, clims=(-Ez_max, Ez_max))
                if save
                    savefig(p, "Ez_"*string(i)*".png")
                end
            end
            display(p)
            
        end
        
    end
end


println("Helpers.jl compiled")
