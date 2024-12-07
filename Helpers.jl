
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

function find_first_e2(A::Vector, tolerance::Real)

    max = maximum(A)

    threshold = max * exp(-2)

    for (index, value) in enumerate(A)
        if abs(value - threshold) <= tolerance
            return index
        end
    end

    # Return nothing if no index is found (should not occur if threshold is reasonable)
    return nothing
end

function find_last_e2(A::Vector, tolerance::Real)

    max = maximum(A)

    threshold = max * exp(-2)

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


function Visualize3D(Pol::String, Comp::String, fn_params::FN_Params, diff_params::Diffract, 
                        zmin::Real, zmax::Real, zsteps::Int, l::Int; slicex=false, slicey=false, 
                            save=false, intensity=true, phase=false)
    # Provide Pol as "P", "S", "LHC", or "RHC"
    # Provide Comp as "t", "x", "y", or "z" for total intensity or x-, y-, or z-component of the intensity respectively
    @unpack N = fn_params
    
    Ef, x, y, z = FullSpatialProfile(fn_params, diff_params, Pol, zmin, zmax, zsteps, l)

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
    if slice == false
        fig = Figure(size=(940,940))
        ax1 = Axis3(fig[1, 1], aspect = :equal, xlabel = "x (μm)", ylabel = "y (μm)", zlabel = "z (μm)")
    end

    if Comp == "t"

        It = I[1]
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

        if slice
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Ix[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "x (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Ix-xz.png", fig)
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

                Ex ./= maximum(Ex[:, :, z0])

                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.8, colormap=:berlin,
                                        isorange = 0.1, isovalue = maximum(Ex[:, :, z0])*exp(-0.5))
                volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.8, colormap=:berlin,
                                        isorange = 0.1, isovalue = minimum(Ex[:, :, z0])*exp(-0.5))
                if save
                    Makie.save("FocusField-Ex-3D.png", fig)
                end    

            end
            
            display(fig)
        end

    elseif Comp == "y"
        
        if slice
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iy[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "x (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Iy-xz.png", fig)
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

                Ey ./= maximum(Ey[:, :, z0])

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
        
        if slice
            fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iz[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
            ax1.xlabel = "z (μm)"
            ax1.ylabel = "x (μm)"
            Colorbar(fig[1, 2], hm1)
            if save
                Makie.save("FocusField-Iz-xz.png", fig)
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

                Ez ./= maximum(Ez[:, :, z0])

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

function Visualize4D(Pol::String, Comp::String, fn_params::FN_Params, diff_params::Diffract, zmin::Real, zmax::Real, zsteps::Int, tmin::Real, tmax::Real, tsteps::Int, l::Int; slice=false, save=false, intensity=true, phase=false)
    @unpack N, t0, τs, ωs, ϕ0 = fn_params

    Ef, x, y, z = FullSpatialProfile(fn_params, diff_params, Pol, zmin, zmax, zsteps, l)

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

        ϕx .= angle.(Ef[1])
        ϕy .= angle.(Ef[2])
        ϕz .= angle.(Ef[3])

    end

    scale = sqrt(1/(-2*log(0.5))) / sqrt(2)
    temp = ComplexEnvelope(1, t0, ϕ0, τs * scale, 0)
    t = collect(LinRange(tmin, tmax, tsteps))
    E_t = real.(temp.(t) .* exp.(1im .* ωs .* t))
    I_t = abs2.(E_t)

    z0 = Int(zsteps/2)
    M = Int((fn_params.N-1) / 2)

    It_max = maximum(It[:, :, z0])
    Ix_max = maximum(Ix[:, :, z0])
    Iy_max = maximum(Iy[:, :, z0])
    Iz_max = maximum(Iz[:, :, z0])


    foreach(eachindex(t)) do i
        if Comp == "t"

            It ./= It_max
    
            if slice
                fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(It[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
                ax1.xlabel = "z (μm)"
                ax1.ylabel = "x (μm)"
                Colorbar(fig[1, 2], hm1)
                if save
                    Makie.save("FocusField-It-xz.png", fig)
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
    
            Ix ./= Ix_max
    
            if slice
                fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Ix[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
                ax1.xlabel = "z (μm)"
                ax1.ylabel = "x (μm)"
                Colorbar(fig[1, 2], hm1)
                if save
                    Makie.save("FocusField-Ix-xz.png", fig)
                end
                display(fig)
            else
                if intensity
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
                else
                    Ex = abs.(E[1])
                    Ex ./= maximum(Ex[:, :, z0])
    
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.8, colormap=:jet,
                                            isorange = 0.1, isovalue = maximum(Ex[:, :, z0])*exp(-0.25))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.5, colormap=:jet,
                                            isorange = 0.05, isovalue = maximum(Ex[:, :, z0])*exp(-0.5))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.3, colormap=:jet,
                                            isorange = 0.01, isovalue = maximum(Ex[:, :, z0])*exp(-1))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ex, algorithm = :iso, alpha = 0.1, colormap=:jet,
                                            isorange = 0.001, isovalue = maximum(Ex[:, :, z0])*exp(-2))
                    if save
                        Makie.save("FocusField-Ex-3D.png", fig)
                    end    
                end
                
                display(fig)
            end
    
        elseif Comp == "y"
    
            Iy ./= Iy_max
            
            if slice
                fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iy[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
                ax1.xlabel = "z (μm)"
                ax1.ylabel = "x (μm)"
                Colorbar(fig[1, 2], hm1)
                if save
                    Makie.save("FocusField-Iy-xz.png", fig)
                end
                display(fig)
            else
                if intensity
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
                    Ey = abs.(E[2])
                    Ey ./= maximum(Ey[:, :, z0])
    
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ey, algorithm = :iso, alpha = 0.8,
                                            isorange = 0.1, isovalue = maximum(Ey[:, :, z0])*exp(-0.25))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ey, algorithm = :iso, alpha = 0.5,
                                            isorange = 0.05, isovalue = maximum(Ey[:, :, z0])*exp(-0.5))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ey, algorithm = :iso, alpha = 0.3,
                                            isorange = 0.01, isovalue = maximum(Ey[:, :, z0])*exp(-1))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ey, algorithm = :iso, alpha = 0.1,
                                            isorange = 0.001, isovalue = maximum(Ey[:, :, z0])*exp(-2))
                    if save
                        Makie.save("FocusField-Ey-3D.png", fig)
                    end    
                end
                display(fig)
            end
    
        elseif Comp == "z"
    
            Iz ./= Iz_max
            
            if slice
                fig, ax1, hm1 = CairoMakie.heatmap(z.*1e6, x.*1e6, transpose(Iz[M, :, :]), colorscale=log10, colormap=:jet, colorrange=(1e-7, 1), interpolate=false)
                ax1.xlabel = "z (μm)"
                ax1.ylabel = "x (μm)"
                Colorbar(fig[1, 2], hm1)
                if save
                    Makie.save("FocusField-Iz-xz.png", fig)
                end
                display(fig)
            else
                if intensity
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
                    Ez = abs.(E[3])
                    Ez ./= maximum(Ez[:, :, z0])
    
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ez, algorithm = :iso, alpha = 0.8,
                                            isorange = 0.1, isovalue = maximum(Ez[:, :, z0])*exp(-0.25))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ez, algorithm = :iso, alpha = 0.5,
                                            isorange = 0.05, isovalue = maximum(Ez[:, :, z0])*exp(-0.5))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ez, algorithm = :iso, alpha = 0.3,
                                            isorange = 0.01, isovalue = maximum(Ez[:, :, z0])*exp(-1))
                    volume!(ax1, x.*1e6, y.*1e6, z.*1e6, Ez, algorithm = :iso, alpha = 0.1,
                                            isorange = 0.001, isovalue = maximum(Ez[:, :, z0])*exp(-2))
                    if save
                        Makie.save("FocusField-Ez-3D.png", fig)
                    end    
                end
                display(fig)
            end
        end
    end
end


function DiffractionMovie(Pol::String, Comp::String, fn_params::FN_Params, diff_params::Diffract, zmin::Real, zmax::Real, zsteps::Int, l::Int; save=false, intensity=true, phase=false)
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
            if intensity == false || phase == false || phase == true
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
