struct P end
struct S end
struct D end
struct RHC end
struct LHC end
struct Radial end
struct Azimuthal end

function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        ::P, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= Gaussian(fn_params, w, w)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, 0)
        Ex .*= exp.(1im.*Φxy_x)
    end

    if hole
        HoleyMirror!(fn_params, -w, 0, 10.25e-3, Ex)
    end

    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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
                        ::S, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ey .= Gaussian(fn_params, w, w)

    if aberration
        Φxy_y = Zernike(fn_params, Ey, Z, 0)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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
                        ::D, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
    Ey .= Gaussian(fn_params, w, w) ./ sqrt(2)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, 0)
        Φxy_y = Zernike(fn_params, Ey, Z, 0)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end

    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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
                        ::RHC, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
    Ey .= -1im .* Gaussian(fn_params, w, w) ./ sqrt(2)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, 0)
        Φxy_y = Zernike(fn_params, Ey, Z, 0)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
        end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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
                        ::LHC, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= Gaussian(fn_params, w, w) ./ sqrt(2)
    Ey .= 1im .* Gaussian(fn_params, w, w) ./ sqrt(2)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, 0)
        Φxy_y = Zernike(fn_params, Ey, Z, 0)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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
                        ::Radial, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= cosϕ.* Gaussian(fn_params, w, w)
    Ey .= sinϕ .* Gaussian(fn_params, w, w)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, 0)
        Φxy_y = Zernike(fn_params, Ey, Z, 0)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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
                        ::Azimuthal, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= -sinϕ .* Gaussian(fn_params, w, w)
    Ey .= cosϕ .* Gaussian(fn_params, w, w)
    
    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, 0)
        Φxy_y = Zernike(fn_params, Ey, Z, 0)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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

function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::P, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= LaguerreGauss(fn_params, 0, l, 1, w)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, l)
        Ex .*= exp.(1im.*Φxy_x)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
    end

    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::S, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ey .= LaguerreGauss(fn_params, 0, l, 1, w)

    if aberration
        Φxy_y = Zernike(fn_params, Ey, Z, l)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::D, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
    Ey .= LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, l)
        Φxy_y = Zernike(fn_params, Ey, Z, l)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end

    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::RHC, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
    Ey .= -1im .* LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, l)
        Φxy_y = Zernike(fn_params, Ey, Z, l)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
        end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::LHC, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)
    Ey .= 1im .* LaguerreGauss(fn_params, 0, l, 1, w) ./ sqrt(2)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, l)
        Φxy_y = Zernike(fn_params, Ey, Z, l)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real,
                        ::Radial, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= cosϕ.* LaguerreGauss(fn_params, 0, l, 1, w)
    Ey .= sinϕ .* LaguerreGauss(fn_params, 0, l, 1, w)

    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, l)
        Φxy_y = Zernike(fn_params, Ey, Z, l)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real,
                        ::Azimuthal, Z::Vector; verbose=false, aberration=false, hole=false)
    @unpack N, x, y, λs = fn_params
    @unpack w, θ, ϕ, Ein = diff_params

    cosθ = cos.(θ)
    sinθ = sin.(θ)
    cosϕ = cos.(ϕ)
    sinϕ = sin.(ϕ)

    Ex = zeros(ComplexF64, N, N)
    Ey = zeros(ComplexF64, N, N)

    Ex .= -sinϕ .* LaguerreGauss(fn_params, 0, l, 1, w)
    Ey .= cosϕ .* LaguerreGauss(fn_params, 0, l, 1, w)
    
    if aberration
        Φxy_x = Zernike(fn_params, Ex, Z, l)
        Φxy_y = Zernike(fn_params, Ey, Z, l)
        Ex .*= exp.(1im.*Φxy_x)
        Ey .*= exp.(1im.*Φxy_y)
    end

    if hole
        Ex = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ex)
        Ey = HoleyMirror(fn_params, -w, 0, 10.25e-3, Ey)
    end


    scaleField!(x, y, Ex, Ey, Ein)
    rp, rs = FresnelCoefficients(0, λs)

    # Print some useful info about initial field
    if verbose
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficiencts are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

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

println("Polarization.jl compiled")