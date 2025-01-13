struct P end
struct S end
struct D end
struct RHC end
struct LHC end
struct Radial end
struct Azimuthal end

function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        ::P, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
        Pol = "p-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    
    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.-2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
        
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        ::S, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        Pol = "s-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        ::D, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        Pol = "45°-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        ::RHC, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        Pol = "right-hand circularly-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        ::LHC, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        Pol = "left-hand circularly-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        ::Radial, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        Pol = "radially-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, 
                        ::Azimuthal, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        Pol = "azimuthally-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end

function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::P, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        OAM = l
        @show OAM
        Pol = "p-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::S, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        OAM = l
        @show OAM
        Pol = "s-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::D, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        OAM = l
        @show OAM
        Pol = "45°-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::RHC, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        OAM = l
        @show OAM
        Pol = "right-hand circularly-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real, 
                        ::LHC, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        OAM = l
        @show OAM
        Pol = "left-hand circularly-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real,
                        ::Radial, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        OAM = l
        @show OAM
        Pol = "radially-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end


function Polarization(fn_params::FN_Params, diff_params::Diffract, l::Real,
                        ::Azimuthal, Z::Vector; verbose=false, aberration=false, hole=false, magnetic=false)

    @unpack N, x, y, λs, Ein, c = fn_params
    @unpack w, θ, ϕ = diff_params

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
    rp, rs = FresnelCoefficients(12.75, λs)

    # Print some useful info about initial field
    if verbose
        OAM = l
        @show OAM
        Pol = "azimuthally-polarized"
        @show Pol
        I_in = abs2.(Ex) + abs2.(Ey)
        E_in = calcEnergy(x, y, I_in)
        println("Energy before parabola = ", round(E_in * 1e3, digits=3), " mJ")
        w0_x, w0_y = e22D(x, y, I_in)
        println("Beam spot size (1/e2) before parabola =", round(2*w0_x*1e3, digits=2), " mm x ", round(2*w0_y*1e3, digits=2), " mm")
        println("Reflection coefficients are : Rp = ", round(abs2.(rp), digits=4), " and Rs = ", round(abs2.(rs), digits=4))
    end

    rp_cosθ = -rp .* cosθ
    rs_cosθ = rs .* cosθ
    sinϕ_cosϕ = sinϕ .* cosϕ

    if magnetic == false
        
        M00 = rp_cosθ .* cosϕ.^2 .- rs .* sinϕ.^2
        M01 = sinϕ_cosϕ .* (rp_cosθ .- rs)

        M10 = M01
        M11 = rs .* cosϕ.^2 .+ rp_cosθ .* sinϕ.^2

        M20 = -rp .* sinθ .* cosϕ
        M21 = -rp .* sinθ .* sinϕ

        Epx = M00 .* Ex .+ M01 .* Ey
        Epy = M10 .* Ex .+ M11 .* Ey
        Epz = M20 .* Ex .+ M21 .* Ey

        return Epx, Epy, Epz
    else

        Hx = Ex / 376
        Hy = Ey / 376

        M00 = sinϕ_cosϕ .* (rs_cosθ .- rp)
        M01 = -rs_cosθ .* cosϕ.^2 .+ rp .* sinϕ.^2

        M10 = -M01
        M11 = -M00

        M20 = -rs .* sinθ .* sinϕ
        M21 = rs .* sinθ .* cosϕ

        Hpx = M00 .* Hx .+ M01 .* Hy
        Hpy = M10 .* Hx .+ M11 .* Hy
        Hpz = M20 .* Hx .+ M21 .* Hy

        return Hpx, Hpy, Hpz

    end
end

println("Polarization.jl compiled")