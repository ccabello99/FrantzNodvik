# (SI units for all quantities)

@with_kw mutable struct FN_Params{T}

    # Physical constants
    h::T = 6.62607015e-34 
    c::T = 2.998e8
    kb::T = 1.38e-23
    λp::T = 527e-9
    λs::T = 785e-9
    ωs::T = 2π * c / λs 
    ωp::T = 2π * c / λp
    frep::T = 1e3

    # Crystal Parameters
    α::T = 1.92e2 
    l::T = 8e-3
    A::T = (1 - exp(-α*l))
    T0::T = 150
    T1rad::T = 3.97e-6
    T1nr::T = 2.93e-9
    ΔE::T = 3.56e-20
    #ηc = (1 + (T1rad / T1nr) * exp(-ΔE / (kb * T0)))
    ηc::T = 1
    ηq::T = (λp / λs)

    # Thermal lens TODO
    #=
    dndT::T = 4.11e-6 - 1.565e-10 * T0 + 6.449e-11 * T0^2
    αth::T = 5.8e-6
    K::T = 150
    Pth::T = (A * Ep0 * frep * (1 - (ηc * ηq)^2))
    r1::T = w_xp
    r2::T = 6e-3
    fth::T = (2π * r1^2 * K) / (Pth * (dndT + ((2 * r2 * αth * (n0-1)) / l)))
    =#

    # Chirped pulse duration
    τs::T = 30e-15
    GDD::T = 208000e-30
    ϕ0::T = 0
    τ::T = τs * √(1+(4log(2)*GDD/(τs^2))^2) # Real pulse duration

    # Grid size
    N::Int = 100
    pass::Int = 13
    xmin::T = 0
    ymin::T = 0
    zmin::T = 0
    xmax::T = 80e-4
    ymax::T = 80e-4
    zmax::T = 8e-3
    dx::T = xmax / N
    dy::T = ymax / N
    dz::T = zmax / N
    Npass::Vector{T} = collect(0:pass-1)

    # Spatial domain
    x::Vector{T} = collect(range(xmin, xmax, N))
    y::Vector{T} = collect(range(ymin, ymax, N))
    z::Vector{T} = collect(range(zmin, zmax, N))
    x0::T = xmax / 2
    y0::T = ymax / 2

    # Time domain
    nt::Int = 2^12
    tmin::T = 0
    tmax::T = 4 * τ
    dt::T = tmax / nt
    t::Vector{T} = collect(LinRange(-nt/2, nt/2, nt) .* dt)
    t0::T = tmax / 2

    # Frequency domain
    fs::T = 1 / dt
    ntot::Int = Int(nt/2)
    f::Vector{T} = (fs/nt) .* range(-nt/2, nt/2, nt)
    f0::T = ωs / 2π

    # Index to start saving results after a single pass
    start = 2
end  

fn_params = FN_Params{Float64}()

println("FN-Params.jl compiled")
