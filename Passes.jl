

function createFluenceArrays(fn_params::FN_Params, w_xs, w_ys, w_xp)

    @unpack x0, y0, N, x, y = fn_params

    type = typeof(x0)
    Jin = abs2.(Gaussian(fn_params, w_xs, w_ys))
    Jout = zeros(type, N, N)

    n_sg = 3
    supergauss = abs2.(SuperGaussian(fn_params, w_xp, n_sg))
    gaussian_p = abs2.(Gaussian(fn_params, w_xp, w_xp))
    scale = NumericalIntegration.integrate((x,y), gaussian_p) / NumericalIntegration.integrate((x,y), supergauss)
    Jsto = scale .* supergauss
    
    return Jin, Jout, Jsto
end


function createIntensityArray(A::Function, p::Int, t::Vector, Jin0::Number)

    I(tt) = abs2.(A(tt))
    n = length(t)
    It = zeros(eltype(Jin0), p, n)
    It[1,:] .= I.(t)
    I0 = Jin0 / NumericalIntegration.integrate(t, It[1, :])
    It[1,:] .*= I0

    return It

end

function createPassArrays(fn_params::FN_Params, I::Matrix, t::Vector, Jin::Matrix, Jsto::Matrix)

    @unpack pass, x, y = fn_params

    Epass = zeros(pass)
    Esto = zeros(pass)
    B = zeros(pass)
    extraction = zeros(pass)
    Epass[1] = NumericalIntegration.integrate((x,y), NumericalIntegration.integrate(t, I[1, :]) .* Jin)
    Esto[1] = NumericalIntegration.integrate((x,y), Jsto)
    B[1] = 0
    extraction[1] = 0

    return Epass, Esto, B, extraction

end

function conserveEnergy(Jout::Matrix, Aeff::Number, tslide::Vector, It::Matrix, p::Int, dgt::Number, J::Matrix, x::Vector, y::Vector)

    Ein0 = NumericalIntegration.integrate((x,y), Jout)
    Jin0 = (2 * Ein0) / Aeff
    I0 = Jin0 / NumericalIntegration.integrate(tslide, It[p,:] ./ maximum(It[p,:]))
    It[p,:] .*= I0 ./ maximum(It[p,:])
    newEin0 = NumericalIntegration.integrate((x,y), sum(It[p, :] * dgt) .* J)
    scaling = Ein0 / newEin0
    J .*= scaling

    return J, scaling

end

# Find index/indices for the maximum value of A
function findMaxIndex(A)
    return findfirst(x -> x == maximum(A), A)
end

# Get instantaneous wavelength
function instWavelength(gabor::Gabor, A::Vector, t::Number, λ::Vector, fn_params::FN_Params, fft_plan)

    @unpack f, c, λs  = fn_params

    # Windowed FFT at specific time t (can be over a specific interval)
    gabor_transform!(A, gabor, fft_plan, t)

    # Zero out half of freq. axis to avoid findMaxIndex() confusion
    gabor.Xgt_spec[:,Int(length(f)/2):end] .= 0

    # Find current max wavelength
    diff = (c / λs) - (1e13)
    index = findMaxIndex(gabor.Xgt_spec[t,:]')[2]
    λ .= c / abs(f[index] .- diff)

    return λ

end

function one_pass(fn_params::FN_Params, fft_plan, gabor::Gabor, Jin::Matrix, Jin0::Number, Jsto::Matrix, At::Vector, 
                    Jsat::Function, Epass::Vector, Esto::Vector, B::Vector, It::Matrix, Aeffs::Vector, extraction::Vector, 
                        profile::String, p::Int, stop::Int, w::Real, w_xs::Real, w_ys::Real, n2::Function; visualize=false)

    type = typeof(Jin0)
    J_in = similar(Jin)
    Jout = similar(Jin)
    Jout .= 0
    J_out = similar(Jin)
    Jinn = similar(Jin)
    Jinn .= 0
    G0 = similar(Jin)
    int = similar(Jin)
    Isav = similar(It[1,:])
    λ = zeros(1)

    @unpack τ, c, λs, x, y, z, N, t, start = fn_params
    @unpack tslide, ng = gabor
    dgt = tslide[2] - tslide[1]
    f = 2π*1e13
    At_g = At.*exp.(-1im*f.*t)

    for t in 1:ng
        λ .= instWavelength(gabor, At_g, t, λ, fn_params, fft_plan)
        J_sat = Jsat(λ[1])
        G0 .= exp.(Jsto ./ J_sat)
        J_in .= It[p-1, t] .* Jin * dgt
        Jinn .+= J_in
        int .= log.(1 .+ (exp.(J_in ./ J_sat) .- 1) .* (G0))
        J_out .= J_sat .* int
        Jout .+= J_out
        Eout = NumericalIntegration.integrate((x,y), J_out)
        Isav[t] = Eout
        Epass[p] += Eout
        Jsto .-= (J_out .- J_in)
        if visualize && t%5 == 0
            display(heatmap(x.*1e3, y.*1e3, Jsto.*1e-4, clims=(0, maximum(Jsto).*1e-4), xlabel="x (mm)", ylabel="y (mm)", clabel="Stored Fluence (J/cm^2)", title="Stored Fluence (J/cm^2) during pass # " * string(p-1) * " at " *  @sprintf("%.1f", t .* dgt*1e12) * " ps"))
            #display(heatmap(x.*1e3, y.*1e3, Jinn.*1e-4, ylims = (2, 6), xlims=(2,6), clims=(0, Jin0.*1e-4), xlabel="x (mm)", ylabel="y (mm)", clabel="Input Fluence (J/cm^2)", title="Input Fluence (J/cm^2) during pass # " * string(p-1) * " at " * @sprintf("%.1f", t .* dgt*1e12) * " ps"))
            #display(heatmap(x.*1e3, y.*1e3, G0, clims=(0, maximum(G0)), xlabel="x (mm)", ylabel="y (mm)", clabel="Gain", title="Gain during pass # " * string(p-1) * " at " * @sprintf("%.1f", t .* dgt*1e12) * " ps"))
            #savefig("Jsto_t"*string(t)*".png")
        end
    end

    Esto0 = NumericalIntegration.integrate((x,y), Jsto)
    Esto[p] = Esto0
    Aeff_s = calcAeff(x, y, Jin)
    Aeffs[p] = Aeff_s

    Isav ./= maximum(Isav) 
    Jin0 = (2 * Epass[p]) / Aeff_s
    I0 = Jin0 / NumericalIntegration.integrate(tslide, Isav)
    It[p,:] .= I0 .* Isav

    extraction[p] = Epass[p] / Esto[1]
    

    Ppeak = (0.94 * Epass[p] / τ) .* Isav
    ϕmax = 2π/λs * NumericalIntegration.integrate(c.*z, n2(λs) * (2 * Ppeak / Aeff_s), SimpsonEven())
    B[p] = B[p-1] + 2 * ϕmax

    if profile == "gauss"
        Jin = ifelse(p < stop, abs2.(Gaussian(fn_params, w_xs, w_ys)), abs2.(Gaussian(fn_params, w, w)))
        Jin, scaling = conserveEnergy(Jout, Aeff_s, tslide, It, p, dgt, Jin, x, y)
    elseif profile == "LG"
        P, L, a = 0, 1, 1
        Jin = abs2.(LaguerreGauss(fn_params, P, L, a, w))
        Jin, scaling = conserveEnergy(Jout, Aeff_s, tslide, It, p, dgt, Jin, x, y)
    end

    Jin .*= 0.95

    return Jin, Jin0, Jout, Jsto, Epass, Esto, B, It, Aeffs, extraction

end



function several_passes(fn_params::FN_Params, fft_plan, gabor::Gabor, start::Int, stop::Int, Jin::Matrix, Jin0::Number, Jout::Matrix, Jsto::Matrix, At::Vector,
                        Jsat::Function, Epass::Vector, Esto::Vector, B::Vector, It::Matrix, Aeffs::Vector, extraction::Vector, 
                                profile::String, w::Real, w_xs::Real, w_ys::Real, n2::Function, visualize::Bool)

    prof = profile == "gauss" ? "gauss" : "LG"

    for p in start:stop
        Jin, Jin0, Jout, Jsto, Epass, Esto, B, It = one_pass(fn_params, fft_plan, gabor, Jin, Jin0, Jsto, At, Jsat, Epass, Esto, B, It, Aeffs, extraction, prof, p, stop, w, w_xs, w_ys, n2; visualize)
    end

    return Jin, Jin0, Jout, Jsto, Epass, Esto, B, It, Aeffs, extraction
end


println("Passes.jl compiled")