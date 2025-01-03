
# Create object for performing Gabor transformations (i.e. STFT)

mutable struct Gabor{T}
    ng::Int64
    nt::Int64
    t_axis::Vector{T}
    Xgt_spec::Matrix{T}
    g::Vector{T}
    gt::Vector{T}
    ft::Vector{ComplexF64}
    tslide::Vector{T}
    window::Float64
end


function Gabor(fn_params::FN_Params)
    @unpack N, t, nt, τ = fn_params
    type = typeof(τ)
    ng = N
    Xgt_spec = zeros(type, ng, nt)
    g = zeros(type, nt)
    gt = zeros(type, nt)
    ft = zeros(ComplexF64, nt)
    tslide = collect(range(t[1], t[end], ng))
    window = 0.1*τ
    Gabor{type}(ng, nt, t, Xgt_spec, g, gt, ft, tslide, window)
end


# Transform over whole signal
function gabor_transform!(A::Vector, gabor::Gabor, fft_plan)
    @unpack ng = gabor
    @simd for t in 1:ng
        gabor_transform!(A, gabor, fft_plan, t)
    end

end
    
# Transform at specific time
function gabor_transform!(A::Vector, gabor::Gabor, fft_plan, t::Int)

    @unpack Xgt_spec, t_axis, tslide, window, g, gt, ft, nt = gabor
    type = typeof(g[1])
    g .= exp.(-(t_axis .- tslide[t]).^2 / (window)^2)
    gt .= g .* real.(A)
    ft .= fftshift(fft_plan * gt)
    @inbounds Xgt_spec[t, :] .= convert.(type, abs.(ft)) / nt

end

# Transform over specific interval
function gabor_transform!(A::Vector, gabor::Gabor, fft_plan, int::Vector)

    @simd for t in int[1]:int[end]
        gabor_transform!(A, gabor, fft_plan, t)
    end

end

println("Gabor.jl compiled")
