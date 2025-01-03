using CSV, Tables, DataFrames, Printf, JLD, Plots
using LinearAlgebra, FFTW, ForwardDiff, NumericalIntegration, DSP
using ThreadsX, Parameters, Dierckx, ClassicalOrthogonalPolynomials, SpecialFunctions

include("FN-Params.jl")
include("Gabor.jl")
include("CrystalProperties.jl")
include("EM-Field.jl")
include("Passes.jl")
include("Diffraction.jl")
include("Polarization.jl")
include("Helpers.jl")

FFTW.set_provider!("fftw")


function ParamScan(fn_param::FN_Params, w_init::Real, w_exp::Real, wp_init::Real, Ep_init::Real, Ein_init::Real, steps::Int, LG::Bool, visualize::Bool)

    @unpack t0, c, ϕ0, τs, GDD, h, ωs, t, A, ηc, ηq, x0, y0, N, x, y, start, pass = fn_param
    

    # Sapphire properties
    sapphire = SapphireSellmeier()
    n0 = SaphRefractiveIndex(sapphire)
    n2 = SaphKerrRefractiveIndex("input_data/n2_sapphire.csv")
    σa, σe = getCrossSections("input_data/absorption_crosssection.csv", "input_data/emission_crosssection.csv")


    # Wavelength-dependent saturation fluence
    Jsat(λ) = h*(ωs/(2π)) / σe(λ)

    # FFT plan to speed up Gabor transforms
    scale = sqrt(1/(-2*log(0.5)))*sqrt(1/2) # scale such that pulse duration @ FWHM = τs
    At = ComplexEnvelope(1, t0, ϕ0, τs*scale, GDD)
    fft_plan = plan_fft(At.(t); flags=FFTW.MEASURE)
    
    # Arrays for post-processing
    Epasses = zeros(steps, steps, steps, steps, steps, steps, pass)
    Epass_max = zeros(steps, steps, steps, steps, steps, steps)
    pass_num = zeros(steps, steps, steps, steps, steps, steps)
    extract_all = zeros(steps, steps, steps, steps, steps, steps, pass)
    B_max = zeros(steps, steps, steps, steps, steps, steps)
    B_all = zeros(steps, steps, steps, steps, steps, steps, pass)
    Aeff_all = zeros(steps, steps, steps, steps, steps, steps, pass)
    Aeff_max = zeros(steps, steps, steps, steps, steps, steps)


    # LG Mode parameters
    if LG
        P = 0
        L = 1
        a = 1
    end

    lk = Threads.SpinLock()

    # Multi-threaded nested for loops over initial pump and seed waists and expanded seed waist
    ThreadsX.foreach(CartesianIndices((steps, steps, steps, steps, steps, steps)); simd=true) do I
        # Initialize all computational variables locally to avoid data races

        # Gabor transform (i.e. windowed FFT) to have instantaneous freq
        gabor = Gabor(fn_param)
        @unpack tslide = gabor
        dgt = tslide[2] - tslide[1]

        # Possible modes to use
        profiles = ["gauss", "LG"]

        # Define parameters for this specific iteration
        if steps != 1
            # Current pump waist size
            w_xp = round(wp_init + (200e-6 * (I[1] - 1)), digits=5)

            # Current initial seed waist size
            w_xs = round(w_init + (200e-6 * (I[2] - 1)), digits=5)
            w_ys = round(w_init + (200e-6 * (I[2] - 1)), digits=5)

            # Current expanded waist size
            w = round(w_exp + (200e-6 * (I[3] - 1)), digits=5)

            # Current pass number to change beam size
            stop = 2 + I[4]

            # Current pump and seed energies
            Ep0 = round(Ep_init + (5e-3 * (I[5] - 1)), digits=5)
            Ein0 = round(Ein_init + (0.6e-3 * (I[6] - 1)), digits=5)
        else
            # Current pump waist size
            w_xp = wp_init

            # Current initial seed waist size
            w_xs = w_init
            w_ys = w_init

            # Current expanded waist size
            w = w_exp

            # Current pass number to change beam size
            stop = 3

            # Current pump and seed energies
            Ep0 = Ep_init
            Ein0 = Ein_init
        end

        # Seed pulse temporal envelope profile
        Ps0 = 0.94 * Ein0 / τs
        At = ComplexEnvelope(sqrt(Ps0), t0, ϕ0, τs, GDD)

        # Initialize fluences
        Jin, Jout, Jsto = createFluenceArrays(fn_param, w_xs, w_ys, w_xp)
        println("w_xs = ", w_xs, "; w_ys = ", w_ys, "; w_xp = ", w_xp)

        # Calculate effective mode areas
        Aeff_s = calcAeff(x, y, Jin)
        Aeff_p = calcAeff(x, y, Jsto)

        Aeffs = zeros(pass)
        Aeffs[1] = Aeff_s

        # Initial fluences
        Jin0 = 2 * Ein0 / Aeff_s
        Eabs = ((Ep0) * (A * ηc * ηq)) * (2 - (A * ηc * ηq))
        Jsto0 = Eabs / Aeff_p
        Jsto .*= Jsto0
        scaleJsto = Eabs / NumericalIntegration.integrate((x,y), Jsto)
        Jsto .*= scaleJsto

        # Create intensity temporal profiles
        It = createIntensityArray(At, pass, tslide, Jin0)

        # Create arrays to store energies and B integral for visualization
        Epass, Esto, B, extraction = createPassArrays(fn_param, It, tslide, Jin, Jsto)

        println(Eabs)
        println(Esto[1])

        # Initial set of passes for seed profile
        prof = profiles[1]
        Jin, Jin0, Jout, Jsto, Epass, Esto, B, It, Aeffs, extraction = several_passes(fn_param, fft_plan, gabor, start, stop, Jin, Jin0, Jout, Jsto, At.(t), Jsat, 
                                                                                    Epass, Esto, B, It, Aeffs, extraction, prof, w, w_xs, w_ys, n2, visualize)

        if LG
        # Transform into expanded Laguerre-Gaussian beam profile
            prof = profiles[2]
            Jin = abs2.(LaguerreGauss(fn_params, P, L, a, w))
        else
        # Transform into expanded Gaussian beam profile
            prof = profiles[1]
            Jin = abs2.(Gaussian(fn_param, w, w))
        end

        # Calculate new fluence
        Aeff_s = calcAeff(x, y, Jin)
        Jin0 = 2 * Epass[stop] / Aeff_s

        # Properly scale to conserve energy according to new beam profile
        Jin, scaling = conserveEnergy(Jout, Aeff_s, tslide, It, stop, dgt, Jin, x, y)

        # Remaining passes for new seed profile
        Jin, Jin0, Jout, Jsto, Epass, Esto, B, It, Aeffs, extraction = several_passes(fn_param, fft_plan, gabor, stop+1, pass, Jin, Jin0, Jout, Jsto, At.(t), Jsat, 
                                                                                Epass, Esto, B, It, Aeffs, extraction, prof, w, w_xs, w_ys, n2, visualize)

        # Record relevant values
        lock(lk) do
        Emax = maximum(Epass)
        index = findfirst(x -> x == Emax, Epass)
        @inbounds pass_num[I] = index
        @inbounds Epass_max[I] = Emax
        @inbounds Epasses[I, :] .= Epass
        @inbounds Aeff_all[I, :] .= Aeffs
        @inbounds Aeff_max[I] = Aeffs[Int(index)]
        @inbounds B_max[I] = B[Int(index)]
        @inbounds extract_all[I, :] .= extraction
        @inbounds B_all[I, :] .= B

        println(I[1],I[2],I[3],I[4],I[5],I[6], " step done")
        println("")
        end
    end

    if steps == 1
        w_init_arr = [w_init]
        w_exp_arr = [w_exp]
        wp_init_arr = [wp_init]
        pass_exp = [2]
        Ep_arr = [Ep_init]
        Ein_arr = [Ein_init]
    else
        w_init_arr = range(w_init, round(w_init + (200e-6 * steps), digits=4), steps)
        w_exp_arr = range(w_exp, round(w_exp + (200e-6 * steps), digits=4), steps)
        wp_init_arr = range(wp_init, round(wp_init + (200e-6 * steps), digits=4), steps)
        pass_exp = range(2, 2+steps; step=1)
        Ep_arr = range(Ep_init, round(Ep_init + (5e-3 * steps), digits=4), steps)
        Ein_arr = range(Ein_init, round(Ein_init + (0.6e-3 * steps), digits=4), steps)
    end

    # Save arrays after all scans

    if isdir("data")
        cd("data")
        CSV.write("w_init_arr.csv", Tables.table(w_init_arr), writeheader=true)
        CSV.write("w_exp_arr.csv", Tables.table(w_exp_arr), writeheader=true)
        CSV.write("wp_init_arr.csv", Tables.table(wp_init_arr), writeheader=true)
        CSV.write("pass_exp.csv", Tables.table(pass_exp), writeheader=true)
        CSV.write("Ep_arr.csv", Tables.table(Ep_arr), writeheader=true)
        CSV.write("Ein_arr.csv", Tables.table(Ein_arr), writeheader=true)

        if LG
            save("LG_E_allpasses_allsteps.jld", "Epasses", Epasses)
            save("LG_extraction_allpasses_allsteps.jld", "extraction", extract_all)
            save("LG_B_allpasses_allsteps.jld", "B", B_all)
            save("LG_Aeff_allpasses_allsteps.jld", "Aeff_all", Aeff_all)
            save("LG_Emax.jld", "Emax", Epass_max)
            save("LG_maxpass.jld", "pass", pass_num)
            save("LG_Bmax.jld", "B", B_max)
            save("LG_Aeffs.jld", "Aeffs", Aeff_max)
        else
            save("Gaussian_E_allpasses_allsteps.jld", "Epasses", Epasses)
            save("Gaussian_extraction_allpasses_allsteps.jld", "extraction", extract_all)
            save("Gaussian_B_allpasses_allsteps.jld", "B", B_all)
            save("Gaussian_Aeff_allpasses_allsteps.jld", "Aeff_all", Aeff_all)
            save("Gaussian_Emax.jld", "Emax", Epass_max)
            save("Gaussian_maxpass.jld", "pass", pass_num)
            save("Gaussian_Bmax.jld", "B", B_max)
            save("Gaussian_Aeffs.jld", "Aeffs", Aeff_max)
        end
        cd("..")
    else
        mkdir("data")
        cd("data")
        CSV.write("w_init_arr.csv", Tables.table(w_init_arr), writeheader=true)
        CSV.write("w_exp_arr.csv", Tables.table(w_exp_arr), writeheader=true)
        CSV.write("wp_init_arr.csv", Tables.table(wp_init_arr), writeheader=true)
        CSV.write("pass_exp.csv", Tables.table(pass_exp), writeheader=true)
        CSV.write("Ep_arr.csv", Tables.table(Ep_arr), writeheader=true)
        CSV.write("Ein_arr.csv", Tables.table(Ein_arr), writeheader=true)

        if LG
            save("LG_E_allpasses_allsteps.jld", "Epasses", Epasses)
            save("LG_extraction_allpasses_allsteps.jld", "extraction", extract_all)
            save("LG_B_allpasses_allsteps.jld", "B", B_all)
            save("LG_Aeff_allpasses_allsteps.jld", "Aeff_all", Aeff_all)
            save("LG_Emax.jld", "Emax", Epass_max)
            save("LG_maxpass.jld", "pass", pass_num)
            save("LG_Bmax.jld", "B", B_max)
            save("LG_Aeffs.jld", "Aeffs", Aeff_max)
        else
            save("Gaussian_E_allpasses_allsteps.jld", "Epasses", Epasses)
            save("Gaussian_extraction_allpasses_allsteps.jld", "extraction", extract_all)
            save("Gaussian_B_allpasses_allsteps.jld", "B", B_all)
            save("Gaussian_Aeff_allpasses_allsteps.jld", "Aeff_all", Aeff_all)
            save("Gaussian_Emax.jld", "Emax", Epass_max)
            save("Gaussian_maxpass.jld", "pass", pass_num)
            save("Gaussian_Bmax.jld", "B", B_max)
            save("Gaussian_Aeffs.jld", "Aeffs", Aeff_max)
        end
        cd("..")
    end

end


println("ParamScan.jl compiled")
