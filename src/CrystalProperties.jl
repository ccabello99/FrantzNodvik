

# Absorption and emission cross sections of Sapphire
function getCrossSections(absorptionCSV::String, emissionCSV::String)
    # Cross section normalization
    σe0 = 1e-23 
    σa0 = 1e-24

    # Read-in data & delete bad points
    data1 = CSV.read(absorptionCSV, DataFrame)
    wl1 = data1[!, 1].*1e-9
    absp = data1[!, 2]
    data2 = CSV.read(emissionCSV, DataFrame)
    wl2 = data2[!, 1].*1e-9
    deleteat!(wl2, 19)
    deleteat!(wl2, 78)
    emiss = data2[!, 2]
    deleteat!(emiss, 19)
    deleteat!(emiss, 78)

    # Fit a spline to raw data
    spl1 = Spline1D(wl1, absp, k=3)
    spl2 = Spline1D(wl2, emiss, k=3)
    σa(λ) = σa0 .* spl1(λ)
    σe(λ) = σe0 .* spl2(λ)
    return σa, σe
end

# Sellmeier Eq coefficients
# Sapphire
@with_kw struct SapphireSellmeier
    B1::Float64 = 1.4313493
    B2::Float64 = 0.65054713
    B3::Float64 = 5.3414021
    C1::Float64 = 0.0726631
    C2::Float64 = 0.1193242
    C3::Float64 = 18.028251
end

sapphire = SapphireSellmeier()

#Ethlyne Glycol
@with_kw struct EthyleneGlycSellmeier
    B1::Float64 = 0.01778
    B2::Float64 = 1.01887
    C1::Float64 = 9.02
    C2::Float64 = 0.01028
    B3::Float64 = 0
    C3::Float64 = 0
end

#=
#BBO stuff (need to wrap this in functions TODO)
no(λ) = sqrt(2.7405 + ((0.0184) / ((λ*1e6)^2 - 0.0179)) - (0.0155 * (λ*1e6)^2))
ne(λ) = sqrt(2.3730 + ((0.0128) / ((λ*1e6)^2 - 0.0156)) - (0.0044 * (λ*1e6)^2))
dnodλ(λ) = ForwardDiff.derivative.(no, λ)
dnedλ(λ) = ForwardDiff.derivative.(ne, λ)
βo1(λ) = (c ./ (no(λ) .- ((λ) .* dnodλ(λ))))^-1 # Group Delay
βe1(λ) = (c ./ (ne(λ) .- ((λ) .* dnedλ(λ))))^-1 # Group Delay
=#

# Refractive index of Sapphire based on Sellmeier equation & tabulated coefficients (refractiveindex.info)
function SaphRefractiveIndex(sapphire::SapphireSellmeier)
    @unpack B1, B2, B3, C1, C2, C3 = sapphire
    n(λ) = sqrt(1 + ((B1 * (λ*1e6)^2) / ((λ*1e6)^2 - C1^2)) + ((B2 * (λ*1e6)^2) / ((λ*1e6)^2 - C2^2)) + ((B3 * (λ*1e6)^2) / ((λ*1e6)^2 - C3^2)))
    return n
end

# TODO if you want to test stuff for propagation in Ethlyne Glycol later
#Ethylene Glycol Refractive Index 
#data4 = CSV.read("EthlyeneGlycol_RefInd.csv",DataFrame, header=false)
#wl4 = data4[!, 1].*1e-6
#n_data = data4[!, 2]
#spl4 = Spline1D(wl4, n_data, k=1)
#global n(λ) = spl4(λ)
#dndλ = Spline1D(wl4, Dierckx.derivative(spl4, wl4), k=1)
#d2ndλ2 = Spline1D(wl4, Dierckx.derivative(dndλ, wl4), k=1)
#global β1(λ) = (c ./ (n(λ) .- ((λ) .* dndλ(λ))))^-1 # Group Delay
#global β2(λ) = (-(λ)^3/(2*π*c^2)) .* d2ndλ2(λ) #GVD

# Spectral Phases from refractive index n
function SpectralPhases(n::Function)
    c = 2.998e8
    dndλ(λ) = ForwardDiff.derivative.(n, λ)
    d2ndλ2(λ) = ForwardDiff.derivative.(dndλ, λ)
    d3ndλ3(λ) = ForwardDiff.derivative.(d2ndλ2, λ)
    β0(λ) = (2π/λ) .* n(λ)# CEP
    β1(λ) = (c ./ (n(λ) .- ((λ) .* dndλ(λ))))^-1 # Group Delay
    β2(λ) = (-(λ)^3/(2*π*c^2)) .* d2ndλ2(λ) #GVD
    β3(λ) = ((λ)^4/(4*π^2*c^3)) .* (3 .* d2ndλ2(λ) + (λ) .* d3ndλ3(λ)) #TOD
    return β0, β1, β2, β3
end

# Wavelength dependent group velocity
function groupVel(β1)
    vg(λ) = 1 / β1(λ)
    return vg
end


# Nonlinear refractive index of Sapphire
function SaphKerrRefractiveIndex(n2CSV::String)
    data3 = CSV.read(n2CSV, DataFrame)
    wl3 = data3[!, 1].*1e-9
    deleteat!(wl3, 24)
    n2_data = data3[!, 2]
    deleteat!(n2_data, 24)
    
    spl3 = Spline1D(wl3, n2_data, k=1)
    n2(λ) = 1e-20 .* spl3(λ) # (m^2/W)
    return n2
end



# Visualization of cross-sections and spectral phases

#=
nlamb = 20000
wl = 1e-6 .* range(0.300, 1.100, nlamb)
σa, σe = getCrossSections("absorption_crosssection.csv", "emission_crosssection.csv")

fig = Figure(fontsize = 42, resolution = (1920, 1080)) 
ax = Axis(fig[1,1], xlabel="Wavelength (nm)", ylabelpadding = 40,
ylabel="Cross-Section (cm²)", xticks = ([300, 500, 700, 900, 1100],[L"300", L"500", L"700", L"900", L"1100"]),
    title = L"Al_{2}O_{3}:Ti^{3+} \text{ Absorption and Emission Cross-Sections}", yticks=([0, 10, 20, 30, 40], [L"0", L"1 \times 10^{-19}", L"2.0 \times 10^{-19}", L"3.0 \times 10^{-19}", L"4.0 \times 10^{-19}"]))

lines!(wl .* 1e9, 1e20 .* σe.(wl) .* 1e4, label="Emission", linewidth=4)
lines!(wl .* 1e9, 1e20 .* σa.(wl) .* 1e4, label="Absorption", linewidth=4)
axislegend(ax, position = :lc)
fig
#Makie.save("Ti-Sapphire_cross_sections.png", fig)


ax = Axis(fig[1,1], ylabel="Normalized values (arb. u.)", 
            title = "Spectral Responses of Sapphire", xlabel="Wavelength (nm)",
            #xticks = ([300, 500, 700, 900, 1100],[L"300", L"500", L"700", L"900", L"1100"]),
            yticks = ([-1.0, -0.5, 0, 0.5, 1.0], [L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"]))

lines!(wl, N.(wl) ./ maximum(abs.(N.(wl))), label="Refractive Index", linewidth=4) 
lines!(wl, n2.(wl) ./ maximum(abs.(n2.(wl))), label="Kerr Refractive Index", linewidth=4) 
lines!(wl, β1.(wl) ./ maximum(abs.(β1.(wl))), label="β1", linewidth=4)
lines!(wl, β2.(wl) ./ maximum(abs.(β2.(wl))), label="β2", linewidth=4)
lines!(wl, β3.(wl) ./ maximum(abs.(β3.(wl))), label="β3", linewidth=4)
axislegend(ax, position = :lc)
#Makie.save("Sapphire_spectral_responses.png", fig)

fig
=#
println("CrystalProperties.jl compiled")
