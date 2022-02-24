module TestPropellers

using Test
using Plots
using UnPack
using StructArrays
using LaTeXStrings
using LinearAlgebra
using Interpolations
using Flight.Propellers
using Flight.Propellers: DefaultAirfoil, cL, cD, cL_α
using Flight.Propellers: PropBlade, PropDataset

function test_propellers()
    @testset verbose = true "Propellers" begin
        @testset verbose = true "DefaultAirfoil" begin test_default_airfoil() end
    end
end

function test_default_airfoil(show = false)
    airfoil = DefaultAirfoil()

    α = range(-π/6, π/3, length = 100)
    M = range(0, 1.5, length = 6)

    cL_data = Array{Float64}(undef, length(α), length(M))
    cL_α_data = similar(cL_data); cD_data = similar(cL_data)

    for (j, M) in enumerate(M)
        for (i, α) in enumerate(α)
            cL_data[i, j] = cL(airfoil, α, M)
            cD_data[i, j] = cD(airfoil, α, M)
            cL_α_data[i, j] = cL_α(airfoil, α, M)
        end
    end

    p = Vector{Plots.Plot}()
    push!(p, plot(α, cL_data, label = "M = ".*string.(M'), title = L"c_L", xlabel = L"\alpha \ (rad)"))
    push!(p, plot(α, cL_α_data, label = "M = ".*string.(M'), title = L"c_{L,\alpha}", xlabel = L"\alpha \ (rad)"))
    push!(p, plot(α, cD_data, label = "M = ".*string.(M'), title = L"cD", xlabel = L"\alpha \ (rad)"))

    if show
        display.(p)
    end

end

function test_prop_coefficients()
    n_blades = 2
    blade = PropBlade()
    # blade = PropBlade(c̃ = Propellers.ConstantFunction(0.0573))

    PropCoefficients(blade, n_blades; J = 0.5, M_t = 1.2, Δβ = 0)
end

function test_fixed_pitch_dataset(show = false)

    #generate dataset
    n_blades = 2
    blade = PropBlade()
    pitch = Propellers.FixedPitch()
    dataset = PropDataset(pitch, blade, n_blades)

    #define evaluation grid
    J_bounds, M_tip_bounds = bounds(dataset)
    J = range(J_bounds[1], J_bounds[2], length = 100)
    M_tip = range(M_tip_bounds[1], stop = M_tip_bounds[2], step = 0.4)
    iter = Iterators.product(J, M_tip)

    coeffs = Array{PropCoefficients{Float64}}(undef, size(iter))

    for (i, (J, M_tip)) in enumerate(iter)
        coeffs[i] = PropCoefficients(dataset, J, M_tip)
    end

    coeffs_sa = coeffs |> StructArray |> StructArrays.components

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = coeffs_sa

    p = Vector{Plots.Plot}()
    push!(p, plot(J, C_Fx, label = "Mtip = " .* string.(M_tip'), title = L"C_{Fx}", xlabel = L"J"))
    push!(p, plot(J, C_Mx, label = "Mtip = " .* string.(M_tip'), title = L"C_{Mx}", xlabel = L"J"))
    push!(p, plot(J, C_Fz_α, label = "Mtip = " .* string.(M_tip'), title = L"C_{Fz, \alpha}", xlabel = L"J"))
    push!(p, plot(J, C_Mz_α, label = "Mtip = " .* string.(M_tip'), title = L"C_{Mz, \alpha}", xlabel = L"J"))
    push!(p, plot(J, C_P, label = "Mtip = " .* string.(M_tip'), title = L"C_{P}", xlabel = L"J"))
    push!(p, plot(J, η_p, label = "Mtip = " .* string.(M_tip'), title = L"\eta_p", xlabel = L"J"))

    if show
        display.(p)
    end

    error("Verify no allocations when evaluating a dataset")
    # alloc =  @allocated PropCoefficients(dataset,0,0)
    # @show alloc

    # return dataset

end

function test_variable_pitch_dataset(show = false)
    #generate dataset, fix M_tip and evaluate at different blade pitch offsets
end

function test_sytem_functions()
    # alloc =  @allocated PropCoefficients(dataset,0,0)
    # @show alloc
    #require no allocations

end


end