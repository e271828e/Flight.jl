Base.@kwdef struct CoefficientData{T}
    C0::Float64
    C1::T #derivatives
end


#requirements:
#0) all derivative terms stored in labelled fields
#1) do not require defining zero derivative terms for every coefficient
#2) do not require manually computing all build-ups
# const C_X_linear_terms = (:α, :α², :α³)
const C_X_data = CoefficientData(C0 = 0.34597, C1 = (α = 0.0, α² = 0.0, α³ = 0.0))
const C_Y_data = CoefficientData(C0 = 0.54597, C1 = (β = 1.0, β² = 1.0, β³ = 1.0))
# const C_X_data = CoefficientData(C0 = -0.03554,
#     C1 = (α = 0.002920, α² = 5.459, α³ = -5.162, q_nd = -0.6748, δr = 0.03412, δf = -0.09447, δf_α = 1.106))

# const C_Y_data = CoefficientData(C0 = -0.002226,
#     C1 = (β = -0.7678, p_nd = -0.1240, r_nd = 0.3666, δa = -0.02956, δr = 0.1158, δr_α = 0.5238, β_dot_nd = -0.1600))


function aero_test()
    α = 0.1
    α² = α^2
    α³ = α^3
    β = 1.0
    β² = α^2
    β³ = α^3
    p̃ = 1
    #require a NamedTuple input to avoid potential field mix ups between different coefficients
    C_X_inputs = (α = α, α² = α², α³ = α³)
    C_Y_inputs = (β = β, β² = β², β³ = β³)
    C_X = evaluate_coefficient(C_X_data, C_X_inputs)
    # @show C_X
    C_Y = evaluate_coefficient(C_Y_data, C_Y_inputs)
    # @show C_Y
end

@inline @generated function evaluate_coefficient(data::CoefficientData{T}, input::T) where {T}
    # Core.println("Generated function called for type $T")
    ex = Expr(:block)
    push!(ex.args, :(coef = data.C0))
    for label in fieldnames(T)
        push!(ex.args,
            :(coef += getproperty(data.C1, $(QuoteNode(label))) * getproperty(input, $(QuoteNode(label)))))
    end
    push!(ex.args, :(return coef))
    return ex
end