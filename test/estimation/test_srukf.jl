module TestSRUKF
using ComponentArrays
using StaticArrays
using LinearAlgebra
using LazyArrays
using Test

using Flight.FlightPhysics

function transform(
    z̄::AbstractVector{<:Real},
    S_δz::LowerTriangular{<:Real},
    P_δuδz::AbstractMatrix{<:Real},
    ū::AbstractVector{<:Real},
    S_δu::LowerTriangular{<:Real},
    S_w::LowerTriangular{<:Real};
    f!::Function = (z, u, w) -> z .= u,
    g_u_plus!::Function = (u1, u0, δu) -> u1 .= u0 .+ δu,
    g_z_plus!::Function = (z1, z0, δz) -> z1 .= z0 .+ δz,
    g_z_minus!::Function = (δz, z1, z0) -> δz .= z1 .- z0,
    )

    N_u = length(ū)
    N_δu = size(S_δu)[1]
    N_w = size(S_w)[1]
    N_δa = N_δu + N_w

    N_z = length(z̄)
    N_δz = size(S_δz)[1]

    #the SRUT operates on the augmented error variable δa, so the
    #input dimension for the weights must be N_δa
    α = 1e-2
    β = 2
    κ = 0
    α² = α^2
    γ² = α² * (N_δa + κ)
    γ = √γ²
    λ = γ² - N_δa

    #mean and covariance weights
    μ_m0 = λ / γ²
    μ_c0 = μ_m0 + (1 - α² + β)
    μ_m1 = 1/(2γ²)
    μ_c1 = μ_m1

    μ_m = [i == 0 ? μ_m0 : μ_m1 for i in 0:2N_δa] #2N_δa+1 sigma points
    μ_c = [i == 0 ? μ_c0 : μ_c1 for i in 0:2N_δa]

    S_δa = zeros(N_δa, N_δa) #BUFFER
    S_δa = LowerTriangular(S_δa)
    S_δa[1:N_δu, 1:N_δu] .= S_δu
    S_δa[N_δu+1:end, N_δu+1:end] .= S_w

    #generate augmented error sigma points
    δA = zeros(N_δa, 2N_δa + 1) #BUFFER
    δA[:, 1] .= 0 #0th column is zero
    δA1 = @view δA[:, 2:N_δa+1] #positive sigma points
    δA2 = @view δA[:, N_δa+2:end] #negative sigma points
    for (s, δa1, δa2) in zip(eachcol(S_δa), eachcol(δA1), eachcol(δA2))
        # _δa= γ .* s
        δa1 .= γ .* s
        δa2 .= .-δa1
    end

    Z = zeros(N_z, 2N_δa + 1) #BUFFER
    _u = zeros(N_u) #BUFFER

    #compute total transformed sigma point distribution Z
    δU = @view δA[1:N_δu, :]
    W = @view δA[N_δu+1:end, :]
    for (δu, w, z) in zip(eachcol(δU) , eachcol(W), eachcol(Z))
        g_u_plus!(_u, ū, δu)
        f!(z, _u, w)
    end

    #compute transformed mean and error sigma point distribution
    δZ = zeros(N_δz, 2N_δa + 1) #BUFFER
    _z  = zeros(N_z) #BUFFER
    _δz  = zeros(N_δz) #BUFFER


    ε_mean = 1e-12
    max_iter = 5
    z̄ .= Z[:, 1] #take the first sigma point as initial guess
    for iter in Iterators.countfrom(0)
        for (z, δz) in zip(eachcol(Z), eachcol(δZ))
            g_z_minus!(δz, z, z̄)
        end
        weighted_mean!(_δz, δZ, μ_m)
        norm_δz̄ = norm(_δz)
        norm_δz̄ <= ε_mean && break
        iter >= max_iter && (@error "Generalized mean error failed to converge ($norm_δz̄ > $ε_mean) after $iter iterations"; break)
        g_z_plus!(_z, z̄, _δz) #compute corrected guess
        z̄ .= _z #assign to z̄ (g_z_plus!(z̄, z̄, _δz̄) directly could be dangerous)
        iter += 1
    end

    δz0 = @view δZ[:, 1]
    δZi = @view δZ[:, 2:end]

    μ_c0 = μ_c[1]
    μ_ci = @view μ_c[2:end]

    Ap_t = zeros(2N_δa, N_δz) #BUFFER
    Ap_t .= sqrt.(μ_ci) .* δZi'

    R = qr(Ap_t).R
    transpose!(S_δz, R)
    S_δz[1,1] < 0 && S_δz .*= -1 #make the main diagonal positive (not required)

    _δz .= √abs(μ_c0) .* δz0
    C_δz = Cholesky(S_δz)
    μ_c0 >= 0 ? lowrankupdate!(C_δz, _δz) : lowrankdowndate!(C_δz, _δz) #mutates S indirectly

    P_δuδz .= 0
    for (δu, δz, μ) in zip(eachcol(δU), eachcol(δZ), μ_c)
        P_δuδz .+= μ .* δu .* δz'
    end

end

function weighted_mean!(z̄::AbstractVector{<:Real}, Z::AbstractMatrix{<:Real},
                        μ) #μ will typically be a vector or an iterator
    N1, N2 = size(Z)
    @assert length(z̄) == N1
    @assert length(μ) == N2
    z̄ .= 0
    for (z, μ) in zip(eachcol(Z), μ)
        z̄ .+= μ .* z
    end
end

function weighted_mean!(z̄::AbstractVector{<:Real}, Z::AbstractMatrix{<:Real})
    N2 = size(Z)[2]
    μ = Iterators.repeated(1/N2, N2)
    weighted_mean!(z̄, Z, μ)
end



################################# Tests ########################################

const δx_template = ComponentVector(ρ = zeros(3), δv = zeros(3));
const x_template = ComponentVector(q = zeros(4), v = zeros(3));
const w_template = ComponentVector(ρ = zeros(3), δv = zeros(3));

const δx_axes = getaxes(δx_template)[1]
const x_axes = getaxes(x_template)[1]
const w_axes = getaxes(w_template)[1]


function g_x_plus!(x1::AbstractVector{<:Real}, x0::AbstractVector{<:Real}, δx::AbstractVector{<:Real})
    x0_ca = ComponentVector(x0, x_axes)
    x1_ca = ComponentVector(x1, x_axes)
    δx_ca = ComponentVector(δx, δx_axes)

    q0 = RQuat(x0_ca.q)
    v0 = SVector{3,Float64}(x0_ca.v)

    ρ = RVec(δx_ca.ρ)
    δv = SVector{3,Float64}(δx_ca.δv)

    q1 = q0 ∘ ρ
    v1 = v0 + δv

    #mutates x1
    x1_ca.q = q1[:]
    x1_ca.v = v1
end

function g_x_minus!(δx::AbstractVector{<:Real}, x1::AbstractVector{<:Real}, x0::AbstractVector{<:Real})
    x0_ca = ComponentVector(x0, x_axes)
    x1_ca = ComponentVector(x1, x_axes)
    δx_ca = ComponentVector(δx, δx_axes)

    q0 = RQuat(x0_ca.q)
    v0 = SVector{3,Float64}(x0_ca.v)

    q1 = RQuat(x1_ca.q)
    v1 = SVector{3,Float64}(x1_ca.v)

    #here, q1 = q0 ∘ q01, so q01 = q0' ∘ q1
    ρ = RVec(q0' ∘ q1)
    δv = v1 - v0

    #mutates δx
    δx_ca.ρ = ρ[:]
    δx_ca.δv = δv
end

function f_x!(x1::AbstractVector{<:Real}, x0::AbstractVector{<:Real}, w::AbstractVector{<:Real})
    #extract q0, construct q01 from w.ρ, obtain q1
    #extract v0, compute v1 = v0 + w.δv
    #it's just noise for now
    x0_ca = ComponentVector(x0, x_axes)
    x1_ca = ComponentVector(x1, x_axes)
    w_ca = ComponentVector(w, w_axes)
    println(w_ca)
    x1 .= x0
end

function test00()

    println("Next: Update f_x! with noise")

    # A = randn(6, 6)
    # P_δx0 = A * A'
    # S_δx0 = cholesky(P_δx0).L

    σ_δx0 = ComponentVector(zeros(6), δx_axes)
    σ_δx0.ρ .= 1e-2
    σ_δx0.δv .= 2e-2
    S_δx0 = LowerTriangular(ComponentMatrix(diagm(σ_δx0), δx_axes, δx_axes))

    # S_w0 = Array{Float64}(undef, 0, 0) |> LowerTriangular

    σ_w0 = ComponentVector(zeros(6), w_axes)
    σ_w0.ρ .= 1e-2
    σ_w0.δv .= 3e-2
    S_w0 = LowerTriangular(ComponentMatrix(diagm(σ_w0), w_axes, w_axes))

    S_δx1 = copy(S_δx0)
    P_δx0δx1 = ComponentMatrix(zeros(size(S_δx0)), δx_axes, δx_axes)

    x̄0 = copy(x_template)
    x̄1 = copy(x_template)
    δx = copy(δx_template)

    x̄0.q .= RQuat([1,2,-2.5,3])[:]
    x̄0.v .= [4, 5, -3]

    δx.ρ = randn(3)
    δx.δv = randn(3)
    δx_copy = copy(δx)

    g_x_plus!(x̄1, x̄0, δx)
    g_x_minus!(δx, x̄1, x̄0)

    #first test the generalized addition and subtraction
    @test all(isapprox.(δx, δx_copy))


    transform(x̄1, S_δx1, P_δx0δx1, x̄0, S_δx0, S_w0;
                f! = f_x!,
                g_u_plus! = g_x_plus!,
                g_z_plus! = g_x_plus!,
                g_z_minus! = g_x_minus!)

    # @show x̄1, S_δx1
    return S_δx1

    #from here on, we are inside the SRUT

    #defaults:
    # g_u_plus! = (u1, u0, δu) -> (u1 .= u0 .+ δu)
    # g_u_minus! = (δx, u1, u0) -> (δx .= u1 .- u0)


    #do i really need to make arrays SizedArrays?. Can't I simply have a
    #non-parametric SRUTWorkspace with generic Vectors and Matrices, and the
    #only place where sizes appear is upon construction. They are not type
    #parameters. However, size type parameters are helpful for inspection. But
    #they make things more difficult to implement, and they are more prone to
    #allocations. Also, if any of the dimensions are incorrect, an error will
    #appear sooner or later when multiplying matrices

    #SizedArrays NOT WORTH IT. They aren't faster for medium-sized matrices, and
    #they are a source of headaches. Now, MMatrix could be worth it, but using
    #it essentially places an upper bound on the filter state size before things
    #go bad. Also, lowrankupdate! and similar mutating methods won't work on
    #SMatrix, and MMatrix for some reason is not as fast as SMatrix

    #conclusion: Non parametrically typed plain old Vectors and Matrices. We can
    #query a workspace for its sizes any time we want.

    #more changes. qr_target AKA A't will no longer be allocated in
    #DGERQFWorkspace, but directly in SRUTWorkspace. in fact, DFERQFWorkspace
    #can be removed, and tau and work moved to SRUTWorkspace

    #then in code we will have
    #Ap_tr = @~ w_i * dZi'
    #dgeqrf!(Ap_tr, tau, work)
    # R = ... see QRFactorization

    #need to see wh
    #if At is 3x5, then R is 3x5

    #we know that our qr factorization target is 2NdA x NdZ. When we initialize
    #the SRUT, we instantiate this as a storage buffer. And we pass it to the
    #get_lwork function

    #the SRUT should NOT store z̄ or S_δz. it should receive them as mutable
    #inputs.
    #it should provide internal storage for any buffer required in the
    #computations

end

end