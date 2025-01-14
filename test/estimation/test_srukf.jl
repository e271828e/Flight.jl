module TestSRUKF

using Test

using ComponentArrays, StaticArrays
using BenchmarkTools
using UnPack

using LinearAlgebra
using LinearAlgebra: libblastrampoline, BlasFloat, BlasInt
using LinearAlgebra.BLAS: @blasfunc
using LinearAlgebra.LAPACK: chklapackerror

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
        iter >= max_iter && (@error "Generalized mean error failed to converge "*"
                            ($norm_δz̄ > $ε_mean) after $iter iterations"; break)
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


############################## SRUTWorkspace ###################################
################################################################################

struct SRUTWorkspace
    S_δa::LowerTriangular{Float64, Matrix{Float64}} #augmented error square-root covariance
    δA::Matrix{Float64} #augmented error sigma points
    Z::Matrix{Float64} #transformed sigma points
    δZ::Matrix{Float64} #transformed error sigma points
    QRt::Matrix{Float64} #QR factorization target matrix
    μ_m::Vector{Float64} #sample mean weights
    μ_c::Vector{Float64} #sample covariance weights
    u_b::Vector{Float64} #buffer for initial variable samples
    z_b::Vector{Float64} #buffer for transformed variable samples
    δz_b::Vector{Float64} #buffer for transformed error variable samples
    work::Vector{Float64} #used by as temporary storage
    τ::Vector{Float64} #used by dgeqrf to store scalar reflector factors for Q
    γ::Float64 #sigma point scaling parameter
    ε_mean::Float64 #tolerance for generalized mean computation
    max_iter::Int #maximum iterations for generalized mean computation
end

function SRUTWorkspace(;
    N_u::Integer, N_δu::Integer, N_z::Integer, N_δz::Integer, N_w::Integer,
    α::AbstractFloat = 1e-3, β::Real = 2, κ::Real = 0,
    ε_mean::AbstractFloat = 1e-11, max_iter::Integer = 5)

    N_δa = N_δu + N_w
    @assert 2N_δa >= N_δz #tall QRt (necessary condition for full-rank S_δz)

    S_δa = LowerTriangular(zeros(N_δa, N_δa))
    δA = zeros(N_δa, 2N_δa + 1)
    Z = zeros(N_z, 2N_δa + 1)
    δZ = zeros(N_δz, 2N_δa + 1)
    QRt = zeros(2N_δa, N_δz)
    u_b = zeros(N_u)
    z_b  = zeros(N_z)
    δz_b  = zeros(N_δz)
    τ = zeros(N_δz) #number of columns of QRt
    work = zeros(lwork_query(QRt, τ))

    α² = α^2
    γ² = α² * (N_δa + κ)
    λ = γ² - N_δa
    γ = √γ²

    #mean and covariance weights
    μ_m0 = λ / γ²
    μ_c0 = μ_m0 + (1 - α² + β)
    μ_m1 = 1/(2γ²)
    μ_c1 = μ_m1

    μ_m = [i == 0 ? μ_m0 : μ_m1 for i in 0:2N_δa] #2N_δa+1 sigma points
    μ_c = [i == 0 ? μ_c0 : μ_c1 for i in 0:2N_δa]


    SRUTWorkspace( S_δa, δA, Z, δZ, QRt, μ_m, μ_c, u_b, z_b, δz_b, work, τ, γ, ε_mean, max_iter)

end

function transform!(
    ws::SRUTWorkspace;
    z̄::AbstractVector{<:Real},
    S_δz::LowerTriangular{<:Real},
    P_δuδz::AbstractMatrix{<:Real},
    ū::AbstractVector{<:Real},
    S_δu::LowerTriangular{<:Real},
    S_w::LowerTriangular{<:Real},
    f!::Function = (z, u, w) -> z .= u,
    g_u_plus!::Function = (u1, u0, δu) -> u1 .= u0 .+ δu,
    g_z_plus!::Function = (z1, z0, δz) -> z1 .= z0 .+ δz,
    g_z_minus!::Function = (δz, z1, z0) -> δz .= z1 .- z0,
    )

    @unpack S_δa, δA, Z, δZ, QRt, μ_m, μ_c, u_b, z_b, δz_b, work, τ, γ, ε_mean, max_iter = ws

    N_δu = size(S_δu)[1]
    N_δz = size(S_δz)[1]
    N_w = size(S_w)[1]
    N_δa = N_δu + N_w

    #initialize augmented error square-root covariance
    S_δa[1:N_δu, 1:N_δu] .= S_δu
    S_δa[N_δu+1:end, N_δu+1:end] .= S_w

    #generate augmented error sigma points
    δA[:, 1] .= 0 #0th column is zero
    δA1 = @view δA[:, 2:N_δa+1] #positive sigma points
    δA2 = @view δA[:, N_δa+2:end] #negative sigma points
    for (s, δa1, δa2) in zip(eachcol(S_δa), eachcol(δA1), eachcol(δA2))
        δa1 .= γ .* s
        δa2 .= .-δa1
    end

    #compute total transformed sigma point distribution
    δU = @view δA[1:N_δu, :]
    W = @view δA[N_δu+1:end, :]
    for (δu, w, z) in zip(eachcol(δU) , eachcol(W), eachcol(Z))
        g_u_plus!(u_b, ū, δu)
        f!(z, u_b, w)
    end

    #compute transformed mean and error sigma point distribution
    z̄ .= @view Z[:, 1] #take the first sigma point as initial guess
    for iter in Iterators.countfrom(0)
        for (z, δz) in zip(eachcol(Z), eachcol(δZ))
            g_z_minus!(δz, z, z̄)
        end
        weighted_mean!(δz_b, δZ, μ_m)
        norm(δz_b) <= ε_mean && break #if threshold is achieved, exit
        g_z_plus!(z_b, z̄, δz_b) #otherwise, compute corrected guess
        z̄ .= z_b #assign it to z̄ (avoid g_z_plus!(z̄, z̄, δz_b), could be dangerous)
        iter >= max_iter && (@error "Generalized mean error failed to converge "*"
                            ($(norm(δz_b)) > $ε_mean) after $iter iterations")
        iter += 1
    end

    δz0 = @view δZ[:, 1]
    δZi = @view δZ[:, 2:end]
    μ_c0 = μ_c[1]
    μ_ci = @view μ_c[2:end]

    #construct QR target matrix, compute its R factor and transpose it into S_δz
    QRt .= sqrt.(μ_ci) .* δZi'
    dgeqrf!(QRt, τ, work)
    R = UpperTriangular(@view QRt[1:N_δz, :])
    transpose!(S_δz, R)

    #force the main diagonal to be positive (not required)
    (S_δz[1,1] < 0) && (S_δz .*= -1)

    #perform the remaining low-rank update or downdate
    δz_b .= √abs(μ_c0) .* δz0
    C_δz = Cholesky(S_δz)
    μ_c0 >= 0 ? lowrankupdate!(C_δz, δz_b) : lowrankdowndate!(C_δz, δz_b) #mutates S indirectly

    #compute cross-covariance
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


function lwork_query(A::Matrix{Float64}, τ::Vector{Float64})

    M, N = size(A)
    @assert length(τ) == min(M, N)
    work  = Vector{Float64}(undef, 1)
    info  = Ref{BlasInt}()
    lwork = BlasInt(-1)
    ccall((@blasfunc(dgeqrf_), libblastrampoline), Cvoid,
            (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
            BlasInt(M), BlasInt(N), A, max(1,stride(A,2)), τ, work, lwork, info)

    chklapackerror(info[])

    return max(BlasInt(1), BlasInt(real(work[1]))) #lwork

end

function dgeqrf!(A::Matrix{Float64}, τ::Vector{Float64}, work::Vector{Float64})

    lwork = length(work)
    M, N = size(A)

    info  = Ref{BlasInt}()
    ccall((@blasfunc(dgeqrf_), libblastrampoline), Cvoid,
            (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
            BlasInt(M), BlasInt(N), A, max(1,stride(A,2)), τ, work, lwork, info)

    chklapackerror(info[])

end




################################# Tests ########################################

function test_dgeqrf()
    A = rand(6, 3)
    τ = zeros(3)
    M, N = size(A)
    @assert M >= N

    A_test = copy(A)
    R_ref = qr!(A_test).R

    A_test = copy(A)
    work = zeros(lwork_query(A, τ))
    dgeqrf!(A_test, τ, work)
    R_dgeqrf = UpperTriangular(@view A_test[1:N, :])

    @test all(isapprox(R_ref, R_dgeqrf))

    # b_qr = @benchmarkable qr!($A_test) setup=(A_test = copy($A))
    # @show run(b_qr)

    b_dgeqrf = @benchmarkable dgeqrf!($A_test, $τ, $work) setup=(A_test = copy($A))
    @show run(b_dgeqrf)

end



function test00()

    #the total state comprises a unit quaternion q and an Euclidean vector v
    x_template = ComponentVector(q = zeros(4), v = zeros(3))
    x_axes = getaxes(x_template)[1]

    #the error state comprises a rotation vector ρ and an Euclidean vector δv
    δx_template = ComponentVector(ρ = zeros(3), δv = zeros(3))
    δx_axes = getaxes(δx_template)[1]

    #the noise vector has the same structure as the error state
    w_template = ComponentVector(ρ = zeros(3), δv = zeros(3))
    w_axes = getaxes(w_template)[1]

    #generalized state addition
    g_x_plus! = let x_axes = x_axes, δx_axes = δx_axes
        function (x1::AbstractVector{<:Real}, x0::AbstractVector{<:Real}, δx::AbstractVector{<:Real})
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
    end

    #generalized state subtraction
    g_x_minus! = let x_axes = x_axes, δx_axes = δx_axes
        function (δx::AbstractVector{<:Real}, x1::AbstractVector{<:Real}, x0::AbstractVector{<:Real})
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
    end

    #####################

    #test consistency between generalized addition and subtraction with
    #arbitrary x and δx values
    x_a = copy(x_template)
    x_b = copy(x_template)
    x_a.q .= RQuat([2, -1, 3, 1])[:] #will be automatically normalized
    x_a.v .= [3, 2, 5]

    δx = copy(δx_template)
    δx.ρ .= [0.3, 0.4, -0.5]
    δx.δv .= [1, 3, 2.0]
    δx_copy = copy(δx)

    g_x_plus!(x_b, x_a, δx)
    g_x_minus!(δx, x_b, x_a)

    @test all(isapprox.(δx, δx_copy))

    ####################

    #set the the initial state mean to the identity rotation and a null v
    x̄0 = copy(x_template)
    x̄0.q .= RQuat()[:]
    x̄0.v .= [0, 0, 0]

    #set initial square-root covariance
    # A = randn(6, 6)
    # P_δx0 = A * A'
    # S_δx0 = cholesky(P_δx0).L
    σ_δx0 = ComponentVector(zeros(6), δx_axes)
    σ_δx0.ρ .= 1e-2
    σ_δx0.δv .= 2e-2
    S_δx0 = LowerTriangular(ComponentMatrix(diagm(σ_δx0), δx_axes, δx_axes))

    #preallocate arrays for the SRUT output
    x̄1 = copy(x̄0)
    S_δx1 = copy(S_δx0)
    P_δx0δx1 = ComponentMatrix(zeros(size(S_δx0)), δx_axes, δx_axes)

    #define the test transformation function, which:
    # 1) composes the initial quaternion with the rotation given by a specified
    # rotation vector ρ01 plus noise w.ρ
    # 2) adds a specified vector Δv01 plus noise w.δv to the initial vector v0

    ρ01 = @SVector[0.0, 0.0, 0.4]
    v01 = @SVector[3.0, 2.0, -0.4]

    f_x! = let x_axes = x_axes, w_axes = w_axes, ρ01 = ρ01, v01 = v01
        function (x1::AbstractVector{<:Real}, x0::AbstractVector{<:Real}, w::AbstractVector{<:Real})
            x0_ca = ComponentVector(x0, x_axes)
            x1_ca = ComponentVector(x1, x_axes)
            w_ca = ComponentVector(w, w_axes)

            q0 = RQuat(x0_ca.q)
            q01w = RQuat(RVec(SVector{3}(ρ01) + SVector{3}(w_ca.ρ)))
            q1 = q0 ∘ q01w

            v0 = SVector{3}(x0_ca.v)
            v01w = SVector{3}(v01) + SVector{3}(w_ca.δv)
            v1 = v0 + v01w

            x1_ca.q = q1[:]
            x1_ca.v = v1
        end
    end

    #perform a direct, noiseless transformation of the initial mean
    x̄1_direct = copy(x̄1)
    f_x!(x̄1_direct, x̄0, zeros(6))

    #now apply the SRUT to the initial random variable with zero noise covariance
    # S_w0 = Array{Float64}(undef, 0, 0) |> LowerTriangular
    S_w0 = zeros(6,6) |> LowerTriangular
    transform(x̄1, S_δx1, P_δx0δx1, x̄0, S_δx0, S_w0;
                f! = f_x!,
                g_u_plus! = g_x_plus!,
                g_z_plus! = g_x_plus!,
                g_z_minus! = g_x_minus!)

    #the mean resulting from the SRUT should be equal to that from the direct
    #noiseless mean transformation
    @test all(isapprox(x̄1_direct, x̄1))

    #also, in the absence of noise, the error state square-root covariance
    #should be preserved by the SRUT (minus a possible sign inversion)
    @test all(isapprox(abs.(S_δx0), abs.(S_δx1)))

    #set noise square-root covariance
    σ_w0 = ComponentVector(zeros(6), w_axes)
    σ_w0.ρ .= 1e-2
    σ_w0.δv .= 3e-2
    S_w0 = LowerTriangular(ComponentMatrix(diagm(σ_w0), w_axes, w_axes))

    #we expect the mean to be the same as before, but covariance should have increased


end

function test01()

    #repeat the tests above for the non-allocating implementation

    srut = SRUTWorkspace(N_u = length(x̄0), N_δu = length(δx), N_z = length(x̄0),
                       N_δz = length(δx), N_w = size(S_w0)[1])


    transform!(srut;
        z̄ = x̄1, S_δz = S_δx1, P_δuδz = P_δx0δx1,
        ū = x̄0, S_δu = S_δx0, S_w = S_w0,
        f! = f_x!,
        g_u_plus! = g_x_plus!,
        g_z_plus! = g_x_plus!,
        g_z_minus! = g_x_minus!)

    return x̄1, S_δx1, P_δx0δx1


    @btime transform!($srut;
        z̄ = $x̄1, S_δz = $S_δx1, P_δuδz = $P_δx0δx1,
        ū = $x̄0, S_δu = $S_δx0, S_w = $S_w0,
        f! = $f_x!,
        g_u_plus! = $g_x_plus!,
        g_z_plus! = $g_x_plus!,
        g_z_minus! = $g_x_minus!)




end

end