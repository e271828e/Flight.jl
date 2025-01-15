module GSRUKF

#efficient, non-allocating generalized square-root unscented Kalman filter
#implementation

using UnPack
using LinearAlgebra
using LinearAlgebra: libblastrampoline, BlasFloat, BlasInt
using LinearAlgebra.BLAS: @blasfunc
using LinearAlgebra.LAPACK: chklapackerror

export GSRUTWorkspace, GSRUKFPredictor, GSRUKFUpdater
export transform!



############################## GSRUTWorkspace ##################################
################################################################################

struct GSRUTWorkspace
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

function GSRUTWorkspace(;
    N_u::Integer, N_z::Integer, N_w::Integer,
    N_δz::Integer = N_z, N_δu::Integer = N_u,
    α::AbstractFloat = 1e-2, β::Real = 2, κ::Real = 0,
    ε_mean::AbstractFloat = 1e-12, max_iter::Integer = 5)

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


    GSRUTWorkspace( S_δa, δA, Z, δZ, QRt, μ_m, μ_c, u_b, z_b, δz_b, work, τ, γ, ε_mean, max_iter)

end

function transform!(
    ws::GSRUTWorkspace;
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
        @. P_δuδz += μ * δu * δz'
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


end #module