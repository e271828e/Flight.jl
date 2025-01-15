module TestGSRUKF

using Test

using LinearAlgebra, ComponentArrays, StaticArrays, BenchmarkTools

using Flight.FlightPhysics.Attitude
using Flight.FlightNavigation.GSRUKF

function test_gsrukf()

    @testset verbose = true "GSRUT" begin test_gsrut() end
    @testset verbose = true "Stages" begin test_stages() end

end

function test_gsrut()

    #non-allocating QR factorization via direct LAPACK function call
    @testset verbose = true "DGEQRF" begin

        A = rand(6, 3)
        τ = zeros(3)
        M, N = size(A)
        @assert M >= N

        A_test = copy(A)
        R_ref = qr!(A_test).R

        A_test = copy(A)
        work = zeros(GSRUKF.lwork_query(A, τ))
        GSRUKF.dgeqrf!(A_test, τ, work)
        R_dgeqrf = UpperTriangular(@view A_test[1:N, :])

        @test all(isapprox(R_ref, R_dgeqrf))

        b_dgeqrf = @benchmarkable GSRUKF.dgeqrf!($A_test, $τ, $work) setup=(A_test = copy($A))
        @test memory(run(b_dgeqrf)) == 0

    end

    #basic tests with a linear transformation on a vector random variable
    @testset verbose = true "Basic" begin

        N = 3
        gsrut = GSRUTWorkspace(; N_u = N, N_z = N, N_w = N)

        f! = (z, u, w) -> @. z = 2*u + 1 + w

        ū = ones(N)
        S_δu = 1.0 * LowerTriangular(Matrix(I, N, N))
        S_w = 1.0 * LowerTriangular(Matrix(I, N, N))

        z̄ = copy(ū)
        S_δz = copy(S_δu)
        P_δuδz = zeros(N, N)

        P_δu = S_δu * S_δu'
        P_w = S_w * S_w'
        transform!(gsrut; z̄, S_δz, P_δuδz, ū, S_δu, S_w, f!)

        #we know this linear transformation must yield the following
        P_δz = S_δz * S_δz'
        @test z̄ ≈ 2ū .+ 1
        @test P_δz ≈ 4*P_δu .+ P_w
        @test P_δuδz ≈ 2I

        #check for allocations
        @test (@ballocated transform!($gsrut;
            z̄ = $z̄, S_δz = $S_δz, P_δuδz = $P_δuδz,
            ū = $ū, S_δu = $S_δu, S_w = $S_w, f! = $f!)) == 0

        #ensure noiseless transforms are correctly handled
        S_w = Array{Float64}(undef, 0, 0) |> LowerTriangular #empty matrix
        gsrut = GSRUTWorkspace(; N_u = N, N_z = N, N_w = 0)

        f_noiseless! = (z, u, _) -> @. z = 2*u + 1

        transform!(gsrut; z̄, S_δz, P_δuδz, ū, S_δu, S_w, f! = f_noiseless!)

        P_δz = S_δz * S_δz'
        @test z̄ ≈ 2ū .+ 1
        @test P_δz ≈ 4*P_δu
        @test P_δuδz ≈ 2I

    end

    #non-linear transformation on a random variable with non-Euclidean manifold
    @testset verbose = true "Manifold" begin

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
            function (x1::AbstractVector{<:Real}, x0::AbstractVector{<:Real},
                    δx::AbstractVector{<:Real})

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
            function (δx::AbstractVector{<:Real}, x1::AbstractVector{<:Real},
                    x0::AbstractVector{<:Real})

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

        #check consistency between generalized addition and subtraction with
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

        #preallocate arrays for the GSRUT output
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
            function (x1::AbstractVector{<:Real}, x0::AbstractVector{<:Real},
                 w::AbstractVector{<:Real})

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

        #apply the nonlinear transformation to the initial mean with zero noise
        x̄1_direct = copy(x̄1)
        f_x!(x̄1_direct, x̄0, zeros(6))

        #create GSRUT workspace
        gsrut = GSRUTWorkspace(N_u = length(x_template), N_δu = length(δx_template),
                            N_z = length(x_template), N_δz = length(δx_template),
                            N_w = length(w_template), ε_mean = 1e-12)

        #apply the GSRUT with zero noise square-root covariance
        S_w0 = zeros(6,6) |> LowerTriangular

        transform!(gsrut;
            z̄ = x̄1, S_δz = S_δx1, P_δuδz = P_δx0δx1,
            ū = x̄0, S_δu = S_δx0, S_w = S_w0,
            f! = f_x!,
            g_u_plus! = g_x_plus!,
            g_z_plus! = g_x_plus!,
            g_z_minus! = g_x_minus!)

        #the transformed mean should be equal to that from the noiseless nonlinear
        #transformation applied directly to the mean
        @test all(isapprox(x̄1_direct, x̄1))

        #also, with zero noise square-root covariance, the error state square-root
        #covariance should be preserved (minus a possible sign inversion)
        @test all(isapprox(abs.(S_δx0), abs.(S_δx1)))

        #now apply the GSRUT with non-zero noise square-root covariance
        σ_w0 = copy(w_template)
        σ_w0.ρ .= 5e-2
        σ_w0.δv .= 3e-2
        S_w0 = LowerTriangular(ComponentMatrix(diagm(σ_w0), w_axes, w_axes))

        transform!(gsrut;
            z̄ = x̄1, S_δz = S_δx1, P_δuδz = P_δx0δx1,
            ū = x̄0, S_δu = S_δx0, S_w = S_w0,
            f! = f_x!,
            g_u_plus! = g_x_plus!,
            g_z_plus! = g_x_plus!,
            g_z_minus! = g_x_minus!)

        #although the noise is zero-mean, because of the non-linear fashion it acts
        #on the quaternion part of the state, we should not expect the transformed
        #mean quaternion to be the same as in the direct noiseless transformation
        @test !all(isapprox(x̄1_direct.q[:], x̄1.q[:]))

        #however, since it does act linearly on the Euclidean part of the state, the
        #transformed mean for that part should indeed be the same
        @test all(isapprox(x̄1_direct.v, x̄1.v))

        #what should be true for both parts is that the main diagonal of S should
        #have increased (allowing for a sign ambiguity in S)
        @test all(abs.(diag(S_δx1)) .>= abs.(diag(S_δx0)))

        #check for allocations
        @test (@ballocated transform!($gsrut;
            z̄ = $x̄1, S_δz = $S_δx1, P_δuδz = $P_δx0δx1,
            ū = $x̄0, S_δu = $S_δx0, S_w = $S_w0,
            f! = $f_x!,
            g_u_plus! = $g_x_plus!,
            g_z_plus! = $g_x_plus!,
            g_z_minus! = $g_x_minus!)) == 0

    end #testset

end #function



function test_stages()

    @testset verbose = true "Predictor" begin

        N = 3
        pred = GSRUKF.Predictor(N_x = N, N_w = N)

        f! = (x1, x0, w0) -> @. x1 = 2*x0 + 1 + w0

        x̂0 = ones(N)
        S_δx0 = 1.0 * LowerTriangular(Matrix(I, N, N))
        S_w = 1.0 * LowerTriangular(Matrix(I, N, N))

        x̂1 = copy(x̂0)
        S_δx1 = copy(S_δx0)

        P_δx0 = S_δx0 * S_δx0'
        P_w = S_w * S_w'
        GSRUKF.predict!(pred; x̂1, S_δx1, x̂0, S_δx0, S_w, f!)

        #we know this linear transformation must yield the following
        P_δx1 = S_δx1 * S_δx1'
        @test x̂1 ≈ 2x̂0 .+ 1
        @test P_δx1 ≈ 4*P_δx0 .+ P_w

        # check for allocations
        @test (@ballocated GSRUKF.predict!( $pred;
            x̂1=$x̂1, S_δx1=$S_δx1, x̂0=$x̂0,
            S_δx0=$S_δx0, S_w=$S_w, f! =$f!)) == 0

        # @btime GSRUKF.predict!( $pred;
        #     x̄1=$x̄1, S_δx1=$S_δx1, x̄0=$x̄0,
        #     S_δx0=$S_δx0, S_wp0=$S_wp0, f! =$f!)
    end

    @testset verbose = true "Updater" begin

        N_x = N_δx = 3
        N_y = 2
        N_w = 2

        updater = GSRUKF.Updater(; N_x, N_y, N_w)

        function h!(y, x, w)
            y[1] = x[1] + w[1]
            y[2] = x[2] + w[2]
        end

        function g_y_minus!(δy, y1, y0)
            δy .= y1 .- y0
        end

        x̂m = ones(N_x)
        S_δxm = 2.0 * LowerTriangular(Matrix(I, N_δx, N_δx))
        S_w = LowerTriangular(diagm([1, 2]))

        x̂p = copy(x̂m)
        S_δxp = copy(S_δxm)

        σ_thr = 3

        #test measurement rejection
        ỹ = [10.0, 2.0]
        @test !GSRUKF.update!(updater; x̂p, S_δxp, x̂m, S_δxm, S_w, ỹ, h!, σ_thr)
        @test any(updater.δη .> 3)

        #test measurement acceptance
        ỹ = @SVector[1.1, 1.1]
        @test GSRUKF.update!(updater; x̂p, S_δxp, x̂m, S_δxm, S_w, ỹ, h!, σ_thr)

        P_δxm = S_δxm * S_δxm
        P_δxp = S_δxp * S_δxp'

        #for the states included in the measurement σ must decrease, for the other
        #one σ must remain unchanged
        @test P_δxp[1,1] < P_δxm[1,1]
        @test P_δxp[2,2] < P_δxm[2,2]
        @test P_δxp[3,3] ≈ P_δxm[3,3]

        #measurement noise is smaller for the first state than for the second,
        #so its σ must decrease more and its measurement residual must be
        #smaller
        ŷp = zeros(N_y)
        h!(ŷp, x̂p, zeros(N_x))
        δỹp = ỹ - ŷp #measurement residual
        @test P_δxp[1,1] < P_δxp[2,2]
        @test abs(δỹp[1]) < abs(δỹp[2])

        @test (@ballocated !GSRUKF.update!($updater; x̂p=$x̂p, S_δxp=$S_δxp, x̂m=$x̂m,
            S_δxm=$S_δxm, S_w=$S_w, ỹ=$ỹ, h! = $h!, σ_thr=$σ_thr)) == 0

        @btime GSRUKF.update!($updater; x̂p=$x̂p, S_δxp=$S_δxp, x̂m=$x̂m,
            S_δxm=$S_δxm, S_w=$S_w, ỹ=$ỹ, h! = $h!, σ_thr=$σ_thr)

    end

end #function

end #module