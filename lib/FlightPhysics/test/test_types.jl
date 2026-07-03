module TestTypes

using Test
using JSON3, StructTypes

using FlightCore.Joysticks: exp_axis_curve
using FlightPhysics.Types

export test_types

#JSON3 target structs must be defined at the top level
@kwdef mutable struct MutableTarget
    a::Ranged{Float64, 0.0, 1.0} = Ranged(0.5, 0.0, 1.0)
    b::Bool = false
end
StructTypes.StructType(::Type{MutableTarget}) = StructTypes.Mutable()

@kwdef struct ImmutableTarget
    a::Ranged{Float64, 0.0, 1.0} = Ranged(0.5, 0.0, 1.0)
    b::Bool = false
end
StructTypes.StructType(::Type{ImmutableTarget}) = StructTypes.Struct()


function test_types()
    @testset verbose = true "Types" begin
        @testset verbose = true "Ranged" begin test_Ranged() end
    end
end

function test_Ranged()

    a = Ranged(0.5, 0., 1.)
    b = Ranged(0.5, 0., 2.) #same value, different bounds
    c = Ranged(0.7, -1., 1.)

    @testset "Constructors" begin

        #within bounds, clamped from above, clamped from below
        @test Ranged(0.5, 0., 1.).val == 0.5
        @test Ranged(1.5, 0., 1.).val == 1.0
        @test Ranged(-0.5, 0., 1.).val == 0.0

        #bounds are converted to the value's type
        @test Ranged(0.5, 0, 1) === a

        #fully parameterized constructor clamps
        R = Ranged{Float64, 0.0, 1.0}
        @test R(0.25) === Ranged(0.25, 0., 1.)
        @test R(3) === Ranged(1.0, 0., 1.)

        @test typemin(R) == 0.0
        @test typemax(R) == 1.0
        @test typemin(a) == 0.0
        @test typemax(a) == 1.0

    end

    @testset "Conversions" begin

        r = Ranged(0.75, 0., 1.)
        @test convert(Float64, r) === 0.75
        @test Float64(r) === 0.75
        @test Float32(r) === 0.75f0

        #conversion to different bounds re-clamps
        @test convert(Ranged{Float64, 0.0, 0.5}, r).val == 0.5

        #if the target does not specify bounds, they carry over from the source
        rf = Ranged{Float32}(r)
        @test rf.val === 0.75f0
        @test typemax(rf) === 1.0f0

    end

    @testset "Comparisons" begin

        #against Reals, both argument orders
        @test a == 0.5 && 0.5 == a
        @test a ≈ 0.5 + 1e-12 && 0.5 + 1e-12 ≈ a
        @test a < 0.7 && 0.3 < a
        @test a <= 0.5 && a >= 0.5

        #comparisons are value-based even across distinct bounds
        @test a == b
        @test a != Ranged(0.6, 0., 2.)
        @test b ≈ a
        @test a < c
        @test isless(a, c)

        #hash must be consistent with value-based equality
        @test hash(a) == hash(0.5)
        @test hash(a) == hash(b)
        @test 0.5 in Set{Any}([a])
        @test length(unique(Any[a, b, 0.5])) == 1

    end

    @testset "Arithmetic" begin

        #Ranged and Real, result clamped to the Ranged operand's bounds
        @test a + 0.3 === Ranged(0.8, 0., 1.)
        @test a + 0.7 === Ranged(1.0, 0., 1.)
        @test a - 1.0 === Ranged(0.0, 0., 1.)

        #negation preserves the bounds, and therefore may clamp
        @test -c === Ranged(-0.7, -1., 1.)
        @test -a === Ranged(0.0, 0., 1.)

        #Ranged pairs require identical bounds
        @test a + Ranged(0.25, 0., 1.) === Ranged(0.75, 0., 1.)
        @test a - Ranged(0.25, 0., 1.) === Ranged(0.25, 0., 1.)
        @test_throws MethodError a + b
        @test_throws MethodError a - b

    end

    @testset "Saturation & Scaling" begin

        @test saturation(Ranged(1.2, 0., 1.)) == 1
        @test saturation(Ranged(-0.2, 0., 1.)) == -1
        @test saturation(a) == 0

        @test linear_scaling(Ranged(0.0, 0., 1.), (10.0, 20.0)) == 10.0
        @test linear_scaling(Ranged(0.5, 0., 1.), (10.0, 20.0)) == 15.0
        @test linear_scaling(Ranged(1.0, 0., 1.), 10:2:20) == 20.0
        @test linear_scaling(Ranged(0.0, -1., 1.), (-30.0, 30.0)) == 0.0

        #Joysticks integration: forwards the underlying value
        r = Ranged(0.4, -1., 1.)
        @test exp_axis_curve(r; strength = 1.0, deadzone = 0.1) ==
              exp_axis_curve(0.4; strength = 1.0, deadzone = 0.1)

    end

    @testset "JSON3" begin

        #numeric value into the Ranged field of a mutable struct
        u = MutableTarget()
        JSON3.read!("""{"a": 0.1209}""", u)
        @test u.a === Ranged(0.1209, 0., 1.)

        #out-of-bounds input clamps
        JSON3.read!("""{"a": 7.5, "b": true}""", u)
        @test u.a === Ranged(1.0, 0., 1.)
        @test u.b

        #construction of an immutable struct with a Ranged field
        y = JSON3.read("""{"a": 0.3, "b": true}""", ImmutableTarget)
        @test y.a === Ranged(0.3, 0., 1.)

        #a Ranged lowers to its wrapped value
        @test JSON3.write(Ranged(0.25, 0., 1.)) == "0.25"

    end

end

end #module
