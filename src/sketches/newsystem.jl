using Flight
using OrdinaryDiffEq
using ComponentArrays


struct NewSystem{C, X, Y, U, P, S}
    dx::X
    x::X
    y::Y
    u::U
    params::P
    subsystems::S
end
# struct NewSystem{C <: NewAbstractComponent}
#     dx::Any
#     x::Any
#     y::Any
#     u::Any
#     params::Any
#     subsystems::Any
#     function NewSystem{C}(dx, x, y, u, params, subsystems) where {C}
#         new{C}(dx, x, y, u, params, subsystems)
#     end
# end

abstract type NewAbstractComponent end

abstract type NewAbstractThruster <: NewAbstractComponent end

Base.@kwdef struct NewEThruster <: NewAbstractThruster
    frame::Frame = Frame()
    battery::Battery = Battery()
    motor::ElectricMotor = ElectricMotor()
    gearbox::Gearbox = Gearbox()
    propeller::SimpleProp = SimpleProp()
end

const NewEThrusterXTemplate = ComponentVector(ω_shaft = 0.0, c_bat = 1.0)
const NewEThrusterYTemplate = ComponentVector(
    throttle = 0.0, ω_shaft = 0.0, ω_prop = 0.0, i = 0.0, c_bat = 1.0,
    wr_Oc_c = ComponentVector(Wrench()), wr_Ob_b = ComponentVector(Wrench()),
    h_Gc_b = zeros(3))

const NewEThrusterX{T, D} = ComponentVector{T, D, typeof(getaxes(NewEThrusterXTemplate))} where {T,D}
const NewEThrusterY{T, D} = ComponentVector{T, D, typeof(getaxes(NewEThrusterYTemplate))} where {T,D}
Base.@kwdef mutable struct NewEThrusterU
    throttle::Float64 = 0.0
end

X(::NewEThruster) = copy(NewEThrusterXTemplate)
Y(::NewEThruster) = copy(NewEThrusterYTemplate)
U(::NewEThruster) = NewEThrusterU()

function NewSystem(thr::NewEThruster, dx::NewEThrusterX = X(thr), x::NewEThrusterX = X(thr),
    y::NewEThrusterY = Y(thr), u::NewEThrusterU = U(thr))
    params = thr
    subsystems = nothing
    NewSystem{typeof(thr), typeof(x), typeof(y), typeof(u), typeof(params), typeof(subsystems)}(dx, x, y, u, params, subsystems)
end

# function f_cont!(sys::NewSystem{NewEThruster, X, Y, U, P, S} where {X,Y,U,P,S})
function f_cont!(sys::NewSystem{NewEThruster})
    sys.dx.ω_shaft = 0.1
    sys.dx.c_bat = -0.1
    sys.y.throttle = sys.u.throttle
    sys.y.ω_shaft = sys.x.ω_shaft
    sys.y.c_bat = sys.x.c_bat
end

get_wr_Ob_b(sys::NewSystem{NewEThruster}) = sys.y.wr_Ob_b
get_h_Gc_b(sys::NewSystem{NewEThruster}) = sys.y.h_Gc_b


#question: we need the C type parameter for dispatch. but do we really need to
#include X, Y, U, P, S as type parameters for type stability? if the compiler
#evaluates all it can at compile time, it may not be required

####################### Dual Thruster #######################

EThrusterGroup{L} = NamedTuple{L, T} where {L, T <: NTuple{N, NewEThruster} where {N}}

X(g::EThrusterGroup) = ComponentVector(NamedTuple{keys(g)}(X.(values(g))))
Y(g::EThrusterGroup) = ComponentVector(NamedTuple{keys(g)}(Y.(values(g))))
U(g::EThrusterGroup) = NamedTuple{keys(g)}(U.(values(g)))

function NewSystem(g::EThrusterGroup{L}, dx = X(g), x = X(g), y = Y(g), u = U(g)) where {L}
    #having L allows us to know the length of g and therefore the number of
    #expressions we need to generate
    subs = Dict{Symbol, NewSystem}()
    for label in L
        comp = getproperty(g, label)
        dx_comp = getproperty(dx, label); x_comp = getproperty(x, label);
        y_comp = getproperty(y, label); u_comp = getproperty(u, label)
        sys_comp = NewSystem(comp, dx_comp, x_comp, y_comp, u_comp)
        push!(subs, label => sys_comp)
    end
    params = nothing #everything is already stored in the subsystem's parameters
    subsystems = NamedTuple{Tuple(keys(subs))}(values(subs))
    NewSystem{typeof(g), typeof(x), typeof(y), typeof(u), typeof(params), typeof(subsystems)}(dx, x, y, u, params, subsystems)
end


#sys::NewSystem{EThrusterGroup} would not work here due to the non-covariance of
#the type system
function get_wr_Ob_b(sys::NewSystem{T}) where {T<:EThrusterGroup{L}} where {L}
    #to generate the unrolled expressions for the @generated function, we only
    #have access to the type, not the values. so we can't do length(g). we must
    #keep the length or the labels as a type parameter
    wr_Ob_b = Wrench()
    for label in L #having L allows us to know the length of g and therefore the number of expressions we need to generate
        wr_Ob_b .+= sys.subsystems[label] |> get_wr_Ob_b
    end
    return wr_Ob_b
end

function get_h_Gc_b(sys::NewSystem{T}) where {T<:EThrusterGroup{L}} where {L}
    h_Gc_b = SVector{3,Float64}(0,0,0)
    for label in L
        h_Gc_b .+= sys.subsystems[label] |> get_h_Gc_b
    end
    return h_Gc_b
end

function f_cont!(sys::NewSystem{T}) where {T<:EThrusterGroup{L}} where {L}
    for label in L #having L allows us to know the length of g and therefore the number of expressions we need to generate
        f_cont!(getproperty(sys.subsystems, label))
    end
end


####################### Dual Thruster #######################

# Base.@kwdef struct DualEThruster <: NewAbstractComponent #could be generated via macro
#     left::NewEThruster = NewEThruster()
#     right::NewEThruster = NewEThruster()
# end

# X(g::DualEThruster) = ComponentVector(left = X(g.left), right)
# Y(g::DualEThruster) = ComponentVector(NamedTuple{keys(g)}(Y.(values(g))))
# U(g::DualEThruster) = NamedTuple{keys(g)}(U.(values(g)))

# function NewSystem(g::EThrusterGroup{L}, dx = X(g), x = X(g), y = Y(g), u = U(g)) where {L}
#     #having L allows us to know the length of g and therefore the number of
#     #expressions we need to generate
#     subs = Dict{Symbol, NewSystem}()
#     for label in L
#         comp = getproperty(g, label)
#         dx_comp = getproperty(dx, label)
#         x_comp = getproperty(x, label)
#         y_comp = getproperty(y, label)
#         u_comp = getproperty(u, label)
#         sys_comp = NewSystem(comp, dx_comp, x_comp, y_comp, u_comp)
#         push!(subs, label => sys_comp)
#     end
#     params = nothing #all that we need is stored in the subsystem's parameters
#     subsystems = NamedTuple{Tuple(keys(subs))}(values(subs))
#     NewSystem{typeof(g), typeof(x), typeof(y), typeof(u), typeof(params), typeof(subsystems)}(dx, x, y, u, params, subsystems)
# end

# function get_wr_Ob_b(sys::NewSystem{EThrusterGroup{L}}) where {L}
#     #to generate the unrolled expressions for the @generated function, we only
#     #have access to the type, not the values. so we can't do length(g). we must
#     #keep the length or the labels as a type parameter
#     wr_Ob_b = Wrench()
#     for label in L #having L allows us to know the length of g and therefore the number of expressions we need to generate
#         wr_Ob_b .+= sys.subsystems[label] |> get_wr_Ob_b
#     end
# end
# function get_h_Gc_b(sys::NewSystem{EThrusterGroup{L}}) where {L}
#     h_Gc_b = SVector{3,Float64}(0,0,0)
#     for label in L
#         h_Gc_b .+= sys.subsystems[label] |> get_h_Gc_b
#     end
# end

# function f_cont!(sys::NewSystem{EThrusterGroup{L}}) where {L}
#     for label in L #having L allows us to know the length of g and therefore the number of expressions we need to generate
#         f_cont!(getproperty(sys.subsystems, label))
#     end
# end