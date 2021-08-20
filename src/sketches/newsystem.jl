using Flight
using OrdinaryDiffEq
using ComponentArrays

#we need the Component type parameter for dispatch, and the rest for type stability
struct NewSystem{C, X, Y, U, P, S}
    dx::X
    x::X
    y::Y
    u::U
    params::P
    subsystems::S
end

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
    params = thr #params is the component itself
    subsystems = nothing #no subsystems to define
    NewSystem{map(typeof, (thr, x, y, u, params, subsystems))...}(dx, x, y, u, params, subsystems)
end

function f_cont!(sys::NewSystem{NewEThruster})
    sys.dx.ω_shaft = 0.1
    sys.dx.c_bat = -0.1
    sys.y.throttle = sys.u.throttle
    sys.y.ω_shaft = sys.x.ω_shaft
    sys.y.c_bat = sys.x.c_bat
end

get_wr_Ob_b(sys::NewSystem{NewEThruster}) = sys.y.wr_Ob_b
get_h_Gc_b(sys::NewSystem{NewEThruster}) = sys.y.h_Gc_b

####################### Dual Thruster #######################

#NEXT STEPS:
#1) Generalize to any NewAbstractComponent. maybe add the Component type as a
#   type parameter (not really necessary)
#2) Change ComponentGroup methods to @generated
#3) comment out everything affected by Component, System and Model, and update
#   their implementations with the new one. replace the propulsion module with
#   the new one. test the new Model and System with ComponentGroup containing
#   the new EThruster
#4) create a new branch new_arch_with_struct_y to change back to nested struct
#   arrays for Y, everything still in newsystem.jl. Remove field y from
#   NewSystem, and modify f_cont!() to return y. Modify NewContinuousModel to
#   account for all this. benchmark and compare with the committed
#   implementation

EThrusterGroup{L} = NamedTuple{L, T} where {L, T <: NTuple{N, NewEThruster} where {N}}
#now generalize this, but restricting it to New AbstractComponent

X(g::EThrusterGroup) = ComponentVector(NamedTuple{keys(g)}(X.(values(g))))
Y(g::EThrusterGroup) = ComponentVector(NamedTuple{keys(g)}(Y.(values(g))))
U(g::EThrusterGroup) = NamedTuple{keys(g)}(U.(values(g)))

function NewSystem(g::EThrusterGroup{L}, dx = X(g), x = X(g), y = Y(g), u = U(g)) where {L}
    #having L allows us to know the length of g and therefore the number of
    #expressions we need to generate
    subs = Dict{Symbol, NewSystem}()
    for label in L
        sys_comp = NewSystem(map((λ)->getproperty(λ, label), (g, dx, x, y, u))...)
        push!(subs, label => sys_comp)
    end
    params = nothing #everything is already stored in the subsystem's parameters
    subsystems = NamedTuple{Tuple(keys(subs))}(values(subs))
    NewSystem{map(typeof, (g, x, y, u, params, subsystems))...}(dx, x, y, u, params, subsystems)
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



################### MORE INSIGHT #################
#=
#if instead of having a ComponentVector Y as the only output, i define a
#discrete Yd and a continuous Y and pack them in a struct to be saved by DiffEq
#on each time step, the resulting log would no longer be a Vector of
#ComponentVectors, but a Vector of structs. this means i would have to extract
#both y and yd for each t and assemble them into two separate VectorOfArrays.
#that is ugly to do by hand. StructArrays to the rescue!!!

Base.@kwdef struct Output{Y, YD} #this is the type we pass to SavedValues
    y::Y = ComponentVector(a = fill(1.0, 2), b = fill(2.0, 3))
    yd::YD = ComponentVector(m = fill(4,3), n = fill(-2, 2))
end
log = collect(Output() for i in 1:5)
sa = StructArray(log) #now all the y fields of log lie in a contiguous array
y = sa.y
y_voa = VectorOfArray(y) #now we have a vector of Y's that indexes like a matrix
y_mat = convert(Array, y_voa)

# with this, Y no longer needs to be an array! we could use nested immutable
#structs or NamedTuples thereof. because StructArrays not only accepts Vectors
#of structs as inputs but also Vectors of NamedTuples. and, if we keep
#everything immutable, it no longer has to be inefficient to allocate it. the
#only problem is doing a deepcopy on each step.

# but NO! this does not have to be. we can simply define y to be the output of
#f_cont!(sys). being an immutable struct freshly generated from scratch on every
#iteration, it can go directly to saveval without deepcopying, because it is NOT the
#same object as the one in the next step. #in fact, with this approach, we no
#longer need to define y as a System field!!

#ok, but what happens if we store things like RQuat or Wrench objects? Vectors,
#and therefore cannot be handled by the VectorOfArray workflow. on the other
#hand, it doesn't have to. because we can apply StructArray recursively. for
#example, if we get an array y_kin of YKin structs, each of them containing a
#RQuat q_lb and another RQuat q_el, first we do sa = StructArray(y_kin), then
#v_q_lb = sa.q_lb. now we have a Vector of RQuats. yep, but the pain does not
#stop here, because then RQuat holds a UnitQuat, which in turn holds a Quat.
#this is not practical. so here is a self-imposed limitation to facilitate
#simulation output post-processing: only use fields of SVectors of Floats, Ints,
#Bools or non-deeply nested structs. A Wrench is acceptable because it is
#handled directly by StructArrays to extract F and M. a WGS84 pos maybe too.

#advantages of the nested struct approach:
#1) we can hold different data types in a single struct, so we no longer need
#   y and yd
#2) for postprocessing we can use StructArrays to extract a Vector of outputs
#   corresponding to a subsystem for plot dispatch, so we're still good there
#3) we don't need to make Wrench a ComponentVector, which recovers some
#   performance
#4) because the overall y no longer needs to be allocated we don't need Y
#   methods, the strct is directly instantiated in the fcont call
#   we only need to define Y structs for Leaf components. in most other cases
#   they will be compose using NTs
#4) we can use ComponentGroups and other forms of composition without worrying
#   about having to deal a particular component returning nothing for y. a NT of
#   nothing is a NT!

#disadvantages:
#1) although they are stack allocated, these nested structs have still some
#   allocation time. but this is probably comparable or even less than the cost
#   of in-place assignment to the heap-preallocated yblocks
#2) we're still restricted to storing SVectors, Floats and such things to make
#   post-processing easier


#nah, it IS better to have y as a preallocated ComponentVector with views stored
#throughout the hierarchy. especially if we can define another y_d for discrete
#outputs if required.


=#