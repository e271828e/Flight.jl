module TestSystem

using Flight.LBV
using Flight.System

export XPwp, UPwp, PwpSystem, XLdg, ULdg, LdgSystem, get_x_type, get_u_type
export @testmacro, test_macro

const SA_Float64 = SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}

@define_node XPwp (ω = LBVLeaf{1}, β = LBVLeaf{1})
@define_node UPwp (throttle = LBVLeaf{1}, prop_pitch = LBVLeaf{1})

struct PwpSystem{X<:XPwp{Float64}, U<:UPwp{Float64}} <: AbstractContinuousSystem{X, U}
	x::X
	u::U
	J::Vector{Float64}
end

get_x_type(::Type{<:PwpSystem}) = XPwp{Float64}
get_u_type(::Type{<:PwpSystem}) = UPwp{Float64}
get_x_type(::PwpSystem) = get_x_type(PwpSystem)
get_u_type(::PwpSystem) = get_u_type(PwpSystem)

@define_node XLdg NamedTuple{(:a2, :b2)}((LBVLeaf{4}, LBVLeaf{2}))
@define_node ULdg NamedTuple{(:u2,)}((LBVLeaf{4},))

struct LdgSystem{X<:XLdg{Float64}, U<:ULdg{Float64}} <: AbstractContinuousSystem{X, U}
	x::X
	u::U
	ξ::Float64
end

get_x_type(::Type{<:LdgSystem}) = XLdg{Float64}
get_u_type(::Type{<:LdgSystem}) = ULdg{Float64}
get_x_type(::LdgSystem) = get_x_type(LdgSystem)
get_u_type(::LdgSystem) = get_u_type(LdgSystem)

macro testmacro(sys_name, subsystems)

	#given Symbol sys_name and NamedTuple subsystems whose values are systems
	#for which methods get_x_type and get_u_type have been defined, construct
	#the LBVNode for the state vector x and the input vector u of the joint
	#system sys_name
	x_name = Symbol(:X, sys_name)
	u_name = Symbol(:U, sys_name)
	println(x_name)
	println(u_name)

    ex = quote
		let
        ss = $(subsystems)
		labels = keys(ss)

        x_types = get_x_type.(values(ss))
		x_descriptor = NamedTuple{labels}(x_types)
		println(x_descriptor)
		@define_node $x_name x_descriptor

        u_types = get_u_type.(values(ss))
		u_descriptor = NamedTuple{labels}(u_types)
		println(u_descriptor)
		@define_node $u_name u_descriptor
		end
    end
	#need to escape the whole macro, otherwise XParent is mangled
	#even if we escape the whole return expression, we are safe, because all the
	#temporary variables used by @define_node are created within a let block!
	#only the intended methods are exported
    return esc(Base.remove_linenums!(ex))
    # return Base.remove_linenums!(ex)

end

# @testmacro XPwp (left = System1(), right = System2())



#create a lbv implementation in which the whole macro is escaped. hygiene is
#guaranteed by the let block


end