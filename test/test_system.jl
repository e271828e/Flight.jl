module TestSystem

using Flight.LBV
# using Flight.System

export XSystem1, XSystem2, USystem1, USystem2, System1, System2, get_x_type, get_u_type
export @testmacro, test_macro

@define_node XSystem1 (a1 = LBVLeaf{2}, b1 = LBVLeaf{3})
@define_node USystem1 (u1 = LBVLeaf{2},)

struct System1
	x::XSystem1
	x_dot::XSystem1
	u::USystem1
end
System1() = System1(XSystem1(), XSystem1(), USystem1())

get_x_type(::Type{<:System1}) = XSystem1
get_u_type(::Type{<:System1}) = USystem1
get_x_type(::System1) = get_x_type(System1)
get_u_type(::System1) = get_u_type(System1)


@define_node XSystem2 NamedTuple{(:a2, :b2)}((LBVLeaf{4}, LBVLeaf{2}))
@define_node USystem2 NamedTuple{(:u2,)}((LBVLeaf{4},))

struct System2
	x::XSystem2
	x_dot::XSystem2
	u::USystem2
end
System2() = System2(XSystem2(), XSystem2(), USystem2())

get_x_type(::Type{<:System2}) = XSystem2
get_u_type(::Type{<:System2}) = USystem2
get_x_type(::System2) = get_x_type(System2)
get_u_type(::System2) = get_u_type(System2)

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