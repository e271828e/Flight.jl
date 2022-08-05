#each of these functions needs its own aircraft instance to compute its outputs,
#since it will mutate it in the process, butchering trim values that other
#functions rely on

#we should define xdot0, x0, u0, y0
#for u0 and y0, we define their corresponding ComponentVector subtypes LinU and
#LinY
#for convenience, we also define AircraftU as typeof(ac.u)
#Then we write functions that map the aircraft's u to this reduced

UTemplate = ComponentVector()
YTemplate = ComponentVector()

function assign!(u::U, ac::System{<:Cessna172R}) end
function assign!(y::Y, ac::System{<:Cessna172R}) end

function assign!(ac::System{<:Cessna172R}, y::Y) end
function assign!(ac::System{<:Cessna172R}, u::U) end #LinU is a ComponentVector

ac = System(Cessna172R())
state = C172R.Trim.State()
params = C172R.Trim.Parameters()
state_opt = C172R.Trim.trim!(ac, env, state, params)
assign!(ac, env, state_opt, params) #set the aircraft to its trim state
ẋ0 = copy(ac.ẋ)
x0 = copy(ac.x)
u0 = copy(UTemplate); assign!(u0, ac) #get the equilibrium value from the trimmed aircraft
y0 = copy(YTemplate); assign!(y0, ac) #idem

# function assign!(ac::System{<:Cessna172R}, x, u)
#     ac.x .= x #this can be assigned directly, because LinX is of the same type as ac.x
#     assign!(ac, u) #assign the input u to the aircraft
#     f_cont!(ac, env) #this also modifies ac.y, but we don't care about that here
# end

# #the same aircraft instance is shared by g and h, so we need to make sure we
# #start with everything at its trimmed condition
# g! = let ac = ac, env = env, state = state, params = params
# 	function (ẋ, x, u)
#         assign!(ac, env, state, params) #restore the aircraft to its trimmed status
#         assign!(ac, x, u) #assign the given state and input to the aircraft and update its ẋ and y
# 		ẋ .= ac.xdot #assign the aircraft's ẋ to the provided ẋ
#     end
# end

# h! = let ac = ac, env = env, state = state, params = params
# 	function (y, x, u)
#         assign!(ac, env, state, params) #restore the aircraft to its trimmed status
#         assign!(ac, x, u) #assign the given state and input to the aircraft and update its ẋ and y
# 		assign!(y, ac) #assign the aircraft's y to the provided y
#     end
# end

#doubly mutating function mirroring f_cont!(), which internally mutates the
#system's ẋ and y
f! = let ac = ac, env = env, state = state, params = params
	function (ẋ, y, x, u)
        assign!(ac, env, state, params) #restore the aircraft to its trimmed status
        ac.x .= x #this can be assigned directly, because LinX is of the same type as ac.x
        assign!(ac, u) #assign the input vector u to the aircraft
        f_cont!(ac, env) #this also modifies ac.y, but we don't care about that here
		ẋ .= ac.xdot #assign the aircraft's ẋ to the provided ẋ
		assign!(y, ac) #assign the aircraft's y to the provided y
    end
end

ẋ_tmp = copy(ẋ0)
y_tmp = copy(y0)

f!(ẋ_tmp, y, x0, u0)
@assert ẋ_tmp ≈ ẋ0 #sanity check
@assert y_tmp ≈ y0 #sanity check

# # g!(ẋ, x0, u0)
# # h!(y, x0, u0)
# g_a! = let u = u0
#     (ẋ, x) = g!(ẋ, x, u)
# end

# g_b! = let x = x0
#     (ẋ, u) = g!(ẋ, x, u)
# end

f_a! = let u = u0, y = y_tmp
    (ẋ, x) = f!(ẋ, y, x, u)
end

f_b! = let x = x0, y = y_tmp
    (ẋ, u) = f!(ẋ, y, x, u)
end

f_c! = let u = u0, y = y_tmp
    (y, x) = f!(ẋ, y, x, u)
end

f_d! = let x = x0, ẋ = ẋ_tmp
    (y, u) = f!(ẋ, y, x, u)
end

f_a!(ẋ_tmp, x0)
@assert ẋ_tmp ≈ ẋ0 #sanity check
f_b!(ẋ_tmp, u0)
@assert ẋ_tmp ≈ ẋ0 #sanity check
f_c!(y_tmp, x0)
@assert y_tmp ≈ y0 #sanity check
f_d!(y_tmp, u0)
@assert y_tmp ≈ y0 #sanity check

B = x0 * u0' #nx x nu
C = y0 * x0'
D = y0 * u0'
A = x0 * x0' #nx x nx
finite_difference_jacobian!(A, f_a!, x0)
finite_difference_jacobian!(B, f_b!, u0)
finite_difference_jacobian!(C, f_c!, x0)
finite_difference_jacobian!(D, f_d!, u0)
